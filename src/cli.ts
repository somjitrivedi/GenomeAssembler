#!/usr/bin/env node
import { resolve } from "node:path";
import { parseArgs } from "node:util";
import { countCommand, evaluateCommand, simulateCommand, statsCommand } from "./commands.js";
import { openDatabase } from "./db/database.js";
import { Repository } from "./db/repository.js";
import { bars, bytes, duration, fixed, int, percent, table } from "./format.js";

const USAGE = `genome-assembler <command> [options]

Commands
  simulate   Generate a random genome with repeats and paired-end reads
             --out <dir> (required) --length 1000000 --gc 0.5
             --repeats 0 --repeat-length 2000 --repeat-copies 2 --repeat-divergence 0
             --coverage 30 --read-length 150 --insert 400 --insert-sd 40
             --error-rate 0.005 --n-rate 0 --seed 42 --gzip
  stats      Read statistics for FASTA/FASTQ files (gzip allowed)
             <files...>
  count      Count canonical k-mers, filter singletons with a Bloom filter,
             estimate coverage and genome size from the k-mer spectrum
             <files...> -k 31 --mode bloom|naive --fpr 0.01 --expected <n>
             --approx --min-count <n> --out solid.kmers --histogram hist.tsv
  evaluate   Assembly statistics; with a reference also genome fraction and
             misassemblies
             --contigs <fasta> (required) --reference <fasta> --min-contig 500 -k 31
  runs       List recent runs        [--limit 20]
  show       Show one run as JSON    <run id>

Global
  --db <path>   SQLite database (default: $GA_DB or ./genome-assembler.db)
  -h, --help`;

function num(value: string | undefined, fallback: number, name: string): number {
  if (value === undefined) return fallback;
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) throw new Error(`--${name} must be a number, got '${value}'`);
  return parsed;
}

function optionalNum(value: string | undefined, name: string): number | undefined {
  return value === undefined ? undefined : num(value, 0, name);
}

async function main(argv: string[]): Promise<void> {
  const [command, ...rest] = argv;
  if (!command || command === "-h" || command === "--help" || command === "help") {
    console.log(USAGE);
    return;
  }
  const { values, positionals } = parseArgs({
    args: rest,
    allowPositionals: true,
    options: {
      db: { type: "string" },
      help: { type: "boolean", short: "h" },
      out: { type: "string" },
      length: { type: "string" },
      gc: { type: "string" },
      repeats: { type: "string" },
      "repeat-length": { type: "string" },
      "repeat-copies": { type: "string" },
      "repeat-divergence": { type: "string" },
      coverage: { type: "string" },
      "read-length": { type: "string" },
      insert: { type: "string" },
      "insert-sd": { type: "string" },
      "error-rate": { type: "string" },
      "n-rate": { type: "string" },
      seed: { type: "string" },
      gzip: { type: "boolean" },
      k: { type: "string", short: "k" },
      mode: { type: "string" },
      fpr: { type: "string" },
      expected: { type: "string" },
      approx: { type: "boolean" },
      "min-count": { type: "string" },
      histogram: { type: "string" },
      contigs: { type: "string" },
      reference: { type: "string" },
      "min-contig": { type: "string" },
      limit: { type: "string" }
    }
  });
  if (values.help) {
    console.log(USAGE);
    return;
  }
  const db = openDatabase(resolve(values.db ?? process.env.GA_DB ?? "genome-assembler.db"));
  const repo = new Repository(db);
  try {
    switch (command) {
      case "simulate": {
        if (!values.out) throw new Error("simulate needs --out <dir>");
        const out = await simulateCommand(repo, {
          outDir: values.out,
          genomeLength: num(values.length, 1_000_000, "length"),
          gc: num(values.gc, 0.5, "gc"),
          repeatFamilies: num(values.repeats, 0, "repeats"),
          repeatLength: num(values["repeat-length"], 2000, "repeat-length"),
          repeatCopies: num(values["repeat-copies"], 2, "repeat-copies"),
          repeatDivergence: num(values["repeat-divergence"], 0, "repeat-divergence"),
          coverage: num(values.coverage, 30, "coverage"),
          readLength: num(values["read-length"], 150, "read-length"),
          insertMean: num(values.insert, 400, "insert"),
          insertSd: num(values["insert-sd"], 40, "insert-sd"),
          errorRate: num(values["error-rate"], 0.005, "error-rate"),
          nRate: num(values["n-rate"], 0, "n-rate"),
          seed: num(values.seed, 42, "seed"),
          gzip: values.gzip ?? false
        });
        console.log(`run ${out.runId}: simulated ${int(out.pairs)} pairs`);
        console.log(
          table([
            ["reference", out.referencePath],
            ["reads 1", out.reads1Path],
            ["reads 2", out.reads2Path],
            ["repeats", out.repeatsPath],
            ["bases", int(out.bases)],
            ["errors", `${int(out.errorsIntroduced)} (${percent(out.observedErrorRate, 3)})`]
          ])
        );
        break;
      }
      case "stats": {
        if (positionals.length === 0) throw new Error("stats needs at least one file");
        const out = await statsCommand(repo, positionals);
        console.log(`run ${out.runId}`);
        for (const file of out.files) {
          console.log(file.path);
          console.log(
            table([
              ["reads", int(file.stats.reads)],
              ["bases", int(file.stats.bases)],
              ["length min/mean/max", `${file.stats.minLength} / ${fixed(file.stats.meanLength, 1)} / ${file.stats.maxLength}`],
              ["GC", percent(file.stats.gc)],
              ["N bases", int(file.stats.nBases)],
              ["mean quality", fixed(file.stats.meanQuality, 1)]
            ])
          );
        }
        break;
      }
      case "count": {
        if (positionals.length === 0) throw new Error("count needs at least one read file");
        const mode = values.mode ?? "bloom";
        if (mode !== "bloom" && mode !== "naive") throw new Error("--mode must be bloom or naive");
        const out = await countCommand(repo, {
          paths: positionals,
          k: num(values.k, 31, "k"),
          mode,
          falsePositiveRate: num(values.fpr, 0.01, "fpr"),
          expectedDistinct: optionalNum(values.expected, "expected"),
          exact: !(values.approx ?? false),
          minCount: optionalNum(values["min-count"], "min-count"),
          outPath: values.out,
          histogramPath: values.histogram,
          onPass: (pass) => console.error(`  ${pass}...`)
        });
        const { count, spectrum } = out;
        console.log(`run ${out.runId}: k=${count.k}, mode=${count.mode}${count.mode === "bloom" ? (count.exact ? " (exact)" : " (approx)") : ""}`);
        const rows: Array<[string, string]> = [
          ["reads", int(count.reads)],
          ["k-mers seen", int(count.totalKmers)],
          ["distinct k-mers", `${int(count.distinctKmers)}${count.distinctEstimated ? " (estimated)" : ""}`],
          ["stored in table", int(count.table.size)],
          ["peak memory (structures)", bytes(count.peakStructureBytes)]
        ];
        if (count.bloom) {
          rows.push(
            ["bloom filter", `${bytes(count.bloom.bytes)}, ${count.bloom.hashes} hashes, target FPR ${percent(count.bloom.targetFalsePositiveRate)}`],
            ["bloom FPR (estimated)", percent(count.bloom.estimatedFalsePositiveRate, 3)],
            ["false positives removed", int(count.bloom.falsePositivesRemoved)]
          );
        }
        rows.push(
          ["solid cutoff", `${spectrum.cutoff} (${spectrum.cutoffSource})`],
          ["spectrum peak", String(spectrum.peak)],
          ["k-mer coverage", fixed(spectrum.kmerCoverage)],
          ["read coverage", fixed(spectrum.readCoverage)],
          ["est. base error rate", percent(spectrum.baseErrorRate, 3)],
          ["genome size estimate", `${int(spectrum.genomeSizeEstimate)} bp`],
          ["solid / weak k-mers", `${int(spectrum.solidKmers)} / ${int(spectrum.weakKmers)}`],
          ["time", `${duration(count.elapsedMs)} over ${count.passes} pass(es)`]
        );
        if (out.outPath) rows.push(["solid k-mers written", `${int(out.solidWritten ?? 0)} -> ${out.outPath}`]);
        console.log(table(rows));
        const limit = Math.max(10, Math.ceil(spectrum.peak * 2));
        const shown = out.histogram.bins.filter((bin) => bin.multiplicity <= limit);
        if (shown.length > 0) {
          console.log("\nk-mer spectrum (multiplicity | distinct k-mers)");
          const scale = Math.max(1, ...shown.filter((bin) => bin.multiplicity >= spectrum.cutoff).map((bin) => bin.kmers));
          console.log(bars(shown.map((bin) => ({ label: String(bin.multiplicity), value: bin.kmers })), 40, scale));
        }
        break;
      }
      case "evaluate": {
        if (!values.contigs) throw new Error("evaluate needs --contigs <fasta>");
        const { runId, report } = await evaluateCommand(repo, {
          contigsPath: values.contigs,
          referencePath: values.reference,
          minContig: num(values["min-contig"], 500, "min-contig"),
          k: num(values.k, 31, "k")
        });
        console.log(`run ${runId}: contigs >= ${report.minContig} bp`);
        const rows: Array<[string, string]> = [
          ["contigs", int(report.contigs)],
          ["total length", int(report.totalLength)],
          ["largest", int(report.largest)],
          ["N50 / L50", `${int(report.n50)} / ${int(report.l50)}`],
          ["GC", percent(report.gc)]
        ];
        if (report.genomeFraction !== null) {
          rows.push(
            ["reference length", int(report.referenceLength)],
            ["NG50 / LG50", `${report.ng50 === null ? "-" : int(report.ng50)} / ${report.lg50 === null ? "-" : int(report.lg50)}`],
            ["genome fraction", percent(report.genomeFraction, 3)],
            ["duplication ratio", fixed(report.duplicationRatio, 3)],
            ["misassemblies", `${int(report.misassemblies ?? 0)} in ${int(report.misassembledContigs ?? 0)} contig(s)`],
            ["unaligned length", int(report.unalignedLength ?? 0)]
          );
        }
        console.log(table(rows));
        for (const m of report.misassemblyDetails.slice(0, 20)) console.log(`  ${m.kind} in ${m.contig} at ${int(m.contigPosition)}`);
        break;
      }
      case "runs": {
        const runs = repo.listRuns(num(values.limit, 20, "limit"));
        if (runs.length === 0) console.log("no runs yet");
        for (const run of runs) {
          const time = run.elapsed_ms === null ? "" : duration(run.elapsed_ms);
          console.log(`${String(run.id).padStart(4)}  ${run.command.padEnd(9)} ${run.status.padEnd(9)} ${run.started_at}  ${time}${run.error ? `  ${run.error}` : ""}`);
        }
        break;
      }
      case "show": {
        const id = Number(positionals[0]);
        if (!Number.isInteger(id)) throw new Error("show needs a run id");
        console.log(JSON.stringify(repo.runDetails(id), null, 2));
        break;
      }
      default:
        throw new Error(`unknown command '${command}'\n\n${USAGE}`);
    }
  } finally {
    db.close();
  }
}

main(process.argv.slice(2)).catch((error: unknown) => {
  console.error(`error: ${error instanceof Error ? error.message : String(error)}`);
  process.exitCode = 1;
});
