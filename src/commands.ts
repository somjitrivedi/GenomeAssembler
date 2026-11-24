import { mkdir, writeFile } from "node:fs/promises";
import { join, resolve } from "node:path";
import type { Repository } from "./db/repository.js";
import { evaluateAssembly, type EvaluationReport } from "./eval/evaluate.js";
import { FastxWriter, forEachRecord, readAllRecords } from "./io/fastx.js";
import { ReadStatsAccumulator, readStatsForFile, type ReadStats } from "./io/readStats.js";
import type { ReadSource } from "./io/source.js";
import { countKmers, type CountMode, type CountResult } from "./kmer/counter.js";
import { writeKmerFile } from "./kmer/kmerFile.js";
import { analyzeSpectrum, histogramFromTable, type Histogram, type SpectrumAnalysis } from "./kmer/spectrum.js";
import { simulateGenome } from "./sim/genome.js";
import { pairName, simulatePairs, validateReadOptions } from "./sim/reads.js";

export interface SimulateInput {
  outDir: string;
  genomeLength: number;
  gc: number;
  repeatFamilies: number;
  repeatLength: number;
  repeatCopies: number;
  repeatDivergence: number;
  coverage: number;
  readLength: number;
  insertMean: number;
  insertSd: number;
  errorRate: number;
  nRate: number;
  seed: number;
  gzip: boolean;
}

export interface SimulateOutput {
  runId: number;
  referencePath: string;
  reads1Path: string;
  reads2Path: string;
  repeatsPath: string;
  pairs: number;
  reads: number;
  bases: number;
  errorsIntroduced: number;
  observedErrorRate: number;
}

export async function simulateCommand(repo: Repository, input: SimulateInput): Promise<SimulateOutput> {
  const { runId, result } = await repo.track("simulate", input, async (id) => {
    const readOptions = {
      coverage: input.coverage,
      readLength: input.readLength,
      insertMean: input.insertMean,
      insertSd: input.insertSd,
      errorRate: input.errorRate,
      nRate: input.nRate,
      seed: input.seed + 1
    };
    validateReadOptions(input.genomeLength, readOptions);
    const outDir = resolve(input.outDir);
    await mkdir(outDir, { recursive: true });
    const genome = simulateGenome({
      length: input.genomeLength,
      gc: input.gc,
      repeatFamilies: input.repeatFamilies,
      repeatLength: input.repeatLength,
      repeatCopies: input.repeatCopies,
      repeatDivergence: input.repeatDivergence,
      seed: input.seed
    });
    const suffix = input.gzip ? ".fq.gz" : ".fq";
    const referencePath = join(outDir, "reference.fa");
    const reads1Path = join(outDir, `reads_1${suffix}`);
    const reads2Path = join(outDir, `reads_2${suffix}`);
    const repeatsPath = join(outDir, "repeats.tsv");
    const reference = new FastxWriter(referencePath);
    await reference.fasta("sim_genome", genome.sequence);
    await reference.close();
    await writeFile(
      repeatsPath,
      ["family\tstart\tend\tstrand", ...genome.repeats.map((r) => `${r.family}\t${r.start}\t${r.start + r.length}\t${r.reverse ? "-" : "+"}`)].join("\n") + "\n"
    );
    const w1 = new FastxWriter(reads1Path);
    const w2 = new FastxWriter(reads2Path);
    let errors = 0;
    let bases = 0;
    let pairs = 0;
    for (const pair of simulatePairs(genome.sequence, readOptions)) {
      await w1.fastq(pairName(pair, 1), pair.read1.seq, pair.read1.qual);
      await w2.fastq(pairName(pair, 2), pair.read2.seq, pair.read2.qual);
      errors += pair.read1.errors + pair.read2.errors;
      bases += pair.read1.seq.length + pair.read2.seq.length;
      pairs++;
    }
    await Promise.all([w1.close(), w2.close()]);
    repo.saveSimulation(id, {
      seed: input.seed,
      genomeLength: input.genomeLength,
      gc: input.gc,
      repeatFamilies: input.repeatFamilies,
      repeatLength: input.repeatLength,
      repeatCopies: input.repeatCopies,
      repeatDivergence: input.repeatDivergence,
      coverage: input.coverage,
      readLength: input.readLength,
      insertMean: input.insertMean,
      insertSd: input.insertSd,
      errorRate: input.errorRate,
      nRate: input.nRate,
      pairs,
      errorsIntroduced: errors,
      referencePath,
      reads1Path,
      reads2Path
    });
    return {
      referencePath,
      reads1Path,
      reads2Path,
      repeatsPath,
      pairs,
      reads: pairs * 2,
      bases,
      errorsIntroduced: errors,
      observedErrorRate: bases === 0 ? 0 : errors / bases
    };
  });
  return { runId, ...result };
}

export interface StatsOutput {
  runId: number;
  files: Array<{ path: string; stats: ReadStats }>;
}

export async function statsCommand(repo: Repository, paths: string[]): Promise<StatsOutput> {
  const { runId, result } = await repo.track("stats", { paths }, async (id) => {
    const files: StatsOutput["files"] = [];
    for (const path of paths) {
      const stats = await readStatsForFile(path);
      repo.saveReadSet(id, resolve(path), stats);
      files.push({ path: resolve(path), stats });
    }
    return files;
  });
  return { runId, files: result };
}

export interface CountInput {
  paths: string[];
  k: number;
  mode: CountMode;
  falsePositiveRate: number;
  expectedDistinct?: number;
  exact: boolean;
  minCount?: number;
  outPath?: string;
  histogramPath?: string;
  onPass?: (pass: string) => void;
}

export interface CountOutput {
  runId: number;
  count: CountResult;
  histogram: Histogram;
  spectrum: SpectrumAnalysis;
  readStats: ReadStats;
  outPath: string | null;
  solidWritten: number | null;
}

export async function countCommand(repo: Repository, input: CountInput): Promise<CountOutput> {
  const { onPass, ...recorded } = input;
  const { runId, result } = await repo.track("count", recorded, async (id) => {
    const paths = input.paths.map((path) => resolve(path));
    const perFile = paths.map(() => new ReadStatsAccumulator());
    const total = new ReadStatsAccumulator();
    let statsTaken = false;
    const trackedSource: ReadSource = {
      description: paths.join(", "),
      async forEach(visit) {
        let records = 0;
        for (let i = 0; i < paths.length; i++) {
          const acc = perFile[i];
          records += await forEachRecord(paths[i], (record) => {
            if (!statsTaken) {
              acc.add(record);
              total.add(record);
            }
            visit(record);
          });
        }
        statsTaken = true;
        return records;
      }
    };
    const count = await countKmers(trackedSource, {
      k: input.k,
      mode: input.mode,
      falsePositiveRate: input.falsePositiveRate,
      expectedDistinct: input.expectedDistinct,
      exact: input.exact,
      onPass
    });
    const readStats = total.result();
    const histogram = histogramFromTable(count.table, count.distinctEstimated ? count.distinctKmers : undefined);
    const spectrum = analyzeSpectrum(histogram, {
      k: input.k,
      meanReadLength: readStats.meanLength,
      minCount: input.minCount,
      totalKmers: count.totalKmers
    });
    let outPath: string | null = null;
    let solidWritten: number | null = null;
    if (input.outPath) {
      outPath = resolve(input.outPath);
      solidWritten = await writeKmerFile(outPath, count.table, input.k, spectrum.cutoff);
    }
    if (input.histogramPath) {
      const lines = ["multiplicity\tkmers", ...histogram.bins.map((bin) => `${bin.multiplicity}\t${bin.kmers}`)];
      await writeFile(resolve(input.histogramPath), lines.join("\n") + "\n");
    }
    paths.forEach((path, i) => repo.saveReadSet(id, path, perFile[i].result()));
    repo.saveKmerCount(id, paths, count, histogram, spectrum, outPath);
    return { count, histogram, spectrum, readStats, outPath, solidWritten };
  });
  return { runId, ...result };
}

export interface EvaluateInput {
  contigsPath: string;
  referencePath?: string;
  minContig: number;
  k: number;
}

export async function evaluateCommand(repo: Repository, input: EvaluateInput): Promise<{ runId: number; report: EvaluationReport }> {
  const { runId, result } = await repo.track("evaluate", input, async (id) => {
    const contigs = await readAllRecords(input.contigsPath);
    const reference = input.referencePath ? await readAllRecords(input.referencePath) : null;
    const report = evaluateAssembly(contigs, reference, { minContig: input.minContig, k: input.k });
    repo.saveEvaluation(id, resolve(input.contigsPath), input.referencePath ? resolve(input.referencePath) : null, report);
    return report;
  });
  return { runId, report: result };
}
