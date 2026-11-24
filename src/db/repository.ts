import type { EvaluationReport } from "../eval/evaluate.js";
import type { ReadStats } from "../io/readStats.js";
import type { CountResult } from "../kmer/counter.js";
import type { Histogram, SpectrumAnalysis } from "../kmer/spectrum.js";
import type { Db } from "./database.js";

export interface RunRow {
  id: number;
  command: string;
  args: string;
  status: "running" | "succeeded" | "failed";
  started_at: string;
  finished_at: string | null;
  elapsed_ms: number | null;
  error: string | null;
}

export interface SimulationRecord {
  seed: number;
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
  pairs: number;
  errorsIntroduced: number;
  referencePath: string;
  reads1Path: string;
  reads2Path: string;
}

export class Repository {
  constructor(private readonly db: Db) {}

  startRun(command: string, args: unknown): number {
    const result = this.db
      .prepare("INSERT INTO runs (command, args, status, started_at) VALUES (?, ?, 'running', ?)")
      .run(command, JSON.stringify(args), new Date().toISOString());
    return Number(result.lastInsertRowid);
  }

  finishRun(runId: number, elapsedMs: number, error?: unknown): void {
    const message = error === undefined ? null : error instanceof Error ? error.message : String(error);
    this.db
      .prepare("UPDATE runs SET status = ?, finished_at = ?, elapsed_ms = ?, error = ? WHERE id = ?")
      .run(message === null ? "succeeded" : "failed", new Date().toISOString(), elapsedMs, message, runId);
  }

  async track<T>(command: string, args: unknown, work: (runId: number) => Promise<T>): Promise<{ runId: number; result: T }> {
    const runId = this.startRun(command, args);
    const started = performance.now();
    try {
      const result = await work(runId);
      this.finishRun(runId, performance.now() - started);
      return { runId, result };
    } catch (error) {
      this.finishRun(runId, performance.now() - started, error);
      throw error;
    }
  }

  saveSimulation(runId: number, s: SimulationRecord): void {
    this.db
      .prepare(
        `INSERT INTO simulations (run_id, seed, genome_length, gc, repeat_families, repeat_length, repeat_copies,
          repeat_divergence, coverage, read_length, insert_mean, insert_sd, error_rate, n_rate, pairs,
          errors_introduced, reference_path, reads1_path, reads2_path)
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)`
      )
      .run(
        runId,
        s.seed,
        s.genomeLength,
        s.gc,
        s.repeatFamilies,
        s.repeatLength,
        s.repeatCopies,
        s.repeatDivergence,
        s.coverage,
        s.readLength,
        s.insertMean,
        s.insertSd,
        s.errorRate,
        s.nRate,
        s.pairs,
        s.errorsIntroduced,
        s.referencePath,
        s.reads1Path,
        s.reads2Path
      );
  }

  saveReadSet(runId: number, path: string, stats: ReadStats): void {
    this.db
      .prepare(
        `INSERT INTO read_sets (run_id, path, reads, bases, min_length, max_length, mean_length, gc, n_bases, mean_quality)
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)`
      )
      .run(
        runId,
        path,
        stats.reads,
        stats.bases,
        stats.minLength,
        stats.maxLength,
        stats.meanLength,
        stats.gc,
        stats.nBases,
        stats.meanQuality
      );
  }

  saveKmerCount(
    runId: number,
    inputs: string[],
    count: CountResult,
    histogram: Histogram,
    spectrum: SpectrumAnalysis,
    outputPath: string | null
  ): void {
    const insertHistogram = this.db.prepare(
      "INSERT INTO kmer_histogram (run_id, multiplicity, kmers, estimated) VALUES (?, ?, ?, ?)"
    );
    this.db.transaction(() => {
      this.db
        .prepare(
          `INSERT INTO kmer_counts (run_id, k, mode, exact, inputs, reads, total_kmers, distinct_kmers, distinct_estimated,
            stored_kmers, bloom_bits, bloom_hashes, bloom_target_fpr, bloom_estimated_fpr, bloom_false_positives_removed,
            peak_structure_bytes, passes, cutoff, cutoff_source, peak_multiplicity, kmer_coverage, read_coverage,
            solid_fraction, base_error_rate, genome_size_estimate, solid_kmers, weak_kmers, output_path, elapsed_ms)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)`
        )
        .run(
          runId,
          count.k,
          count.mode,
          count.exact ? 1 : 0,
          JSON.stringify(inputs),
          count.reads,
          count.totalKmers,
          count.distinctKmers,
          count.distinctEstimated ? 1 : 0,
          count.table.size,
          count.bloom?.bits ?? null,
          count.bloom?.hashes ?? null,
          count.bloom?.targetFalsePositiveRate ?? null,
          count.bloom?.estimatedFalsePositiveRate ?? null,
          count.bloom?.falsePositivesRemoved ?? null,
          count.peakStructureBytes,
          count.passes,
          spectrum.cutoff,
          spectrum.cutoffSource,
          spectrum.peak,
          spectrum.kmerCoverage,
          spectrum.readCoverage,
          spectrum.solidFraction,
          spectrum.baseErrorRate,
          spectrum.genomeSizeEstimate,
          spectrum.solidKmers,
          spectrum.weakKmers,
          outputPath,
          count.elapsedMs
        );
      for (const bin of histogram.bins) {
        const estimated = histogram.singletonsEstimated && bin.multiplicity === 1 ? 1 : 0;
        insertHistogram.run(runId, bin.multiplicity, bin.kmers, estimated);
      }
    })();
  }

  saveEvaluation(runId: number, contigsPath: string, referencePath: string | null, report: EvaluationReport): void {
    this.db
      .prepare(
        `INSERT INTO evaluations (run_id, contigs_path, reference_path, min_contig, contigs, total_length, largest, n50, l50,
          ng50, lg50, gc, reference_length, genome_fraction, duplication_ratio, misassemblies, misassembled_contigs,
          unaligned_length)
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)`
      )
      .run(
        runId,
        contigsPath,
        referencePath,
        report.minContig,
        report.contigs,
        report.totalLength,
        report.largest,
        report.n50,
        report.l50,
        report.ng50,
        report.lg50,
        report.gc,
        report.referenceLength,
        report.genomeFraction,
        report.duplicationRatio,
        report.misassemblies,
        report.misassembledContigs,
        report.unalignedLength
      );
  }

  listRuns(limit = 20): RunRow[] {
    return this.db.prepare("SELECT * FROM runs ORDER BY id DESC LIMIT ?").all(limit) as RunRow[];
  }

  getRun(runId: number): RunRow | undefined {
    return this.db.prepare("SELECT * FROM runs WHERE id = ?").get(runId) as RunRow | undefined;
  }

  runDetails(runId: number): Record<string, unknown> {
    const run = this.getRun(runId);
    if (!run) throw new Error(`run ${runId} not found`);
    return {
      run,
      simulation: this.db.prepare("SELECT * FROM simulations WHERE run_id = ?").get(runId) ?? null,
      readSets: this.db.prepare("SELECT * FROM read_sets WHERE run_id = ? ORDER BY id").all(runId),
      kmerCount: this.db.prepare("SELECT * FROM kmer_counts WHERE run_id = ?").get(runId) ?? null,
      histogram: this.db
        .prepare("SELECT multiplicity, kmers, estimated FROM kmer_histogram WHERE run_id = ? ORDER BY multiplicity")
        .all(runId),
      evaluation: this.db.prepare("SELECT * FROM evaluations WHERE run_id = ?").get(runId) ?? null
    };
  }
}
