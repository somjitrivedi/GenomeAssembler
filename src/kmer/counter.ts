import { assertCanonicalK, forEachCanonicalKmer, kmerShape } from "../dna/kmer.js";
import type { ReadSource } from "../io/source.js";
import { BloomFilter } from "./bloom.js";
import { HyperLogLog } from "./hyperloglog.js";
import { KmerTable } from "./table.js";

export type CountMode = "bloom" | "naive";

export interface CountOptions {
  k: number;
  mode?: CountMode;
  falsePositiveRate?: number;
  expectedDistinct?: number;
  exact?: boolean;
  onPass?: (pass: string) => void;
}

export interface BloomSummary {
  bits: number;
  hashes: number;
  bytes: number;
  targetFalsePositiveRate: number;
  estimatedFalsePositiveRate: number;
  falsePositivesRemoved: number;
}

export interface CountResult {
  k: number;
  mode: CountMode;
  exact: boolean;
  table: KmerTable;
  reads: number;
  totalKmers: number;
  distinctKmers: number;
  distinctEstimated: boolean;
  bloom: BloomSummary | null;
  passes: number;
  peakStructureBytes: number;
  elapsedMs: number;
}

export async function countKmers(source: ReadSource, options: CountOptions): Promise<CountResult> {
  const started = performance.now();
  assertCanonicalK(options.k);
  const shape = kmerShape(options.k);
  const mode = options.mode ?? "bloom";
  const report = options.onPass ?? (() => undefined);
  if (mode === "naive") {
    const table = new KmerTable(1 << 16);
    let totalKmers = 0;
    report("counting all k-mers");
    const reads = await source.forEach((record) => {
      totalKmers += forEachCanonicalKmer(record.seq, shape, (hi, lo) => {
        table.add(hi, lo, 1);
      });
    });
    return {
      k: options.k,
      mode,
      exact: true,
      table,
      reads,
      totalKmers,
      distinctKmers: table.size,
      distinctEstimated: false,
      bloom: null,
      passes: 1,
      peakStructureBytes: table.bytes,
      elapsedMs: performance.now() - started
    };
  }
  const fpr = options.falsePositiveRate ?? 0.01;
  const exact = options.exact ?? true;
  let passes = 0;
  let distinct = options.expectedDistinct;
  let estimatedFromSketch = false;
  if (distinct === undefined) {
    report("estimating distinct k-mers (HyperLogLog)");
    const sketch = new HyperLogLog(14);
    await source.forEach((record) => {
      forEachCanonicalKmer(record.seq, shape, (hi, lo) => sketch.add(hi, lo));
    });
    passes++;
    distinct = sketch.estimate();
    estimatedFromSketch = true;
  }
  const filter = BloomFilter.forCapacity(Math.max(1, distinct), fpr);
  const table = new KmerTable(Math.max(1 << 12, Math.ceil(distinct / 8)));
  let totalKmers = 0;
  report("filtering singletons through the Bloom filter");
  const reads = await source.forEach((record) => {
    totalKmers += forEachCanonicalKmer(record.seq, shape, (hi, lo) => {
      if (table.addIfPresent(hi, lo, 1)) return;
      if (filter.testAndAdd(hi, lo)) table.set(hi, lo, 2);
    });
  });
  passes++;
  const peakStructureBytes = filter.bytes + table.bytes;
  const bloomBits = filter.bits;
  const bloomHashes = filter.hashes;
  const bloomBytes = filter.bytes;
  const estimatedFalsePositiveRate = filter.estimatedFalsePositiveRate();
  let falsePositivesRemoved = 0;
  if (exact) {
    report("recounting candidate k-mers exactly");
    table.fillValues(0);
    await source.forEach((record) => {
      forEachCanonicalKmer(record.seq, shape, (hi, lo) => {
        table.addIfPresent(hi, lo, 1);
      });
    });
    passes++;
    falsePositivesRemoved = table.removeWhere((_hi, _lo, count) => count < 2);
  }
  return {
    k: options.k,
    mode,
    exact,
    table,
    reads,
    totalKmers,
    distinctKmers: estimatedFromSketch ? Math.max(distinct, table.size) : table.size,
    distinctEstimated: estimatedFromSketch,
    bloom: {
      bits: bloomBits,
      hashes: bloomHashes,
      bytes: bloomBytes,
      targetFalsePositiveRate: fpr,
      estimatedFalsePositiveRate,
      falsePositivesRemoved
    },
    passes,
    peakStructureBytes,
    elapsedMs: performance.now() - started
  };
}
