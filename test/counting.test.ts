import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { decodeKmer, encodeKmer, forEachCanonicalKmer, kmerShape } from "../src/dna/kmer.js";
import { memorySource } from "../src/io/source.js";
import { countKmers } from "../src/kmer/counter.js";
import { readKmerFile, writeKmerFile } from "../src/kmer/kmerFile.js";
import { analyzeSpectrum, histogramFromTable, type Histogram } from "../src/kmer/spectrum.js";
import { KmerTable } from "../src/kmer/table.js";
import { simulateGenome } from "../src/sim/genome.js";
import { simulatePairs } from "../src/sim/reads.js";
import { tempDir } from "./helpers.js";

function simulatedReads(options: { length: number; coverage: number; errorRate: number; seed: number; repeats?: number }): string[] {
  const genome = simulateGenome({
    length: options.length,
    repeatFamilies: options.repeats ?? 0,
    repeatLength: 1000,
    repeatCopies: 3,
    seed: options.seed
  }).sequence;
  const reads: string[] = [];
  for (const pair of simulatePairs(genome, {
    coverage: options.coverage,
    readLength: 150,
    insertMean: 400,
    insertSd: 40,
    errorRate: options.errorRate,
    seed: options.seed + 1
  })) {
    reads.push(pair.read1.seq, pair.read2.seq);
  }
  return reads;
}

function tableToMap(table: KmerTable): Map<string, number> {
  const out = new Map<string, number>();
  table.forEach((hi, lo, value) => out.set(`${hi}:${lo}`, value));
  return out;
}

describe("k-mer counting", () => {
  const reads = simulatedReads({ length: 100000, coverage: 30, errorRate: 0.01, seed: 5 });

  it("counts exactly like a brute force map in naive mode", async () => {
    const small = reads.slice(0, 500);
    const shape = kmerShape(21);
    const model = new Map<string, number>();
    for (const read of small) {
      forEachCanonicalKmer(read, shape, (hi, lo) => {
        const key = decodeKmer(hi, lo, 21);
        model.set(key, (model.get(key) ?? 0) + 1);
      });
    }
    const result = await countKmers(memorySource(small), { k: 21, mode: "naive" });
    expect(result.distinctKmers).toBe(model.size);
    expect(result.totalKmers).toBe(small.length * 130);
    for (const [kmer, count] of model) {
      const { hi, lo } = encodeKmer(kmer);
      expect(result.table.get(hi, lo)).toBe(count);
    }
  });

  it("keeps every repeated k-mer with its exact count in bloom mode", async () => {
    const naive = await countKmers(memorySource(reads), { k: 31, mode: "naive" });
    const bloom = await countKmers(memorySource(reads), { k: 31, mode: "bloom", falsePositiveRate: 0.01 });
    const expected = new Map([...tableToMap(naive.table)].filter(([, count]) => count >= 2));
    expect(tableToMap(bloom.table)).toEqual(expected);
    expect(bloom.passes).toBe(3);
    expect(bloom.bloom!.falsePositivesRemoved).toBeGreaterThan(0);
    expect(Math.abs(bloom.distinctKmers - naive.distinctKmers) / naive.distinctKmers).toBeLessThan(0.03);
    expect(bloom.totalKmers).toBe(naive.totalKmers);
  });

  it("is nearly exact in approximate mode and skips the sketch pass when told the size", async () => {
    const naive = await countKmers(memorySource(reads), { k: 31, mode: "naive" });
    const approx = await countKmers(memorySource(reads), { k: 31, exact: false, expectedDistinct: naive.distinctKmers });
    expect(approx.passes).toBe(1);
    let wrong = 0;
    approx.table.forEach((hi, lo, value) => {
      if (naive.table.get(hi, lo) !== value) wrong++;
    });
    expect(wrong / approx.table.size).toBeLessThan(0.02);
    expect(approx.distinctEstimated).toBe(false);
  });

  it("uses at least 3x less memory than naive counting on error-prone reads", async () => {
    const naive = await countKmers(memorySource(reads), { k: 31, mode: "naive" });
    const bloom = await countKmers(memorySource(reads), { k: 31 });
    expect(naive.peakStructureBytes / bloom.peakStructureBytes).toBeGreaterThanOrEqual(3);
  });

  it("rejects even k", async () => {
    await expect(countKmers(memorySource(reads.slice(0, 2)), { k: 30 })).rejects.toThrow(/odd/);
  });
});

describe("k-mer spectrum", () => {
  it("finds the valley, peak and coverage of a synthetic spectrum", () => {
    const histogram: Histogram = {
      bins: [
        [1, 100000],
        [2, 9000],
        [3, 800],
        [4, 200],
        [5, 400],
        [6, 900],
        [7, 1200],
        [8, 1300],
        [9, 1200],
        [10, 900],
        [11, 400],
        [16, 60]
      ].map(([multiplicity, kmers]) => ({ multiplicity, kmers })),
      singletonsEstimated: false
    };
    const spectrum = analyzeSpectrum(histogram, { k: 31, meanReadLength: 150 });
    expect(spectrum).toMatchObject({ cutoff: 4, cutoffSource: "valley", peak: 8 });
    expect(spectrum.kmerCoverage).toBeCloseTo(8, 0);
    expect(spectrum.readCoverage).toBeCloseTo((spectrum.kmerCoverage * 150) / 120, 6);
    expect(spectrum.weakKmers).toBe(109800);
    expect(spectrum.baseErrorRate).toBeNull();
    const withTotal = analyzeSpectrum(histogram, { k: 31, totalKmers: spectrum.solidOccurrences * 2 });
    expect(withTotal.solidFraction).toBeCloseTo(0.5);
    expect(withTotal.baseErrorRate).toBeCloseTo(1 - Math.pow(0.5, 1 / 31), 6);
    expect(analyzeSpectrum(histogram, { k: 31, minCount: 6 })).toMatchObject({ cutoff: 6, cutoffSource: "manual" });
  });

  it("fills in the singleton bin from the distinct estimate", () => {
    const table = new KmerTable();
    table.set(1, 1, 2);
    table.set(2, 2, 3);
    expect(histogramFromTable(table, 10)).toEqual({
      bins: [
        { multiplicity: 1, kmers: 8 },
        { multiplicity: 2, kmers: 1 },
        { multiplicity: 3, kmers: 1 }
      ],
      singletonsEstimated: true
    });
  });

  for (const scenario of [
    { length: 150000, coverage: 25, errorRate: 0.005, seed: 21, repeats: 0 },
    { length: 150000, coverage: 40, errorRate: 0.01, seed: 22, repeats: 5 },
    { length: 80000, coverage: 60, errorRate: 0.002, seed: 23, repeats: 3 }
  ]) {
    it(`estimates genome size within 5% (${scenario.length} bp, ${scenario.coverage}x, ${scenario.errorRate * 100}% errors)`, async () => {
      const reads = simulatedReads(scenario);
      const count = await countKmers(memorySource(reads), { k: 31 });
      const spectrum = analyzeSpectrum(histogramFromTable(count.table, count.distinctKmers), {
        k: 31,
        meanReadLength: 150,
        totalKmers: count.totalKmers
      });
      expect(spectrum.cutoffSource).toBe("valley");
      expect(Math.abs(spectrum.genomeSizeEstimate - scenario.length) / scenario.length).toBeLessThan(0.05);
      expect(Math.abs((spectrum.readCoverage ?? 0) - scenario.coverage) / scenario.coverage).toBeLessThan(0.05);
      expect(Math.abs((spectrum.baseErrorRate ?? 0) - scenario.errorRate) / scenario.errorRate).toBeLessThan(0.2);
    });
  }
});

describe("solid k-mer file", () => {
  it("round trips entries at or above the cutoff", async () => {
    const dir = tempDir();
    const table = new KmerTable();
    for (let i = 0; i < 100000; i++) table.set(i & 0xfff, i * 7919, (i % 5) + 1);
    const path = join(dir, "solid.kmers");
    const written = await writeKmerFile(path, table, 31, 3);
    expect(written).toBe(60000);
    const { header, table: loaded } = await readKmerFile(path);
    expect(header).toEqual({ k: 31, minCount: 3, entries: 60000 });
    table.forEach((hi, lo, value) => expect(loaded.get(hi, lo)).toBe(value >= 3 ? value : -1));
  });
});
