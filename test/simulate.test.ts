import { describe, expect, it } from "vitest";
import { reverseComplement } from "../src/dna/alphabet.js";
import { simulateGenome } from "../src/sim/genome.js";
import { Rng } from "../src/sim/random.js";
import { pairCount, pairName, simulatePairs } from "../src/sim/reads.js";

describe("Rng", () => {
  it("is deterministic per seed and roughly uniform", () => {
    const a = new Rng(5);
    const b = new Rng(5);
    const c = new Rng(6);
    const xs = Array.from({ length: 10 }, () => a.nextUint32());
    expect(Array.from({ length: 10 }, () => b.nextUint32())).toEqual(xs);
    expect(Array.from({ length: 10 }, () => c.nextUint32())).not.toEqual(xs);
    const rng = new Rng(1);
    let sum = 0;
    let sq = 0;
    for (let i = 0; i < 100000; i++) {
      const v = rng.normal(10, 2);
      sum += v;
      sq += v * v;
    }
    const mean = sum / 100000;
    expect(mean).toBeCloseTo(10, 1);
    expect(Math.sqrt(sq / 100000 - mean * mean)).toBeCloseTo(2, 1);
  });
});

describe("genome simulation", () => {
  it("produces the requested length and GC content deterministically", () => {
    const g = simulateGenome({ length: 100000, gc: 0.6, seed: 3 });
    expect(g.sequence).toHaveLength(100000);
    expect(simulateGenome({ length: 100000, gc: 0.6, seed: 3 }).sequence).toBe(g.sequence);
    const gc = [...g.sequence].filter((c) => c === "G" || c === "C").length / g.sequence.length;
    expect(gc).toBeCloseTo(0.6, 1);
  });

  it("plants identical non-overlapping repeat copies on either strand", () => {
    const g = simulateGenome({ length: 200000, repeatFamilies: 4, repeatLength: 1000, repeatCopies: 3, seed: 9 });
    expect(g.repeats).toHaveLength(12);
    for (let i = 1; i < g.repeats.length; i++) expect(g.repeats[i].start).toBeGreaterThan(g.repeats[i - 1].start + 1000);
    for (let family = 0; family < 4; family++) {
      const copies = g.repeats
        .filter((r) => r.family === family)
        .map((r) => {
          const s = g.sequence.slice(r.start, r.start + r.length);
          return r.reverse ? reverseComplement(s) : s;
        });
      expect(new Set(copies).size).toBe(1);
    }
  });

  it("rejects impossible repeat layouts", () => {
    expect(() => simulateGenome({ length: 1000, repeatFamilies: 2, repeatLength: 400, repeatCopies: 2 })).toThrow();
    expect(() => simulateGenome({ length: 0 })).toThrow();
  });
});

describe("read simulation", () => {
  const genome = simulateGenome({ length: 50000, seed: 4 }).sequence;

  it("places error-free pairs exactly where their names say, facing each other", () => {
    let pairs = 0;
    for (const pair of simulatePairs(genome, { coverage: 10, readLength: 100, insertMean: 300, insertSd: 30, errorRate: 0, seed: 1 })) {
      let fragment = genome.slice(pair.start, pair.start + pair.fragmentLength);
      if (pair.reverse) fragment = reverseComplement(fragment);
      expect(pair.read1.seq).toBe(fragment.slice(0, 100));
      expect(reverseComplement(pair.read2.seq)).toBe(fragment.slice(fragment.length - 100));
      expect(pair.read1.qual).toHaveLength(100);
      expect(pairName(pair, 2)).toMatch(new RegExp(`^sim${pair.index}_${pair.start}_${pair.fragmentLength}_[fr]/2$`));
      pairs++;
    }
    expect(pairs).toBe(pairCount(50000, { coverage: 10, readLength: 100, insertMean: 0, insertSd: 0, errorRate: 0 }));
    expect(pairs).toBe(2500);
  });

  it("matches the requested error rate, insert size and N rate", () => {
    let errors = 0;
    let bases = 0;
    let ns = 0;
    let insertSum = 0;
    let pairs = 0;
    let lowQualityErrors = 0;
    for (const pair of simulatePairs(genome, { coverage: 40, readLength: 150, insertMean: 400, insertSd: 40, errorRate: 0.01, nRate: 0.001, seed: 2 })) {
      for (const read of [pair.read1, pair.read2]) {
        errors += read.errors;
        bases += read.seq.length;
        for (const c of read.seq) if (c === "N") ns++;
      }
      let fragment = genome.slice(pair.start, pair.start + pair.fragmentLength);
      if (pair.reverse) fragment = reverseComplement(fragment);
      for (let i = 0; i < 150; i++) {
        if (pair.read1.seq[i] !== fragment[i] && pair.read1.seq[i] !== "N" && pair.read1.qual.charCodeAt(i) - 33 < 15) lowQualityErrors++;
      }
      insertSum += pair.fragmentLength;
      pairs++;
    }
    expect(errors / bases).toBeGreaterThan(0.009);
    expect(errors / bases).toBeLessThan(0.011);
    expect(ns / bases).toBeCloseTo(0.001, 3);
    expect(insertSum / pairs).toBeCloseTo(400, -1);
    expect(lowQualityErrors / (errors / 2)).toBeGreaterThan(0.8);
  });

  it("validates options", () => {
    const run = (overrides: Record<string, number>) => [
      ...simulatePairs(genome, { coverage: 1, readLength: 100, insertMean: 300, insertSd: 30, errorRate: 0, ...overrides })
    ];
    expect(() => run({ readLength: 0 })).toThrow();
    expect(() => run({ insertMean: 50 })).toThrow();
    expect(() => run({ errorRate: 0.6 })).toThrow();
    expect(() => run({ coverage: 0 })).toThrow();
  });
});
