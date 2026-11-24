import { describe, expect, it } from "vitest";
import { addComposition, emptyComposition, gcFraction, reverseComplement } from "../src/dna/alphabet.js";
import {
  assertCanonicalK,
  compareKmers,
  decodeKmer,
  encodeKmer,
  forEachCanonicalKmer,
  hashKmer,
  kmerShape
} from "../src/dna/kmer.js";
import { bruteCanonicalKmers, randomDna } from "./helpers.js";

describe("alphabet", () => {
  it("reverse complements including N and lower case", () => {
    expect(reverseComplement("AACGTN")).toBe("NACGTT");
    expect(reverseComplement("acgt")).toBe("acgt");
    expect(reverseComplement("")).toBe("");
  });

  it("computes GC over ACGT bases only", () => {
    const comp = addComposition("GGCCAATTNN", emptyComposition());
    expect(comp.n).toBe(2);
    expect(gcFraction(comp)).toBe(0.5);
  });
});

describe("k-mer encoding", () => {
  it("round trips every supported k", () => {
    for (let k = 1; k <= 31; k++) {
      const kmer = randomDna(k, k);
      const { hi, lo } = encodeKmer(kmer);
      expect(decodeKmer(hi, lo, k)).toBe(kmer);
    }
  });

  it("orders encoded k-mers lexicographically", () => {
    for (let n = 0; n < 200; n++) {
      const a = randomDna(31, n * 2 + 1);
      const b = randomDna(31, n * 2 + 2);
      const ea = encodeKmer(a);
      const eb = encodeKmer(b);
      expect(Math.sign(compareKmers(ea.hi, ea.lo, eb.hi, eb.lo))).toBe(Math.sign(a < b ? -1 : a > b ? 1 : 0));
    }
  });

  it("rejects invalid k and even k for canonical counting", () => {
    expect(() => kmerShape(0)).toThrow();
    expect(() => kmerShape(32)).toThrow();
    expect(() => assertCanonicalK(30)).toThrow(/odd/);
    expect(() => assertCanonicalK(31)).not.toThrow();
  });
});

describe("canonical k-mer iteration", () => {
  for (const k of [1, 5, 15, 16, 17, 21, 31]) {
    it(`matches a brute force scan for k=${k}`, () => {
      const seq = randomDna(400, 100 + k).replace(/(.{37})./g, "$1N") + "acgtacgtacgtacgtacgtacgtacgtacgtacgt";
      const shape = kmerShape(k);
      const seen: Array<{ kmer: string; pos: number; forward: boolean }> = [];
      const emitted = forEachCanonicalKmer(seq, shape, (hi, lo, pos, forward) => {
        seen.push({ kmer: decodeKmer(hi, lo, k), pos, forward });
      });
      const expected = bruteCanonicalKmers(seq, k);
      expect(emitted).toBe(expected.length);
      if (k % 2 === 1) expect(seen).toEqual(expected);
      else expect(seen.map((s) => s.kmer)).toEqual(expected.map((e) => e.kmer));
    });
  }

  it("gives a sequence and its reverse complement the same canonical k-mers", () => {
    const seq = randomDna(1000, 7);
    const shape = kmerShape(31);
    const collect = (s: string): string[] => {
      const out: string[] = [];
      forEachCanonicalKmer(s, shape, (hi, lo) => out.push(decodeKmer(hi, lo, 31)));
      return out.sort();
    };
    expect(collect(reverseComplement(seq))).toEqual(collect(seq));
  });

  it("emits nothing for reads shorter than k or full of N", () => {
    const shape = kmerShape(21);
    expect(forEachCanonicalKmer("ACGT", shape, () => undefined)).toBe(0);
    expect(forEachCanonicalKmer("N".repeat(100), shape, () => undefined)).toBe(0);
  });

  it("hashes with different seeds independently", () => {
    const { hi, lo } = encodeKmer(randomDna(31, 3));
    expect(hashKmer(hi, lo, 1)).not.toBe(hashKmer(hi, lo, 2));
    expect(hashKmer(hi, lo, 1)).toBe(hashKmer(hi, lo, 1));
  });
});
