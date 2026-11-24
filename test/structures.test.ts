import { describe, expect, it } from "vitest";
import { BloomFilter, bloomParameters } from "../src/kmer/bloom.js";
import { HyperLogLog } from "../src/kmer/hyperloglog.js";
import { KmerTable } from "../src/kmer/table.js";
import { Rng } from "../src/sim/random.js";

function randomKeys(count: number, seed: number): Array<[number, number]> {
  const rng = new Rng(seed);
  const keys = new Map<string, [number, number]>();
  while (keys.size < count) {
    const hi = rng.nextUint32() & 0x3fffffff;
    const lo = rng.nextUint32();
    keys.set(`${hi}:${lo}`, [hi, lo]);
  }
  return [...keys.values()];
}

describe("KmerTable", () => {
  it("behaves like a map through growth, updates and deletions", () => {
    const table = new KmerTable(4);
    const model = new Map<string, number>();
    const keys = randomKeys(5000, 11);
    const rng = new Rng(12);
    for (let step = 0; step < 40000; step++) {
      const [hi, lo] = keys[rng.int(keys.length)];
      const id = `${hi}:${lo}`;
      const action = rng.int(10);
      if (action < 6) {
        table.add(hi, lo, 1);
        model.set(id, (model.get(id) ?? 0) + 1);
      } else if (action < 8) {
        expect(table.delete(hi, lo)).toBe(model.delete(id));
      } else {
        expect(table.get(hi, lo)).toBe(model.get(id) ?? -1);
      }
    }
    expect(table.size).toBe(model.size);
    const dumped = new Map<string, number>();
    table.forEach((hi, lo, value) => dumped.set(`${hi}:${lo}`, value));
    expect(dumped).toEqual(model);
  });

  it("removes matching entries in place without losing the rest", () => {
    const table = new KmerTable(16);
    const keys = randomKeys(20000, 21);
    keys.forEach(([hi, lo], i) => table.set(hi, lo, i % 3));
    const capacity = table.capacity;
    const removed = table.removeWhere((_hi, _lo, value) => value === 0);
    expect(removed).toBe(Math.ceil(keys.length / 3));
    expect(table.capacity).toBe(capacity);
    keys.forEach(([hi, lo], i) => expect(table.get(hi, lo)).toBe(i % 3 === 0 ? -1 : i % 3));
  });

  it("supports addIfPresent, fillValues and filter", () => {
    const table = new KmerTable();
    table.set(1, 2, 5);
    expect(table.addIfPresent(1, 2, 3)).toBe(true);
    expect(table.addIfPresent(9, 9, 3)).toBe(false);
    expect(table.get(1, 2)).toBe(8);
    table.set(3, 4, 1);
    table.fillValues(0);
    expect(table.get(1, 2)).toBe(0);
    table.set(1, 2, 7);
    const kept = table.filter((_hi, _lo, value) => value > 0);
    expect(kept.size).toBe(1);
    expect(kept.get(1, 2)).toBe(7);
    expect(table.bytes).toBe(table.capacity * 12);
  });
});

describe("BloomFilter", () => {
  it("never reports a false negative", () => {
    const bloom = BloomFilter.forCapacity(10000, 0.01);
    const keys = randomKeys(10000, 31);
    let collisions = 0;
    for (const [hi, lo] of keys) if (bloom.testAndAdd(hi, lo)) collisions++;
    expect(collisions / keys.length).toBeLessThan(0.01);
    for (const [hi, lo] of keys) expect(bloom.has(hi, lo)).toBe(true);
    for (const [hi, lo] of keys) expect(bloom.testAndAdd(hi, lo)).toBe(true);
  });

  it("keeps the false positive rate near the target", () => {
    for (const target of [0.05, 0.01, 0.001]) {
      const bloom = BloomFilter.forCapacity(20000, target);
      const keys = randomKeys(120000, 41);
      for (const [hi, lo] of keys.slice(0, 20000)) bloom.add(hi, lo);
      let hits = 0;
      for (const [hi, lo] of keys.slice(20000)) if (bloom.has(hi, lo)) hits++;
      const observed = hits / 100000;
      expect(observed).toBeLessThan(target * 1.6);
      expect(Math.abs(bloom.estimatedFalsePositiveRate() - target) / target).toBeLessThan(0.3);
    }
  });

  it("sizes itself with the textbook formula", () => {
    const { bits, hashes } = bloomParameters(1_000_000, 0.01);
    expect(bits / 1_000_000).toBeCloseTo(9.585, 1);
    expect(hashes).toBe(7);
    expect(() => bloomParameters(10, 0)).toThrow();
    expect(() => new BloomFilter(8, 1)).toThrow();
  });
});

describe("HyperLogLog", () => {
  it("estimates cardinality within a few percent", () => {
    for (const n of [100, 5000, 200000]) {
      const hll = new HyperLogLog(14);
      for (const [hi, lo] of randomKeys(n, n)) {
        hll.add(hi, lo);
        hll.add(hi, lo);
      }
      expect(Math.abs(hll.estimate() - n) / n).toBeLessThan(0.03);
    }
  });

  it("returns zero when empty", () => {
    expect(new HyperLogLog().estimate()).toBe(0);
  });
});
