import { describe, expect, it } from "vitest";
import { reverseComplement } from "../src/dna/alphabet.js";
import { assemblyStats } from "../src/eval/assemblyStats.js";
import { evaluateAssembly } from "../src/eval/evaluate.js";
import { simulateGenome } from "../src/sim/genome.js";
import { randomDna } from "./helpers.js";

describe("assembly statistics", () => {
  it("computes N50, L50, NG50 and LG50", () => {
    const stats = assemblyStats([2, 3, 4, 5, 6, 7, 8, 9, 10], 100);
    expect(stats).toEqual({ contigs: 9, totalLength: 54, largest: 10, n50: 8, l50: 3, ng50: 3, lg50: 8 });
    expect(assemblyStats([2, 3, 4, 5, 6, 7, 8, 9, 10], 200)).toMatchObject({ ng50: null, lg50: null });
    expect(assemblyStats([80, 70, 50, 40, 30, 20], 200)).toMatchObject({ n50: 70, l50: 2, ng50: 70, lg50: 2 });
    expect(assemblyStats([])).toMatchObject({ contigs: 0, n50: 0, largest: 0 });
  });
});

describe("reference-based evaluation", () => {
  const genome = simulateGenome({ length: 60000, repeatFamilies: 2, repeatLength: 800, repeatCopies: 2, seed: 17 });
  const reference = [{ name: "ref", seq: genome.sequence }];

  it("scores the reference against itself as perfect", () => {
    const report = evaluateAssembly(reference, reference);
    expect(report).toMatchObject({ contigs: 1, n50: 60000, ng50: 60000, misassemblies: 0, unalignedLength: 0 });
    expect(report.genomeFraction).toBe(1);
    expect(report.duplicationRatio).toBe(1);
  });

  it("handles shredded and reverse-complemented pieces", () => {
    const pieces = [];
    for (let start = 0; start < 60000; start += 10000) {
      const seq = genome.sequence.slice(start, start + 10000);
      pieces.push({ name: `c${start}`, seq: start % 20000 === 0 ? reverseComplement(seq) : seq });
    }
    const report = evaluateAssembly(pieces, reference);
    expect(report).toMatchObject({ contigs: 6, n50: 10000, l50: 3, misassemblies: 0 });
    expect(report.genomeFraction).toBeGreaterThan(0.99);
  });

  it("detects relocations and inversions in chimeric contigs", () => {
    const s = genome.sequence;
    const relocation = { name: "reloc", seq: s.slice(1000, 6000) + s.slice(30000, 35000) };
    const inversion = { name: "inv", seq: s.slice(40000, 45000) + reverseComplement(s.slice(45000, 50000)) };
    const report = evaluateAssembly([relocation, inversion], reference);
    expect(report.misassemblies).toBe(2);
    expect(report.misassembledContigs).toBe(2);
    expect(report.misassemblyDetails.map((m) => [m.contig, m.kind])).toEqual([
      ["reloc", "relocation"],
      ["inv", "inversion"]
    ]);
  });

  it("tolerates small indels and counts foreign sequence as unaligned", () => {
    const s = genome.sequence;
    const withIndel = { name: "indel", seq: s.slice(10000, 15000) + s.slice(15005, 20000) };
    const foreign = { name: "foreign", seq: randomDna(3000, 99) };
    const report = evaluateAssembly([withIndel, foreign], reference);
    expect(report.misassemblies).toBe(0);
    expect(report.unalignedLength).toBeGreaterThanOrEqual(3000);
    expect(report.unalignedLength).toBeLessThan(3100);
  });

  it("filters short contigs and works without a reference", () => {
    const report = evaluateAssembly(
      [
        { name: "a", seq: "A".repeat(400) },
        { name: "b", seq: genome.sequence.slice(0, 1000) }
      ],
      null
    );
    expect(report).toMatchObject({ contigs: 1, totalLength: 1000, genomeFraction: null, misassemblies: null });
  });
});
