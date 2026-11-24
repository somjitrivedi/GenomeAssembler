import { writeFileSync } from "node:fs";
import { join } from "node:path";
import { gzipSync } from "node:zlib";
import { describe, expect, it } from "vitest";
import { FastxWriter, forEachRecord, readAllRecords, writeFasta } from "../src/io/fastx.js";
import { readStatsForFile } from "../src/io/readStats.js";
import { fileSource, memorySource } from "../src/io/source.js";
import { tempDir } from "./helpers.js";

describe("FASTA/FASTQ reading", () => {
  it("joins multi-line FASTA records and keeps names", async () => {
    const dir = tempDir();
    const path = join(dir, "a.fa");
    writeFileSync(path, ">chr1 first\nACGT\nAC\n\n>chr2\nGG\n");
    expect(await readAllRecords(path)).toEqual([
      { name: "chr1 first", seq: "ACGTAC" },
      { name: "chr2", seq: "GG" }
    ]);
  });

  it("parses FASTQ with CRLF line endings and '@' in quality lines", async () => {
    const dir = tempDir();
    const path = join(dir, "r.fq");
    writeFileSync(path, "@r1\r\nACGT\r\n+\r\n@@II\r\n@r2\r\nGG\r\n+r2\r\nII\r\n");
    expect(await readAllRecords(path)).toEqual([
      { name: "r1", seq: "ACGT", qual: "@@II" },
      { name: "r2", seq: "GG", qual: "II" }
    ]);
  });

  it("reads gzip input transparently", async () => {
    const dir = tempDir();
    const path = join(dir, "r.fq.gz");
    writeFileSync(path, gzipSync("@r1\nACGT\n+\nIIII\n"));
    expect(await readAllRecords(path)).toEqual([{ name: "r1", seq: "ACGT", qual: "IIII" }]);
  });

  it("reports malformed input with a line number", async () => {
    const dir = tempDir();
    const bad = join(dir, "bad.fq");
    writeFileSync(bad, "@r1\nACGT\n+\nII\n");
    await expect(readAllRecords(bad)).rejects.toThrow(/bad\.fq:4: quality length/);
    const truncated = join(dir, "trunc.fq");
    writeFileSync(truncated, "@r1\nACGT\n");
    await expect(readAllRecords(truncated)).rejects.toThrow(/truncated/);
    const other = join(dir, "x.txt");
    writeFileSync(other, "hello\n");
    await expect(readAllRecords(other)).rejects.toThrow(/not a FASTA or FASTQ/);
  });

  it("returns no records for an empty file", async () => {
    const dir = tempDir();
    const path = join(dir, "empty.fa");
    writeFileSync(path, "");
    expect(await forEachRecord(path, () => undefined)).toBe(0);
  });
});

describe("FASTA/FASTQ writing", () => {
  it("round trips plain and gzip output", async () => {
    const dir = tempDir();
    for (const name of ["out.fq", "out.fq.gz"]) {
      const path = join(dir, name);
      const writer = new FastxWriter(path);
      for (let i = 0; i < 3000; i++) await writer.fastq(`read${i}`, "ACGT".repeat(25), "I".repeat(100));
      await writer.close();
      const records = await readAllRecords(path);
      expect(records).toHaveLength(3000);
      expect(records[2999]).toEqual({ name: "read2999", seq: "ACGT".repeat(25), qual: "I".repeat(100) });
    }
  });

  it("wraps FASTA lines", async () => {
    const dir = tempDir();
    const path = join(dir, "g.fa");
    await writeFasta(path, [{ name: "g", seq: "A".repeat(170) }], 80);
    expect(await readAllRecords(path)).toEqual([{ name: "g", seq: "A".repeat(170) }]);
  });
});

describe("read sources and stats", () => {
  it("computes read statistics", async () => {
    const dir = tempDir();
    const path = join(dir, "r.fq");
    writeFileSync(path, "@a\nACGN\n+\n!!II\n@b\nGG\n+\nII\n");
    const stats = await readStatsForFile(path);
    expect(stats).toMatchObject({ reads: 2, bases: 6, minLength: 2, maxLength: 4, meanLength: 3, nBases: 1 });
    expect(stats.gc).toBeCloseTo(4 / 5);
    expect(stats.meanQuality).toBeCloseTo((0 + 0 + 40 + 40 + 40 + 40) / 6);
  });

  it("iterates several files and in-memory reads", async () => {
    const dir = tempDir();
    const a = join(dir, "a.fa");
    const b = join(dir, "b.fa");
    writeFileSync(a, ">x\nAC\n");
    writeFileSync(b, ">y\nGT\n>z\nTT\n");
    const names: string[] = [];
    expect(await fileSource([a, b]).forEach((r) => names.push(r.name))).toBe(3);
    expect(names).toEqual(["x", "y", "z"]);
    expect(() => fileSource([])).toThrow();
    const seqs: string[] = [];
    await memorySource(["AC", { name: "n", seq: "GG" }]).forEach((r) => seqs.push(`${r.name}:${r.seq}`));
    expect(seqs).toEqual(["r0:AC", "n:GG"]);
  });
});
