import { existsSync, statSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { countCommand, evaluateCommand, simulateCommand, statsCommand, type SimulateInput } from "../src/commands.js";
import { migrate, openDatabase } from "../src/db/database.js";
import { Repository } from "../src/db/repository.js";
import { readAllRecords, writeFasta } from "../src/io/fastx.js";
import { readKmerFile } from "../src/kmer/kmerFile.js";
import { tempDir } from "./helpers.js";

function simulation(outDir: string, overrides: Partial<SimulateInput> = {}): SimulateInput {
  return {
    outDir,
    genomeLength: 80000,
    gc: 0.5,
    repeatFamilies: 2,
    repeatLength: 1000,
    repeatCopies: 2,
    repeatDivergence: 0,
    coverage: 30,
    readLength: 150,
    insertMean: 400,
    insertSd: 40,
    errorRate: 0.005,
    nRate: 0,
    seed: 7,
    gzip: false,
    ...overrides
  };
}

describe("database", () => {
  it("migrates idempotently", () => {
    const db = openDatabase(":memory:");
    expect(db.pragma("user_version", { simple: true })).toBe(1);
    migrate(db);
    expect(db.pragma("user_version", { simple: true })).toBe(1);
    const tables = db.prepare("SELECT name FROM sqlite_master WHERE type = 'table' ORDER BY name").all() as Array<{ name: string }>;
    expect(tables.map((t) => t.name)).toEqual(
      expect.arrayContaining(["evaluations", "kmer_counts", "kmer_histogram", "read_sets", "runs", "simulations"])
    );
    db.close();
  });

  it("records failed runs with their error", async () => {
    const db = openDatabase(":memory:");
    const repo = new Repository(db);
    await expect(statsCommand(repo, [join(tempDir(), "missing.fq")])).rejects.toThrow();
    const [run] = repo.listRuns();
    expect(run).toMatchObject({ command: "stats", status: "failed" });
    expect(run.error).toMatch(/ENOENT/);
    db.close();
  });
});

describe("pipeline commands", () => {
  it("simulates, counts, evaluates and records everything in SQLite", async () => {
    const dir = tempDir();
    const db = openDatabase(join(dir, "ga.db"));
    const repo = new Repository(db);

    const sim = await simulateCommand(repo, simulation(join(dir, "sim"), { gzip: true }));
    expect(sim.pairs).toBe(8000);
    expect(sim.observedErrorRate).toBeCloseTo(0.005, 2);
    for (const path of [sim.referencePath, sim.reads1Path, sim.reads2Path, sim.repeatsPath]) expect(existsSync(path)).toBe(true);
    expect((await readAllRecords(sim.reads1Path)).length).toBe(8000);

    const stats = await statsCommand(repo, [sim.reads1Path, sim.reads2Path]);
    expect(stats.files.map((f) => f.stats.reads)).toEqual([8000, 8000]);
    expect(stats.files[0].stats.meanLength).toBe(150);

    const kmersPath = join(dir, "solid.kmers");
    const histogramPath = join(dir, "hist.tsv");
    const count = await countCommand(repo, {
      paths: [sim.reads1Path, sim.reads2Path],
      k: 31,
      mode: "bloom",
      falsePositiveRate: 0.01,
      exact: true,
      outPath: kmersPath,
      histogramPath
    });
    expect(count.readStats.reads).toBe(16000);
    expect(Math.abs(count.spectrum.genomeSizeEstimate - 80000) / 80000).toBeLessThan(0.05);
    expect(statSync(histogramPath).size).toBeGreaterThan(0);
    const solid = await readKmerFile(kmersPath);
    expect(solid.header.entries).toBe(count.solidWritten);
    expect(solid.header.minCount).toBe(count.spectrum.cutoff);

    const contigsPath = join(dir, "contigs.fa");
    const reference = await readAllRecords(sim.referencePath);
    await writeFasta(contigsPath, [
      { name: "a", seq: reference[0].seq.slice(0, 40000) },
      { name: "b", seq: reference[0].seq.slice(40000) }
    ]);
    const evaluation = await evaluateCommand(repo, { contigsPath, referencePath: sim.referencePath, minContig: 500, k: 31 });
    expect(evaluation.report).toMatchObject({ contigs: 2, n50: 40000, misassemblies: 0 });

    const runs = repo.listRuns();
    expect(runs.map((r) => [r.command, r.status])).toEqual([
      ["evaluate", "succeeded"],
      ["count", "succeeded"],
      ["stats", "succeeded"],
      ["simulate", "succeeded"]
    ]);
    const details = repo.runDetails(count.runId) as {
      kmerCount: { genome_size_estimate: number; mode: string; stored_kmers: number };
      histogram: Array<{ multiplicity: number; kmers: number; estimated: number }>;
      readSets: unknown[];
    };
    expect(details.kmerCount).toMatchObject({ mode: "bloom", genome_size_estimate: count.spectrum.genomeSizeEstimate });
    expect(details.readSets).toHaveLength(2);
    expect(details.histogram[0]).toMatchObject({ multiplicity: 1, estimated: 1 });
    expect(details.histogram.reduce((sum, bin) => sum + bin.kmers, 0)).toBe(count.count.distinctKmers);
    expect((repo.runDetails(sim.runId) as { simulation: { pairs: number } }).simulation.pairs).toBe(8000);
    expect(() => repo.runDetails(999)).toThrow(/not found/);
    db.close();
  });
});
