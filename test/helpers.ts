import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { afterEach } from "vitest";
import { reverseComplement } from "../src/dna/alphabet.js";
import { Rng } from "../src/sim/random.js";

const created: string[] = [];

export function tempDir(): string {
  const dir = mkdtempSync(join(tmpdir(), "genome-assembler-"));
  created.push(dir);
  return dir;
}

afterEach(() => {
  while (created.length > 0) rmSync(created.pop()!, { recursive: true, force: true });
});

export function randomDna(length: number, seed: number): string {
  const rng = new Rng(seed);
  let out = "";
  for (let i = 0; i < length; i++) out += "ACGT"[rng.int(4)];
  return out;
}

export function bruteCanonicalKmers(seq: string, k: number): Array<{ kmer: string; pos: number; forward: boolean }> {
  const out: Array<{ kmer: string; pos: number; forward: boolean }> = [];
  for (let i = 0; i + k <= seq.length; i++) {
    const fwd = seq.slice(i, i + k).toUpperCase();
    if (!/^[ACGT]+$/.test(fwd)) continue;
    const rev = reverseComplement(fwd);
    out.push(fwd <= rev ? { kmer: fwd, pos: i, forward: true } : { kmer: rev, pos: i, forward: false });
  }
  return out;
}
