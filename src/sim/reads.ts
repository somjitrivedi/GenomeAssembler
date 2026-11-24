import { reverseComplement } from "../dna/alphabet.js";
import { mutateBase } from "./genome.js";
import { Rng } from "./random.js";

export interface ReadOptions {
  coverage: number;
  readLength: number;
  insertMean: number;
  insertSd: number;
  errorRate: number;
  nRate?: number;
  seed?: number;
}

export interface SimulatedRead {
  seq: string;
  qual: string;
  errors: number;
}

export interface SimulatedPair {
  index: number;
  start: number;
  fragmentLength: number;
  reverse: boolean;
  read1: SimulatedRead;
  read2: SimulatedRead;
}

export function pairCount(genomeLength: number, options: ReadOptions): number {
  return Math.max(1, Math.round((options.coverage * genomeLength) / (2 * options.readLength)));
}

export function validateReadOptions(genomeLength: number, options: ReadOptions): void {
  if (!(options.coverage > 0)) throw new Error("coverage must be positive");
  if (!Number.isInteger(options.readLength) || options.readLength < 1) throw new Error("read length must be a positive integer");
  if (options.readLength > genomeLength) throw new Error("read length cannot exceed genome length");
  if (!(options.insertMean >= options.readLength)) throw new Error("insert size mean must be at least the read length");
  if (!(options.insertSd >= 0)) throw new Error("insert size sd cannot be negative");
  if (!(options.errorRate >= 0 && options.errorRate < 0.5)) throw new Error("error rate must be in [0, 0.5)");
  if (!((options.nRate ?? 0) >= 0 && (options.nRate ?? 0) < 0.5)) throw new Error("N rate must be in [0, 0.5)");
}

function phred(probability: number): number {
  if (probability <= 0) return 41;
  return Math.max(2, Math.min(41, Math.round(-10 * Math.log10(probability))));
}

function sequenceRead(rng: Rng, template: string, options: ReadOptions): SimulatedRead {
  const length = template.length;
  const seq = Buffer.allocUnsafe(length);
  const qual = Buffer.allocUnsafe(length);
  const nRate = options.nRate ?? 0;
  let errors = 0;
  for (let i = 0; i < length; i++) {
    const position = length === 1 ? 0.5 : i / (length - 1);
    const p = options.errorRate * (0.5 + position);
    let base = template.charCodeAt(i);
    let q = phred(p) + Math.round(rng.normal(0, 1.5));
    if (nRate > 0 && rng.next() < nRate) {
      base = 78;
      q = 2;
    } else if (rng.next() < p) {
      base = mutateBase(rng, base);
      q -= 12 + rng.int(10);
      errors++;
    }
    seq[i] = base;
    qual[i] = 33 + Math.max(2, Math.min(41, q));
  }
  return { seq: seq.toString("latin1"), qual: qual.toString("latin1"), errors };
}

export function* simulatePairs(genome: string, options: ReadOptions): Generator<SimulatedPair> {
  validateReadOptions(genome.length, options);
  const rng = new Rng(options.seed ?? 2);
  const pairs = pairCount(genome.length, options);
  const L = options.readLength;
  for (let index = 0; index < pairs; index++) {
    const drawn = Math.round(rng.normal(options.insertMean, options.insertSd));
    const fragmentLength = Math.min(genome.length, Math.max(L, drawn));
    const start = rng.int(genome.length - fragmentLength + 1);
    const reverse = rng.next() < 0.5;
    let fragment = genome.slice(start, start + fragmentLength);
    if (reverse) fragment = reverseComplement(fragment);
    const read1 = sequenceRead(rng, fragment.slice(0, L), options);
    const read2 = sequenceRead(rng, reverseComplement(fragment.slice(fragmentLength - L)), options);
    yield { index, start, fragmentLength, reverse, read1, read2 };
  }
}

export function pairName(pair: SimulatedPair, mate: 1 | 2): string {
  return `sim${pair.index}_${pair.start}_${pair.fragmentLength}_${pair.reverse ? "r" : "f"}/${mate}`;
}
