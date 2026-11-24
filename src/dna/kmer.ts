import { BASE_CODE, BASES } from "./alphabet.js";

export const MAX_K = 31;

export interface KmerShape {
  readonly k: number;
  readonly hiMask: number;
  readonly loMask: number;
  readonly topShift: number;
  readonly topInHi: boolean;
}

export function kmerShape(k: number): KmerShape {
  if (!Number.isInteger(k) || k < 1 || k > MAX_K) {
    throw new Error(`k must be an integer between 1 and ${MAX_K}, got ${k}`);
  }
  const bits = 2 * k;
  const hiBits = Math.max(0, bits - 32);
  const loBits = Math.min(32, bits);
  const top = 2 * (k - 1);
  return {
    k,
    hiMask: hiBits === 0 ? 0 : (1 << hiBits) - 1,
    loMask: loBits === 32 ? 0xffffffff : (1 << loBits) - 1,
    topShift: top >= 32 ? top - 32 : top,
    topInHi: top >= 32
  };
}

export function assertCanonicalK(k: number): void {
  kmerShape(k);
  if (k % 2 === 0) throw new Error(`k must be odd so no k-mer equals its own reverse complement, got ${k}`);
}

export type KmerVisitor = (hi: number, lo: number, pos: number, forward: boolean) => void;

export function compareKmers(aHi: number, aLo: number, bHi: number, bLo: number): number {
  if (aHi !== bHi) return aHi < bHi ? -1 : 1;
  if (aLo !== bLo) return aLo < bLo ? -1 : 1;
  return 0;
}

export function forEachCanonicalKmer(seq: string, shape: KmerShape, visit: KmerVisitor): number {
  const { k, hiMask, loMask, topShift, topInHi } = shape;
  let hi = 0;
  let lo = 0;
  let rHi = 0;
  let rLo = 0;
  let run = 0;
  let emitted = 0;
  const n = seq.length;
  for (let i = 0; i < n; i++) {
    const b = BASE_CODE[seq.charCodeAt(i) & 0xff];
    if (b < 0) {
      run = 0;
      hi = lo = rHi = rLo = 0;
      continue;
    }
    hi = ((hi << 2) | (lo >>> 30)) & hiMask;
    lo = (((lo << 2) | b) & loMask) >>> 0;
    rLo = ((rLo >>> 2) | ((rHi & 3) << 30)) >>> 0;
    rHi = rHi >>> 2;
    const c = 3 - b;
    if (topInHi) rHi = (rHi | (c << topShift)) >>> 0;
    else rLo = (rLo | (c << topShift)) >>> 0;
    if (++run >= k) {
      const pos = i - k + 1;
      if (hi < rHi || (hi === rHi && lo <= rLo)) visit(hi, lo, pos, true);
      else visit(rHi, rLo, pos, false);
      emitted++;
    }
  }
  return emitted;
}

export function encodeKmer(kmer: string): { hi: number; lo: number } {
  const shape = kmerShape(kmer.length);
  let hi = 0;
  let lo = 0;
  for (let i = 0; i < kmer.length; i++) {
    const b = BASE_CODE[kmer.charCodeAt(i) & 0xff];
    if (b < 0) throw new Error(`invalid base '${kmer[i]}' in k-mer ${kmer}`);
    hi = ((hi << 2) | (lo >>> 30)) & shape.hiMask;
    lo = (((lo << 2) | b) & shape.loMask) >>> 0;
  }
  return { hi: hi >>> 0, lo };
}

export function decodeKmer(hi: number, lo: number, k: number): string {
  const chars: string[] = new Array(k);
  for (let i = k - 1; i >= 0; i--) {
    chars[i] = BASES[lo & 3];
    lo = ((lo >>> 2) | ((hi & 3) << 30)) >>> 0;
    hi = hi >>> 2;
  }
  return chars.join("");
}

function fmix32(h: number): number {
  h ^= h >>> 16;
  h = Math.imul(h, 0x85ebca6b);
  h ^= h >>> 13;
  h = Math.imul(h, 0xc2b2ae35);
  h ^= h >>> 16;
  return h >>> 0;
}

export function hashKmer(hi: number, lo: number, seed: number): number {
  return fmix32(lo ^ fmix32(hi ^ seed) ^ Math.imul(seed, 0x9e3779b1));
}
