import { hashKmer } from "../dna/kmer.js";

const SEED_A = 0x9747b28c;
const SEED_B = 0x5bd1e995;

export function bloomParameters(expectedEntries: number, falsePositiveRate: number): { bits: number; hashes: number } {
  if (!(falsePositiveRate > 0 && falsePositiveRate < 1)) {
    throw new Error(`false positive rate must be between 0 and 1, got ${falsePositiveRate}`);
  }
  const n = Math.max(1, Math.ceil(expectedEntries));
  const raw = Math.ceil((-n * Math.log(falsePositiveRate)) / (Math.LN2 * Math.LN2));
  const bits = Math.max(64, Math.ceil(raw / 32) * 32);
  const hashes = Math.max(1, Math.round((bits / n) * Math.LN2));
  return { bits, hashes };
}

export class BloomFilter {
  readonly bits: number;
  readonly hashes: number;
  private readonly words: Uint32Array;

  constructor(bits: number, hashes: number) {
    if (!Number.isInteger(bits) || bits < 32) throw new Error(`bloom filter needs at least 32 bits, got ${bits}`);
    if (!Number.isInteger(hashes) || hashes < 1) throw new Error(`bloom filter needs at least one hash, got ${hashes}`);
    this.bits = Math.ceil(bits / 32) * 32;
    this.hashes = hashes;
    this.words = new Uint32Array(this.bits / 32);
  }

  static forCapacity(expectedEntries: number, falsePositiveRate: number): BloomFilter {
    const { bits, hashes } = bloomParameters(expectedEntries, falsePositiveRate);
    return new BloomFilter(bits, hashes);
  }

  get bytes(): number {
    return this.words.byteLength;
  }

  testAndAdd(hi: number, lo: number): boolean {
    const h1 = hashKmer(hi, lo, SEED_A);
    const h2 = hashKmer(hi, lo, SEED_B) | 1;
    const m = this.bits;
    const words = this.words;
    let present = true;
    for (let i = 0; i < this.hashes; i++) {
      const bit = ((h1 + Math.imul(i, h2)) >>> 0) % m;
      const word = bit >>> 5;
      const flag = 1 << (bit & 31);
      if ((words[word] & flag) === 0) {
        present = false;
        words[word] |= flag;
      }
    }
    return present;
  }

  add(hi: number, lo: number): void {
    this.testAndAdd(hi, lo);
  }

  has(hi: number, lo: number): boolean {
    const h1 = hashKmer(hi, lo, SEED_A);
    const h2 = hashKmer(hi, lo, SEED_B) | 1;
    const m = this.bits;
    for (let i = 0; i < this.hashes; i++) {
      const bit = ((h1 + Math.imul(i, h2)) >>> 0) % m;
      if ((this.words[bit >>> 5] & (1 << (bit & 31))) === 0) return false;
    }
    return true;
  }

  fillRatio(): number {
    let set = 0;
    for (let i = 0; i < this.words.length; i++) {
      let w = this.words[i];
      w = w - ((w >>> 1) & 0x55555555);
      w = (w & 0x33333333) + ((w >>> 2) & 0x33333333);
      set += (Math.imul((w + (w >>> 4)) & 0x0f0f0f0f, 0x01010101) >>> 24);
    }
    return set / this.bits;
  }

  estimatedFalsePositiveRate(): number {
    return Math.pow(this.fillRatio(), this.hashes);
  }
}
