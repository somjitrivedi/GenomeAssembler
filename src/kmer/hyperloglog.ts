import { hashKmer } from "../dna/kmer.js";

const INDEX_SEED = 0x1b873593;
const RANK_SEED = 0xe6546b64;

export class HyperLogLog {
  readonly precision: number;
  private readonly registers: Uint8Array;

  constructor(precision = 14) {
    if (!Number.isInteger(precision) || precision < 4 || precision > 18) {
      throw new Error(`HyperLogLog precision must be between 4 and 18, got ${precision}`);
    }
    this.precision = precision;
    this.registers = new Uint8Array(1 << precision);
  }

  add(hi: number, lo: number): void {
    const index = hashKmer(hi, lo, INDEX_SEED) >>> (32 - this.precision);
    const rank = Math.clz32(hashKmer(hi, lo, RANK_SEED)) + 1;
    if (rank > this.registers[index]) this.registers[index] = rank;
  }

  estimate(): number {
    const m = this.registers.length;
    let sum = 0;
    let zeros = 0;
    for (let i = 0; i < m; i++) {
      const r = this.registers[i];
      sum += Math.pow(2, -r);
      if (r === 0) zeros++;
    }
    const alpha = 0.7213 / (1 + 1.079 / m);
    const raw = (alpha * m * m) / sum;
    if (raw <= 2.5 * m && zeros > 0) return Math.round(m * Math.log(m / zeros));
    return Math.round(raw);
  }
}
