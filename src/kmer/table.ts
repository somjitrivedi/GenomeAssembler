import { hashKmer } from "../dna/kmer.js";

const EMPTY = 0xffffffff;
const SLOT_SEED = 0x2545f491;
const MAX_LOAD = 0.7;
export const BYTES_PER_SLOT = 12;

function capacityFor(entries: number): number {
  let capacity = 16;
  while (capacity * MAX_LOAD < entries) capacity *= 2;
  return capacity;
}

export class KmerTable {
  private his: Uint32Array;
  private los: Uint32Array;
  private vals: Uint32Array;
  private mask: number;
  private count = 0;

  constructor(expectedEntries = 1 << 12) {
    const capacity = capacityFor(expectedEntries);
    this.his = new Uint32Array(capacity).fill(EMPTY);
    this.los = new Uint32Array(capacity);
    this.vals = new Uint32Array(capacity);
    this.mask = capacity - 1;
  }

  get size(): number {
    return this.count;
  }

  get capacity(): number {
    return this.mask + 1;
  }

  get bytes(): number {
    return this.capacity * BYTES_PER_SLOT;
  }

  private slot(hi: number, lo: number): number {
    const his = this.his;
    const los = this.los;
    const mask = this.mask;
    let i = hashKmer(hi, lo, SLOT_SEED) & mask;
    for (;;) {
      const h = his[i];
      if (h === EMPTY || (h === hi && los[i] === lo)) return i;
      i = (i + 1) & mask;
    }
  }

  private grow(): void {
    const oldHis = this.his;
    const oldLos = this.los;
    const oldVals = this.vals;
    const capacity = oldHis.length * 2;
    this.his = new Uint32Array(capacity).fill(EMPTY);
    this.los = new Uint32Array(capacity);
    this.vals = new Uint32Array(capacity);
    this.mask = capacity - 1;
    for (let i = 0; i < oldHis.length; i++) {
      const hi = oldHis[i];
      if (hi === EMPTY) continue;
      const j = this.slot(hi, oldLos[i]);
      this.his[j] = hi;
      this.los[j] = oldLos[i];
      this.vals[j] = oldVals[i];
    }
  }

  private insertAt(i: number, hi: number, lo: number, value: number): number {
    if (this.count + 1 > this.capacity * MAX_LOAD) {
      this.grow();
      i = this.slot(hi, lo);
    }
    this.his[i] = hi;
    this.los[i] = lo;
    this.vals[i] = value;
    this.count++;
    return value;
  }

  get(hi: number, lo: number): number {
    const i = this.slot(hi, lo);
    return this.his[i] === EMPTY ? -1 : this.vals[i];
  }

  has(hi: number, lo: number): boolean {
    return this.his[this.slot(hi, lo)] !== EMPTY;
  }

  set(hi: number, lo: number, value: number): void {
    const i = this.slot(hi, lo);
    if (this.his[i] === EMPTY) this.insertAt(i, hi, lo, value);
    else this.vals[i] = value;
  }

  add(hi: number, lo: number, delta: number): number {
    const i = this.slot(hi, lo);
    if (this.his[i] === EMPTY) return this.insertAt(i, hi, lo, delta);
    return (this.vals[i] += delta);
  }

  addIfPresent(hi: number, lo: number, delta: number): boolean {
    const i = this.slot(hi, lo);
    if (this.his[i] === EMPTY) return false;
    this.vals[i] += delta;
    return true;
  }

  fillValues(value: number): void {
    for (let i = 0; i < this.his.length; i++) if (this.his[i] !== EMPTY) this.vals[i] = value;
  }

  forEach(visit: (hi: number, lo: number, value: number) => void): void {
    const his = this.his;
    for (let i = 0; i < his.length; i++) if (his[i] !== EMPTY) visit(his[i], this.los[i], this.vals[i]);
  }

  removeWhere(remove: (hi: number, lo: number, value: number) => boolean): number {
    const doomedHi: number[] = [];
    const doomedLo: number[] = [];
    this.forEach((hi, lo, value) => {
      if (remove(hi, lo, value)) {
        doomedHi.push(hi);
        doomedLo.push(lo);
      }
    });
    for (let n = 0; n < doomedHi.length; n++) this.deleteSlot(this.slot(doomedHi[n], doomedLo[n]));
    return doomedHi.length;
  }

  delete(hi: number, lo: number): boolean {
    const i = this.slot(hi, lo);
    if (this.his[i] === EMPTY) return false;
    this.deleteSlot(i);
    return true;
  }

  private deleteSlot(hole: number): void {
    const his = this.his;
    const los = this.los;
    const vals = this.vals;
    const mask = this.mask;
    let j = hole;
    for (;;) {
      j = (j + 1) & mask;
      const hi = his[j];
      if (hi === EMPTY) break;
      const home = hashKmer(hi, los[j], SLOT_SEED) & mask;
      if (((j - home) & mask) >= ((j - hole) & mask)) {
        his[hole] = hi;
        los[hole] = los[j];
        vals[hole] = vals[j];
        hole = j;
      }
    }
    his[hole] = EMPTY;
    los[hole] = 0;
    vals[hole] = 0;
    this.count--;
  }

  filter(keep: (hi: number, lo: number, value: number) => boolean): KmerTable {
    let kept = 0;
    this.forEach((hi, lo, value) => {
      if (keep(hi, lo, value)) kept++;
    });
    const next = new KmerTable(kept);
    this.forEach((hi, lo, value) => {
      if (keep(hi, lo, value)) next.set(hi, lo, value);
    });
    return next;
  }
}
