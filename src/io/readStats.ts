import { addComposition, emptyComposition, gcFraction, type BaseComposition } from "../dna/alphabet.js";
import { forEachRecord, type FastxRecord } from "./fastx.js";

export interface ReadStats {
  reads: number;
  bases: number;
  minLength: number;
  maxLength: number;
  meanLength: number;
  gc: number;
  nBases: number;
  meanQuality: number | null;
}

export class ReadStatsAccumulator {
  private reads = 0;
  private bases = 0;
  private minLength = Number.POSITIVE_INFINITY;
  private maxLength = 0;
  private qualitySum = 0;
  private qualityBases = 0;
  private readonly composition: BaseComposition = emptyComposition();

  add(record: FastxRecord): void {
    const length = record.seq.length;
    this.reads++;
    this.bases += length;
    if (length < this.minLength) this.minLength = length;
    if (length > this.maxLength) this.maxLength = length;
    addComposition(record.seq, this.composition);
    if (record.qual !== undefined) {
      const qual = record.qual;
      for (let i = 0; i < qual.length; i++) this.qualitySum += qual.charCodeAt(i) - 33;
      this.qualityBases += qual.length;
    }
  }

  result(): ReadStats {
    return {
      reads: this.reads,
      bases: this.bases,
      minLength: this.reads === 0 ? 0 : this.minLength,
      maxLength: this.maxLength,
      meanLength: this.reads === 0 ? 0 : this.bases / this.reads,
      gc: gcFraction(this.composition),
      nBases: this.composition.n,
      meanQuality: this.qualityBases === 0 ? null : this.qualitySum / this.qualityBases
    };
  }
}

export async function readStatsForFile(path: string): Promise<ReadStats> {
  const acc = new ReadStatsAccumulator();
  await forEachRecord(path, (record) => acc.add(record));
  return acc.result();
}
