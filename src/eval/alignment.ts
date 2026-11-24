import { forEachCanonicalKmer, kmerShape, type KmerShape } from "../dna/kmer.js";
import type { FastxRecord } from "../io/fastx.js";
import { KmerTable } from "../kmer/table.js";

const MULTI = 0xffffffff;

export interface ReferenceIndex {
  shape: KmerShape;
  table: KmerTable;
  totalLength: number;
  offsets: Array<{ name: string; offset: number; length: number }>;
  uniqueKmers: number;
}

export interface AlignedBlock {
  contigStart: number;
  contigEnd: number;
  refStart: number;
  refEnd: number;
  reverse: boolean;
  diagonal: number;
  anchors: number;
}

export function buildReferenceIndex(reference: FastxRecord[], k: number): ReferenceIndex {
  const shape = kmerShape(k);
  const totalLength = reference.reduce((sum, record) => sum + record.seq.length, 0);
  const table = new KmerTable(totalLength);
  const offsets: ReferenceIndex["offsets"] = [];
  let offset = 0;
  for (const record of reference) {
    offsets.push({ name: record.name, offset, length: record.seq.length });
    const base = offset;
    forEachCanonicalKmer(record.seq, shape, (hi, lo, pos, forward) => {
      const value = (base + pos) * 2 + (forward ? 0 : 1);
      if (table.get(hi, lo) === -1) table.set(hi, lo, value);
      else table.set(hi, lo, MULTI);
    });
    offset += record.seq.length;
  }
  let uniqueKmers = 0;
  table.forEach((_hi, _lo, value) => {
    if (value !== MULTI) uniqueKmers++;
  });
  return { shape, table, totalLength, offsets, uniqueKmers };
}

export function alignContig(index: ReferenceIndex, contig: string, diagonalTolerance = 50): AlignedBlock[] {
  const blocks: AlignedBlock[] = [];
  const k = index.shape.k;
  let current: AlignedBlock | null = null;
  forEachCanonicalKmer(contig, index.shape, (hi, lo, pos, forward) => {
    const value = index.table.get(hi, lo);
    if (value === -1 || value === MULTI) return;
    const refPos = Math.floor(value / 2);
    const refForward = value % 2 === 0;
    const reverse = refForward !== forward;
    const diagonal = reverse ? refPos + pos : refPos - pos;
    if (current && current.reverse === reverse && Math.abs(diagonal - current.diagonal) <= diagonalTolerance) {
      current.contigEnd = pos + k;
      current.refStart = Math.min(current.refStart, refPos);
      current.refEnd = Math.max(current.refEnd, refPos + k);
      current.diagonal = diagonal;
      current.anchors++;
      return;
    }
    current = { contigStart: pos, contigEnd: pos + k, refStart: refPos, refEnd: refPos + k, reverse, diagonal, anchors: 1 };
    blocks.push(current);
  });
  return blocks;
}
