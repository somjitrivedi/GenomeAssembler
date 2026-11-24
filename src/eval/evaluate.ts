import { addComposition, emptyComposition, gcFraction } from "../dna/alphabet.js";
import type { FastxRecord } from "../io/fastx.js";
import { alignContig, buildReferenceIndex, type AlignedBlock } from "./alignment.js";
import { assemblyStats, type AssemblyStats } from "./assemblyStats.js";

export interface EvaluateOptions {
  k?: number;
  minContig?: number;
  minBlock?: number;
  misassemblyBlock?: number;
  relocationDistance?: number;
}

export interface Misassembly {
  contig: string;
  kind: "relocation" | "inversion";
  contigPosition: number;
}

export interface EvaluationReport extends AssemblyStats {
  referenceLength: number;
  minContig: number;
  gc: number;
  genomeFraction: number | null;
  duplicationRatio: number | null;
  misassemblies: number | null;
  misassembledContigs: number | null;
  unalignedLength: number | null;
  misassemblyDetails: Misassembly[];
}

function blockLength(block: AlignedBlock): number {
  return block.contigEnd - block.contigStart;
}

export function evaluateAssembly(
  contigs: FastxRecord[],
  reference: FastxRecord[] | null,
  options: EvaluateOptions = {}
): EvaluationReport {
  const k = options.k ?? 31;
  const minContig = options.minContig ?? 500;
  const minBlock = options.minBlock ?? 100;
  const misassemblyBlock = options.misassemblyBlock ?? 500;
  const relocationDistance = options.relocationDistance ?? 1000;
  const kept = contigs.filter((contig) => contig.seq.length >= minContig);
  const composition = emptyComposition();
  for (const contig of kept) addComposition(contig.seq, composition);
  const referenceLength = reference ? reference.reduce((sum, record) => sum + record.seq.length, 0) : 0;
  const stats = assemblyStats(
    kept.map((contig) => contig.seq.length),
    reference ? referenceLength : undefined
  );
  const base: EvaluationReport = {
    ...stats,
    referenceLength,
    minContig,
    gc: gcFraction(composition),
    genomeFraction: null,
    duplicationRatio: null,
    misassemblies: null,
    misassembledContigs: null,
    unalignedLength: null,
    misassemblyDetails: []
  };
  if (!reference || referenceLength === 0) return base;
  const index = buildReferenceIndex(reference, k);
  const covered = new Uint8Array(index.totalLength);
  let alignedReferenceBases = 0;
  let unaligned = 0;
  const details: Misassembly[] = [];
  const misassembled = new Set<string>();
  for (const contig of kept) {
    const blocks = alignContig(index, contig.seq).filter((block) => blockLength(block) >= minBlock);
    let alignedContigBases = 0;
    for (const block of blocks) {
      alignedContigBases += blockLength(block);
      alignedReferenceBases += block.refEnd - block.refStart;
      covered.fill(1, block.refStart, block.refEnd);
    }
    unaligned += Math.max(0, contig.seq.length - alignedContigBases);
    const major = blocks.filter((block) => blockLength(block) >= misassemblyBlock);
    for (let i = 1; i < major.length; i++) {
      const previous = major[i - 1];
      const next = major[i];
      let kind: Misassembly["kind"] | null = null;
      if (previous.reverse !== next.reverse) kind = "inversion";
      else if (Math.abs(next.diagonal - previous.diagonal) > relocationDistance) kind = "relocation";
      if (kind) {
        details.push({ contig: contig.name, kind, contigPosition: next.contigStart });
        misassembled.add(contig.name);
      }
    }
  }
  let coveredBases = 0;
  for (let i = 0; i < covered.length; i++) coveredBases += covered[i];
  return {
    ...base,
    genomeFraction: coveredBases / index.totalLength,
    duplicationRatio: coveredBases === 0 ? null : alignedReferenceBases / coveredBases,
    misassemblies: details.length,
    misassembledContigs: misassembled.size,
    unalignedLength: unaligned,
    misassemblyDetails: details
  };
}
