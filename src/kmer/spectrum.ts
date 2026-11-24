import type { KmerTable } from "./table.js";

export interface HistogramBin {
  multiplicity: number;
  kmers: number;
}

export interface Histogram {
  bins: HistogramBin[];
  singletonsEstimated: boolean;
}

export interface SpectrumAnalysis {
  cutoff: number;
  cutoffSource: "valley" | "default" | "manual";
  peak: number;
  kmerCoverage: number;
  genomeSizeEstimate: number;
  solidKmers: number;
  solidOccurrences: number;
  weakKmers: number;
  solidFraction: number | null;
  baseErrorRate: number | null;
  readCoverage: number | null;
}

export function histogramFromTable(table: KmerTable, distinctKmers?: number): Histogram {
  const counts = new Map<number, number>();
  table.forEach((_hi, _lo, value) => counts.set(value, (counts.get(value) ?? 0) + 1));
  let singletonsEstimated = false;
  if (distinctKmers !== undefined && !counts.has(1)) {
    const missing = distinctKmers - table.size;
    if (missing > 0) {
      counts.set(1, missing);
      singletonsEstimated = true;
    }
  }
  const bins = [...counts.entries()]
    .map(([multiplicity, kmers]) => ({ multiplicity, kmers }))
    .sort((a, b) => a.multiplicity - b.multiplicity);
  return { bins, singletonsEstimated };
}

export interface AnalyzeOptions {
  k: number;
  meanReadLength?: number;
  minCount?: number;
  totalKmers?: number;
}

export function analyzeSpectrum(histogram: Histogram, options: AnalyzeOptions): SpectrumAnalysis {
  const lookup = new Map(histogram.bins.map((bin) => [bin.multiplicity, bin.kmers]));
  const at = (m: number): number => lookup.get(m) ?? 0;
  const maxMultiplicity = histogram.bins.length === 0 ? 0 : histogram.bins[histogram.bins.length - 1].multiplicity;
  let cutoff = 2;
  let cutoffSource: SpectrumAnalysis["cutoffSource"] = "default";
  if (options.minCount !== undefined) {
    cutoff = Math.max(1, Math.floor(options.minCount));
    cutoffSource = "manual";
  } else {
    for (let m = 1; m < maxMultiplicity; m++) {
      if (at(m + 1) > at(m)) {
        cutoff = Math.max(1, m);
        cutoffSource = "valley";
        break;
      }
    }
  }
  let peak = cutoff;
  let peakCount = -1;
  for (const bin of histogram.bins) {
    if (bin.multiplicity >= cutoff && bin.kmers > peakCount) {
      peak = bin.multiplicity;
      peakCount = bin.kmers;
    }
  }
  const spread = Math.max(1, peak - cutoff);
  let windowWeight = 0;
  let windowSum = 0;
  for (let m = Math.max(cutoff, peak - spread); m <= peak + spread; m++) {
    windowWeight += at(m);
    windowSum += m * at(m);
  }
  const kmerCoverage = windowWeight === 0 ? peak : windowSum / windowWeight;
  let solidKmers = 0;
  let solidOccurrences = 0;
  let weakKmers = 0;
  for (const bin of histogram.bins) {
    if (bin.multiplicity >= cutoff) {
      solidKmers += bin.kmers;
      solidOccurrences += bin.kmers * bin.multiplicity;
    } else {
      weakKmers += bin.kmers;
    }
  }
  const genomeSizeEstimate = kmerCoverage > 0 ? Math.round(solidOccurrences / kmerCoverage) : 0;
  const total = options.totalKmers;
  const solidFraction = total !== undefined && total > 0 ? Math.min(1, solidOccurrences / total) : null;
  const baseErrorRate = solidFraction !== null && solidFraction > 0 ? 1 - Math.pow(solidFraction, 1 / options.k) : null;
  const L = options.meanReadLength;
  const readCoverage =
    L !== undefined && L > options.k ? (kmerCoverage * L) / (L - options.k + 1) / (solidFraction !== null && solidFraction > 0 ? solidFraction : 1) : null;
  return {
    cutoff,
    cutoffSource,
    peak,
    kmerCoverage,
    genomeSizeEstimate,
    solidKmers,
    solidOccurrences,
    weakKmers,
    solidFraction,
    baseErrorRate,
    readCoverage
  };
}
