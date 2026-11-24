export interface AssemblyStats {
  contigs: number;
  totalLength: number;
  largest: number;
  n50: number;
  l50: number;
  ng50: number | null;
  lg50: number | null;
}

function nxx(sorted: number[], target: number): { length: number; count: number } | null {
  let running = 0;
  for (let i = 0; i < sorted.length; i++) {
    running += sorted[i];
    if (running >= target) return { length: sorted[i], count: i + 1 };
  }
  return null;
}

export function assemblyStats(lengths: number[], referenceLength?: number): AssemblyStats {
  const sorted = [...lengths].sort((a, b) => b - a);
  const total = sorted.reduce((sum, value) => sum + value, 0);
  const half = nxx(sorted, total / 2);
  const ngHalf = referenceLength !== undefined && referenceLength > 0 ? nxx(sorted, referenceLength / 2) : null;
  return {
    contigs: sorted.length,
    totalLength: total,
    largest: sorted[0] ?? 0,
    n50: half?.length ?? 0,
    l50: half?.count ?? 0,
    ng50: ngHalf?.length ?? null,
    lg50: ngHalf?.count ?? null
  };
}
