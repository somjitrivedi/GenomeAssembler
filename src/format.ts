export function int(value: number): string {
  return Math.round(value).toLocaleString("en-US");
}

export function fixed(value: number | null, digits = 2): string {
  return value === null || !Number.isFinite(value) ? "-" : value.toFixed(digits);
}

export function percent(value: number | null, digits = 2): string {
  return value === null ? "-" : `${(value * 100).toFixed(digits)}%`;
}

export function bytes(value: number): string {
  const units = ["B", "KiB", "MiB", "GiB"];
  let v = value;
  let unit = 0;
  while (v >= 1024 && unit < units.length - 1) {
    v /= 1024;
    unit++;
  }
  return `${v.toFixed(unit === 0 ? 0 : 1)} ${units[unit]}`;
}

export function duration(ms: number): string {
  return ms < 1000 ? `${Math.round(ms)} ms` : `${(ms / 1000).toFixed(1)} s`;
}

export function table(rows: Array<[string, string]>): string {
  const width = Math.max(...rows.map(([label]) => label.length));
  return rows.map(([label, value]) => `  ${label.padEnd(width)}  ${value}`).join("\n");
}

export function bars(points: Array<{ label: string; value: number }>, width = 50, scaleMax?: number): string {
  const max = Math.max(1, scaleMax ?? Math.max(...points.map((p) => p.value)));
  const labelWidth = Math.max(...points.map((p) => p.label.length));
  return points
    .map((p) => `  ${p.label.padStart(labelWidth)} | ${"#".repeat(Math.min(width, Math.round((p.value / max) * width))).padEnd(width)} ${int(p.value)}`)
    .join("\n");
}
