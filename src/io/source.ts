import { forEachRecord, type FastxRecord } from "./fastx.js";

export interface ReadSource {
  readonly description: string;
  forEach(visit: (record: FastxRecord) => void): Promise<number>;
}

export function fileSource(paths: string[]): ReadSource {
  if (paths.length === 0) throw new Error("at least one read file is required");
  return {
    description: paths.join(", "),
    async forEach(visit) {
      let total = 0;
      for (const path of paths) total += await forEachRecord(path, visit);
      return total;
    }
  };
}

export function memorySource(records: Array<FastxRecord | string>): ReadSource {
  const normalized = records.map((record, i) => (typeof record === "string" ? { name: `r${i}`, seq: record } : record));
  return {
    description: `${normalized.length} in-memory reads`,
    async forEach(visit) {
      for (const record of normalized) visit(record);
      return normalized.length;
    }
  };
}
