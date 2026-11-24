import { open, readFile } from "node:fs/promises";
import { KmerTable } from "./table.js";

const MAGIC = 0x4d4b4147;
const VERSION = 1;
const HEADER_BYTES = 24;
const RECORD_BYTES = 12;
const CHUNK_RECORDS = 1 << 16;

export interface KmerFileHeader {
  k: number;
  minCount: number;
  entries: number;
}

export async function writeKmerFile(path: string, table: KmerTable, k: number, minCount: number): Promise<number> {
  const handle = await open(path, "w");
  try {
    let entries = 0;
    table.forEach((_hi, _lo, value) => {
      if (value >= minCount) entries++;
    });
    const header = Buffer.alloc(HEADER_BYTES);
    header.writeUInt32LE(MAGIC, 0);
    header.writeUInt32LE(VERSION, 4);
    header.writeUInt32LE(k, 8);
    header.writeUInt32LE(minCount, 12);
    header.writeUInt32LE(entries % 0x100000000, 16);
    header.writeUInt32LE(Math.floor(entries / 0x100000000), 20);
    await handle.write(header);
    const chunk = Buffer.alloc(CHUNK_RECORDS * RECORD_BYTES);
    let filled = 0;
    const pending: Buffer[] = [];
    table.forEach((hi, lo, value) => {
      if (value < minCount) return;
      const offset = filled * RECORD_BYTES;
      chunk.writeUInt32LE(hi, offset);
      chunk.writeUInt32LE(lo, offset + 4);
      chunk.writeUInt32LE(value, offset + 8);
      filled++;
      if (filled === CHUNK_RECORDS) {
        pending.push(Buffer.from(chunk));
        filled = 0;
      }
    });
    for (const part of pending) await handle.write(part);
    if (filled > 0) await handle.write(chunk.subarray(0, filled * RECORD_BYTES));
    return entries;
  } finally {
    await handle.close();
  }
}

export async function readKmerFile(path: string): Promise<{ header: KmerFileHeader; table: KmerTable }> {
  const data = await readFile(path);
  if (data.length < HEADER_BYTES || data.readUInt32LE(0) !== MAGIC) throw new Error(`${path} is not a k-mer file`);
  const version = data.readUInt32LE(4);
  if (version !== VERSION) throw new Error(`${path} has unsupported k-mer file version ${version}`);
  const header: KmerFileHeader = {
    k: data.readUInt32LE(8),
    minCount: data.readUInt32LE(12),
    entries: data.readUInt32LE(16) + data.readUInt32LE(20) * 0x100000000
  };
  if (data.length !== HEADER_BYTES + header.entries * RECORD_BYTES) throw new Error(`${path} is truncated`);
  const table = new KmerTable(header.entries);
  for (let i = 0; i < header.entries; i++) {
    const offset = HEADER_BYTES + i * RECORD_BYTES;
    table.set(data.readUInt32LE(offset), data.readUInt32LE(offset + 4), data.readUInt32LE(offset + 8));
  }
  return { header, table };
}
