import { createReadStream, createWriteStream, type WriteStream } from "node:fs";
import { open } from "node:fs/promises";
import { createInterface } from "node:readline";
import { createGunzip, createGzip, type Gzip } from "node:zlib";
import type { Readable } from "node:stream";
import { once } from "node:events";

export interface FastxRecord {
  name: string;
  seq: string;
  qual?: string;
}

export type FastxFormat = "fasta" | "fastq";

async function isGzipped(path: string): Promise<boolean> {
  const handle = await open(path, "r");
  try {
    const header = Buffer.alloc(2);
    const { bytesRead } = await handle.read(header, 0, 2, 0);
    return bytesRead === 2 && header[0] === 0x1f && header[1] === 0x8b;
  } finally {
    await handle.close();
  }
}

export async function openTextStream(path: string): Promise<Readable> {
  const gzipped = await isGzipped(path);
  const raw = createReadStream(path, { highWaterMark: 1 << 20 });
  if (!gzipped) return raw;
  const gunzip = createGunzip();
  raw.on("error", (error) => gunzip.destroy(error));
  return raw.pipe(gunzip);
}

export async function forEachRecord(path: string, visit: (record: FastxRecord) => void): Promise<number> {
  const stream = await openTextStream(path);
  const lines = createInterface({ input: stream, crlfDelay: Infinity });
  let format: FastxFormat | undefined;
  let lineNo = 0;
  let records = 0;
  let state = 0;
  let name = "";
  let seq = "";
  let seqParts: string[] = [];
  const fail = (message: string): never => {
    throw new Error(`${path}:${lineNo}: ${message}`);
  };
  for await (const rawLine of lines) {
    lineNo++;
    const line = rawLine.endsWith("\r") ? rawLine.slice(0, -1) : rawLine;
    if (format === undefined) {
      if (line.length === 0) continue;
      if (line[0] === ">") format = "fasta";
      else if (line[0] === "@") format = "fastq";
      else fail("not a FASTA or FASTQ file (first record must start with '>' or '@')");
    }
    if (format === "fasta") {
      if (line.length > 0 && line[0] === ">") {
        if (state === 1) {
          visit({ name, seq: seqParts.join("") });
          records++;
        }
        name = line.slice(1).trim();
        seqParts = [];
        state = 1;
      } else if (line.length > 0) {
        seqParts.push(line.trim());
      }
      continue;
    }
    switch (state) {
      case 0:
        if (line.length === 0) continue;
        if (line[0] !== "@") fail("expected FASTQ header starting with '@'");
        name = line.slice(1).trim();
        state = 1;
        break;
      case 1:
        seq = line.trim();
        state = 2;
        break;
      case 2:
        if (line[0] !== "+") fail("expected FASTQ separator line starting with '+'");
        state = 3;
        break;
      default: {
        const qual = line.trim();
        if (qual.length !== seq.length) fail(`quality length ${qual.length} does not match sequence length ${seq.length}`);
        visit({ name, seq, qual });
        records++;
        state = 0;
      }
    }
  }
  if (format === "fasta" && state === 1) {
    visit({ name, seq: seqParts.join("") });
    records++;
  }
  if (format === "fastq" && state !== 0) fail("truncated FASTQ record at end of file");
  return records;
}

export async function readAllRecords(path: string): Promise<FastxRecord[]> {
  const out: FastxRecord[] = [];
  await forEachRecord(path, (record) => out.push(record));
  return out;
}

const FLUSH_BYTES = 1 << 22;

export class FastxWriter {
  private readonly file: WriteStream;
  private readonly sink: WriteStream | Gzip;
  private buffer: string[] = [];
  private buffered = 0;

  constructor(readonly path: string) {
    this.file = createWriteStream(path);
    if (path.endsWith(".gz")) {
      const gzip = createGzip({ level: 4 });
      gzip.pipe(this.file);
      this.sink = gzip;
    } else {
      this.sink = this.file;
    }
  }

  private async push(text: string): Promise<void> {
    this.buffer.push(text);
    this.buffered += text.length;
    if (this.buffered >= FLUSH_BYTES) await this.flush();
  }

  private async flush(): Promise<void> {
    if (this.buffered === 0) return;
    const chunk = this.buffer.join("");
    this.buffer = [];
    this.buffered = 0;
    if (!this.sink.write(chunk)) await once(this.sink, "drain");
  }

  async fasta(name: string, seq: string, width = 80): Promise<void> {
    const lines = [`>${name}\n`];
    for (let i = 0; i < seq.length; i += width) lines.push(seq.slice(i, i + width), "\n");
    await this.push(lines.join(""));
  }

  async fastq(name: string, seq: string, qual: string): Promise<void> {
    await this.push(`@${name}\n${seq}\n+\n${qual}\n`);
  }

  async close(): Promise<void> {
    await this.flush();
    const finished = once(this.file, "finish");
    this.sink.end();
    await finished;
  }
}

export async function writeFasta(path: string, records: FastxRecord[], width = 80): Promise<void> {
  const writer = new FastxWriter(path);
  for (const record of records) await writer.fasta(record.name, record.seq, width);
  await writer.close();
}
