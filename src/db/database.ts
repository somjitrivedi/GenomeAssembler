import Database from "better-sqlite3";
import { mkdirSync } from "node:fs";
import { dirname } from "node:path";

export type Db = Database.Database;

const MIGRATIONS: string[] = [
  `
  CREATE TABLE runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    command TEXT NOT NULL,
    args TEXT NOT NULL,
    status TEXT NOT NULL CHECK (status IN ('running', 'succeeded', 'failed')),
    started_at TEXT NOT NULL,
    finished_at TEXT,
    elapsed_ms REAL,
    error TEXT
  );

  CREATE TABLE simulations (
    run_id INTEGER PRIMARY KEY REFERENCES runs(id) ON DELETE CASCADE,
    seed INTEGER NOT NULL,
    genome_length INTEGER NOT NULL,
    gc REAL NOT NULL,
    repeat_families INTEGER NOT NULL,
    repeat_length INTEGER NOT NULL,
    repeat_copies INTEGER NOT NULL,
    repeat_divergence REAL NOT NULL,
    coverage REAL NOT NULL,
    read_length INTEGER NOT NULL,
    insert_mean REAL NOT NULL,
    insert_sd REAL NOT NULL,
    error_rate REAL NOT NULL,
    n_rate REAL NOT NULL,
    pairs INTEGER NOT NULL,
    errors_introduced INTEGER NOT NULL,
    reference_path TEXT NOT NULL,
    reads1_path TEXT NOT NULL,
    reads2_path TEXT NOT NULL
  );

  CREATE TABLE read_sets (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    path TEXT NOT NULL,
    reads INTEGER NOT NULL,
    bases INTEGER NOT NULL,
    min_length INTEGER NOT NULL,
    max_length INTEGER NOT NULL,
    mean_length REAL NOT NULL,
    gc REAL NOT NULL,
    n_bases INTEGER NOT NULL,
    mean_quality REAL
  );

  CREATE TABLE kmer_counts (
    run_id INTEGER PRIMARY KEY REFERENCES runs(id) ON DELETE CASCADE,
    k INTEGER NOT NULL,
    mode TEXT NOT NULL CHECK (mode IN ('bloom', 'naive')),
    exact INTEGER NOT NULL,
    inputs TEXT NOT NULL,
    reads INTEGER NOT NULL,
    total_kmers INTEGER NOT NULL,
    distinct_kmers INTEGER NOT NULL,
    distinct_estimated INTEGER NOT NULL,
    stored_kmers INTEGER NOT NULL,
    bloom_bits INTEGER,
    bloom_hashes INTEGER,
    bloom_target_fpr REAL,
    bloom_estimated_fpr REAL,
    bloom_false_positives_removed INTEGER,
    peak_structure_bytes INTEGER NOT NULL,
    passes INTEGER NOT NULL,
    cutoff INTEGER NOT NULL,
    cutoff_source TEXT NOT NULL,
    peak_multiplicity INTEGER NOT NULL,
    kmer_coverage REAL NOT NULL,
    read_coverage REAL,
    solid_fraction REAL,
    base_error_rate REAL,
    genome_size_estimate INTEGER NOT NULL,
    solid_kmers INTEGER NOT NULL,
    weak_kmers INTEGER NOT NULL,
    output_path TEXT,
    elapsed_ms REAL NOT NULL
  );

  CREATE TABLE kmer_histogram (
    run_id INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    multiplicity INTEGER NOT NULL,
    kmers INTEGER NOT NULL,
    estimated INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (run_id, multiplicity)
  );

  CREATE TABLE evaluations (
    run_id INTEGER PRIMARY KEY REFERENCES runs(id) ON DELETE CASCADE,
    contigs_path TEXT NOT NULL,
    reference_path TEXT,
    min_contig INTEGER NOT NULL,
    contigs INTEGER NOT NULL,
    total_length INTEGER NOT NULL,
    largest INTEGER NOT NULL,
    n50 INTEGER NOT NULL,
    l50 INTEGER NOT NULL,
    ng50 INTEGER,
    lg50 INTEGER,
    gc REAL NOT NULL,
    reference_length INTEGER NOT NULL,
    genome_fraction REAL,
    duplication_ratio REAL,
    misassemblies INTEGER,
    misassembled_contigs INTEGER,
    unaligned_length INTEGER
  );

  CREATE INDEX idx_runs_command ON runs(command);
  CREATE INDEX idx_read_sets_run ON read_sets(run_id);
  `
];

export function migrate(db: Db): number {
  const current = db.pragma("user_version", { simple: true }) as number;
  for (let version = current; version < MIGRATIONS.length; version++) {
    db.transaction(() => {
      db.exec(MIGRATIONS[version]);
      db.pragma(`user_version = ${version + 1}`);
    })();
  }
  return MIGRATIONS.length;
}

export function openDatabase(path: string): Db {
  if (path !== ":memory:") mkdirSync(dirname(path), { recursive: true });
  const db = new Database(path);
  db.pragma("journal_mode = WAL");
  db.pragma("foreign_keys = ON");
  migrate(db);
  return db;
}
