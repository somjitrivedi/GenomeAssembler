# GenomeAssembler

A de novo short-read genome assembler in TypeScript. It rebuilds a genome from millions of short DNA reads, like putting a shredded book back together.

```
reads (FASTQ) -> k-mer counting -> error correction -> de Bruijn graph
                 (Bloom filter)    (fix bad bases)     (k-mers = nodes)
                                                              |
scaffolds (FASTA) <- scaffolding <- contigs <- graph cleaning -+
                    (paired reads)            (tips, bubbles)
```

Run metadata, read statistics, k-mer spectra and evaluation results are stored in SQLite. The k-mers themselves stay in typed-array hash tables in memory, and solid k-mers are written to a compact binary file for the later graph phases.

Target scale: bacterial genomes (around 5 Mb, a few million reads).

## Status

| # | Phase | Status |
|---|---|---|
| 0 | Foundations and test harness | Done |
| 1 | k-mer counting with a Bloom filter | Done |
| 2 | Sequencing error correction | Not started |
| 3 | de Bruijn graph construction | Not started |
| 4 | Graph cleaning, then contigs | Not started |
| 5 | Scaffolding with paired reads | Not started |
| 6 | Scale and real data | Not started |
| 7 | Reporting and visualization (optional) | Not started |

## Setup

Requires Node 22 or newer (developed on Node 24).

```
npm install
npm run build
```

If npm blocks install scripts, approve and rebuild the native packages:

```
npm approve-scripts better-sqlite3 esbuild
npm rebuild better-sqlite3 esbuild
```

## Usage

```
node dist/cli.js <command> [options]
npm run ga -- <command> [options]
```

Every command records a run in the database: `--db <path>`, else `$GA_DB`, else `./genome-assembler.db`.

### simulate

Generates a random genome with planted repeats, plus paired-end reads (forward/reverse orientation) with a 3'-rising substitution error profile and matching Phred qualities.

```
node dist/cli.js simulate --out data/sim --length 5000000 --repeats 20 --repeat-length 2000 \
  --repeat-copies 3 --coverage 30 --read-length 150 --insert 400 --insert-sd 40 --error-rate 0.005 --seed 42
```

Writes `reference.fa`, `reads_1.fq`, `reads_2.fq` (`.fq.gz` with `--gzip`) and `repeats.tsv`. Read names encode the truth: `sim<index>_<start>_<fragment length>_<f|r>/<mate>`.

### stats

```
node dist/cli.js stats data/sim/reads_1.fq data/sim/reads_2.fq
```

Reports read count, bases, length range, GC, N bases and mean quality. Accepts FASTA or FASTQ, plain or gzip.

### count

```
node dist/cli.js count data/sim/reads_1.fq data/sim/reads_2.fq -k 31 --out data/sim/solid.kmers --histogram data/sim/hist.tsv
```

| Option | Meaning |
|---|---|
| `-k` | Odd k from 1 to 31 (default 31) |
| `--mode bloom` | Default. Three passes: HyperLogLog sizes the Bloom filter; the Bloom filter keeps first sightings out of the hash table; an exact recount removes Bloom false positives |
| `--mode naive` | Every k-mer goes into the hash table (baseline for comparison) |
| `--fpr` | Bloom filter false positive rate (default 0.01) |
| `--expected` | Known distinct k-mer count; skips the HyperLogLog pass |
| `--approx` | Skips the exact recount. Counts are exact except for k-mers whose first sighting was a Bloom false positive (+1) |
| `--min-count` | Override the automatic solid-k-mer cutoff |
| `--out` | Write solid k-mers (count >= cutoff) to a binary file |
| `--histogram` | Write the k-mer spectrum as TSV |

From the k-mer spectrum it derives:
- the solid cutoff (first valley)
- the peak and k-mer coverage (mean over a window around the peak)
- genome size (solid occurrences / k-mer coverage)
- per-base error rate (from the fraction of k-mer occurrences that are solid)
- read coverage, corrected for k-mers destroyed by errors

### evaluate

```
node dist/cli.js evaluate --contigs contigs.fa --reference data/sim/reference.fa --min-contig 500
```

Reports contigs, total length, largest, N50/L50, GC, and, with a reference, NG50/LG50, genome fraction, duplication ratio, misassemblies (relocations and inversions) and unaligned length.

Alignment chains unique reference k-mers along diagonals and tolerates small indels. A contig made entirely of repeat sequence has no unique anchors, so it counts as unaligned.

### runs / show

```
node dist/cli.js runs --limit 20
node dist/cli.js show 3
```

## Acceptance results

Simulated E. coli-sized data: 5 Mb genome, 20 repeat families × 3 copies × 2 kb, 30× coverage, 150 bp paired reads, 0.5% errors (1,000,000 reads).

| Check | Target | Result |
|---|---|---|
| 1M reads round-trip | works | 1,000,000 reads simulated, parsed and counted; observed error rate 0.501% |
| Evaluator, reference against itself | 100%, 0 misassemblies | 100.000%, 0, N50 = 5,000,000 |
| Bloom vs naive peak memory | at least 3× lower | 121.3 MiB vs 384.0 MiB (3.17×) |
| Genome size estimate | within 5% | 5,005,969 bp (0.12% off) |
| Read coverage / error rate | | 29.96× / 0.498% (true values 30× / 0.501%) |
| Solid k-mers, bloom vs naive | identical | 4,921,071 in both modes |

Peak memory counts the k-mer data structures (Bloom filter bits plus hash table slots), not the whole process.

## Database

| Table | Contents |
|---|---|
| `runs` | Every command, its arguments, status, timing and error |
| `simulations` | Simulation parameters, output paths and errors introduced |
| `read_sets` | Per-file read statistics |
| `kmer_counts` | Counting mode, Bloom parameters, memory, spectrum analysis |
| `kmer_histogram` | Full k-mer spectrum per run (singletons flagged as estimated in bloom mode) |
| `evaluations` | Assembly metrics |

## Solid k-mer file format

Little-endian.
- **Header (24 bytes):** magic `GAKM`, version, k, minimum count, entry count (two 32-bit words).
- **Records (12 bytes each):** k-mer high bits, k-mer low bits, count.

## Development

```
npm run typecheck
npm test
npm run build
```

On Windows, run `npm test` from PowerShell or cmd.

## Layout

```
src/dna      alphabet, 2-bit k-mer encoding, rolling canonical k-mers, hashing
src/io       streaming FASTA/FASTQ reader and writer (gzip), read statistics, read sources
src/kmer     hash table, Bloom filter, HyperLogLog, counter, spectrum analysis, solid k-mer file
src/sim      seeded RNG, genome simulator with repeats, paired-end read simulator
src/eval     N50-style statistics, unique-k-mer alignment, reference-based evaluation
src/db       SQLite migrations and repository
src/commands.ts, src/cli.ts
```
