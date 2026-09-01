# ITSME

**Iterative Targeted Sequence Mining and Extension**

ITSME is a targeted assembly pipeline for recovering eukaryotic nuclear ribosomal DNA from paired-end genomic reads. It recruits rRNA-associated reads, performs post-recruitment quality control, assembles the enriched read set, validates the expected rDNA structure, and reports taxonomic matches against local NCBI databases.

```text
18S (SSU) → ITS1 → 5.8S → ITS2 → 28S (LSU)
```

> **Version:** 1.0.0  
> **Primary application:** assembly and classification of complete or partial eukaryotic rDNA loci from short-read WGS data

## Key features

- Maps raw reads before QC, so trimming and merging are limited to the recruited subset.
- Supports multiple rRNA seed databases or separate 18S and 28S references.
- Retains strong seed hits and their read mates.
- Performs a localized SPAdes assembly.
- Validates the order, orientation, and compatibility of 18S and 28S anchors.
- Automatically orients accepted contigs from 18S to 28S.
- Uses ITSx to annotate SSU, ITS1, 5.8S, ITS2, and LSU coordinates.
- Maps reads back competitively to measure coverage and depth.
- Classifies SSU, ITS, and LSU sequences against targeted NCBI BLAST databases.
- Reports taxonomic ranks from domain through species.
- Produces an integrated master summary for every assembled locus.

## How the pipeline works

### 1. Seed recruitment

Bowtie2 maps the original reads against the supplied rRNA seed sequences. Alignments are filtered by aligned-read fraction and approximate identity. Accepted reads and their mates are retained as candidate rDNA templates.

### 2. Post-mapping QC

Only recruited reads are processed by Trimmomatic and FLASH. Trimmomatic removes adapters and low-quality sequence, while FLASH merges overlapping paired reads.

### 3. Inward sequence extension

One guarded exact-word extension round is enabled by default. When both `--left-db` and `--right-db` are supplied, ITSME automatically uses **inward extension**.

The left database represents the 18S/SSU side of the locus, and the right database represents the 28S/LSU side:

```text
left database                                               right database
18S / SSU  ────────────────→  ITS1 → 5.8S → ITS2  ────────────────→  28S / LSU
                 inward extension →       ← inward extension
```

ITSME begins with reads associated with the 3′ end of the left/18S references and the 5′ end of the right/28S references. These are the seed ends facing the internal ITS region. Their exact words recruit overlapping reads from the original library. Newly accepted reads form the next frontier, allowing recruitment to move from both conserved rRNA genes toward ITS1, 5.8S, and ITS2.

Each round uses only the newest accepted frontier rather than every previously accepted read. Excessive growth or an excessive total accepted fraction causes the entire proposed round to be rejected before assembly.

### 4. Outward or bidirectional extension

`--outward` and `--bidirectional` switch to a broader mode. Instead of starting only from the inward-facing terminal portions of the left and right references, this mode uses the complete seed-hit read pool as the initial frontier. It can therefore extend:

- inward toward the ITS region;
- outward from the 5′ side of 18S; and
- outward from the 3′ side of 28S.

This can recover additional flanking rDNA sequence, but it requires more care. Full-length SSU and LSU contain highly conserved sequence, so unrestricted frontier words are more likely to recruit reads from other rRNA copies, organisms, repeats, or unrelated genomic regions. Outward extension should normally use one round, 31-base words, multiple required word hits, and conservative growth limits.

### 5. Assembly

SPAdes assembles the final recruited read set once. This avoids repeatedly assembling progressively larger datasets.

### 6. rDNA validation and orientation

BLASTN compares assembled contigs with the supplied 18S and 28S references. A validated dual-anchor contig must contain compatible, non-overlapping anchors in the expected biological order. Reverse-strand contigs are automatically reverse-complemented.

### 7. ITS annotation

ITSx annotates the structural components of each oriented locus:

- SSU/18S
- ITS1
- 5.8S
- ITS2
- LSU/28S

Complete ITS regions and their individual components are written to separate FASTA files.

### 8. Read-back validation

The original reads are mapped competitively against all validated candidate contigs. Samtools reports coverage breadth, mean depth, mean mapping quality, and pairing statistics.

### 9. Taxonomic assessment

Extracted SSU, ITS, and LSU sequences are searched against their corresponding NCBI databases. Taxids are resolved using the NCBI taxonomy dump, producing classifications from domain through species.

## Dependencies

| Dependency | Purpose | Requirement |
|---|---|---|
| Bash | Pipeline execution and file handling | Required |
| Bowtie2 | Seed recruitment and read-back mapping | Required |
| Samtools | BAM processing, read extraction, coverage, and mapping statistics | Required |
| SPAdes | Localized assembly of recruited reads | Required |
| NCBI BLAST+ | Anchor validation and taxonomic searches | Required |
| Python 3 | Sequence processing, orientation, taxonomy, and report generation | Required |
| Trimmomatic | Adapter and quality trimming after recruitment | Required unless post-map QC is skipped |
| FLASH or FLASH2 | Merging overlapping read pairs | Required unless post-map QC is skipped |
| ITSx | rDNA structure annotation and ITS extraction | Required for structural and NCBI locus reporting |
| BBMap/BBDuk | Default exact-word extension | Required unless extension is disabled with `-R 0` |
| Java | Runs Trimmomatic and BBDuk in typical installations | Conditional |
| pigz | Parallel compression | Optional |
| seqkit | Convenient FASTA/FASTQ inspection | Optional |

## Installation

```bash
mamba create -n itsme \
    -c conda-forge \
    -c bioconda \
    bowtie2 \
    samtools \
    spades \
    blast \
    itsx \
    bbmap \
    trimmomatic \
    flash2 \
    pigz \
    seqkit \
    --yes

conda activate itsme
chmod +x itsme.sh
```

## Required databases

The standard input is a single database root supplied with `--db-dir`. For a flat installation, place all required data directly inside that directory:

- `silva-euk-18s-id95.fasta` — left/18S/SSU database
- `silva-euk-28s-id98.fasta` — right/28S/LSU database
- `SSU_eukaryote_rRNA`
- `ITS_eukaryote_sequences`
- `LSU_eukaryote_rRNA`
- `nodes.dmp`
- `names.dmp`
- `merged.dmp` — optional but recommended

ITSME also recognizes the NCBI BLAST databases under `NCBI_rRNA_BLAST/` and taxonomy files under `taxdump/` within the database root. Explicit database options remain available for backward compatibility and override automatic discovery. The seed FASTAs must remain directly under the database root.

## Example command

```bash
./itsme.sh \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    --db-dir /home/ark/databases/itsme_db \
    -o itsme_sample
```

Optional unpaired reads may be supplied with `-s` or `--single`.

The output directory must be new or empty.

## Batch controller

`itsme_controller.sh` runs ITSME over a directory of paired FASTQ libraries. It discovers filenames containing `_R1`, derives each mate by replacing the first `_R1` with `_R2`, and uses everything before `_R1` as the library prefix and output-subdirectory name.

```bash
./itsme_controller.sh \
    -i /path/to/fastqs \
    -o /path/to/itsme_results \
    --resume \
    -- \
    --db-dir /home/ark/databases/itsme_db
```

Arguments after `--` are passed directly to `itsme.sh`. The controller supplies `-1`, `-2`, and `-o` for each library and runs libraries sequentially to prevent CPU and memory oversubscription.

The controller creates:

- one result subdirectory per library;
- `logs/LIBRARY.log` for each run;
- `batch_status.tsv` with success, failure, and runtime information; and
- `batch_master_summary.csv`, combining all available per-library master summaries with a leading `library` column.

Use `--dry-run` to inspect commands without running them. Use `--resume` to skip output directories containing both `run_summary.txt` and `master_summary.csv`.

## Important options

### Input and databases

| Option | Description |
|---|---|
| `-1`, `--reads1 FILE` | Forward paired-end reads |
| `-2`, `--reads2 FILE` | Reverse paired-end reads |
| `-s`, `--single FILE` | Optional single or unpaired reads |
| `-o`, `--output-dir DIR` | New or empty output directory |
| `--db-dir DIR` | Unified ITSME database root containing seeds, NCBI BLAST databases, and taxonomy files |
| `-d`, `--database FILE` | General rRNA seed database; repeatable |
| `--left-db FILE` | 18S seed database; repeatable |
| `--right-db FILE` | 28S seed database; repeatable |
| `--inward` | Force inward-only initialization from the 18S 3′ and 28S 5′ ends |
| `--outward`, `--bidirectional` | Use the complete seed-hit frontier and extend in all directions |

### Recruitment and QC

| Option | Default | Description |
|---|---:|---|
| `-R`, `--max-rounds` | `1` | Maximum extension rounds; zero disables extension |
| `-w`, `--word-size` | `31` | Exact word length for BBDuk extension |
| `--min-word-hits` | `3` | Minimum shared words required per read |
| `--max-round-growth` | `0.25` | Reject excessive recruitment growth |
| `--max-accepted-fraction` | `0.03` | Maximum accepted fraction of the raw library |
| `--min-read-aligned` | `0.75` | Minimum aligned fraction for a seed-hit read |
| `--min-read-identity` | `0.85` | Minimum approximate seed-hit identity |
| `--skip-post-qc` | Off | Skip Trimmomatic and FLASH |
| `--qc-quality` | `20` | Sliding-window quality threshold |
| `--min-read-length` | `50` | Minimum retained read length |

### Assembly and validation

| Option | Default | Description |
|---|---:|---|
| `-t`, `--threads` | `8` | Number of threads |
| `-m`, `--memory` | `32` | SPAdes memory limit in GB |
| `-k`, `--kmers` | `auto` | SPAdes k-mer list |
| `--spades-mode` | `standard` | SPAdes mode: `standard` or `meta` |
| `--anchor-min-aligned` | `150` | Minimum anchor alignment length |
| `--anchor-min-identity` | `0.85` | Minimum anchor identity |
| `--itsx-taxa` | `M` | ITSx taxonomic profile; `M` is Metazoa |
| `--skip-readback` | Off | Skip competitive read-back mapping |
| `--skip-ncbi-blast` | Off | Skip NCBI taxonomic assessment |
| `--keep-bam` | Off | Retain the large initial mapping BAM |

Run the following for the complete option list:

```bash
./itsme.sh --help
```

## Output overview

### Primary reports

| File | Description |
|---|---|
| `run_summary.txt` | Run-level recruitment, assembly, validation, and elapsed-time statistics |
| `recruitment.tsv` | Recruitment history and extension decisions |
| `master_summary.csv` | Spreadsheet-friendly master summary written directly to the root output directory |
| `final/master_summary.tsv` | Integrated per-locus structure, coverage, and taxonomy report |
| `final/contig_validation.tsv` | 18S/28S anchor coordinates, identity, orientation, and layout status |
| `final/itsx_validation.tsv` | SSU, ITS1, 5.8S, ITS2, and LSU coordinates and structural status |
| `final/readback_coverage.tsv` | Coverage breadth, mean depth, and mapping quality |
| `final/ncbi_blast_top_hits.tsv` | Ranked SSU, ITS, and LSU hits with domain-to-species taxonomy |

### Primary FASTA files

| File | Description |
|---|---|
| `final/oriented_dual_anchor_contigs.fasta` | Validated rDNA contigs oriented 18S to 28S |
| `final/complete_ITS.fasta` | Complete ITS1–5.8S–ITS2 regions |
| `final/ITS1.fasta` | Extracted ITS1 sequences |
| `final/5_8S.fasta` | Extracted 5.8S sequences |
| `final/ITS2.fasta` | Extracted ITS2 sequences |
| `final/rrna_candidate_contigs.fasta` | All contigs with qualifying seed evidence |
| `final/all_assembled_contigs.fasta` | Complete SPAdes assembly output |

### Supporting directories

| Directory | Contents |
|---|---|
| `seed/` | Combined seed references, Bowtie2 indexes, mapping logs, recruited reads, and QC results |
| `assembly/` | SPAdes output and assembly logs |
| `validation/` | Anchor searches, ITSx output, read-back BAM files, and NCBI database information |
| `final/` | Main FASTA files and tabular reports |

## Interpreting results

A strong complete locus should generally show:

1. `PASS_dual_anchor` layout validation.
2. `PASS_complete_ITS` structural validation.
3. Ordered SSU, ITS1, 5.8S, ITS2, and LSU coordinates.
4. Broad read coverage across the assembled locus.
5. Compatible SSU, ITS, and LSU taxonomic matches.

The best BLAST match should not be treated automatically as a species identification. Consider percent identity, aligned length, query coverage, E-value, and agreement among all three markers. ITS is often more discriminatory than SSU or LSU but is less comprehensively represented in public databases.

Multiple validated loci are allowed. They may represent different organisms, symbionts, environmental DNA, divergent rDNA copies, or assembly alternatives.

## Extension guidance

The default is one inward extension round when both left and right databases are present. This is the recommended mode for recovering the 18S–ITS1–5.8S–ITS2–28S locus.

Use `-R 0` to perform seed recruitment and assembly without extension. Use additional rounds only when the locus remains incomplete and accepted recruitment remains controlled.

For outward extension, begin conservatively:

```bash
./itsme.sh \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    --db-dir /home/ark/databases/itsme_db \
    --outward \
    -R 1 \
    -w 31 \
    --min-word-hits 3 \
    --max-round-growth 0.10 \
    --max-accepted-fraction 0.01 \
    -o itsme_sample_outward
```

Inspect `recruitment.tsv` after every outward run. A rejected round is a successful safety stop: none of the proposed new reads are admitted to the final assembly.

## Limitations

- ITSME is inspired by GetOrganelle's bait-and-extension strategy but does not reproduce its graph-disentangling or multiplicity inference.
- Short reads may not resolve nearly identical rDNA copies or complex mixtures.
- rDNA depth is not a direct measure of organism abundance because copy number varies among taxa.
- Classification depends on the completeness and accuracy of the reference databases.
- A complete ITS structure does not by itself establish species identity.
- The default ITSx profile is Metazoa; other eukaryotic groups require an appropriate `--itsx-taxa` setting.
- ITSME targets nuclear rDNA. Mitochondrial and COI analyses should be handled separately.

## Reproducibility

For each analysis, retain:

- the ITSME script and version;
- the exact command line;
- `run_summary.txt` and `recruitment.tsv`;
- `master_summary.csv`, `final/master_summary.tsv`, and the primary FASTA files;
- software versions;
- seed database names and versions;
- NCBI BLAST database update dates; and
- the NCBI taxonomy dump date.

When reporting results, cite the underlying software and databases used, including Bowtie2, Samtools, SPAdes, BLAST+, ITSx, Trimmomatic, FLASH, BBMap/BBDuk when applicable, the seed database source, and NCBI taxonomy resources.
