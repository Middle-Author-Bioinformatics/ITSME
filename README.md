# ITSME

**Iterative Targeted Sequence Mining and Extension**

```text
18S (SSU) ── ITS1 ── 5.8S ── ITS2 ── 28S (LSU)
```

ITSME recovers eukaryotic nuclear rDNA loci from paired short reads, including
mixed environmental samples. It maps raw reads to 18S and 28S seeds, performs
graph-aware inward recruitment and one conservative outward pass, quality
filters only recruited reads, and assembles the enriched pool with metaSPAdes.

```mermaid
flowchart LR
    A[Raw reads] --> B[18S and 28S baiting]
    B --> C[Guarded extension]
    C --> D[metaSPAdes graph]
    D --> E[Validated rDNA loci]
```

Native SPAdes `NODE_...` contigs are retained when they meet the default 150-bp
reporting minimum. ITSME also enumerates bounded assembly-graph paths and
promotes them as `ITSME_LOCUS_...` sequences when they have the expected
SSU–ITS1–5.8S–ITS2–LSU structure, supported graph junctions, and concordant SSU
and LSU taxonomy. When several reconstructed paths have the same most-resolved
consensus taxonomy, only the longest is reported as the representative. Every
original path remains in `validation/graph_paths/` for audit. SSU/LSU conflicts
at domain, kingdom, or phylum are excluded from the master summary and reported
FASTA, but retained with both assignments in dedicated chimera files.

## Install

```bash
mamba create -n itsme -c conda-forge -c bioconda \
    bowtie2 samtools bcftools spades blast itsx bbmap \
    trimmomatic flash2 pigz seqkit --yes
mamba activate itsme
chmod +x itsme.sh itsme_controller.sh setup_db.sh
```

The database directory supplied with `--db-dir` must contain:

```text
itsme_db/
├── silva-euk-18s-id95.fasta
├── silva-euk-28s-id98.fasta
├── NCBI_rRNA_BLAST/
│   ├── SSU_eukaryote_rRNA.*
│   ├── ITS_eukaryote_sequences.*
│   ├── LSU_eukaryote_rRNA.*
│   └── taxdb.*
└── taxdump/
    ├── nodes.dmp
    ├── names.dmp
    └── merged.dmp
```

Create or populate it with `./setup_db.sh`, optionally using
`--db-dir /path/to/itsme_db`.

## Quick start

Run one library with the balanced sensitivity preset:

```bash
./itsme.sh \
    -1 sample_R1_001.fastq.gz \
    -2 sample_R2_001.fastq.gz \
    -o itsme_sample \
    --db-dir /home/ark/databases/itsme_db \
    --sensitivity 2
```

`--sensitivity` accepts `1` (specific), `2` (balanced; default), or `3`
(sensitive). The minimum reported full contig or reconstructed locus is 150 bp;
change it with `--min-report-length`. Individual ITS subregions may be shorter.

Run all paired libraries in a directory:

```bash
./itsme_controller.sh \
    -i /path/to/fastqs \
    -o /path/to/results \
    --db-dir /home/ark/databases/itsme_db \
    --sensitivity 2 \
    --resume
```

No `--` separator is required before ITSME options.

## Key outputs

| File | Contents |
| --- | --- |
| `master_summary.csv` | Complete nonchimeric ITS-containing loci with region coordinates and consensus taxonomy |
| `partials.csv` | Incomplete 18S/28S and other seed-anchored candidates without a validated complete ITS cassette |
| `final/reported_loci.fasta` | Complete sequences corresponding one-for-one to the master-summary rows, with annotated region coordinates |
| `final/complete_rDNA_loci.fasta` | Complete nonchimeric native and reconstructed rDNA loci |
| `final/reconstructed_graph_loci.fasta` | Longest representative of each concordant reconstructed taxonomic group |
| `final/taxonomic_chimeras.fasta` | Structurally recovered loci excluded for SSU/LSU disagreement at phylum or above |
| `final/taxonomic_chimeras.tsv` | Conflict rank and separate SSU/LSU assignments for excluded chimeras |
| `final/graph_locus_validation.tsv` | Decision, SSU/LSU taxonomy, and representative mapping for every bounded graph path |
| `final/partial_locus_contigs.fasta` | Native partial 18S/28S contigs retained from the assembly |
| `validation/graph_paths/graph_candidate_paths.fasta` | Every bounded path before promotion filtering |
| `validation/graph_paths/collapsed_same_taxonomy_paths.fasta` | Redundant reconstructed paths omitted in favor of the longest taxonomic representative |
| `run_summary.txt` | Run settings, counts, stopping reason, and elapsed time |

Use `./itsme.sh --help` and `./itsme_controller.sh --help` for all options.
