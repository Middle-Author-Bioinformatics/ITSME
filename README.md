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

Native SPAdes `NODE_...` contigs are always retained. ITSME also enumerates
bounded paths through the assembly graph, but promotes a path as an
`ITSME_LOCUS_...` only when it has the expected SSU–ITS1–5.8S–ITS2–LSU
structure, every graph junction has read-template support, and independent SSU
and LSU classifications agree at least at phylum. ITS classification is
secondary evidence: agreement strengthens a result, absence does not reject it,
and conflict is reported as a warning. All unpromoted paths remain available in
`validation/graph_paths/` for audit.

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
    --expected-taxonomy Echinodermata \
    --sensitivity 2
```

`--sensitivity` accepts `1` (specific), `2` (balanced; default), or `3`
(sensitive). `--expected-taxonomy` labels results but never removes non-target
sequences.

Run all paired libraries in a directory:

```bash
./itsme_controller.sh \
    -i /path/to/fastqs \
    -o /path/to/results \
    --db-dir /home/ark/databases/itsme_db \
    --expected-taxonomy Echinodermata \
    --sensitivity 2 \
    --resume
```

No `--` separator is required before ITSME options.

## Key outputs

| File | Contents |
| --- | --- |
| `master_summary.csv` | Concise native-contig and promoted-locus summary with region coordinates and consensus taxonomy |
| `final/complete_rDNA_loci.fasta` | Complete validated native and reconstructed rDNA loci |
| `final/reconstructed_graph_loci.fasta` | Graph paths promoted after structural, junction-support, and taxonomy checks |
| `final/graph_locus_validation.tsv` | Decision and reason for every bounded graph path |
| `final/partial_locus_contigs.fasta` | Native partial 18S/28S contigs retained from the assembly |
| `validation/graph_paths/graph_candidate_paths.fasta` | Every bounded path before promotion filtering |
| `run_summary.txt` | Run settings, counts, stopping reason, and elapsed time |

Use `./itsme.sh --help` and `./itsme_controller.sh --help` for all options.
