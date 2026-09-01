#!/usr/bin/env bash

# ITSME v1: variant-aware eukaryotic rDNA recruitment, assembly, and classification
#
# Recruit eukaryotic rRNA reads with a GetOrganelle-like workflow:
#   1. Map raw reads once to one or more full-length rRNA seed databases.
#   2. Retain strong seed hits and their mates.
#   3. QC only recruited reads.
#   4. Extend inward and conservatively outward with frontier-only exact k-mers.
#   5. Reject abnormal recruitment before it can enter the assembly.
#   6. Assemble the final accepted pool once with SPAdes.
#   7. Enumerate bounded 18S-to-28S paths through the SPAdes graph.
#   8. Report complete paths, partial loci, read support, and residual variants.
#
# This is not GetOrganelle. Graph paths are read-supported candidate sequences;
# short reads cannot prove long-range phase across every rDNA repeat variant.

set -Eeuo pipefail
IFS=$'\n\t'
SECONDS=0

VERSION="1.3.0"

READS_1=""
READS_2=""
READS_SINGLE=""
OUTDIR=""
DB_DIR=""
THREADS=8
SPADES_MEMORY_GB=32
BBTOOLS_MEMORY="16g"
MAX_ROUNDS=3
WORD_SIZE=31
MIN_WORD_HITS=3
BAIT_MIN_ENTROPY=0.45
MAX_ROUND_GROWTH=3.0
MAX_ACCEPTED_FRACTION=0.07
MAX_ACCEPTED_TEMPLATES=250000
OUTWARD_MIN_WORD_HITS=5
MAX_OUTWARD_TEMPLATES=25000
MAX_OUTWARD_GROWTH=0.25
MIN_NEW_TEMPLATES=10
MIN_GROWTH=0.0001
MAX_INSERT=1000
SEED_SCORE_MIN="G,20,8"
MIN_READ_ALIGNED_FRACTION=0.75
MIN_READ_IDENTITY=0.85
ANCHOR_MIN_ALIGNED=150
ANCHOR_MIN_IDENTITY=0.85
INWARD=false
OUTWARD=false
EXTENSION_DIRECTION="auto"
INWARD_END_LENGTH=400
SPADES_MODE="meta"
KMERS="auto"
POST_MAP_QC=true
QC_QUALITY=20
MIN_READ_LENGTH=50
FLASH_MIN_OVERLAP=10
FLASH_MAX_OVERLAP=150
FLASH_MISMATCH_DENSITY=0.25
ADAPTER_CLIP="2:30:10"
ADAPTERS=""
TRIMMOMATIC_CMD="trimmomatic"
TRIMMOMATIC_JAR=""
FLASH_CMD="flash"
BBDUK_CMD="bbduk.sh"
KEEP_INTERMEDIATES=false
KEEP_BAM=false
RUN_ITSX=true
ITSX_TAXA="M"
READBACK=true
LAYOUT_OVERLAP_TOLERANCE=25
RUN_NCBI_BLAST=true
NCBI_DB_DIR="${BLASTDB:-}"
NCBI_DB_DIR_EXPLICIT=false
NCBI_MAX_HITS=50
NCBI_REPORT_HITS=10
TAXDUMP_DIR="${NCBI_TAXDUMP_DIR:-}"
TAXDUMP_DIR_EXPLICIT=false
EXPECTED_TAXONOMY=""
TAXONOMY_NEAR_TOP_FRACTION=0.95
GRAPH_PATHS=true
MAX_GRAPH_PATHS=500
MAX_GRAPH_NODES=40
MAX_GRAPH_BP=25000
MIN_GRAPH_DEPTH=1.0
ADAPTIVE_GRAPH_STOP=true
GRAPH_CHECK_K=31
GRAPH_STOP_BRANCHES=2
GRAPH_STOP_CYCLES=1
GRAPH_STOP_NODES=500
CALL_VARIANTS=true
VARIANT_MIN_MAPQ=10
VARIANT_MIN_BASEQ=20
VARIANT_MAX_DEPTH=100000

DATABASES=()
LEFT_DATABASES=()
RIGHT_DATABASES=()
ASSEMBLY_R1_FILES=()
ASSEMBLY_R2_FILES=()
ASSEMBLY_SINGLE_FILES=()

usage() {
    cat <<'EOF'
ITSME v1.3.0 - variant-aware eukaryotic rDNA recruitment and classification

Usage:
  itsme.sh -1 R1.fastq.gz -2 R2.fastq.gz \
      --db-dir RRNA_DB_DIR -o OUTPUT [options]

Inward 18S-to-28S mode:
  itsme.sh -1 R1.fastq.gz -2 R2.fastq.gz \
      --left-db 18S.fasta --right-db 28S.fasta \
      -o OUTPUT [options]

Required:
  -1, --reads1 FILE                Raw forward reads
  -2, --reads2 FILE                Raw reverse reads
  -o, --output-dir DIR             New or empty output directory

Seed databases:
      --db-dir DIR                  ITSME database root containing:
                                      silva-euk-18s-id95.fasta (left/18S)
                                      silva-euk-28s-id98.fasta (right/28S)
                                      NCBI targeted BLAST databases
                                      nodes.dmp and names.dmp taxonomy files
  -d, --database FILE              General rRNA seed FASTA; repeatable
      --left-db FILE               18S seed FASTA; repeatable
      --right-db FILE              28S seed FASTA; repeatable
      --inward                     Extend only toward ITS from the 3' end of
                                   left/18S and 5' end of right/28S seeds
      --outward                    Extend only from the outward-facing 18S 5'
                                   and 28S 3' seed ends
      --bidirectional              Extend from all four terminal seed ends
                                   [default with left and right databases]
      --inward-end-length INT      Terminal seed bases used to start inward
                                   extension [400]

Recruitment:
  -R, --max-rounds INT             Maximum word-extension rounds [3]
  -w, --word-size INT              Exact k-mer/word length [31]
      --min-word-hits INT          Minimum bait words required per read [3]
      --bait-min-entropy FLOAT     Exclude low-complexity bait sequences [0.45]
      --max-round-growth FLOAT     Reject if new/accepted exceeds this [3.0]
      --max-accepted-fraction FLOAT
                                   Reject if accepted/all exceeds this [0.07]
      --max-accepted-templates INT Absolute recruitment ceiling; 0 disables
                                   the count ceiling [250000]
      --outward-min-word-hits INT  Exact word hits required for the one-pass
                                   outward extension [5]
      --max-outward-templates INT  Maximum new templates accepted outward;
                                   0 disables this count cap [25000]
      --max-outward-growth FLOAT   Reject outward pass if new/accepted exceeds
                                   this value [0.25]
      --min-new-templates INT      Stop below this number of new templates [10]
      --min-growth FLOAT           Stop below this fractional growth [0.0001]
      --bbtools-memory SIZE        Java heap for BBDuk, e.g. 16g [16g]
      --bbduk CMD                  BBDuk executable [bbduk.sh]

Seed alignment:
      --seed-score-min FUNC        Bowtie2 local score minimum [G,20,8]
      --min-read-aligned FLOAT     Minimum aligned fraction of a seed-hit read [0.75]
      --min-read-identity FLOAT    Minimum approximate seed-hit identity [0.85]
      --max-insert INT             Maximum Bowtie2 fragment length [1000]

Post-mapping QC (enabled by default):
      --skip-post-qc               Do not trim or merge recruited reads
      --qc-quality INT             Trimmomatic sliding-window quality [20]
      --min-read-length INT        Minimum retained read length [50]
      --adapters FILE              Adapter FASTA; otherwise auto-detected
      --adapter-clip PARAMS        ILLUMINACLIP parameters [2:30:10]
      --trimmomatic CMD            Trimmomatic executable [trimmomatic]
      --trimmomatic-jar FILE       Use java -jar FILE instead
      --flash CMD                  FLASH executable [flash; fallback flash2]
      --flash-min-overlap INT      Minimum overlap [10]
      --flash-max-overlap INT      Maximum overlap [150]
      --flash-mismatch FLOAT       Maximum mismatch density [0.25]

Assembly and output:
  -t, --threads INT                Threads [8]
  -m, --memory INT                 SPAdes memory limit in GB [32]
  -k, --kmers LIST                 SPAdes k-mers or auto [auto]
      --spades-mode MODE           standard or meta [meta]
      --anchor-min-aligned INT     Minimum seed alignment on a contig [150]
      --anchor-min-identity FLOAT  Minimum contig seed-hit identity [0.85]
      --layout-overlap INT         Maximum allowed 18S/28S overlap [25]
      --skip-graph-paths           Do not enumerate alternative GFA paths
      --max-graph-paths INT        Maximum complete graph paths reported [500]
      --max-graph-nodes INT        Maximum segments in one graph path [40]
      --max-graph-bp INT           Maximum candidate path length [25000]
      --min-graph-depth FLOAT      Minimum segment k-mer depth in a path [1.0]
      --skip-adaptive-graph-stop   Skip fast graph checks between inward rounds
      --graph-check-k INT          Single k-mer used by checkpoint SPAdes [31]
      --graph-stop-branches INT    Stop at this many branch nodes in a connected
                                   18S-to-28S component [2]
      --graph-stop-cycles INT      Stop at this cycle rank in a connected
                                   18S-to-28S component [1]
      --graph-stop-nodes INT       Emergency stop when a seed-anchored component
                                   reaches this many segments [500]
      --itsx-taxa CODE             ITSx taxon profile [M, Metazoa]
      --skip-itsx                  Skip ITSx structural annotation
      --skip-readback              Skip competitive raw-read mapping to candidates
      --skip-variant-calling       Do not call residual variants after read-back
      --variant-min-mapq INT       Minimum mapping quality for variants [10]
      --variant-min-baseq INT      Minimum base quality for variants [20]
      --variant-max-depth INT      Maximum pileup depth per position [100000]
      --ncbi-db-dir DIR           Directory containing the NCBI targeted
                                   BLAST databases [auto-detected from --db-dir]
      --ncbi-max-hits INT         Maximum BLAST hits retained per query [50]
      --ncbi-report-hits INT      Top hits per locus in summary table [10]
      --expected-taxonomy NAME    Expected clade at any rank, e.g.
                                   Echinodermata; flags but never removes loci
      --ncbi-taxdump-dir DIR      NCBI taxdump directory containing nodes.dmp,
                                   names.dmp and optionally merged.dmp
                                   [auto-detected from --db-dir]
      --skip-ncbi-blast           Skip local NCBI SSU/ITS/LSU classification
      --keep-bam                   Keep the large all-read seed BAM
      --keep-intermediates         Keep raw recruited and temporary FASTQs
  -h, --help                       Show this help
      --version                    Show the version

Principal outputs:
  OUTPUT/final/rrna_candidate_contigs.fasta
  OUTPUT/final/graph_candidate_paths.fasta
  OUTPUT/final/all_locus_candidates.fasta
  OUTPUT/final/partial_18S_contigs.fasta
  OUTPUT/final/partial_28S_contigs.fasta
  OUTPUT/final/rrna_dual_anchor_contigs.fasta
  OUTPUT/final/oriented_dual_anchor_contigs.fasta
  OUTPUT/final/complete_ITS.fasta
  OUTPUT/final/contig_validation.tsv
  OUTPUT/final/readback_coverage.tsv
  OUTPUT/final/path_support.tsv
  OUTPUT/final/residual_variants.vcf.gz
  OUTPUT/final/residual_variants.tsv
  OUTPUT/final/ncbi_SSU_hits.tsv
  OUTPUT/final/ncbi_ITS_hits.tsv
  OUTPUT/final/ncbi_LSU_hits.tsv
  OUTPUT/final/ncbi_blast_top_hits.tsv
  OUTPUT/final/taxonomy_validation.tsv
  OUTPUT/final/rrna_locus_summary.tsv
  OUTPUT/final/master_summary.tsv
  OUTPUT/master_summary.csv
  OUTPUT/final/all_assembled_contigs.fasta
  OUTPUT/final/accepted_R1.fastq.gz
  OUTPUT/final/accepted_R2.fastq.gz
  OUTPUT/final/accepted_single.fastq.gz
  OUTPUT/recruitment.tsv
  OUTPUT/recruitment_graph.tsv
  OUTPUT/run_summary.txt

Dependencies:
  bowtie2, samtools, bcftools, SPAdes, BLAST+, Python 3, BBMap/BBDuk,
  Trimmomatic and FLASH; ITSx, seqkit and pigz optional

Conda installation:
  mamba install -c conda-forge -c bioconda \
      bowtie2 samtools bcftools spades blast itsx bbmap trimmomatic flash2 pigz

Notes:
  * Raw reads are mapped before QC. Trimmomatic and FLASH see only accepted
    rRNA-associated reads, not the full input library.
  * Up to three guarded inward rounds are enabled by default. Each round scans
    the original reads with exact words from only the newest inward frontier.
    A fast single-k graph check stops more inward recruitment once the target
    connects or a connected 18S-to-28S component branches.
  * With left and right databases, the default is bidirectional: inward from
    the 18S 3' and 28S 5' ends, plus conservative outward extension from the
    18S 5' and 28S 3' ends. Outward extension is a stricter, nonrecursive
    one-pass operation capped at 25,000 new templates. Seed records must be
    oriented 5'-to-3'.
  * Recruitment is bounded by both a library fraction and an absolute template
    ceiling. Checkpoint SPAdes uses one k-mer; the full multi-k SPAdes assembly
    is run once after recruitment stops.
  * ITSME enumerates bounded graph paths between compatible 18S and 28S
    anchors. It stops at the anchor segments rather than choosing branches
    beyond the locus.
  * Dual-anchor contigs are validated for non-overlapping, consistently
    oriented 18S and 28S anchors and automatically oriented 18S-to-28S.
  * ITSx runs only on oriented contigs, avoiding reverse-orientation artifacts.
  * NCBI classification uses the local SSU_eukaryote_rRNA,
    ITS_eukaryote_sequences and LSU_eukaryote_rRNA databases.
  * Complete paths and partial 18S/28S loci are retained regardless of taxon.
    Residual VCF records are unphased; graph paths are candidate haplotypes.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

log() {
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

need_value() {
    [[ $# -ge 2 && -n "${2:-}" ]] || die "Option $1 requires a value."
}

is_positive_int() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

is_nonnegative_int() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

is_number() {
    awk -v value="$1" 'BEGIN {
        exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/)
    }'
}

fraction_in_range() {
    is_number "$1" && awk -v value="$1" 'BEGIN { exit !(value >= 0 && value <= 1) }'
}

greater_than() {
    awk -v left="$1" -v right="$2" 'BEGIN { exit !(left > right) }'
}

less_than() {
    awk -v left="$1" -v right="$2" 'BEGIN { exit !(left < right) }'
}

fraction() {
    awk -v numerator="$1" -v denominator="$2" 'BEGIN {
        if (denominator == 0) print "0.00000000"
        else printf "%.8f", numerator / denominator
    }'
}

format_duration() {
    local total_seconds="$1"
    printf '%02d:%02d:%02d\n' \
        "$((total_seconds / 3600))" \
        "$(((total_seconds % 3600) / 60))" \
        "$((total_seconds % 60))"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -1|--reads1)
            need_value "$@"; READS_1="$2"; shift 2 ;;
        -2|--reads2)
            need_value "$@"; READS_2="$2"; shift 2 ;;
        -s|--single)
            need_value "$@"; READS_SINGLE="$2"; shift 2 ;;
        -o|--output-dir)
            need_value "$@"; OUTDIR="$2"; shift 2 ;;
        --db-dir)
            need_value "$@"; DB_DIR="$2"; shift 2 ;;
        -d|--database)
            need_value "$@"; DATABASES+=("$2"); shift 2 ;;
        --left-db)
            need_value "$@"; LEFT_DATABASES+=("$2"); shift 2 ;;
        --right-db)
            need_value "$@"; RIGHT_DATABASES+=("$2"); shift 2 ;;
        --inward)
            EXTENSION_DIRECTION="inward"; shift ;;
        --outward)
            EXTENSION_DIRECTION="outward"; shift ;;
        --bidirectional)
            EXTENSION_DIRECTION="bidirectional"; shift ;;
        --inward-end-length)
            need_value "$@"; INWARD_END_LENGTH="$2"; shift 2 ;;
        -R|--max-rounds)
            need_value "$@"; MAX_ROUNDS="$2"; shift 2 ;;
        -w|--word-size)
            need_value "$@"; WORD_SIZE="$2"; shift 2 ;;
        --min-word-hits)
            need_value "$@"; MIN_WORD_HITS="$2"; shift 2 ;;
        --bait-min-entropy)
            need_value "$@"; BAIT_MIN_ENTROPY="$2"; shift 2 ;;
        --max-round-growth)
            need_value "$@"; MAX_ROUND_GROWTH="$2"; shift 2 ;;
        --max-accepted-fraction)
            need_value "$@"; MAX_ACCEPTED_FRACTION="$2"; shift 2 ;;
        --max-accepted-templates)
            need_value "$@"; MAX_ACCEPTED_TEMPLATES="$2"; shift 2 ;;
        --outward-min-word-hits)
            need_value "$@"; OUTWARD_MIN_WORD_HITS="$2"; shift 2 ;;
        --max-outward-templates)
            need_value "$@"; MAX_OUTWARD_TEMPLATES="$2"; shift 2 ;;
        --max-outward-growth)
            need_value "$@"; MAX_OUTWARD_GROWTH="$2"; shift 2 ;;
        --min-new-templates)
            need_value "$@"; MIN_NEW_TEMPLATES="$2"; shift 2 ;;
        --min-growth)
            need_value "$@"; MIN_GROWTH="$2"; shift 2 ;;
        --bbtools-memory)
            need_value "$@"; BBTOOLS_MEMORY="$2"; shift 2 ;;
        --bbduk)
            need_value "$@"; BBDUK_CMD="$2"; shift 2 ;;
        --seed-score-min)
            need_value "$@"; SEED_SCORE_MIN="$2"; shift 2 ;;
        --min-read-aligned)
            need_value "$@"; MIN_READ_ALIGNED_FRACTION="$2"; shift 2 ;;
        --min-read-identity)
            need_value "$@"; MIN_READ_IDENTITY="$2"; shift 2 ;;
        --max-insert)
            need_value "$@"; MAX_INSERT="$2"; shift 2 ;;
        --skip-post-qc)
            POST_MAP_QC=false; shift ;;
        --qc-quality)
            need_value "$@"; QC_QUALITY="$2"; shift 2 ;;
        --min-read-length)
            need_value "$@"; MIN_READ_LENGTH="$2"; shift 2 ;;
        --adapters)
            need_value "$@"; ADAPTERS="$2"; shift 2 ;;
        --adapter-clip)
            need_value "$@"; ADAPTER_CLIP="$2"; shift 2 ;;
        --trimmomatic)
            need_value "$@"; TRIMMOMATIC_CMD="$2"; shift 2 ;;
        --trimmomatic-jar)
            need_value "$@"; TRIMMOMATIC_JAR="$2"; shift 2 ;;
        --flash)
            need_value "$@"; FLASH_CMD="$2"; shift 2 ;;
        --flash-min-overlap)
            need_value "$@"; FLASH_MIN_OVERLAP="$2"; shift 2 ;;
        --flash-max-overlap)
            need_value "$@"; FLASH_MAX_OVERLAP="$2"; shift 2 ;;
        --flash-mismatch)
            need_value "$@"; FLASH_MISMATCH_DENSITY="$2"; shift 2 ;;
        -t|--threads)
            need_value "$@"; THREADS="$2"; shift 2 ;;
        -m|--memory)
            need_value "$@"; SPADES_MEMORY_GB="$2"; shift 2 ;;
        -k|--kmers)
            need_value "$@"; KMERS="$2"; shift 2 ;;
        --spades-mode)
            need_value "$@"; SPADES_MODE="$2"; shift 2 ;;
        --anchor-min-aligned)
            need_value "$@"; ANCHOR_MIN_ALIGNED="$2"; shift 2 ;;
        --anchor-min-identity)
            need_value "$@"; ANCHOR_MIN_IDENTITY="$2"; shift 2 ;;
        --layout-overlap)
            need_value "$@"; LAYOUT_OVERLAP_TOLERANCE="$2"; shift 2 ;;
        --skip-graph-paths)
            GRAPH_PATHS=false; shift ;;
        --max-graph-paths)
            need_value "$@"; MAX_GRAPH_PATHS="$2"; shift 2 ;;
        --max-graph-nodes)
            need_value "$@"; MAX_GRAPH_NODES="$2"; shift 2 ;;
        --max-graph-bp)
            need_value "$@"; MAX_GRAPH_BP="$2"; shift 2 ;;
        --min-graph-depth)
            need_value "$@"; MIN_GRAPH_DEPTH="$2"; shift 2 ;;
        --skip-adaptive-graph-stop)
            ADAPTIVE_GRAPH_STOP=false; shift ;;
        --graph-check-k)
            need_value "$@"; GRAPH_CHECK_K="$2"; shift 2 ;;
        --graph-stop-branches)
            need_value "$@"; GRAPH_STOP_BRANCHES="$2"; shift 2 ;;
        --graph-stop-cycles)
            need_value "$@"; GRAPH_STOP_CYCLES="$2"; shift 2 ;;
        --graph-stop-nodes)
            need_value "$@"; GRAPH_STOP_NODES="$2"; shift 2 ;;
        --itsx-taxa)
            need_value "$@"; ITSX_TAXA="$2"; shift 2 ;;
        --skip-itsx)
            RUN_ITSX=false; shift ;;
        --skip-readback)
            READBACK=false; shift ;;
        --skip-variant-calling)
            CALL_VARIANTS=false; shift ;;
        --variant-min-mapq)
            need_value "$@"; VARIANT_MIN_MAPQ="$2"; shift 2 ;;
        --variant-min-baseq)
            need_value "$@"; VARIANT_MIN_BASEQ="$2"; shift 2 ;;
        --variant-max-depth)
            need_value "$@"; VARIANT_MAX_DEPTH="$2"; shift 2 ;;
        --ncbi-db-dir)
            need_value "$@"; NCBI_DB_DIR="$2"; NCBI_DB_DIR_EXPLICIT=true; shift 2 ;;
        --ncbi-max-hits)
            need_value "$@"; NCBI_MAX_HITS="$2"; shift 2 ;;
        --ncbi-report-hits)
            need_value "$@"; NCBI_REPORT_HITS="$2"; shift 2 ;;
        --expected-taxonomy)
            need_value "$@"; EXPECTED_TAXONOMY="$2"; shift 2 ;;
        --ncbi-taxdump-dir)
            need_value "$@"; TAXDUMP_DIR="$2"; TAXDUMP_DIR_EXPLICIT=true; shift 2 ;;
        --skip-ncbi-blast)
            RUN_NCBI_BLAST=false; shift ;;
        --keep-bam)
            KEEP_BAM=true; shift ;;
        --keep-intermediates)
            KEEP_INTERMEDIATES=true; shift ;;
        -h|--help)
            usage; exit 0 ;;
        --version)
            printf '%s\n' "$VERSION"; exit 0 ;;
        --)
            shift; break ;;
        *)
            die "Unknown option: $1 (use --help)" ;;
    esac
done

[[ -n "$READS_1" ]] || die "Missing -1/--reads1."
[[ -n "$READS_2" ]] || die "Missing -2/--reads2."
[[ -n "$OUTDIR" ]] || die "Missing -o/--output-dir."

if [[ -n "$DB_DIR" ]]; then
    [[ -d "$DB_DIR" ]] || die "rRNA database directory not found: $DB_DIR"
    DB_DIR=$(cd "$DB_DIR" && pwd)
    DB_LEFT="$DB_DIR/silva-euk-18s-id95.fasta"
    DB_RIGHT="$DB_DIR/silva-euk-28s-id98.fasta"
    [[ -r "$DB_LEFT" && -s "$DB_LEFT" ]] || \
        die "Required left/18S database is missing or empty: $DB_LEFT"
    [[ -r "$DB_RIGHT" && -s "$DB_RIGHT" ]] || \
        die "Required right/28S database is missing or empty: $DB_RIGHT"
    LEFT_DATABASES+=("$DB_LEFT")
    RIGHT_DATABASES+=("$DB_RIGHT")

    # Accept either a flat database root or conventional NCBI/taxdump
    # subdirectories. Explicit command-line paths retain precedence.
    if [[ "$NCBI_DB_DIR_EXPLICIT" != true ]]; then
        if [[ -d "$DB_DIR/NCBI_rRNA_BLAST" ]]; then
            NCBI_DB_DIR="$DB_DIR/NCBI_rRNA_BLAST"
        else
            NCBI_DB_DIR="$DB_DIR"
        fi
    fi
    if [[ "$TAXDUMP_DIR_EXPLICIT" != true ]]; then
        if [[ -r "$DB_DIR/nodes.dmp" && -r "$DB_DIR/names.dmp" ]]; then
            TAXDUMP_DIR="$DB_DIR"
        elif [[ -d "$DB_DIR/taxdump" ]]; then
            TAXDUMP_DIR="$DB_DIR/taxdump"
        else
            TAXDUMP_DIR="$DB_DIR"
        fi
    fi
fi

(( ${#DATABASES[@]} + ${#LEFT_DATABASES[@]} + ${#RIGHT_DATABASES[@]} > 0 )) || \
    die "Supply --db-dir or at least one explicit seed database."

if [[ "$EXTENSION_DIRECTION" == "auto" ]]; then
    if (( ${#LEFT_DATABASES[@]} > 0 && ${#RIGHT_DATABASES[@]} > 0 )); then
        INWARD=true
        OUTWARD=true
        EXTENSION_DIRECTION="bidirectional"
    else
        INWARD=false
        OUTWARD=true
        EXTENSION_DIRECTION="outward"
    fi
elif [[ "$EXTENSION_DIRECTION" == "inward" ]]; then
    INWARD=true
    OUTWARD=false
elif [[ "$EXTENSION_DIRECTION" == "bidirectional" ]]; then
    INWARD=true
    OUTWARD=true
else
    INWARD=false
    OUTWARD=true
fi

[[ -r "$READS_1" && -s "$READS_1" ]] || die "Cannot read forward reads: $READS_1"
[[ -r "$READS_2" && -s "$READS_2" ]] || die "Cannot read reverse reads: $READS_2"
if [[ -n "$READS_SINGLE" ]]; then
    [[ -r "$READS_SINGLE" && -s "$READS_SINGLE" ]] || die "Cannot read single reads: $READS_SINGLE"
fi
for database in "${DATABASES[@]}" "${LEFT_DATABASES[@]}" "${RIGHT_DATABASES[@]}"; do
    [[ -z "$database" ]] && continue
    [[ -r "$database" && -s "$database" ]] || die "Cannot read seed database: $database"
done

if [[ "$INWARD" == true || "$EXTENSION_DIRECTION" == "outward" ]]; then
    if [[ "$EXTENSION_DIRECTION" != "outward" || \
          ${#LEFT_DATABASES[@]} -gt 0 || ${#RIGHT_DATABASES[@]} -gt 0 ]]; then
        (( ${#LEFT_DATABASES[@]} > 0 )) || die "Directional extension requires at least one --left-db."
        (( ${#RIGHT_DATABASES[@]} > 0 )) || die "Directional extension requires at least one --right-db."
    fi
fi
is_positive_int "$THREADS" || die "--threads must be a positive integer."
is_positive_int "$SPADES_MEMORY_GB" || die "--memory must be a positive integer."
is_nonnegative_int "$MAX_ROUNDS" || die "--max-rounds must be a nonnegative integer."
is_positive_int "$WORD_SIZE" || die "--word-size must be a positive integer."
(( WORD_SIZE >= 15 && WORD_SIZE <= 31 )) || die "--word-size must be between 15 and 31 for exact BBDuk words."
is_positive_int "$MIN_WORD_HITS" || die "--min-word-hits must be a positive integer."
is_positive_int "$OUTWARD_MIN_WORD_HITS" || die "--outward-min-word-hits must be a positive integer."
is_positive_int "$INWARD_END_LENGTH" || die "--inward-end-length must be a positive integer."
is_nonnegative_int "$MIN_NEW_TEMPLATES" || die "--min-new-templates must be nonnegative."
is_positive_int "$MAX_INSERT" || die "--max-insert must be a positive integer."
is_positive_int "$MIN_READ_LENGTH" || die "--min-read-length must be a positive integer."
is_nonnegative_int "$QC_QUALITY" || die "--qc-quality must be nonnegative."
is_positive_int "$FLASH_MIN_OVERLAP" || die "--flash-min-overlap must be positive."
is_positive_int "$FLASH_MAX_OVERLAP" || die "--flash-max-overlap must be positive."
(( FLASH_MAX_OVERLAP >= FLASH_MIN_OVERLAP )) || die "FLASH maximum overlap must be >= minimum overlap."
is_positive_int "$ANCHOR_MIN_ALIGNED" || die "--anchor-min-aligned must be positive."
is_nonnegative_int "$LAYOUT_OVERLAP_TOLERANCE" || die "--layout-overlap must be nonnegative."
is_positive_int "$NCBI_MAX_HITS" || die "--ncbi-max-hits must be positive."
is_positive_int "$NCBI_REPORT_HITS" || die "--ncbi-report-hits must be positive."
is_nonnegative_int "$MAX_ACCEPTED_TEMPLATES" || \
    die "--max-accepted-templates must be a nonnegative integer."
is_nonnegative_int "$MAX_OUTWARD_TEMPLATES" || \
    die "--max-outward-templates must be a nonnegative integer."
is_positive_int "$MAX_GRAPH_PATHS" || die "--max-graph-paths must be positive."
is_positive_int "$MAX_GRAPH_NODES" || die "--max-graph-nodes must be positive."
is_positive_int "$MAX_GRAPH_BP" || die "--max-graph-bp must be positive."
is_positive_int "$GRAPH_CHECK_K" || die "--graph-check-k must be positive."
(( GRAPH_CHECK_K >= 21 && GRAPH_CHECK_K <= 127 && GRAPH_CHECK_K % 2 == 1 )) || \
    die "--graph-check-k must be an odd integer from 21 through 127."
is_nonnegative_int "$GRAPH_STOP_BRANCHES" || die "--graph-stop-branches must be nonnegative."
is_nonnegative_int "$GRAPH_STOP_CYCLES" || die "--graph-stop-cycles must be nonnegative."
is_positive_int "$GRAPH_STOP_NODES" || die "--graph-stop-nodes must be positive."
is_number "$MIN_GRAPH_DEPTH" || die "--min-graph-depth must be nonnegative."
greater_than 0 "$MIN_GRAPH_DEPTH" && die "--min-graph-depth must be nonnegative."
is_nonnegative_int "$VARIANT_MIN_MAPQ" || die "--variant-min-mapq must be nonnegative."
is_nonnegative_int "$VARIANT_MIN_BASEQ" || die "--variant-min-baseq must be nonnegative."
is_positive_int "$VARIANT_MAX_DEPTH" || die "--variant-max-depth must be positive."
(( NCBI_REPORT_HITS <= NCBI_MAX_HITS )) || \
    die "--ncbi-report-hits cannot exceed --ncbi-max-hits."

for value in "$BAIT_MIN_ENTROPY" "$MAX_ACCEPTED_FRACTION" \
             "$MIN_READ_ALIGNED_FRACTION" "$MIN_READ_IDENTITY" \
             "$ANCHOR_MIN_IDENTITY" "$FLASH_MISMATCH_DENSITY"; do
    fraction_in_range "$value" || die "Fractional parameters must be numbers from 0 through 1."
done
if [[ "$RUN_NCBI_BLAST" == true ]]; then
    command -v blastdbcmd >/dev/null 2>&1 || \
        die "NCBI classification requested but blastdbcmd was not found."
    [[ -n "$NCBI_DB_DIR" ]] || die \
        "NCBI classification is enabled. Supply --ncbi-db-dir or set BLASTDB."
    [[ -d "$NCBI_DB_DIR" ]] || die "NCBI BLAST database directory not found: $NCBI_DB_DIR"
    [[ -n "$TAXDUMP_DIR" ]] || die \
        "NCBI classification requires --ncbi-taxdump-dir or NCBI_TAXDUMP_DIR."
    [[ -r "$TAXDUMP_DIR/nodes.dmp" && -r "$TAXDUMP_DIR/names.dmp" ]] || die \
        "Taxdump directory must contain readable nodes.dmp and names.dmp: $TAXDUMP_DIR"
    for database in ITS_eukaryote_sequences SSU_eukaryote_rRNA LSU_eukaryote_rRNA; do
        blastdbcmd -db "$NCBI_DB_DIR/$database" -info >/dev/null 2>&1 || \
            die "Required NCBI BLAST database is unavailable: $NCBI_DB_DIR/$database"
    done
fi
is_number "$MAX_ROUND_GROWTH" || die "--max-round-growth must be nonnegative."
is_number "$MAX_OUTWARD_GROWTH" || die "--max-outward-growth must be nonnegative."
is_number "$MIN_GROWTH" || die "--min-growth must be nonnegative."
greater_than 0 "$MAX_ROUND_GROWTH" && die "--max-round-growth must be nonnegative."
greater_than 0 "$MAX_OUTWARD_GROWTH" && die "--max-outward-growth must be nonnegative."
greater_than 0 "$MIN_GROWTH" && die "--min-growth must be nonnegative."
[[ "$BBTOOLS_MEMORY" =~ ^[1-9][0-9]*[mMgG]$ ]] || \
    die "--bbtools-memory must look like 8000m or 16g."
[[ "$ADAPTER_CLIP" =~ ^[0-9]+:[0-9]+:[0-9]+(:[0-9]+:[A-Za-z]+)?$ ]] || \
    die "--adapter-clip is not a recognized ILLUMINACLIP parameter string."
[[ "$SPADES_MODE" == "standard" || "$SPADES_MODE" == "meta" ]] || \
    die "--spades-mode must be standard or meta."
[[ "$KMERS" == "auto" || "$KMERS" =~ ^[0-9]+(,[0-9]+)*$ ]] || \
    die "--kmers must be auto or a comma-separated integer list."

for program in bowtie2 bowtie2-build samtools spades.py blastn python3 awk sort gzip; do
    command -v "$program" >/dev/null 2>&1 || die "Required program not found in PATH: $program"
done
if [[ "$CALL_VARIANTS" == true ]]; then
    command -v bcftools >/dev/null 2>&1 || \
        die "Variant calling is enabled but bcftools was not found; install bcftools or use --skip-variant-calling."
fi
if [[ "$CALL_VARIANTS" == true && "$READBACK" != true ]]; then
    die "Variant calling requires read-back mapping; remove --skip-readback or add --skip-variant-calling."
fi
if (( MAX_ROUNDS > 0 )); then
    command -v "$BBDUK_CMD" >/dev/null 2>&1 || \
        die "Extension requested but BBDuk was not found: $BBDUK_CMD"
fi
if [[ "$RUN_ITSX" == true ]] && ! command -v ITSx >/dev/null 2>&1; then
    log "WARNING: ITSx was not found; structural ITS annotation will be skipped."
    RUN_ITSX=false
fi
if [[ "$RUN_NCBI_BLAST" == true && "$RUN_ITSX" != true ]]; then
    die "NCBI locus classification requires ITSx region extraction; install ITSx or use both --skip-itsx and --skip-ncbi-blast."
fi
samtools view --help 2>&1 | grep -Eq -- '(^|[[:space:]])-N([[:space:],]|$)' || \
    die "This script requires a samtools view command that supports -N."

find_default_adapters() {
    local candidate=""
    local roots=()
    [[ -n "${CONDA_PREFIX:-}" ]] && roots+=("$CONDA_PREFIX/share")
    if [[ -n "$TRIMMOMATIC_JAR" ]]; then
        roots+=("$(dirname "$TRIMMOMATIC_JAR")")
        roots+=("$(dirname "$(dirname "$TRIMMOMATIC_JAR")")")
    fi
    roots+=("/usr/share/trimmomatic" "/usr/local/share/trimmomatic")
    local root
    for root in "${roots[@]}"; do
        [[ -d "$root" ]] || continue
        candidate=$(find "$root" -maxdepth 4 -type f -name 'TruSeq3-PE-2.fa' -print -quit 2>/dev/null)
        if [[ -n "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

if [[ "$POST_MAP_QC" == true ]]; then
    if ! command -v "$FLASH_CMD" >/dev/null 2>&1; then
        if [[ "$FLASH_CMD" == "flash" ]] && command -v flash2 >/dev/null 2>&1; then
            FLASH_CMD="flash2"
        else
            die "FLASH executable not found: $FLASH_CMD"
        fi
    fi
    if [[ -n "$TRIMMOMATIC_JAR" ]]; then
        [[ -r "$TRIMMOMATIC_JAR" && -s "$TRIMMOMATIC_JAR" ]] || die "Cannot read Trimmomatic JAR."
        command -v java >/dev/null 2>&1 || die "java is required with --trimmomatic-jar."
    else
        command -v "$TRIMMOMATIC_CMD" >/dev/null 2>&1 || die "Trimmomatic executable not found."
    fi
    if [[ -z "$ADAPTERS" ]]; then
        ADAPTERS=$(find_default_adapters) || \
            die "Could not find TruSeq3-PE-2.fa; provide --adapters FILE."
    fi
    [[ -r "$ADAPTERS" && -s "$ADAPTERS" ]] || die "Cannot read adapter FASTA: $ADAPTERS"
fi

if [[ -e "$OUTDIR" ]]; then
    [[ -d "$OUTDIR" ]] || die "Output path exists and is not a directory: $OUTDIR"
    [[ -z "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]] || \
        die "Output directory is not empty: $OUTDIR"
else
    mkdir -p "$OUTDIR"
fi
mkdir -p "$OUTDIR/seed" "$OUTDIR/recruitment" "$OUTDIR/assembly" "$OUTDIR/final"

COMBINED_DB="$OUTDIR/seed/combined_rrna_seeds.fasta"
LEFT_DB="$OUTDIR/seed/left_18S_seeds.fasta"
RIGHT_DB="$OUTDIR/seed/right_28S_seeds.fasta"
SEED_INDEX="$OUTDIR/seed/combined_rrna_seeds"
LEFT_INDEX="$OUTDIR/seed/left_18S_seeds"
RIGHT_INDEX="$OUTDIR/seed/right_28S_seeds"
SEED_BAM="$OUTDIR/seed/all_reads_to_seeds.bam"
ACCEPTED_NAMES="$OUTDIR/accepted_read_names.txt"
METRICS="$OUTDIR/recruitment.tsv"
GRAPH_METRICS="$OUTDIR/recruitment_graph.tsv"

read_fasta() {
    case "$1" in
        *.gz|*.bgz) gzip -cd -- "$1" ;;
        *) awk '{ sub(/\r$/, ""); print }' "$1" ;;
    esac
}

stream_fastq() {
    local input="${1:-}"
    [[ -n "$input" ]] || return 0
    case "$input" in
        *.gz|*.bgz) gzip -cd -- "$input" ;;
        *) awk '{ sub(/\r$/, ""); print }' "$input" ;;
    esac
}

fastq_count() {
    local input="${1:-}"
    if [[ -z "$input" ]]; then
        printf '0\n'
        return 0
    fi
    [[ -r "$input" ]] || die "Cannot read FASTQ while counting: $input"
    stream_fastq "$input" | awk -v input_name="$input" '
        END {
            if (NR % 4 != 0) {
                printf "Malformed FASTQ: %s\n", input_name > "/dev/stderr"
                exit 2
            }
            print NR / 4
        }
    '
}

extract_fastq_names() {
    awk '
        NR % 4 == 1 {
            name = substr($0, 2)
            sub(/[[:space:]].*$/, "", name)
            sub(/\/[12]$/, "", name)
            print name
        }
    '
}

normalize_database() {
    local database="$1"
    local role="$2"
    local number="$3"
    local output="$4"
    local prefix
    printf -v prefix '%s|db%04d' "$role" "$number"
    read_fasta "$database" | awk -v prefix="$prefix" '
        BEGIN { record = 0 }
        /^>/ {
            record++
            header = substr($0, 2)
            gsub(/[[:space:]]+/, "_", header)
            printf ">%s|seq%06d|%s\n", prefix, record, header
            next
        }
        { sequence = toupper($0); gsub(/[^ACGTN]/, "N", sequence); print sequence }
    ' > "$output"
}

combine_databases() {
    : > "$COMBINED_DB"
    : > "$LEFT_DB"
    : > "$RIGHT_DB"
    local database number=0 normalized

    for database in "${DATABASES[@]}"; do
        number=$((number + 1))
        normalized="$OUTDIR/seed/general_${number}.fasta"
        normalize_database "$database" general "$number" "$normalized"
        awk '{ print }' "$normalized" >> "$COMBINED_DB"
        rm -f -- "$normalized"
    done
    number=0
    for database in "${LEFT_DATABASES[@]}"; do
        number=$((number + 1))
        normalized="$OUTDIR/seed/left_${number}.fasta"
        normalize_database "$database" left "$number" "$normalized"
        awk '{ print }' "$normalized" >> "$COMBINED_DB"
        awk '{ print }' "$normalized" >> "$LEFT_DB"
        rm -f -- "$normalized"
    done
    number=0
    for database in "${RIGHT_DATABASES[@]}"; do
        number=$((number + 1))
        normalized="$OUTDIR/seed/right_${number}.fasta"
        normalize_database "$database" right "$number" "$normalized"
        awk '{ print }' "$normalized" >> "$COMBINED_DB"
        awk '{ print }' "$normalized" >> "$RIGHT_DB"
        rm -f -- "$normalized"
    done

    local count
    count=$(awk '/^>/ { n++ } END { print n + 0 }' "$COMBINED_DB")
    (( count > 0 )) || die "No FASTA records were found in the seed databases."
    log "Combined seed databases into $count uniquely named sequences."
}

make_directional_seed_baits() {
    local output="$1"
    local direction="$2"
    : > "$output"
    if [[ "$direction" == "inward" || "$direction" == "bidirectional" ]]; then
    awk -v end_len="$INWARD_END_LENGTH" '
        function emit(    start) {
            if (name == "" || sequence == "") return
            start = length(sequence) > end_len ? length(sequence) - end_len + 1 : 1
            print ">" name "|inward_3prime"
            print substr(sequence, start)
        }
        /^>/ { emit(); name = substr($0, 2); sequence = ""; next }
        { sequence = sequence $0 }
        END { emit() }
    ' "$LEFT_DB" >> "$output"
    awk -v end_len="$INWARD_END_LENGTH" '
        function emit() {
            if (name == "" || sequence == "") return
            print ">" name "|inward_5prime"
            print substr(sequence, 1, end_len)
        }
        /^>/ { emit(); name = substr($0, 2); sequence = ""; next }
        { sequence = sequence $0 }
        END { emit() }
    ' "$RIGHT_DB" >> "$output"
    fi
    if [[ "$direction" == "outward" || "$direction" == "bidirectional" ]]; then
    awk -v end_len="$INWARD_END_LENGTH" '
        function emit() {
            if (name == "" || sequence == "") return
            print ">" name "|outward_5prime"
            print substr(sequence, 1, end_len)
        }
        /^>/ { emit(); name = substr($0, 2); sequence = ""; next }
        { sequence = sequence $0 }
        END { emit() }
    ' "$LEFT_DB" >> "$output"
    awk -v end_len="$INWARD_END_LENGTH" '
        function emit(    start) {
            if (name == "" || sequence == "") return
            start = length(sequence) > end_len ? length(sequence) - end_len + 1 : 1
            print ">" name "|outward_3prime"
            print substr(sequence, start)
        }
        /^>/ { emit(); name = substr($0, 2); sequence = ""; next }
        { sequence = sequence $0 }
        END { emit() }
    ' "$RIGHT_DB" >> "$output"
    fi
    [[ -s "$output" ]] || die "Failed to generate directional terminal seed baits."
}

map_raw_reads_to_seeds() {
    local command=(
        bowtie2 --very-sensitive-local --score-min "$SEED_SCORE_MIN"
        -X "$MAX_INSERT" -p "$THREADS" -x "$SEED_INDEX"
        -1 "$READS_1" -2 "$READS_2"
    )
    [[ -n "$READS_SINGLE" ]] && command+=(-U "$READS_SINGLE")
    "${command[@]}" 2> "$OUTDIR/seed/bowtie2.log" | \
        samtools view -@ "$THREADS" -b -o "$SEED_BAM" -
}

mapped_primary_names() {
    local bam="$1"
    local output="$2"
    samtools view -@ "$THREADS" -F 2308 "$bam" | \
        awk -v min_fraction="$MIN_READ_ALIGNED_FRACTION" -v min_identity="$MIN_READ_IDENTITY" '
            function aligned_bases(cigar,    rest, token, total) {
                rest = cigar
                total = 0
                while (match(rest, /[0-9]+[MI=X]/)) {
                    token = substr(rest, RSTART, RLENGTH)
                    total += token + 0
                    rest = substr(rest, RSTART + RLENGTH)
                }
                return total
            }
            {
                query_length = length($10)
                aligned = aligned_bases($6)
                nm = 0
                for (i = 12; i <= NF; i++) {
                    if ($i ~ /^NM:i:/) {
                        split($i, value, ":")
                        nm = value[3] + 0
                        break
                    }
                }
                aligned_fraction = query_length > 0 ? aligned / query_length : 0
                identity = aligned > 0 ? (aligned - nm) / aligned : 0
                if (aligned_fraction >= min_fraction && identity >= min_identity)
                    print $1
            }
        ' | LC_ALL=C sort -u > "$output"
}

extract_templates_from_bam() {
    local names="$1"
    local output_r1="$2"
    local output_r2="$3"
    local output_single="$4"
    local prefix="$5"
    local selected_bam="${prefix}.selected.bam"
    local other_fastq="${prefix}.other.fastq.gz"
    local singleton_fastq="${prefix}.singleton.fastq.gz"
    local log_file="${prefix}.fastq.log"

    samtools view -@ "$THREADS" -b -F 2304 -N "$names" \
        -o "$selected_bam" "$SEED_BAM"
    samtools collate -@ "$THREADS" -u -O "$selected_bam" 2> "$log_file" | \
        samtools fastq -@ "$THREADS" -c 1 -n \
            -1 "$output_r1" -2 "$output_r2" \
            -0 "$other_fastq" -s "$singleton_fastq" \
            - 2>> "$log_file" >/dev/null

    if command -v pigz >/dev/null 2>&1; then
        gzip -cd -- "$other_fastq" "$singleton_fastq" | \
            pigz -1 -p "$THREADS" > "$output_single"
    else
        gzip -cd -- "$other_fastq" "$singleton_fastq" | gzip -1 > "$output_single"
    fi
    rm -f -- "$selected_bam" "$other_fastq" "$singleton_fastq"

    local pairs singles expected
    pairs=$(fastq_count "$output_r1")
    singles=$(fastq_count "$output_single")
    expected=$(awk 'END { print NR + 0 }' "$names")
    (( pairs + singles == expected )) || \
        die "Extracted $((pairs + singles)) templates but expected $expected; read names may not be unique. Inspect $log_file"
}

run_trimmomatic() {
    if [[ -n "$TRIMMOMATIC_JAR" ]]; then
        java -jar "$TRIMMOMATIC_JAR" "$@"
    else
        "$TRIMMOMATIC_CMD" "$@"
    fi
}

qc_recruited_batch() {
    local input_r1="$1"
    local input_r2="$2"
    local input_single="$3"
    local qc_dir="$4"
    local output_r1="$5"
    local output_r2="$6"
    local output_single="$7"
    local qc_log="$qc_dir/qc_run.log"
    local trim_r1p="$qc_dir/preflash_R1P.fastq.gz"
    local trim_r1u="$qc_dir/preflash_R1U.fastq.gz"
    local trim_r2p="$qc_dir/preflash_R2P.fastq.gz"
    local trim_r2u="$qc_dir/preflash_R2U.fastq.gz"
    local trim_single="$qc_dir/trimmed_single.fastq.gz"
    local flash_prefix="flashed"
    local flash_merged="$qc_dir/${flash_prefix}.extendedFrags.fastq.gz"
    local flash_r1="$qc_dir/${flash_prefix}.notCombined_1.fastq.gz"
    local flash_r2="$qc_dir/${flash_prefix}.notCombined_2.fastq.gz"
    local input_pairs input_reverse input_singles

    mkdir -p "$qc_dir"
    input_pairs=$(fastq_count "$input_r1")
    input_reverse=$(fastq_count "$input_r2")
    input_singles=$(fastq_count "$input_single")
    (( input_pairs == input_reverse )) || die "Recruited R1/R2 counts differ before QC."

    printf '%s\n' \
        "Post-mapping QC settings" \
        "Input paired templates: $input_pairs" \
        "Input single reads: $input_singles" \
        "Adapter FASTA: $ADAPTERS" \
        "ILLUMINACLIP: $ADAPTER_CLIP" \
        "SLIDINGWINDOW: 4:$QC_QUALITY" \
        "MINLEN: $MIN_READ_LENGTH" \
        "FLASH minimum overlap: $FLASH_MIN_OVERLAP" \
        "FLASH maximum overlap: $FLASH_MAX_OVERLAP" \
        "FLASH maximum mismatch density: $FLASH_MISMATCH_DENSITY" \
        > "$qc_log"

    if (( input_pairs > 0 )); then
        run_trimmomatic PE -threads "$THREADS" -phred33 \
            "$input_r1" "$input_r2" \
            "$trim_r1p" "$trim_r1u" "$trim_r2p" "$trim_r2u" \
            "ILLUMINACLIP:${ADAPTERS}:${ADAPTER_CLIP}" \
            "SLIDINGWINDOW:4:${QC_QUALITY}" "MINLEN:${MIN_READ_LENGTH}" \
            >> "$qc_log" 2>&1 || die "Post-map paired QC failed; inspect $qc_log"
    else
        gzip -c </dev/null > "$trim_r1p"
        gzip -c </dev/null > "$trim_r1u"
        gzip -c </dev/null > "$trim_r2p"
        gzip -c </dev/null > "$trim_r2u"
    fi

    local surviving_pairs surviving_reverse
    surviving_pairs=$(fastq_count "$trim_r1p")
    surviving_reverse=$(fastq_count "$trim_r2p")
    (( surviving_pairs == surviving_reverse )) || die "Trimmomatic paired outputs differ."

    if (( surviving_pairs > 0 )); then
        "$FLASH_CMD" --threads "$THREADS" \
            --min-overlap "$FLASH_MIN_OVERLAP" \
            --max-overlap "$FLASH_MAX_OVERLAP" \
            --max-mismatch-density "$FLASH_MISMATCH_DENSITY" \
            --output-prefix "$flash_prefix" --output-directory "$qc_dir" --compress \
            "$trim_r1p" "$trim_r2p" >> "$qc_log" 2>&1 || \
            die "FLASH failed; inspect $qc_log"
    else
        gzip -c </dev/null > "$flash_merged"
        gzip -c </dev/null > "$flash_r1"
        gzip -c </dev/null > "$flash_r2"
    fi

    if (( input_singles > 0 )); then
        run_trimmomatic SE -threads "$THREADS" -phred33 \
            "$input_single" "$trim_single" \
            "ILLUMINACLIP:${ADAPTERS}:${ADAPTER_CLIP}" \
            "SLIDINGWINDOW:4:${QC_QUALITY}" "MINLEN:${MIN_READ_LENGTH}" \
            >> "$qc_log" 2>&1 || die "Post-map single-read QC failed; inspect $qc_log"
    else
        gzip -c </dev/null > "$trim_single"
    fi

    mv -- "$flash_r1" "$output_r1"
    mv -- "$flash_r2" "$output_r2"
    if command -v pigz >/dev/null 2>&1; then
        gzip -cd -- "$flash_merged" "$trim_r1u" "$trim_r2u" "$trim_single" | \
            pigz -1 -p "$THREADS" > "$output_single"
    else
        gzip -cd -- "$flash_merged" "$trim_r1u" "$trim_r2u" "$trim_single" | \
            gzip -1 > "$output_single"
    fi

    local output_pairs output_singles merged
    output_pairs=$(fastq_count "$output_r1")
    output_singles=$(fastq_count "$output_single")
    merged=$(fastq_count "$flash_merged")
    printf 'metric\tcount\n' > "$qc_dir/qc_summary.tsv"
    printf '%s\t%d\n' \
        recruited_input_pairs "$input_pairs" \
        recruited_input_singles "$input_singles" \
        pairs_surviving_trimmomatic "$surviving_pairs" \
        flash_merged_fragments "$merged" \
        output_paired_templates "$output_pairs" \
        output_single_reads "$output_singles" \
        >> "$qc_dir/qc_summary.tsv"

    if [[ "$KEEP_INTERMEDIATES" != true ]]; then
        rm -f -- "$trim_r1p" "$trim_r1u" "$trim_r2p" "$trim_r2u" \
            "$trim_single" "$flash_merged"
    fi
}

prepare_recruited_batch() {
    local raw_r1="$1"
    local raw_r2="$2"
    local raw_single="$3"
    local batch_dir="$4"
    local clean_r1="$5"
    local clean_r2="$6"
    local clean_single="$7"
    if [[ "$POST_MAP_QC" == true ]]; then
        qc_recruited_batch "$raw_r1" "$raw_r2" "$raw_single" \
            "$batch_dir/postmap_qc" "$clean_r1" "$clean_r2" "$clean_single"
    else
        cp -- "$raw_r1" "$clean_r1"
        cp -- "$raw_r2" "$clean_r2"
        cp -- "$raw_single" "$clean_single"
    fi
}

make_frontier_baits() {
    local output="$1"
    local log_file="$2"
    shift 2
    local raw_fasta="${output%.fasta}.unfiltered.fasta"
    : > "$raw_fasta"
    local input tag=0
    for input in "$@"; do
        [[ -r "$input" ]] || continue
        tag=$((tag + 1))
        stream_fastq "$input" | awk -v tag="$tag" '
            NR % 4 == 1 {
                record++
                name = substr($0, 2)
                sub(/[[:space:]].*$/, "", name)
                printf ">bait%02d|%09d|%s\n", tag, record, name
            }
            NR % 4 == 2 { print toupper($0) }
        ' >> "$raw_fasta"
    done
    "$BBDUK_CMD" "-Xmx${BBTOOLS_MEMORY}" \
        in="$raw_fasta" out="$output" overwrite=t \
        minlen="$WORD_SIZE" entropy="$BAIT_MIN_ENTROPY" \
        entropywindow=50 entropyk=5 \
        2> "$log_file"
    rm -f -- "$raw_fasta"
    [[ -s "$output" ]] || die "No frontier sequences survived bait preparation; inspect $log_file"
}

scan_raw_with_baits() {
    local bait_fasta="$1"
    local output_names="$2"
    local log_prefix="$3"
    local word_size="${4:-$WORD_SIZE}"
    local min_word_hits="${5:-$MIN_WORD_HITS}"
    local paired_names="${output_names%.txt}.paired.txt"
    local single_names="${output_names%.txt}.single.txt"

    "$BBDUK_CMD" "-Xmx${BBTOOLS_MEMORY}" \
        in1="$READS_1" in2="$READS_2" outm=stdout.fq \
        ref="$bait_fasta" k="$word_size" hdist=0 \
        minkmerhits="$min_word_hits" mm=f rcomp=t \
        t="$THREADS" overwrite=t ordered=f \
        2> "${log_prefix}.paired.log" | extract_fastq_names > "$paired_names"

    : > "$single_names"
    if [[ -n "$READS_SINGLE" ]]; then
        "$BBDUK_CMD" "-Xmx${BBTOOLS_MEMORY}" \
            in="$READS_SINGLE" outm=stdout.fq \
            ref="$bait_fasta" k="$word_size" hdist=0 \
            minkmerhits="$min_word_hits" mm=f rcomp=t \
            t="$THREADS" overwrite=t ordered=f \
            2> "${log_prefix}.single.log" | extract_fastq_names > "$single_names"
    fi

    LC_ALL=C sort -u "$paired_names" "$single_names" > "$output_names"
    rm -f -- "$paired_names" "$single_names"
}

combine_fastqs() {
    local output="$1"
    shift
    if command -v pigz >/dev/null 2>&1; then
        gzip -cd -- "$@" | pigz -1 -p "$THREADS" > "$output"
    else
        gzip -cd -- "$@" | gzip -1 > "$output"
    fi
}

assemble_final_reads() {
    local input_r1="$1"
    local input_r2="$2"
    local input_single="$3"
    local assembly_dir="$OUTDIR/assembly/spades"
    local command=(
        spades.py --only-assembler -o "$assembly_dir"
        -t "$THREADS" -m "$SPADES_MEMORY_GB"
    )
    local pairs reverse singles
    pairs=$(fastq_count "$input_r1")
    reverse=$(fastq_count "$input_r2")
    singles=$(fastq_count "$input_single")
    (( pairs == reverse )) || die "Final paired FASTQs have different counts."
    (( pairs + singles > 0 )) || die "No recruited reads remain for assembly."
    (( pairs > 0 )) && command+=(-1 "$input_r1" -2 "$input_r2")
    (( singles > 0 )) && command+=(-s "$input_single")
    if [[ "$SPADES_MODE" == "meta" ]]; then
        if (( pairs > 0 )); then
            command+=(--meta)
        else
            log "WARNING: no unmerged paired reads remained; falling back from metaSPAdes to standard SPAdes."
        fi
    fi
    [[ "$KMERS" != "auto" ]] && command+=(-k "$KMERS")
    "${command[@]}" > "$OUTDIR/assembly/spades.stdout.log" \
        2> "$OUTDIR/assembly/spades.stderr.log" || \
        die "SPAdes failed; inspect $OUTDIR/assembly/spades.stderr.log"
    [[ -s "$assembly_dir/contigs.fasta" ]] || die "SPAdes produced no contigs."
}

contig_anchor_hits() {
    local assembly="$1"
    local index="$2"
    local hits="$3"
    local log_file="$4"
    bowtie2 --very-sensitive-local --score-min "$SEED_SCORE_MIN" \
        -p "$THREADS" -x "$index" -f -U "$assembly" \
        2> "$log_file" | samtools view -@ "$THREADS" -F 2308 - | \
        awk -v min_aligned="$ANCHOR_MIN_ALIGNED" -v min_identity="$ANCHOR_MIN_IDENTITY" '
            function aligned_bases(cigar,    rest, token, total) {
                rest = cigar
                total = 0
                while (match(rest, /[0-9]+[MI=X]/)) {
                    token = substr(rest, RSTART, RLENGTH)
                    total += token + 0
                    rest = substr(rest, RSTART + RLENGTH)
                }
                return total
            }
            {
                aligned = aligned_bases($6)
                nm = 0
                for (i = 12; i <= NF; i++) {
                    if ($i ~ /^NM:i:/) {
                        split($i, value, ":")
                        nm = value[3] + 0
                        break
                    }
                }
                identity = aligned > 0 ? (aligned - nm) / aligned : 0
                if (aligned >= min_aligned && identity >= min_identity)
                    printf "%s\t%s\t%d\t%d\t%.6f\t%s\n", $1, $3, aligned, nm, identity, $5
            }
        ' > "$hits"
}

run_graph_checkpoint() {
    local round="$1"
    local round_dir="$2"
    local checkpoint_tmp checkpoint_r1 checkpoint_r2 checkpoint_single
    local checkpoint_assembly graph segment_fasta left_hits right_hits stats
    checkpoint_tmp=$(mktemp -d "$round_dir/graph_checkpoint.tmp.XXXXXX")
    checkpoint_r1="$checkpoint_tmp/accepted_R1.fastq.gz"
    checkpoint_r2="$checkpoint_tmp/accepted_R2.fastq.gz"
    checkpoint_single="$checkpoint_tmp/accepted_single.fastq.gz"
    checkpoint_assembly="$checkpoint_tmp/spades"
    graph="$checkpoint_assembly/assembly_graph_with_scaffolds.gfa"
    segment_fasta="$checkpoint_tmp/segments.fasta"
    left_hits="$checkpoint_tmp/left_hits.tsv"
    right_hits="$checkpoint_tmp/right_hits.tsv"
    stats="$checkpoint_tmp/stats.tsv"

    combine_fastqs "$checkpoint_r1" "${ASSEMBLY_R1_FILES[@]}"
    combine_fastqs "$checkpoint_r2" "${ASSEMBLY_R2_FILES[@]}"
    combine_fastqs "$checkpoint_single" "${ASSEMBLY_SINGLE_FILES[@]}"

    local pairs singles
    pairs=$(fastq_count "$checkpoint_r1")
    singles=$(fastq_count "$checkpoint_single")
    local command=(
        spades.py --only-assembler -o "$checkpoint_assembly"
        -t "$THREADS" -m "$SPADES_MEMORY_GB" -k "$GRAPH_CHECK_K"
    )
    (( pairs > 0 )) && command+=(-1 "$checkpoint_r1" -2 "$checkpoint_r2")
    (( singles > 0 )) && command+=(-s "$checkpoint_single")
    if [[ "$SPADES_MODE" == "meta" && $pairs -gt 0 ]]; then
        command+=(--meta)
    fi

    if ! "${command[@]}" > "$round_dir/checkpoint_spades.stdout.log" \
            2> "$round_dir/checkpoint_spades.stderr.log"; then
        printf '%s\t%s\t0\t0\t0\t0\t0\t0\t0\tcheckpoint_failed\n' \
            "$round" "$ACCEPTED_COUNT" >> "$GRAPH_METRICS"
        log "WARNING: round-$round checkpoint SPAdes failed; recruitment remains protected by the read-count and growth caps."
        rm -rf -- "$checkpoint_tmp"
        CHECKPOINT_DECISION="continue"
        return 0
    fi
    if [[ ! -s "$graph" ]]; then
        printf '%s\t%s\t0\t0\t0\t0\t0\t0\t0\tgraph_unavailable\n' \
            "$round" "$ACCEPTED_COUNT" >> "$GRAPH_METRICS"
        log "WARNING: round-$round checkpoint produced no GFA graph; continuing with count-based guards."
        rm -rf -- "$checkpoint_tmp"
        CHECKPOINT_DECISION="continue"
        return 0
    fi

    awk -F '\t' '$1 == "S" && $3 != "*" { print ">" $2; print $3 }' \
        "$graph" > "$segment_fasta"
    if [[ ! -s "$segment_fasta" ]]; then
        printf '%s\t%s\t0\t0\t0\t0\t0\t0\t0\tsegments_unavailable\n' \
            "$round" "$ACCEPTED_COUNT" >> "$GRAPH_METRICS"
        rm -rf -- "$checkpoint_tmp"
        CHECKPOINT_DECISION="continue"
        return 0
    fi

    contig_anchor_hits "$segment_fasta" "$LEFT_INDEX" "$left_hits" \
        "$round_dir/checkpoint_left_anchor.log"
    contig_anchor_hits "$segment_fasta" "$RIGHT_INDEX" "$right_hits" \
        "$round_dir/checkpoint_right_anchor.log"

    python3 - "$graph" "$left_hits" "$right_hits" > "$stats" <<'PY'
import collections
import sys

graph_path, left_path, right_path = sys.argv[1:]
nodes = set()
adj = collections.defaultdict(set)
with open(graph_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if not fields:
            continue
        if fields[0] == "S" and len(fields) >= 2:
            nodes.add(fields[1])
        elif fields[0] == "L" and len(fields) >= 5:
            a, b = fields[1], fields[3]
            nodes.update((a, b))
            if a != b:
                adj[a].add(b)
                adj[b].add(a)

def names(path):
    result = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t", 1)
            if fields and fields[0]:
                result.add(fields[0])
    return result

left = names(left_path)
right = names(right_path)
anchors = (left | right) & nodes
seen = set()
anchored_components = []
for start in anchors:
    if start in seen:
        continue
    stack = [start]
    component = set()
    while stack:
        node = stack.pop()
        if node in component:
            continue
        component.add(node)
        stack.extend(adj[node] - component)
    seen.update(component)
    anchored_components.append(component)

anchored_nodes = 0
anchored_edges = 0
dual_components = 0
dual_branch_nodes = 0
dual_cycle_rank = 0
max_component_nodes = 0
for component in anchored_components:
    edge_count = sum(len(adj[node] & component) for node in component) // 2
    branch_count = sum(1 for node in component if len(adj[node] & component) > 2)
    cycle_rank = max(0, edge_count - len(component) + 1)
    anchored_nodes += len(component)
    anchored_edges += edge_count
    max_component_nodes = max(max_component_nodes, len(component))
    if component & left and component & right:
        dual_components += 1
        dual_branch_nodes += branch_count
        dual_cycle_rank += cycle_rank

print("\t".join(map(str, (
    len(anchored_components), dual_components, anchored_nodes, anchored_edges,
    dual_branch_nodes, dual_cycle_rank, max_component_nodes,
))))
PY

    local anchored_components dual_components anchored_nodes anchored_edges
    local dual_branches dual_cycles max_component_nodes status
    IFS=$'\t' read -r anchored_components dual_components anchored_nodes \
        anchored_edges dual_branches dual_cycles max_component_nodes < "$stats"
    status="continue"
    CHECKPOINT_DECISION="continue"
    if (( max_component_nodes >= GRAPH_STOP_NODES )); then
        status="stop_large_anchored_graph"
        CHECKPOINT_DECISION="stop"
    elif (( dual_components > 0 )); then
        if (( (GRAPH_STOP_BRANCHES > 0 && dual_branches >= GRAPH_STOP_BRANCHES) || \
              (GRAPH_STOP_CYCLES > 0 && dual_cycles >= GRAPH_STOP_CYCLES) )); then
            status="stop_connected_complex_graph"
        else
            status="stop_target_connected"
        fi
        CHECKPOINT_DECISION="stop"
    elif (( anchored_components == 0 )); then
        status="no_anchor_segments"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$round" "$ACCEPTED_COUNT" "$anchored_components" "$dual_components" \
        "$anchored_nodes" "$anchored_edges" "$dual_branches" "$dual_cycles" \
        "$max_component_nodes" "$status" >> "$GRAPH_METRICS"

    if [[ "$KEEP_INTERMEDIATES" == true ]]; then
        cp -- "$graph" "$round_dir/checkpoint_graph.gfa"
        cp -- "$segment_fasta" "$round_dir/checkpoint_segments.fasta"
    fi
    rm -rf -- "$checkpoint_tmp"
    CHECKPOINT_STATUS="$status"
    CHECKPOINT_DUAL_COMPONENTS="$dual_components"
    CHECKPOINT_BRANCHES="$dual_branches"
    CHECKPOINT_CYCLES="$dual_cycles"
    CHECKPOINT_MAX_NODES="$max_component_nodes"
}

extract_fasta_by_names() {
    local names="$1"
    local input="$2"
    local output="$3"
    if [[ ! -s "$names" ]]; then
        : > "$output"
        return 0
    fi
    awk '
        NR == FNR { keep[$1] = 1; next }
        /^>/ {
            name = substr($0, 2)
            sub(/[[:space:]].*$/, "", name)
            printing = (name in keep)
        }
        printing { print }
    ' "$names" "$input" > "$output"
}

fasta_stats() {
    awk '
        /^>/ { sequences++; next }
        { bases += length($0) }
        END { printf "%d\t%d\n", sequences + 0, bases + 0 }
    ' "$1"
}

validate_and_orient_dual_anchors() {
    local contigs="$1"
    local validation_dir="$OUTDIR/validation"
    mkdir -p "$validation_dir"

    local blast_fields='6 qseqid qlen qstart qend sseqid slen sstart send sstrand length pident evalue bitscore'
    blastn -query "$contigs" -subject "$LEFT_DB" -task blastn \
        -evalue 1e-20 -max_target_seqs 25 -max_hsps 3 \
        -outfmt "$blast_fields" > "$validation_dir/contigs_vs_18S.tsv"
    blastn -query "$contigs" -subject "$RIGHT_DB" -task blastn \
        -evalue 1e-20 -max_target_seqs 25 -max_hsps 3 \
        -outfmt "$blast_fields" > "$validation_dir/contigs_vs_28S.tsv"

    python3 - "$contigs" "$validation_dir/contigs_vs_18S.tsv" \
        "$validation_dir/contigs_vs_28S.tsv" "$OUTDIR/final" \
        "$ANCHOR_MIN_ALIGNED" "$ANCHOR_MIN_IDENTITY" \
        "$LAYOUT_OVERLAP_TOLERANCE" <<'PY'
import sys
from pathlib import Path

fasta, left_tsv, right_tsv, outdir, min_aln, min_ident, tolerance = sys.argv[1:]
min_aln, min_ident, tolerance = int(min_aln), float(min_ident), int(tolerance)
outdir = Path(outdir)

seqs = {}
name = None
with open(fasta) as handle:
    for line in handle:
        line = line.rstrip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = ''
        elif name:
            seqs[name] += line

def best_hits(path):
    hits = {}
    with open(path) as handle:
        for line in handle:
            f = line.rstrip('\n').split('\t')
            if len(f) < 13:
                continue
            aln, ident, bits = int(f[9]), float(f[10]), float(f[12])
            if aln < min_aln or ident / 100.0 < min_ident:
                continue
            hit = dict(ref=f[4], qstart=int(f[2]), qend=int(f[3]),
                       strand=f[8], aln=aln, identity=ident, bitscore=bits)
            if f[0] not in hits or bits > hits[f[0]]['bitscore']:
                hits[f[0]] = hit
    return hits

left, right = best_hits(left_tsv), best_hits(right_tsv)
comp = str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn',
                     'TGCAYRKMVHDBNtgcayrkmvhdbn')
def rc(seq):
    return seq.translate(comp)[::-1]
def emit(handle, header, seq):
    handle.write('>' + header + '\n')
    for i in range(0, len(seq), 80):
        handle.write(seq[i:i+80] + '\n')

oriented = open(outdir / 'oriented_dual_anchor_contigs.fasta', 'w')
ambiguous = open(outdir / 'ambiguous_dual_anchor_contigs.fasta', 'w')
report = open(outdir / 'contig_validation.tsv', 'w')
report.write('contig\tlength\t18S_qstart\t18S_qend\t18S_strand\t18S_aligned_bp\t18S_identity\t'
             '28S_qstart\t28S_qend\t28S_strand\t28S_aligned_bp\t28S_identity\t'
             'inter_anchor_bp\torientation\tlayout_status\n')

for contig, seq in seqs.items():
    l, r = left.get(contig), right.get(contig)
    if not l or not r:
        continue
    gap = max(l['qstart'], r['qstart']) - min(l['qend'], r['qend']) - 1
    orientation, status = 'unresolved', 'FAIL_inconsistent_anchors'
    if (l['strand'] == r['strand'] == 'plus' and
            l['qstart'] < r['qstart'] and l['qend'] <= r['qstart'] + tolerance):
        orientation, status = 'forward', 'PASS_dual_anchor'
        emit(oriented, contig, seq)
    elif (l['strand'] == r['strand'] == 'minus' and
            r['qstart'] < l['qstart'] and r['qend'] <= l['qstart'] + tolerance):
        orientation, status = 'reverse_complemented', 'PASS_dual_anchor'
        emit(oriented, contig + ' orientation=reverse_complemented', rc(seq))
    else:
        emit(ambiguous, contig + ' layout=ambiguous', seq)
    report.write(
        f"{contig}\t{len(seq)}\t{l['qstart']}\t{l['qend']}\t{l['strand']}\t{l['aln']}\t{l['identity']:.3f}\t"
        f"{r['qstart']}\t{r['qend']}\t{r['strand']}\t{r['aln']}\t{r['identity']:.3f}\t"
        f"{gap}\t{orientation}\t{status}\n")

oriented.close(); ambiguous.close(); report.close()
PY
}

enumerate_graph_paths() {
    local graph="$OUTDIR/assembly/spades/assembly_graph_with_scaffolds.gfa"
    local graph_dir="$OUTDIR/validation/graph_paths"
    local nodes="$graph_dir/graph_segments.fasta"
    local segment_map="$graph_dir/segment_map.tsv"
    local left_hits="$graph_dir/segments_vs_18S.tsv"
    local right_hits="$graph_dir/segments_vs_28S.tsv"
    local paths="$OUTDIR/final/graph_candidate_paths.fasta"
    local report="$OUTDIR/final/graph_path_summary.tsv"

    : > "$paths"
    printf 'path_id\tlength_bp\tnode_count\tminimum_graph_depth\tmean_graph_depth\tleft_segment\tright_segment\tnode_path\tjunction_positions\tpath_status\n' \
        > "$report"
    [[ "$GRAPH_PATHS" == true && -s "$graph" ]] || return 0
    mkdir -p "$graph_dir"

    python3 - "$graph" "$nodes" "$segment_map" <<'PY'
import sys
from pathlib import Path

graph, fasta_out, map_out = map(Path, sys.argv[1:])
segments = []
with graph.open() as handle:
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 3 or fields[0] != 'S' or fields[2] == '*':
            continue
        depth = 0.0
        for tag in fields[3:]:
            if tag.startswith('DP:f:'):
                try: depth = float(tag[5:])
                except ValueError: pass
        segments.append((fields[1], fields[2].upper(), depth))

with fasta_out.open('w') as fasta, map_out.open('w') as mapping:
    mapping.write('safe_id\tsegment_id\tlength_bp\tgraph_depth\n')
    for number, (segment, sequence, depth) in enumerate(segments, 1):
        safe = f'GRAPHSEG{number:07d}'
        mapping.write(f'{safe}\t{segment}\t{len(sequence)}\t{depth:.6f}\n')
        fasta.write(f'>{safe}\n')
        for start in range(0, len(sequence), 80):
            fasta.write(sequence[start:start + 80] + '\n')
PY

    [[ -s "$nodes" ]] || return 0
    local fields='6 qseqid qstart qend sstrand length pident bitscore'
    blastn -query "$nodes" -subject "$LEFT_DB" -task blastn \
        -evalue 1e-20 -max_target_seqs 25 -max_hsps 3 \
        -outfmt "$fields" > "$left_hits"
    blastn -query "$nodes" -subject "$RIGHT_DB" -task blastn \
        -evalue 1e-20 -max_target_seqs 25 -max_hsps 3 \
        -outfmt "$fields" > "$right_hits"

    python3 - "$graph" "$segment_map" "$left_hits" "$right_hits" \
        "$paths" "$report" "$ANCHOR_MIN_ALIGNED" "$ANCHOR_MIN_IDENTITY" \
        "$MAX_GRAPH_PATHS" "$MAX_GRAPH_NODES" "$MAX_GRAPH_BP" \
        "$MIN_GRAPH_DEPTH" "$graph_dir/path_limit_reached.txt" <<'PY'
import hashlib, re, sys
from collections import defaultdict
from pathlib import Path

(graph_path, map_path, left_path, right_path, fasta_path, report_path,
 min_aln, min_identity, max_paths, max_nodes, max_bp, min_depth,
 limit_path) = sys.argv[1:]
min_aln, min_identity = int(min_aln), float(min_identity)
max_paths, max_nodes, max_bp = int(max_paths), int(max_nodes), int(max_bp)
min_depth = float(min_depth)

safe_to_segment = {}
with open(map_path) as handle:
    next(handle, None)
    for line in handle:
        safe, segment, *_ = line.rstrip('\n').split('\t')
        safe_to_segment[safe] = segment

seqs, depths = {}, {}
raw_links = []
with open(graph_path) as handle:
    for line in handle:
        f = line.rstrip('\n').split('\t')
        if not f:
            continue
        if f[0] == 'S' and len(f) >= 3 and f[2] != '*':
            seqs[f[1]] = f[2].upper()
            depth = 0.0
            for tag in f[3:]:
                if tag.startswith('DP:f:'):
                    try: depth = float(tag[5:])
                    except ValueError: pass
            depths[f[1]] = depth
        elif f[0] == 'L' and len(f) >= 6:
            match = re.fullmatch(r'(\d+)[M=]', f[5])
            overlap = int(match.group(1)) if match else 0
            raw_links.append((f[1], f[2], f[3], f[4], overlap))

def flip(orientation):
    return '-' if orientation == '+' else '+'

adjacency = defaultdict(dict)
for left, left_o, right, right_o, overlap in raw_links:
    if left not in seqs or right not in seqs:
        continue
    adjacency[(left, left_o)][(right, right_o)] = max(
        overlap, adjacency[(left, left_o)].get((right, right_o), 0))
    adjacency[(right, flip(right_o))][(left, flip(left_o))] = max(
        overlap, adjacency[(right, flip(right_o))].get((left, flip(left_o)), 0))

def read_hits(path):
    best = {}
    with open(path) as handle:
        for line in handle:
            f = line.rstrip('\n').split('\t')
            if len(f) < 7 or f[0] not in safe_to_segment:
                continue
            aligned, identity, bitscore = int(f[4]), float(f[5]), float(f[6])
            if aligned < min_aln or identity / 100.0 < min_identity:
                continue
            segment = safe_to_segment[f[0]]
            hit = {'state': (segment, '+' if f[3] == 'plus' else '-'),
                   'qstart': int(f[1]), 'qend': int(f[2]), 'bitscore': bitscore}
            if segment not in best or bitscore > best[segment]['bitscore']:
                best[segment] = hit
    return best

left_hits, right_hits = read_hits(left_path), read_hits(right_path)
left_states = {hit['state'] for hit in left_hits.values()}
right_states = {hit['state'] for hit in right_hits.values()}

def oriented_position(hit):
    segment, orientation = hit['state']
    midpoint = (hit['qstart'] + hit['qend']) / 2.0
    return midpoint if orientation == '+' else len(seqs[segment]) - midpoint + 1

same_segment_valid = set()
for segment in set(left_hits) & set(right_hits):
    left, right = left_hits[segment], right_hits[segment]
    if left['state'] == right['state'] and oriented_position(left) < oriented_position(right):
        same_segment_valid.add(left['state'])

comp = str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN')
def sequence_for(state):
    segment, orientation = state
    sequence = seqs[segment]
    return sequence if orientation == '+' else sequence.translate(comp)[::-1]

def usable(state):
    depth = depths.get(state[0], 0.0)
    return depth <= 0 or depth >= min_depth

def assembled_sequence(path):
    sequence = sequence_for(path[0])
    junctions = []
    for previous, current in zip(path, path[1:]):
        overlap = adjacency[previous].get(current, 0)
        next_sequence = sequence_for(current)
        if overlap <= 0:
            maximum = min(len(sequence), len(next_sequence), 255)
            for size in range(maximum, 0, -1):
                if sequence[-size:] == next_sequence[:size]:
                    overlap = size
                    break
        junctions.append(len(sequence))
        sequence += next_sequence[overlap:]
    return sequence, junctions

found, truncated = [], False
for start in sorted(left_states):
    if not usable(start):
        continue
    stack = [([start], len(sequence_for(start)))]
    while stack:
        path, path_bp = stack.pop()
        current = path[-1]
        if current in right_states and (len(path) > 1 or current in same_segment_valid):
            found.append(tuple(path))
            if len(found) >= max_paths:
                truncated = True
                break
            continue
        if len(path) >= max_nodes or path_bp >= max_bp:
            continue
        used_segments = {state[0] for state in path}
        for nxt, overlap in sorted(adjacency.get(current, {}).items(), reverse=True):
            if nxt[0] in used_segments or not usable(nxt):
                continue
            proposed_bp = path_bp + max(0, len(sequence_for(nxt)) - overlap)
            if proposed_bp <= max_bp:
                stack.append((path + [nxt], proposed_bp))
    if truncated:
        break

unique = []
seen_sequences = set()
for path in found:
    sequence, junctions = assembled_sequence(path)
    digest = hashlib.sha256(sequence.encode()).digest()
    if digest in seen_sequences:
        continue
    seen_sequences.add(digest)
    unique.append((path, sequence, junctions))

with open(fasta_path, 'w') as fasta, open(report_path, 'w') as report:
    report.write('path_id\tlength_bp\tnode_count\tminimum_graph_depth\tmean_graph_depth\t'
                 'left_segment\tright_segment\tnode_path\tjunction_positions\tpath_status\n')
    for number, (path, sequence, junctions) in enumerate(unique, 1):
        path_id = f'ITSME_PATH_{number:06d}'
        observed_depths = [depths.get(state[0], 0.0) for state in path
                           if depths.get(state[0], 0.0) > 0]
        minimum = min(observed_depths) if observed_depths else 0.0
        weights = [len(seqs[state[0]]) for state in path]
        weighted = sum(depths.get(state[0], 0.0) * weight
                       for state, weight in zip(path, weights))
        mean = weighted / sum(weights) if weights else 0.0
        node_text = ','.join(f'{segment}{orientation}' for segment, orientation in path)
        junction_text = ','.join(map(str, junctions)) if junctions else 'NA'
        status = 'CANDIDATE_COMPLETE_PATH'
        fasta.write(f'>{path_id} nodes={node_text} min_graph_depth={minimum:.3f}\n')
        for start in range(0, len(sequence), 80):
            fasta.write(sequence[start:start + 80] + '\n')
        report.write(f'{path_id}\t{len(sequence)}\t{len(path)}\t{minimum:.6f}\t{mean:.6f}\t'
                     f'{path[0][0]}\t{path[-1][0]}\t{node_text}\t{junction_text}\t{status}\n')

Path(limit_path).write_text('true\n' if truncated else 'false\n')
PY

    local count
    count=$(awk '/^>/ { n++ } END { print n + 0 }' "$paths")
    log "Graph-path analysis recovered $count distinct bounded 18S-to-28S candidate path(s)."
    if [[ -s "$graph_dir/path_limit_reached.txt" ]] && \
       grep -qx 'true' "$graph_dir/path_limit_reached.txt"; then
        log "WARNING: graph-path enumeration reached --max-graph-paths=$MAX_GRAPH_PATHS; all emitted paths are retained, but additional graph alternatives may exist."
    fi
}

merge_locus_candidates() {
    local graph_paths="$OUTDIR/final/graph_candidate_paths.fasta"
    local contigs="$OUTDIR/final/rrna_candidate_contigs.fasta"
    local output="$OUTDIR/final/all_locus_candidates.fasta"
    local provenance="$OUTDIR/final/locus_candidate_provenance.tsv"
    python3 - "$graph_paths" "$contigs" "$output" "$provenance" <<'PY'
import hashlib, sys
from pathlib import Path

graph, contigs, output, provenance = map(Path, sys.argv[1:])
comp = str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN')

def records(path):
    if not path.exists():
        return
    name, description, sequence = None, '', []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, description, ''.join(sequence).upper()
                description = line[1:]
                name = description.split()[0]
                sequence = []
            elif name is not None:
                sequence.append(line.strip())
        if name is not None:
            yield name, description, ''.join(sequence).upper()

seen = set()
with output.open('w') as fasta, provenance.open('w') as table:
    table.write('candidate\tsource\tsource_header\tlength_bp\n')
    for source, path in (('graph_path', graph), ('spades_contig', contigs)):
        for name, description, sequence in records(path):
            reverse = sequence.translate(comp)[::-1]
            canonical = min(sequence, reverse)
            digest = hashlib.sha256(canonical.encode()).digest()
            if digest in seen:
                continue
            seen.add(digest)
            fasta.write(f'>{description}\n')
            for start in range(0, len(sequence), 80):
                fasta.write(sequence[start:start + 80] + '\n')
            table.write(f'{name}\t{source}\t{description}\t{len(sequence)}\n')
PY
}

catalog_locus_candidates() {
    local candidates="$OUTDIR/final/all_locus_candidates.fasta"
    local left_names="$OUTDIR/final/left_anchor_contig_names.txt"
    local right_names="$OUTDIR/final/right_anchor_contig_names.txt"
    python3 - "$candidates" "$left_names" "$right_names" \
        "$OUTDIR/final/contig_validation.tsv" \
        "$OUTDIR/final/locus_candidate_provenance.tsv" "$OUTDIR/final" <<'PY'
import csv, sys
from pathlib import Path

candidates, left_path, right_path, validation_path, provenance_path, outdir = sys.argv[1:]
outdir = Path(outdir)

def names(path):
    with open(path) as handle:
        return {line.strip().split()[0] for line in handle if line.strip()}

def fasta(path):
    records, name, description, sequence = {}, None, '', []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    records[name] = (description, ''.join(sequence))
                description = line[1:]
                name = description.split()[0]
                sequence = []
            elif name is not None:
                sequence.append(line.strip())
        if name is not None:
            records[name] = (description, ''.join(sequence))
    return records

records = fasta(candidates)
left, right = names(left_path), names(right_path)
validation = {}
with open(validation_path) as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        validation[row['contig']] = row.get('layout_status', 'NA')
provenance = {}
with open(provenance_path) as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        provenance[row['candidate']] = row.get('source', 'NA')

outputs = {
    'PARTIAL_18S': open(outdir / 'partial_18S_contigs.fasta', 'w'),
    'PARTIAL_28S': open(outdir / 'partial_28S_contigs.fasta', 'w'),
    'AMBIGUOUS_DUAL_ANCHOR': open(outdir / 'ambiguous_locus_candidates.fasta', 'w'),
    'SEED_ANCHORED_OTHER': open(outdir / 'other_seed_anchored_contigs.fasta', 'w'),
}
catalog = open(outdir / 'locus_catalog.tsv', 'w')
catalog.write('candidate\tlength_bp\tlocus_type\tsource\tlayout_status\n')

for name in sorted(records):
    description, sequence = records[name]
    layout = validation.get(name, 'NA')
    if layout == 'PASS_dual_anchor':
        locus_type = 'COMPLETE_DUAL_ANCHOR'
    elif name in left and name in right:
        locus_type = 'AMBIGUOUS_DUAL_ANCHOR'
    elif name in left:
        locus_type = 'PARTIAL_18S'
    elif name in right:
        locus_type = 'PARTIAL_28S'
    else:
        locus_type = 'SEED_ANCHORED_OTHER'
    catalog.write(f'{name}\t{len(sequence)}\t{locus_type}\t{provenance.get(name, "NA")}\t{layout}\n')
    if locus_type in outputs:
        outputs[locus_type].write(f'>{description} locus_type={locus_type}\n')
        for start in range(0, len(sequence), 80):
            outputs[locus_type].write(sequence[start:start + 80] + '\n')

catalog.close()
for handle in outputs.values():
    handle.close()

with open(outdir / 'partial_locus_contigs.fasta', 'w') as combined:
    for filename in ('partial_18S_contigs.fasta', 'partial_28S_contigs.fasta',
                     'ambiguous_locus_candidates.fasta', 'other_seed_anchored_contigs.fasta'):
        with open(outdir / filename) as source:
            combined.write(source.read())
PY
}

run_itsx_validation() {
    local oriented="$OUTDIR/final/oriented_dual_anchor_contigs.fasta"
    local itsx_dir="$OUTDIR/validation/itsx"
    : > "$OUTDIR/final/complete_ITS.fasta"
    : > "$OUTDIR/final/ITS1.fasta"
    : > "$OUTDIR/final/5_8S.fasta"
    : > "$OUTDIR/final/ITS2.fasta"
    printf 'contig\tlength\tSSU\tITS1\t5.8S\tITS2\tLSU\tstatus\n' \
        > "$OUTDIR/final/itsx_validation.tsv"
    [[ "$RUN_ITSX" == true && -s "$oriented" ]] || return 0
    mkdir -p "$itsx_dir"
    ITSx -i "$oriented" -o "$itsx_dir/rrna" -t "$ITSX_TAXA" \
        --nhmmer T --save_regions all --positions T \
        --detailed_results T --preserve T --cpu "$THREADS" \
        > "$itsx_dir/itsx.stdout.log" 2> "$itsx_dir/itsx.stderr.log" || \
        die "ITSx validation failed; inspect $itsx_dir/itsx.stderr.log"
    [[ -s "$itsx_dir/rrna.full.fasta" ]] && \
        cp -- "$itsx_dir/rrna.full.fasta" "$OUTDIR/final/complete_ITS.fasta"
    [[ -s "$itsx_dir/rrna.ITS1.fasta" ]] && \
        cp -- "$itsx_dir/rrna.ITS1.fasta" "$OUTDIR/final/ITS1.fasta"
    [[ -s "$itsx_dir/rrna.5_8S.fasta" ]] && \
        cp -- "$itsx_dir/rrna.5_8S.fasta" "$OUTDIR/final/5_8S.fasta"
    [[ -s "$itsx_dir/rrna.ITS2.fasta" ]] && \
        cp -- "$itsx_dir/rrna.ITS2.fasta" "$OUTDIR/final/ITS2.fasta"
    cp -- "$itsx_dir/rrna.positions.txt" "$OUTDIR/final/itsx_positions.txt"
    python3 - "$itsx_dir/rrna.positions.txt" "$OUTDIR/final/itsx_validation.tsv" <<'PY'
import re, sys
with open(sys.argv[2], 'a') as out, open(sys.argv[1]) as src:
    for line in src:
        fields = line.rstrip('\n').split('\t')
        if not fields:
            continue
        text = '\t'.join(fields)
        def region(label):
            m = re.search(rf'{re.escape(label)}:\s*([^\t]+)', text)
            return m.group(1).strip() if m else 'NA'
        length_match = re.search(r'\b(\d+) bp\.', text)
        values = [region(x) for x in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')]
        complete = all(re.fullmatch(r'\d+-\d+', x) for x in values)
        chimeric = 'Chimeric!' in text or any('Overlap' in x for x in values)
        status = 'FAIL_chimeric' if chimeric else ('PASS_complete_ITS' if complete else 'PARTIAL_ITS')
        out.write('\t'.join([fields[0], length_match.group(1) if length_match else 'NA',
                             *values, status]) + '\n')
PY
}

run_competitive_readback() {
    local targets="$OUTDIR/final/oriented_dual_anchor_contigs.fasta"
    local readback_dir="$OUTDIR/validation/readback"
    : > "$OUTDIR/final/readback_coverage.tsv"
    printf 'candidate\tlength_bp\tall_mapped_reads\tcoverage_percent\tmean_depth\tmean_mapq\tunique_mapped_reads\tunique_coverage_percent\tunique_mean_depth\tgraph_junctions\tjunctions_with_read_support\tminimum_junction_spanning_reads\n' \
        > "$OUTDIR/final/path_support.tsv"
    printf 'candidate\tposition\treference\talternate\tquality\ttotal_depth\tallele_depths\n' \
        > "$OUTDIR/final/residual_variants.tsv"
    [[ "$READBACK" == true && -s "$targets" ]] || return 0
    mkdir -p "$readback_dir"
    bowtie2-build --threads "$THREADS" "$targets" "$readback_dir/candidates" \
        > "$readback_dir/build.stdout.log" 2> "$readback_dir/build.stderr.log"
    local command=(bowtie2 --very-sensitive --no-unal -p "$THREADS"
        -x "$readback_dir/candidates" -X "$MAX_INSERT" -1 "$READS_1" -2 "$READS_2")
    [[ -n "$READS_SINGLE" ]] && command+=(-U "$READS_SINGLE")
    "${command[@]}" 2> "$readback_dir/bowtie2.log" | \
        samtools view -@ "$THREADS" -b -F 4 - | \
        samtools sort -@ "$THREADS" -o "$readback_dir/candidates.bam" -
    samtools index "$readback_dir/candidates.bam"
    samtools coverage "$readback_dir/candidates.bam" \
        > "$OUTDIR/final/readback_coverage.tsv"
    samtools flagstat "$readback_dir/candidates.bam" \
        > "$OUTDIR/final/readback_flagstat.txt"

    samtools view -@ "$THREADS" -b -q 20 \
        -o "$readback_dir/candidates.mapq20.bam" "$readback_dir/candidates.bam"
    samtools index "$readback_dir/candidates.mapq20.bam"
    samtools coverage "$readback_dir/candidates.mapq20.bam" \
        > "$readback_dir/coverage.mapq20.tsv"
    samtools view -@ "$THREADS" -F 2308 \
        -o "$readback_dir/primary_alignments.sam" "$readback_dir/candidates.bam"
    python3 - "$OUTDIR/final/readback_coverage.tsv" \
        "$readback_dir/coverage.mapq20.tsv" \
        "$OUTDIR/final/graph_path_summary.tsv" \
        "$readback_dir/primary_alignments.sam" \
        "$OUTDIR/final/path_support.tsv" <<'PY'
import csv, re, sys

def read(path):
    result = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            name = row.get('#rname') or row.get('rname')
            if name:
                result[name] = row
    return result

all_reads, unique = read(sys.argv[1]), read(sys.argv[2])
junctions = {}
with open(sys.argv[3]) as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        text = row.get('junction_positions', 'NA')
        junctions[row['path_id']] = [] if text == 'NA' else [int(x) for x in text.split(',')]

spanning = {name: [0] * len(positions) for name, positions in junctions.items()}
cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
with open(sys.argv[4]) as handle:
    for line in handle:
        if not line or line.startswith('@'):
            continue
        fields_row = line.rstrip('\n').split('\t')
        if len(fields_row) < 6 or fields_row[2] not in junctions:
            continue
        start = int(fields_row[3])
        reference_bases = sum(int(length) for length, operation in cigar_pattern.findall(fields_row[5])
                              if operation in {'M','D','N','=','X'})
        end = start + reference_bases - 1
        for index, junction in enumerate(junctions[fields_row[2]]):
            if start <= junction - 10 and end >= junction + 10:
                spanning[fields_row[2]][index] += 1

fields = ['candidate','length_bp','all_mapped_reads','coverage_percent','mean_depth',
          'mean_mapq','unique_mapped_reads','unique_coverage_percent','unique_mean_depth',
          'graph_junctions','junctions_with_read_support','minimum_junction_spanning_reads']
with open(sys.argv[5], 'w', newline='') as handle:
    out = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
    out.writeheader()
    for name in sorted(all_reads):
        overall, high = all_reads[name], unique.get(name, {})
        counts = spanning.get(name, [])
        out.writerow({
            'candidate': name,
            'length_bp': overall.get('endpos', 'NA'),
            'all_mapped_reads': overall.get('numreads', '0'),
            'coverage_percent': overall.get('coverage', '0'),
            'mean_depth': overall.get('meandepth', '0'),
            'mean_mapq': overall.get('meanmapq', '0'),
            'unique_mapped_reads': high.get('numreads', '0'),
            'unique_coverage_percent': high.get('coverage', '0'),
            'unique_mean_depth': high.get('meandepth', '0'),
            'graph_junctions': len(counts) if name in junctions else 'NA',
            'junctions_with_read_support': sum(value > 0 for value in counts) if name in junctions else 'NA',
            'minimum_junction_spanning_reads': min(counts) if counts else ('NA' if name not in junctions else 0)})
PY
    [[ "$KEEP_INTERMEDIATES" == true ]] || rm -f -- "$readback_dir/primary_alignments.sam"

    if [[ "$CALL_VARIANTS" == true ]]; then
        samtools faidx "$targets"
        bcftools mpileup -Ou -f "$targets" \
            -q "$VARIANT_MIN_MAPQ" -Q "$VARIANT_MIN_BASEQ" \
            -d "$VARIANT_MAX_DEPTH" -a FORMAT/DP,FORMAT/AD \
            "$readback_dir/candidates.bam" 2> "$readback_dir/bcftools_mpileup.log" | \
            bcftools call -m -v -Oz \
                -o "$OUTDIR/final/residual_variants.vcf.gz" \
                2> "$readback_dir/bcftools_call.log"
        bcftools index -f "$OUTDIR/final/residual_variants.vcf.gz"
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP[\t%AD]\n' \
            "$OUTDIR/final/residual_variants.vcf.gz" | \
            awk 'BEGIN { OFS="\t"; print "candidate","position","reference","alternate","quality","total_depth","allele_depths" } { print }' \
            > "$OUTDIR/final/residual_variants.tsv"
    fi
}

run_ncbi_classification() {
    local blast_dir="$OUTDIR/validation/ncbi_blast"
    mkdir -p "$blast_dir"
    local fields='6 qseqid qlen saccver pident length qcovhsp evalue bitscore staxids stitle'

    : > "$OUTDIR/final/ncbi_SSU_hits.tsv"
    : > "$OUTDIR/final/ncbi_ITS_hits.tsv"
    : > "$OUTDIR/final/ncbi_LSU_hits.tsv"
    printf 'marker\tquery\trank\tquery_length\taccession\tpercent_identity\taligned_bp\tquery_coverage_hsp\tevalue\tbitscore\ttaxids\tscientific_names\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\ttitle\n' \
        > "$OUTDIR/final/ncbi_blast_top_hits.tsv"
    [[ "$RUN_NCBI_BLAST" == true ]] || return 0

    log "NCBI classification: searching SSU, ITS, and LSU targeted databases."
    : > "$blast_dir/database_info.txt"
    for database in SSU_eukaryote_rRNA ITS_eukaryote_sequences LSU_eukaryote_rRNA; do
        printf '===== %s =====\n' "$database" >> "$blast_dir/database_info.txt"
        blastdbcmd -db "$NCBI_DB_DIR/$database" -info >> "$blast_dir/database_info.txt"
    done
    local ssu_queries="$blast_dir/SSU_queries.fasta"
    local lsu_queries="$blast_dir/LSU_queries.fasta"
    : > "$ssu_queries"
    : > "$lsu_queries"
    for query in "$OUTDIR/validation/itsx/rrna.SSU.fasta" \
                 "$OUTDIR/final/partial_18S_contigs.fasta" \
                 "$OUTDIR/final/ambiguous_locus_candidates.fasta"; do
        [[ -s "$query" ]] && awk '{ print }' "$query" >> "$ssu_queries"
    done
    for query in "$OUTDIR/validation/itsx/rrna.LSU.fasta" \
                 "$OUTDIR/final/partial_28S_contigs.fasta" \
                 "$OUTDIR/final/ambiguous_locus_candidates.fasta"; do
        [[ -s "$query" ]] && awk '{ print }' "$query" >> "$lsu_queries"
    done
    if [[ -s "$ssu_queries" ]]; then
        blastn -query "$ssu_queries" \
            -db "$NCBI_DB_DIR/SSU_eukaryote_rRNA" -task blastn \
            -evalue 1e-20 -max_target_seqs "$NCBI_MAX_HITS" -max_hsps 1 \
            -num_threads "$THREADS" -outfmt "$fields" \
            -out "$OUTDIR/final/ncbi_SSU_hits.tsv"
    fi
    if [[ -s "$OUTDIR/final/complete_ITS.fasta" ]]; then
        blastn -query "$OUTDIR/final/complete_ITS.fasta" \
            -db "$NCBI_DB_DIR/ITS_eukaryote_sequences" -task blastn \
            -word_size 7 -dust no -evalue 1e-5 \
            -max_target_seqs "$NCBI_MAX_HITS" -max_hsps 1 \
            -num_threads "$THREADS" -outfmt "$fields" \
            -out "$OUTDIR/final/ncbi_ITS_hits.tsv"
    fi
    if [[ -s "$lsu_queries" ]]; then
        blastn -query "$lsu_queries" \
            -db "$NCBI_DB_DIR/LSU_eukaryote_rRNA" -task blastn \
            -evalue 1e-20 -max_target_seqs "$NCBI_MAX_HITS" -max_hsps 1 \
            -num_threads "$THREADS" -outfmt "$fields" \
            -out "$OUTDIR/final/ncbi_LSU_hits.tsv"
    fi

    python3 - "$OUTDIR/final/ncbi_blast_top_hits.tsv" "$NCBI_REPORT_HITS" "$TAXDUMP_DIR" \
        "SSU=$OUTDIR/final/ncbi_SSU_hits.tsv" \
        "ITS=$OUTDIR/final/ncbi_ITS_hits.tsv" \
        "LSU=$OUTDIR/final/ncbi_LSU_hits.tsv" <<'PY'
import sys
from pathlib import Path
out_path, report_n, taxdump, *specs = sys.argv[1:]
report_n = int(report_n)
taxdump = Path(taxdump)

parents, ranks, names, merged = {}, {}, {}, {}
with (taxdump / 'nodes.dmp').open() as handle:
    for line in handle:
        f = [x.strip() for x in line.split('|')]
        parents[f[0]], ranks[f[0]] = f[1], f[2]
with (taxdump / 'names.dmp').open() as handle:
    for line in handle:
        f = [x.strip() for x in line.split('|')]
        if len(f) > 3 and f[3] == 'scientific name':
            names[f[0]] = f[1]
if (taxdump / 'merged.dmp').exists():
    with (taxdump / 'merged.dmp').open() as handle:
        for line in handle:
            f = [x.strip() for x in line.split('|')]
            merged[f[0]] = f[1]

wanted = ('superkingdom','kingdom','phylum','class','order','family','genus','species')
def lineage(taxids):
    taxid = next((x for x in taxids.replace(',', ';').split(';') if x.isdigit()), '')
    while taxid in merged:
        taxid = merged[taxid]
    found, seen = {}, set()
    current = taxid
    while current and current not in seen and current in parents:
        seen.add(current)
        rank = ranks.get(current, '')
        if rank in wanted and rank not in found:
            found[rank] = names.get(current, 'NA')
        parent = parents[current]
        if parent == current:
            break
        current = parent
    return [found.get(rank, 'NA') for rank in wanted], names.get(taxid, 'NA')

with open(out_path, 'a') as out:
    for spec in specs:
        marker, path = spec.split('=', 1)
        grouped = {}
        with open(path) as handle:
            for line in handle:
                f = line.rstrip('\n').split('\t', 9)
                if len(f) != 10:
                    continue
                grouped.setdefault(f[0], []).append(f)
        for query in sorted(grouped):
            hits = sorted(grouped[query], key=lambda x: float(x[7]), reverse=True)
            for rank, f in enumerate(hits[:report_n], 1):
                taxonomy, resolved_name = lineage(f[8])
                out.write('\t'.join([marker, query, str(rank), f[1], f[2], f[3],
                                     f[4], f[5], f[6], f[7], f[8], resolved_name,
                                     *taxonomy, f[9]]) + '\n')
PY
}

build_locus_summary() {
    python3 - "$OUTDIR/final" "$EXPECTED_TAXONOMY" \
        "$TAXONOMY_NEAR_TOP_FRACTION" <<'PY'
import csv, re, sys
from pathlib import Path

root = Path(sys.argv[1])
expected = sys.argv[2].strip()
near_top_fraction = float(sys.argv[3])
ranks = ('domain','kingdom','phylum','class','order','family','genus','species')

def rows(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def region_length(value):
    m = re.fullmatch(r'(\d+)-(\d+)', value or '')
    return str(abs(int(m.group(2)) - int(m.group(1))) + 1) if m else 'NA'

def resolved(value):
    return bool(value and value not in {'NA', 'N/A', 'Unclassified', 'unclassified'})

layout = {r['contig']: r for r in rows(root / 'contig_validation.tsv')}
itsx = {r['contig']: r for r in rows(root / 'itsx_validation.tsv')}
catalog = {r['candidate']: r for r in rows(root / 'locus_catalog.tsv')}
coverage = {}
for r in rows(root / 'readback_coverage.tsv'):
    key = r.get('#rname') or r.get('rname')
    if key:
        coverage[key] = r

all_hits = rows(root / 'ncbi_blast_top_hits.tsv')
best, grouped = {}, {}
for row in all_hits:
    grouped.setdefault((row['query'], row['marker']), []).append(row)
    if row['rank'] == '1':
        best[(row['query'], row['marker'])] = row

def hit(contig, marker, field):
    return best.get((contig, marker), {}).get(field, 'NA') or 'NA'

def marker_consensus(contig, marker):
    marker_hits = grouped.get((contig, marker), [])
    if not marker_hits:
        return {rank: 'NA' for rank in ranks}, 0
    top_bitscore = max(float(row['bitscore']) for row in marker_hits)
    selected = [row for row in marker_hits
                if float(row['bitscore']) >= top_bitscore * near_top_fraction]
    consensus = {rank: 'NA' for rank in ranks}
    for rank in ranks:
        values = {row[rank].strip() for row in selected if resolved(row.get(rank, ''))}
        if len(values) != 1:
            break
        consensus[rank] = values.pop()
    return consensus, len(selected)

def lineage_text(consensus):
    values = [consensus[rank] for rank in ranks if resolved(consensus.get(rank, ''))]
    return '; '.join(values) if values else 'NA'

def locus_consensus(ssu, lsu):
    consensus = {rank: 'NA' for rank in ranks}
    for rank in ranks:
        left, right = ssu.get(rank, 'NA'), lsu.get(rank, 'NA')
        if not resolved(left) or not resolved(right) or left.casefold() != right.casefold():
            break
        consensus[rank] = left
    return consensus

def taxonomy_result(contig, locus_type):
    ssu, ssu_n = marker_consensus(contig, 'SSU')
    its, its_n = marker_consensus(contig, 'ITS')
    lsu, lsu_n = marker_consensus(contig, 'LSU')
    if locus_type == 'PARTIAL_18S':
        combined = ssu
    elif locus_type == 'PARTIAL_28S':
        combined = lsu
    else:
        combined = locus_consensus(ssu, lsu)
    ssu_phylum, lsu_phylum = ssu.get('phylum', 'NA'), lsu.get('phylum', 'NA')

    note = ''
    if locus_type in {'PARTIAL_18S', 'PARTIAL_28S'}:
        marker = ssu if locus_type == 'PARTIAL_18S' else lsu
        marker_name = 'SSU' if locus_type == 'PARTIAL_18S' else 'LSU'
        if not resolved(marker.get('phylum', 'NA')):
            status = 'PARTIAL_unresolved'
            note = f'{marker_name}-only locus was unresolved at phylum'
        elif expected:
            lineage_names = {value.casefold() for value in marker.values() if resolved(value)}
            if expected.casefold() in lineage_names:
                status = 'PARTIAL_expected_taxon'
            else:
                status = 'PARTIAL_non_target_taxon'
                note = f'{marker_name}-only lineage does not contain expected taxon {expected}'
        else:
            status = 'PARTIAL_taxonomically_assigned'
    elif (resolved(ssu_phylum) and resolved(lsu_phylum) and
            ssu_phylum.casefold() != lsu_phylum.casefold()):
        status = 'FAIL_taxonomic_chimera'
        note = f'SSU phylum {ssu_phylum} conflicts with LSU phylum {lsu_phylum}'
    elif not resolved(ssu_phylum) or not resolved(lsu_phylum):
        status = 'UNRESOLVED'
        note = 'SSU or LSU consensus was unresolved at phylum'
    elif expected:
        lineage_names = {value.casefold() for value in combined.values() if resolved(value)}
        if expected.casefold() in lineage_names:
            status = 'PASS_expected_taxon'
        else:
            status = 'PASS_non_target_taxon'
            note = f'Consensus lineage does not contain expected taxon {expected}'
    else:
        status = 'PASS_taxonomically_coherent'

    return {
        'taxonomy_status': status,
        'expected_taxonomy': expected or 'NA',
        **{f'consensus_{rank}': combined[rank] for rank in ranks},
        'consensus_taxonomy': lineage_text(combined),
        'SSU_consensus_taxonomy': lineage_text(ssu),
        'ITS_consensus_taxonomy': lineage_text(its),
        'LSU_consensus_taxonomy': lineage_text(lsu),
        'SSU_near_top_hits': str(ssu_n),
        'ITS_near_top_hits': str(its_n),
        'LSU_near_top_hits': str(lsu_n),
        'validation_note': note or 'NA'
    }

taxonomy_fields = ['contig','locus_type','taxonomy_status','expected_taxonomy',
                   *[f'consensus_{rank}' for rank in ranks],
                   'consensus_taxonomy','SSU_consensus_taxonomy',
                   'ITS_consensus_taxonomy','LSU_consensus_taxonomy',
                   'SSU_near_top_hits','ITS_near_top_hits','LSU_near_top_hits',
                   'validation_note']
taxonomy = {contig: taxonomy_result(contig, row.get('locus_type', 'NA'))
            for contig, row in catalog.items()}
with (root / 'taxonomy_validation.tsv').open('w', newline='') as handle:
    out = csv.DictWriter(handle, fieldnames=taxonomy_fields, delimiter='\t', lineterminator='\n')
    out.writeheader()
    for contig in sorted(taxonomy):
        out.writerow({'contig': contig,
                      'locus_type': catalog[contig].get('locus_type', 'NA'),
                      **taxonomy[contig]})

detailed_fields = ['contig','contig_length','locus_type','source','orientation','layout_status','structure_status',
          'SSU_coordinates','SSU_bp','ITS1_coordinates','ITS1_bp','5.8S_coordinates','5.8S_bp',
          'ITS2_coordinates','ITS2_bp','LSU_coordinates','LSU_bp',
          '18S_anchor_identity','28S_anchor_identity','readback_coverage_percent','readback_mean_depth',
          'readback_mean_mapq','unique_mapped_reads','unique_coverage_percent',
          'unique_mean_depth','graph_junctions','junctions_with_read_support',
          'minimum_junction_spanning_reads','taxonomy_status','expected_taxonomy',
          *[f'consensus_{rank}' for rank in ranks],
          'consensus_taxonomy','SSU_consensus_taxonomy','ITS_consensus_taxonomy',
          'LSU_consensus_taxonomy','SSU_near_top_hits','ITS_near_top_hits',
          'LSU_near_top_hits','validation_note']
for marker in ('SSU','ITS','LSU'):
    detailed_fields += [f'{marker}_top_accession',f'{marker}_top_scientific_name',
               f'{marker}_top_domain',f'{marker}_top_kingdom',f'{marker}_top_phylum',
               f'{marker}_top_class',f'{marker}_top_order',f'{marker}_top_family',
               f'{marker}_top_genus',f'{marker}_top_species',
               f'{marker}_top_identity',f'{marker}_top_aligned_bp',
               f'{marker}_top_query_coverage',f'{marker}_top_evalue',
               f'{marker}_top_bitscore',f'{marker}_top_title']

with (root / 'rrna_locus_summary.tsv').open('w', newline='') as handle:
    out = csv.DictWriter(handle, fieldnames=detailed_fields, delimiter='\t', lineterminator='\n')
    out.writeheader()
    support = {r['candidate']: r for r in rows(root / 'path_support.tsv')}
    for contig in sorted(catalog):
        x, l, c = itsx.get(contig, {}), layout.get(contig, {}), coverage.get(contig, {})
        candidate, p = catalog[contig], support.get(contig, {})
        row = {'contig':contig, 'contig_length':candidate.get('length_bp', x.get('length','NA')),
               'locus_type':candidate.get('locus_type','NA'),
               'source':candidate.get('source','NA'),
               'orientation':l.get('orientation','NA'), 'layout_status':l.get('layout_status','NA'),
               'structure_status':x.get('status',candidate.get('locus_type','NA')),
               '18S_anchor_identity':l.get('18S_identity','NA'),
               '28S_anchor_identity':l.get('28S_identity','NA'),
               'readback_coverage_percent':c.get('coverage','NA'),
               'readback_mean_depth':c.get('meandepth','NA'),
               'readback_mean_mapq':c.get('meanmapq','NA'),
               'unique_mapped_reads':p.get('unique_mapped_reads','NA'),
               'unique_coverage_percent':p.get('unique_coverage_percent','NA'),
               'unique_mean_depth':p.get('unique_mean_depth','NA'),
               'graph_junctions':p.get('graph_junctions','NA'),
               'junctions_with_read_support':p.get('junctions_with_read_support','NA'),
               'minimum_junction_spanning_reads':p.get('minimum_junction_spanning_reads','NA')}
        for marker in ('SSU','ITS1','5.8S','ITS2','LSU'):
            row[f'{marker}_coordinates'] = x.get(marker, 'NA')
            row[f'{marker}_bp'] = region_length(x.get(marker, ''))
        row.update(taxonomy[contig])
        for marker in ('SSU','ITS','LSU'):
            mapping = {'accession':'accession','scientific_name':'scientific_names',
                       'domain':'domain','kingdom':'kingdom','phylum':'phylum',
                       'class':'class','order':'order','family':'family',
                       'genus':'genus','species':'species',
                       'identity':'percent_identity','aligned_bp':'aligned_bp',
                       'query_coverage':'query_coverage_hsp','evalue':'evalue',
                       'bitscore':'bitscore','title':'title'}
            for suffix, source in mapping.items():
                row[f'{marker}_top_{suffix}'] = hit(contig, marker, source)
        out.writerow(row)

# Keep the master summaries intentionally concise. The complete audit table
# remains available as rrna_locus_summary.tsv.
master_fields = ['contig','locus_type','length_bp','regions','mean_depth','taxonomy_status',
                 'consensus_taxonomy','expected_taxonomy',
                 'SSU_top_hit','SSU_identity','SSU_aligned_bp',
                 'ITS_top_hit','ITS_identity','ITS_aligned_bp',
                 'LSU_top_hit','LSU_identity','LSU_aligned_bp']
master_rows = []
for contig in sorted(catalog):
    x, c, tax = itsx.get(contig, {}), coverage.get(contig, {}), taxonomy[contig]
    candidate = catalog[contig]
    if x:
        regions = '; '.join(f'{marker}={x.get(marker, "NA")}'
                            for marker in ('SSU','ITS1','5.8S','ITS2','LSU'))
    elif candidate.get('locus_type') == 'PARTIAL_18S':
        regions = 'SSU=partial'
    elif candidate.get('locus_type') == 'PARTIAL_28S':
        regions = 'LSU=partial'
    elif candidate.get('locus_type') == 'AMBIGUOUS_DUAL_ANCHOR':
        regions = 'SSU=present; LSU=present; layout=ambiguous'
    else:
        regions = 'rRNA=seed_anchored'
    header_depth = 'NA'
    match = re.search(r'_cov_([0-9.]+)', contig)
    if match:
        header_depth = match.group(1)
    master_rows.append({
        'contig': contig, 'locus_type':candidate.get('locus_type','NA'),
        'length_bp': candidate.get('length_bp',x.get('length','NA')), 'regions': regions,
        'mean_depth': c.get('meandepth',header_depth),
        'taxonomy_status': tax['taxonomy_status'],
        'consensus_taxonomy': tax['consensus_taxonomy'],
        'expected_taxonomy': tax['expected_taxonomy'],
        'SSU_top_hit': hit(contig,'SSU','scientific_names'),
        'SSU_identity': hit(contig,'SSU','percent_identity'),
        'SSU_aligned_bp': hit(contig,'SSU','aligned_bp'),
        'ITS_top_hit': hit(contig,'ITS','scientific_names'),
        'ITS_identity': hit(contig,'ITS','percent_identity'),
        'ITS_aligned_bp': hit(contig,'ITS','aligned_bp'),
        'LSU_top_hit': hit(contig,'LSU','scientific_names'),
        'LSU_identity': hit(contig,'LSU','percent_identity'),
        'LSU_aligned_bp': hit(contig,'LSU','aligned_bp')
    })

with (root / 'master_summary.tsv').open('w', newline='') as destination:
    writer = csv.DictWriter(destination, fieldnames=master_fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(master_rows)
with (root.parent / 'master_summary.csv').open('w', newline='') as destination:
    writer = csv.DictWriter(destination, fieldnames=master_fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows(master_rows)
PY
}

# Count raw input templates before any mapping or QC.
RAW_PAIR_COUNT=$(fastq_count "$READS_1")
RAW_REVERSE_COUNT=$(fastq_count "$READS_2")
RAW_SINGLE_COUNT=$(fastq_count "$READS_SINGLE")
(( RAW_PAIR_COUNT == RAW_REVERSE_COUNT )) || \
    die "Raw R1/R2 files contain different read counts."
TOTAL_TEMPLATES=$((RAW_PAIR_COUNT + RAW_SINGLE_COUNT))
(( TOTAL_TEMPLATES > 0 )) || die "Input FASTQs contain no reads."
log "Raw input: $RAW_PAIR_COUNT paired templates and $RAW_SINGLE_COUNT single reads."

combine_databases
log "Building Bowtie2 seed index."
bowtie2-build --threads "$THREADS" "$COMBINED_DB" "$SEED_INDEX" \
    > "$OUTDIR/seed/bowtie2-build.stdout.log" \
    2> "$OUTDIR/seed/bowtie2-build.stderr.log"

if (( ${#LEFT_DATABASES[@]} > 0 )); then
    bowtie2-build --threads "$THREADS" "$LEFT_DB" "$LEFT_INDEX" \
        > "$OUTDIR/seed/left-build.stdout.log" 2> "$OUTDIR/seed/left-build.stderr.log"
fi
if (( ${#RIGHT_DATABASES[@]} > 0 )); then
    bowtie2-build --threads "$THREADS" "$RIGHT_DB" "$RIGHT_INDEX" \
        > "$OUTDIR/seed/right-build.stdout.log" 2> "$OUTDIR/seed/right-build.stderr.log"
fi

log "Seed mapping: aligning the original raw reads once."
map_raw_reads_to_seeds
mapped_primary_names "$SEED_BAM" "$OUTDIR/seed/strict_seed_hit_names.txt"
cp -- "$OUTDIR/seed/strict_seed_hit_names.txt" "$ACCEPTED_NAMES"

ACCEPTED_COUNT=$(awk 'END { print NR + 0 }' "$ACCEPTED_NAMES")
(( ACCEPTED_COUNT > 0 )) || die "No reads passed the strict seed-alignment filters."
ACCEPTED_FRACTION=$(fraction "$ACCEPTED_COUNT" "$TOTAL_TEMPLATES")
log "Seed mapping accepted $ACCEPTED_COUNT templates ($ACCEPTED_FRACTION of input)."
if greater_than "$ACCEPTED_FRACTION" "$MAX_ACCEPTED_FRACTION"; then
    log "WARNING: seed hits already exceed --max-accepted-fraction; verify database specificity."
fi
if (( MAX_ACCEPTED_TEMPLATES > 0 && ACCEPTED_COUNT > MAX_ACCEPTED_TEMPLATES )); then
    log "WARNING: seed hits already exceed --max-accepted-templates; extension will remain guarded but seed hits are retained."
fi

printf 'round\tfrontier_templates\tcandidate_templates\tnew_templates\taccepted_total\tgrowth_fraction\taccepted_fraction\tdecision\n' \
    > "$METRICS"
printf '0\t%d\t%d\t%d\t%d\tNA\t%s\tseed\n' \
    "$ACCEPTED_COUNT" "$ACCEPTED_COUNT" "$ACCEPTED_COUNT" \
    "$ACCEPTED_COUNT" "$ACCEPTED_FRACTION" >> "$METRICS"
printf 'round\taccepted_templates\tanchored_components\tdual_anchor_components\tanchored_nodes\tanchored_edges\tdual_branch_nodes\tdual_cycle_rank\tmax_component_nodes\tdecision\n' \
    > "$GRAPH_METRICS"

# Extract and QC the initial full-seed accepted pool for final assembly.
SEED_RAW_R1="$OUTDIR/seed/accepted_raw_R1.fastq.gz"
SEED_RAW_R2="$OUTDIR/seed/accepted_raw_R2.fastq.gz"
SEED_RAW_SINGLE="$OUTDIR/seed/accepted_raw_single.fastq.gz"
SEED_CLEAN_R1="$OUTDIR/seed/accepted_R1.fastq.gz"
SEED_CLEAN_R2="$OUTDIR/seed/accepted_R2.fastq.gz"
SEED_CLEAN_SINGLE="$OUTDIR/seed/accepted_single.fastq.gz"
extract_templates_from_bam "$ACCEPTED_NAMES" \
    "$SEED_RAW_R1" "$SEED_RAW_R2" "$SEED_RAW_SINGLE" "$OUTDIR/seed/extract"
log "Post-map QC: processing the initial recruited subset."
prepare_recruited_batch "$SEED_RAW_R1" "$SEED_RAW_R2" "$SEED_RAW_SINGLE" \
    "$OUTDIR/seed" "$SEED_CLEAN_R1" "$SEED_CLEAN_R2" "$SEED_CLEAN_SINGLE"
ASSEMBLY_R1_FILES+=("$SEED_CLEAN_R1")
ASSEMBLY_R2_FILES+=("$SEED_CLEAN_R2")
ASSEMBLY_SINGLE_FILES+=("$SEED_CLEAN_SINGLE")

STOP_REASON="extension disabled (--max-rounds 0)"
ROUNDS_ACCEPTED=0
INWARD_ROUNDS_ACCEPTED=0
OUTWARD_ACCEPTED=false
OUTWARD_REASON="not requested"
HAS_DIRECTIONAL_DBS=false
if (( ${#LEFT_DATABASES[@]} > 0 && ${#RIGHT_DATABASES[@]} > 0 )); then
    HAS_DIRECTIONAL_DBS=true
fi

if (( MAX_ROUNDS > 0 )); then
    RUN_RECURSIVE=false
    if [[ "$INWARD" == true || "$HAS_DIRECTIONAL_DBS" != true ]]; then
        RUN_RECURSIVE=true
    fi
    STOP_REASON="no recursive inward extension requested"

    if [[ "$RUN_RECURSIVE" == true ]]; then
        FRONTIER_NAMES="$OUTDIR/recruitment/round_00_frontier_names.txt"
        FRONTIER_R1="$SEED_CLEAN_R1"
        FRONTIER_R2="$SEED_CLEAN_R2"
        FRONTIER_SINGLE="$SEED_CLEAN_SINGLE"

        if [[ "$HAS_DIRECTIONAL_DBS" == true ]]; then
            TERMINAL_SEEDS="$OUTDIR/seed/inward_terminal_seeds.fasta"
            TERMINAL_CANDIDATES="$OUTDIR/seed/inward_terminal_candidate_names.txt"
            make_directional_seed_baits "$TERMINAL_SEEDS" "inward"
            log "Inward mode: identifying templates at the 18S 3' and 28S 5' seed ends."
            scan_raw_with_baits "$TERMINAL_SEEDS" "$TERMINAL_CANDIDATES" \
                "$OUTDIR/seed/inward_terminal_scan"
            LC_ALL=C comm -12 "$ACCEPTED_NAMES" "$TERMINAL_CANDIDATES" > "$FRONTIER_NAMES"
            FRONTIER_COUNT=$(awk 'END { print NR + 0 }' "$FRONTIER_NAMES")
            (( FRONTIER_COUNT > 0 )) || \
                die "Inward mode found no terminal seed templates; increase --inward-end-length or use --outward."

            FRONTIER_RAW_R1="$OUTDIR/seed/frontier_raw_R1.fastq.gz"
            FRONTIER_RAW_R2="$OUTDIR/seed/frontier_raw_R2.fastq.gz"
            FRONTIER_RAW_SINGLE="$OUTDIR/seed/frontier_raw_single.fastq.gz"
            FRONTIER_R1="$OUTDIR/seed/frontier_R1.fastq.gz"
            FRONTIER_R2="$OUTDIR/seed/frontier_R2.fastq.gz"
            FRONTIER_SINGLE="$OUTDIR/seed/frontier_single.fastq.gz"
            extract_templates_from_bam "$FRONTIER_NAMES" \
                "$FRONTIER_RAW_R1" "$FRONTIER_RAW_R2" "$FRONTIER_RAW_SINGLE" \
                "$OUTDIR/seed/frontier_extract"
            prepare_recruited_batch "$FRONTIER_RAW_R1" "$FRONTIER_RAW_R2" "$FRONTIER_RAW_SINGLE" \
                "$OUTDIR/seed/frontier_qc" "$FRONTIER_R1" "$FRONTIER_R2" "$FRONTIER_SINGLE"
        else
            cp -- "$ACCEPTED_NAMES" "$FRONTIER_NAMES"
            FRONTIER_COUNT="$ACCEPTED_COUNT"
        fi

        if [[ "$KEEP_INTERMEDIATES" != true && "$HAS_DIRECTIONAL_DBS" == true ]]; then
            rm -f -- "$FRONTIER_RAW_R1" "$FRONTIER_RAW_R2" "$FRONTIER_RAW_SINGLE"
        fi

        FRONTIER_BAITS="$OUTDIR/recruitment/round_00_frontier_baits.fasta"
        make_frontier_baits "$FRONTIER_BAITS" \
            "$OUTDIR/recruitment/round_00_bait_filter.log" \
            "$FRONTIER_R1" "$FRONTIER_R2" "$FRONTIER_SINGLE"
        STOP_REASON="maximum inward recruitment rounds reached"

        for (( round = 1; round <= MAX_ROUNDS; round++ )); do
            printf -v round_name 'round_%02d' "$round"
            ROUND_DIR="$OUTDIR/recruitment/$round_name"
            mkdir -p "$ROUND_DIR"
            CANDIDATE_NAMES="$ROUND_DIR/candidate_names.txt"
            NEW_NAMES="$ROUND_DIR/new_names.txt"

            FRONTIER_COUNT=$(awk 'END { print NR + 0 }' "$FRONTIER_NAMES")
            log "Inward recruitment round $round: scanning raw reads from a frontier of $FRONTIER_COUNT templates."
            scan_raw_with_baits "$FRONTIER_BAITS" "$CANDIDATE_NAMES" "$ROUND_DIR/bbduk"
            LC_ALL=C comm -23 "$CANDIDATE_NAMES" "$ACCEPTED_NAMES" > "$NEW_NAMES"

            CANDIDATE_COUNT=$(awk 'END { print NR + 0 }' "$CANDIDATE_NAMES")
            NEW_COUNT=$(awk 'END { print NR + 0 }' "$NEW_NAMES")
            GROWTH=$(fraction "$NEW_COUNT" "$ACCEPTED_COUNT")
            PROPOSED_COUNT=$((ACCEPTED_COUNT + NEW_COUNT))
            PROPOSED_FRACTION=$(fraction "$PROPOSED_COUNT" "$TOTAL_TEMPLATES")

            if greater_than "$GROWTH" "$MAX_ROUND_GROWTH"; then
                printf '%d\t%d\t%d\t%d\t%d\t%s\t%s\trejected_growth\n' \
                    "$round" "$FRONTIER_COUNT" "$CANDIDATE_COUNT" "$NEW_COUNT" \
                    "$ACCEPTED_COUNT" "$GROWTH" "$ACCEPTED_FRACTION" >> "$METRICS"
                STOP_REASON="inward round-$round growth $GROWTH exceeded $MAX_ROUND_GROWTH; recruitment rejected"
                log "WARNING: $STOP_REASON."
                break
            fi
            if greater_than "$PROPOSED_FRACTION" "$MAX_ACCEPTED_FRACTION"; then
                printf '%d\t%d\t%d\t%d\t%d\t%s\t%s\trejected_total_fraction\n' \
                    "$round" "$FRONTIER_COUNT" "$CANDIDATE_COUNT" "$NEW_COUNT" \
                    "$ACCEPTED_COUNT" "$GROWTH" "$ACCEPTED_FRACTION" >> "$METRICS"
                STOP_REASON="inward round-$round accepted fraction $PROPOSED_FRACTION exceeded $MAX_ACCEPTED_FRACTION; recruitment rejected"
                log "WARNING: $STOP_REASON."
                break
            fi
            if (( MAX_ACCEPTED_TEMPLATES > 0 && PROPOSED_COUNT > MAX_ACCEPTED_TEMPLATES )); then
                printf '%d\t%d\t%d\t%d\t%d\t%s\t%s\trejected_total_templates\n' \
                    "$round" "$FRONTIER_COUNT" "$CANDIDATE_COUNT" "$NEW_COUNT" \
                    "$ACCEPTED_COUNT" "$GROWTH" "$ACCEPTED_FRACTION" >> "$METRICS"
                STOP_REASON="inward round-$round proposed total $PROPOSED_COUNT exceeded $MAX_ACCEPTED_TEMPLATES templates; recruitment rejected"
                log "WARNING: $STOP_REASON."
                break
            fi
            if (( NEW_COUNT == 0 )); then
                printf '%d\t%d\t%d\t0\t%d\t0.00000000\t%s\tconverged\n' \
                    "$round" "$FRONTIER_COUNT" "$CANDIDATE_COUNT" \
                    "$ACCEPTED_COUNT" "$ACCEPTED_FRACTION" >> "$METRICS"
                STOP_REASON="no new inward templates were recruited"
                break
            fi

            LC_ALL=C sort -u "$ACCEPTED_NAMES" "$NEW_NAMES" > "$ROUND_DIR/accepted_names.next.txt"
            mv -- "$ROUND_DIR/accepted_names.next.txt" "$ACCEPTED_NAMES"
            ACCEPTED_COUNT="$PROPOSED_COUNT"
            ACCEPTED_FRACTION="$PROPOSED_FRACTION"
            ROUNDS_ACCEPTED=$((ROUNDS_ACCEPTED + 1))
            INWARD_ROUNDS_ACCEPTED=$((INWARD_ROUNDS_ACCEPTED + 1))

            RAW_R1="$ROUND_DIR/new_raw_R1.fastq.gz"
            RAW_R2="$ROUND_DIR/new_raw_R2.fastq.gz"
            RAW_SINGLE="$ROUND_DIR/new_raw_single.fastq.gz"
            CLEAN_R1="$ROUND_DIR/new_R1.fastq.gz"
            CLEAN_R2="$ROUND_DIR/new_R2.fastq.gz"
            CLEAN_SINGLE="$ROUND_DIR/new_single.fastq.gz"
            extract_templates_from_bam "$NEW_NAMES" "$RAW_R1" "$RAW_R2" "$RAW_SINGLE" \
                "$ROUND_DIR/extract"
            log "Inward recruitment round $round: QC of $NEW_COUNT newly accepted templates."
            prepare_recruited_batch "$RAW_R1" "$RAW_R2" "$RAW_SINGLE" \
                "$ROUND_DIR" "$CLEAN_R1" "$CLEAN_R2" "$CLEAN_SINGLE"
            ASSEMBLY_R1_FILES+=("$CLEAN_R1")
            ASSEMBLY_R2_FILES+=("$CLEAN_R2")
            ASSEMBLY_SINGLE_FILES+=("$CLEAN_SINGLE")

            printf '%d\t%d\t%d\t%d\t%d\t%s\t%s\taccepted_inward\n' \
                "$round" "$FRONTIER_COUNT" "$CANDIDATE_COUNT" "$NEW_COUNT" \
                "$ACCEPTED_COUNT" "$GROWTH" "$ACCEPTED_FRACTION" >> "$METRICS"
            log "Inward recruitment round $round: accepted $NEW_COUNT new templates; total=$ACCEPTED_COUNT, growth=$GROWTH."
            [[ "$KEEP_INTERMEDIATES" == true ]] || rm -f -- "$RAW_R1" "$RAW_R2" "$RAW_SINGLE"

            if [[ "$ADAPTIVE_GRAPH_STOP" == true && "$HAS_DIRECTIONAL_DBS" == true ]]; then
                log "Inward recruitment round $round: checking anchored assembly-graph complexity (k=$GRAPH_CHECK_K)."
                run_graph_checkpoint "$round" "$ROUND_DIR"
                if [[ "$CHECKPOINT_DECISION" == "stop" ]]; then
                    STOP_REASON="adaptive graph stop after inward round-$round: $CHECKPOINT_STATUS (dual_components=$CHECKPOINT_DUAL_COMPONENTS, branches=$CHECKPOINT_BRANCHES, cycles=$CHECKPOINT_CYCLES, max_nodes=$CHECKPOINT_MAX_NODES)"
                    log "$STOP_REASON."
                    break
                fi
            fi
            if (( NEW_COUNT < MIN_NEW_TEMPLATES )); then
                STOP_REASON="new inward-template count $NEW_COUNT fell below $MIN_NEW_TEMPLATES"
                break
            fi
            if less_than "$GROWTH" "$MIN_GROWTH"; then
                STOP_REASON="inward growth $GROWTH fell below $MIN_GROWTH"
                break
            fi

            FRONTIER_NAMES="$NEW_NAMES"
            FRONTIER_BAITS="$ROUND_DIR/frontier_baits.fasta"
            make_frontier_baits "$FRONTIER_BAITS" "$ROUND_DIR/bait_filter.log" \
                "$CLEAN_R1" "$CLEAN_R2" "$CLEAN_SINGLE"
        done
    fi

    # Outward growth is intentionally separate from the recursive inward
    # frontier. It is one-pass and can never seed another outward scan.
    if [[ "$OUTWARD" == true && "$HAS_DIRECTIONAL_DBS" == true ]]; then
        OUTWARD_DIR="$OUTDIR/recruitment/outward_once"
        mkdir -p "$OUTWARD_DIR"
        OUTWARD_SEEDS="$OUTWARD_DIR/outward_terminal_seeds.fasta"
        OUTWARD_TERMINAL_CANDIDATES="$OUTWARD_DIR/terminal_candidate_names.txt"
        OUTWARD_FRONTIER_NAMES="$OUTWARD_DIR/frontier_names.txt"
        OUTWARD_FRONTIER_BAITS="$OUTWARD_DIR/frontier_baits.fasta"
        OUTWARD_CANDIDATES="$OUTWARD_DIR/candidate_names.txt"
        OUTWARD_NEW_NAMES="$OUTWARD_DIR/new_names.txt"
        make_directional_seed_baits "$OUTWARD_SEEDS" "outward"
        log "Conservative outward extension: locating seed-hit templates at the 18S 5' and 28S 3' ends."
        scan_raw_with_baits "$OUTWARD_SEEDS" "$OUTWARD_TERMINAL_CANDIDATES" \
            "$OUTWARD_DIR/terminal_scan"
        LC_ALL=C comm -12 "$ACCEPTED_NAMES" "$OUTWARD_TERMINAL_CANDIDATES" \
            > "$OUTWARD_FRONTIER_NAMES"
        OUTWARD_FRONTIER_COUNT=$(awk 'END { print NR + 0 }' "$OUTWARD_FRONTIER_NAMES")
        if (( OUTWARD_FRONTIER_COUNT == 0 )); then
            OUTWARD_REASON="no outward terminal frontier"
            printf 'outward\t0\t0\t0\t%d\t0.00000000\t%s\tno_outward_frontier\n' \
                "$ACCEPTED_COUNT" "$ACCEPTED_FRACTION" >> "$METRICS"
            log "WARNING: conservative outward extension found no terminal frontier and was skipped."
            STOP_REASON="$STOP_REASON; outward: $OUTWARD_REASON"
        else
        OUTWARD_FRONTIER_RAW_R1="$OUTWARD_DIR/frontier_raw_R1.fastq.gz"
        OUTWARD_FRONTIER_RAW_R2="$OUTWARD_DIR/frontier_raw_R2.fastq.gz"
        OUTWARD_FRONTIER_RAW_SINGLE="$OUTWARD_DIR/frontier_raw_single.fastq.gz"
        OUTWARD_FRONTIER_R1="$OUTWARD_DIR/frontier_R1.fastq.gz"
        OUTWARD_FRONTIER_R2="$OUTWARD_DIR/frontier_R2.fastq.gz"
        OUTWARD_FRONTIER_SINGLE="$OUTWARD_DIR/frontier_single.fastq.gz"
        extract_templates_from_bam "$OUTWARD_FRONTIER_NAMES" \
            "$OUTWARD_FRONTIER_RAW_R1" "$OUTWARD_FRONTIER_RAW_R2" \
            "$OUTWARD_FRONTIER_RAW_SINGLE" "$OUTWARD_DIR/frontier_extract"
        prepare_recruited_batch "$OUTWARD_FRONTIER_RAW_R1" "$OUTWARD_FRONTIER_RAW_R2" \
            "$OUTWARD_FRONTIER_RAW_SINGLE" "$OUTWARD_DIR/frontier_qc" \
            "$OUTWARD_FRONTIER_R1" "$OUTWARD_FRONTIER_R2" "$OUTWARD_FRONTIER_SINGLE"
        make_frontier_baits "$OUTWARD_FRONTIER_BAITS" "$OUTWARD_DIR/frontier_bait_filter.log" \
            "$OUTWARD_FRONTIER_R1" "$OUTWARD_FRONTIER_R2" "$OUTWARD_FRONTIER_SINGLE"
        log "Conservative outward extension: one nonrecursive frontier scan requiring $OUTWARD_MIN_WORD_HITS exact word hits."
        scan_raw_with_baits "$OUTWARD_FRONTIER_BAITS" "$OUTWARD_CANDIDATES" \
            "$OUTWARD_DIR/bbduk" "$WORD_SIZE" "$OUTWARD_MIN_WORD_HITS"
        LC_ALL=C comm -23 "$OUTWARD_CANDIDATES" "$ACCEPTED_NAMES" > "$OUTWARD_NEW_NAMES"
        OUTWARD_CANDIDATE_COUNT=$(awk 'END { print NR + 0 }' "$OUTWARD_CANDIDATES")
        OUTWARD_NEW_COUNT=$(awk 'END { print NR + 0 }' "$OUTWARD_NEW_NAMES")
        OUTWARD_GROWTH=$(fraction "$OUTWARD_NEW_COUNT" "$ACCEPTED_COUNT")
        OUTWARD_PROPOSED_COUNT=$((ACCEPTED_COUNT + OUTWARD_NEW_COUNT))
        OUTWARD_PROPOSED_FRACTION=$(fraction "$OUTWARD_PROPOSED_COUNT" "$TOTAL_TEMPLATES")
        OUTWARD_REASON="accepted"
        if (( OUTWARD_NEW_COUNT == 0 )); then
            OUTWARD_REASON="no new outward templates"
        elif (( MAX_OUTWARD_TEMPLATES > 0 && OUTWARD_NEW_COUNT > MAX_OUTWARD_TEMPLATES )); then
            OUTWARD_REASON="new outward templates $OUTWARD_NEW_COUNT exceeded $MAX_OUTWARD_TEMPLATES"
        elif greater_than "$OUTWARD_GROWTH" "$MAX_OUTWARD_GROWTH"; then
            OUTWARD_REASON="outward growth $OUTWARD_GROWTH exceeded $MAX_OUTWARD_GROWTH"
        elif greater_than "$OUTWARD_PROPOSED_FRACTION" "$MAX_ACCEPTED_FRACTION"; then
            OUTWARD_REASON="outward proposed fraction $OUTWARD_PROPOSED_FRACTION exceeded $MAX_ACCEPTED_FRACTION"
        elif (( MAX_ACCEPTED_TEMPLATES > 0 && OUTWARD_PROPOSED_COUNT > MAX_ACCEPTED_TEMPLATES )); then
            OUTWARD_REASON="outward proposed total $OUTWARD_PROPOSED_COUNT exceeded $MAX_ACCEPTED_TEMPLATES"
        fi

        if [[ "$OUTWARD_REASON" == "no new outward templates" ]]; then
            printf 'outward\t%d\t%d\t0\t%d\t0.00000000\t%s\tconverged_outward\n' \
                "$OUTWARD_FRONTIER_COUNT" "$OUTWARD_CANDIDATE_COUNT" \
                "$ACCEPTED_COUNT" "$ACCEPTED_FRACTION" >> "$METRICS"
            log "Conservative outward extension converged without new templates."
        elif [[ "$OUTWARD_REASON" == "accepted" ]]; then
            LC_ALL=C sort -u "$ACCEPTED_NAMES" "$OUTWARD_NEW_NAMES" > "$OUTWARD_DIR/accepted_names.next.txt"
            mv -- "$OUTWARD_DIR/accepted_names.next.txt" "$ACCEPTED_NAMES"
            ACCEPTED_COUNT="$OUTWARD_PROPOSED_COUNT"
            ACCEPTED_FRACTION="$OUTWARD_PROPOSED_FRACTION"
            ROUNDS_ACCEPTED=$((ROUNDS_ACCEPTED + 1))
            OUTWARD_ACCEPTED=true
            OUTWARD_RAW_R1="$OUTWARD_DIR/new_raw_R1.fastq.gz"
            OUTWARD_RAW_R2="$OUTWARD_DIR/new_raw_R2.fastq.gz"
            OUTWARD_RAW_SINGLE="$OUTWARD_DIR/new_raw_single.fastq.gz"
            OUTWARD_CLEAN_R1="$OUTWARD_DIR/new_R1.fastq.gz"
            OUTWARD_CLEAN_R2="$OUTWARD_DIR/new_R2.fastq.gz"
            OUTWARD_CLEAN_SINGLE="$OUTWARD_DIR/new_single.fastq.gz"
            extract_templates_from_bam "$OUTWARD_NEW_NAMES" \
                "$OUTWARD_RAW_R1" "$OUTWARD_RAW_R2" "$OUTWARD_RAW_SINGLE" \
                "$OUTWARD_DIR/extract"
            log "Conservative outward extension: QC of $OUTWARD_NEW_COUNT newly accepted templates."
            prepare_recruited_batch "$OUTWARD_RAW_R1" "$OUTWARD_RAW_R2" "$OUTWARD_RAW_SINGLE" \
                "$OUTWARD_DIR" "$OUTWARD_CLEAN_R1" "$OUTWARD_CLEAN_R2" "$OUTWARD_CLEAN_SINGLE"
            ASSEMBLY_R1_FILES+=("$OUTWARD_CLEAN_R1")
            ASSEMBLY_R2_FILES+=("$OUTWARD_CLEAN_R2")
            ASSEMBLY_SINGLE_FILES+=("$OUTWARD_CLEAN_SINGLE")
            printf 'outward\t%d\t%d\t%d\t%d\t%s\t%s\taccepted_outward_once\n' \
                "$OUTWARD_FRONTIER_COUNT" "$OUTWARD_CANDIDATE_COUNT" "$OUTWARD_NEW_COUNT" "$ACCEPTED_COUNT" \
                "$OUTWARD_GROWTH" "$ACCEPTED_FRACTION" >> "$METRICS"
            log "Conservative outward extension accepted $OUTWARD_NEW_COUNT templates; total=$ACCEPTED_COUNT."
            [[ "$KEEP_INTERMEDIATES" == true ]] || \
                rm -f -- "$OUTWARD_RAW_R1" "$OUTWARD_RAW_R2" "$OUTWARD_RAW_SINGLE"
        else
            printf 'outward\t%d\t%d\t%d\t%d\t%s\t%s\trejected_outward_guard\n' \
                "$OUTWARD_FRONTIER_COUNT" "$OUTWARD_CANDIDATE_COUNT" "$OUTWARD_NEW_COUNT" "$ACCEPTED_COUNT" \
                "$OUTWARD_GROWTH" "$ACCEPTED_FRACTION" >> "$METRICS"
            log "WARNING: conservative outward extension rejected: $OUTWARD_REASON."
        fi
        if [[ "$KEEP_INTERMEDIATES" != true ]]; then
            rm -f -- "$OUTWARD_FRONTIER_RAW_R1" "$OUTWARD_FRONTIER_RAW_R2" \
                "$OUTWARD_FRONTIER_RAW_SINGLE" "$OUTWARD_FRONTIER_R1" \
                "$OUTWARD_FRONTIER_R2" "$OUTWARD_FRONTIER_SINGLE"
        fi
        STOP_REASON="$STOP_REASON; outward: $OUTWARD_REASON"
        fi
    fi

    [[ "$KEEP_INTERMEDIATES" == true ]] || \
        rm -f -- "$SEED_RAW_R1" "$SEED_RAW_R2" "$SEED_RAW_SINGLE"
else
    log "Word extension disabled; proceeding directly from seed recruitment to assembly."
    [[ "$KEEP_INTERMEDIATES" == true ]] || \
        rm -f -- "$SEED_RAW_R1" "$SEED_RAW_R2" "$SEED_RAW_SINGLE"
fi

log "Recruitment stopped: $STOP_REASON."

FINAL_R1="$OUTDIR/final/accepted_R1.fastq.gz"
FINAL_R2="$OUTDIR/final/accepted_R2.fastq.gz"
FINAL_SINGLE="$OUTDIR/final/accepted_single.fastq.gz"
combine_fastqs "$FINAL_R1" "${ASSEMBLY_R1_FILES[@]}"
combine_fastqs "$FINAL_R2" "${ASSEMBLY_R2_FILES[@]}"
combine_fastqs "$FINAL_SINGLE" "${ASSEMBLY_SINGLE_FILES[@]}"

log "Final assembly: running SPAdes once with all accepted, QC-filtered reads."
assemble_final_reads "$FINAL_R1" "$FINAL_R2" "$FINAL_SINGLE"

ALL_CONTIGS="$OUTDIR/final/all_assembled_contigs.fasta"
CANDIDATE_CONTIGS="$OUTDIR/final/rrna_candidate_contigs.fasta"
DUAL_CONTIGS="$OUTDIR/final/rrna_dual_anchor_contigs.fasta"
cp -- "$OUTDIR/assembly/spades/contigs.fasta" "$ALL_CONTIGS"

contig_anchor_hits "$ALL_CONTIGS" "$SEED_INDEX" \
    "$OUTDIR/final/all_seed_anchor_hits.tsv" "$OUTDIR/final/all_seed_anchor.log"
awk '{ print $1 }' "$OUTDIR/final/all_seed_anchor_hits.tsv" | LC_ALL=C sort -u \
    > "$OUTDIR/final/rrna_candidate_contig_names.txt"
[[ -s "$OUTDIR/final/rrna_candidate_contig_names.txt" ]] || \
    die "No final contig retained a direct seed hit."
extract_fasta_by_names "$OUTDIR/final/rrna_candidate_contig_names.txt" \
    "$ALL_CONTIGS" "$CANDIDATE_CONTIGS"

: > "$DUAL_CONTIGS"
DUAL_COUNT=0
if (( ${#LEFT_DATABASES[@]} > 0 && ${#RIGHT_DATABASES[@]} > 0 )); then
    log "Variant-aware graph analysis: enumerating bounded paths between compatible anchors."
    enumerate_graph_paths
    merge_locus_candidates
    ALL_LOCUS_CANDIDATES="$OUTDIR/final/all_locus_candidates.fasta"

    contig_anchor_hits "$ALL_LOCUS_CANDIDATES" "$LEFT_INDEX" \
        "$OUTDIR/final/left_anchor_hits.tsv" "$OUTDIR/final/left_anchor.log"
    contig_anchor_hits "$ALL_LOCUS_CANDIDATES" "$RIGHT_INDEX" \
        "$OUTDIR/final/right_anchor_hits.tsv" "$OUTDIR/final/right_anchor.log"
    awk '{ print $1 }' "$OUTDIR/final/left_anchor_hits.tsv" | LC_ALL=C sort -u \
        > "$OUTDIR/final/left_anchor_contig_names.txt"
    awk '{ print $1 }' "$OUTDIR/final/right_anchor_hits.tsv" | LC_ALL=C sort -u \
        > "$OUTDIR/final/right_anchor_contig_names.txt"
    LC_ALL=C comm -12 "$OUTDIR/final/left_anchor_contig_names.txt" \
        "$OUTDIR/final/right_anchor_contig_names.txt" \
        > "$OUTDIR/final/dual_anchor_contig_names.txt"
    extract_fasta_by_names "$OUTDIR/final/dual_anchor_contig_names.txt" \
        "$ALL_LOCUS_CANDIDATES" "$DUAL_CONTIGS"
    DUAL_COUNT=$(awk '/^>/ { n++ } END { print n + 0 }' "$DUAL_CONTIGS")

    log "Validating dual-anchor layout and orienting accepted contigs 18S-to-28S."
    validate_and_orient_dual_anchors "$ALL_LOCUS_CANDIDATES"
    ORIENTED_COUNT=$(awk '/^>/ { n++ } END { print n + 0 }' \
        "$OUTDIR/final/oriented_dual_anchor_contigs.fasta")
    AMBIGUOUS_COUNT=$(awk '/^>/ { n++ } END { print n + 0 }' \
        "$OUTDIR/final/ambiguous_dual_anchor_contigs.fasta")
    log "Layout validation retained $ORIENTED_COUNT oriented dual-anchor contig(s); $AMBIGUOUS_COUNT ambiguous."
    cp -- "$OUTDIR/final/oriented_dual_anchor_contigs.fasta" "$DUAL_CONTIGS"
    awk '/^>/ { name=substr($0,2); sub(/[[:space:]].*$/, "", name); print name }' \
        "$DUAL_CONTIGS" > "$OUTDIR/final/dual_anchor_contig_names.txt"
    DUAL_COUNT="$ORIENTED_COUNT"
    catalog_locus_candidates
    run_itsx_validation
    log "Competitive read-back validation against all oriented candidates."
    run_competitive_readback
    run_ncbi_classification
    build_locus_summary
else
    ORIENTED_COUNT=0
    AMBIGUOUS_COUNT=0
    : > "$OUTDIR/final/oriented_dual_anchor_contigs.fasta"
    : > "$OUTDIR/final/ambiguous_dual_anchor_contigs.fasta"
    : > "$OUTDIR/final/graph_candidate_paths.fasta"
    : > "$OUTDIR/final/graph_path_summary.tsv"
    cp -- "$CANDIDATE_CONTIGS" "$OUTDIR/final/all_locus_candidates.fasta"
    : > "$OUTDIR/final/partial_18S_contigs.fasta"
    : > "$OUTDIR/final/partial_28S_contigs.fasta"
    : > "$OUTDIR/final/partial_locus_contigs.fasta"
    : > "$OUTDIR/final/ambiguous_locus_candidates.fasta"
    : > "$OUTDIR/final/other_seed_anchored_contigs.fasta"
    : > "$OUTDIR/final/locus_catalog.tsv"
    : > "$OUTDIR/final/contig_validation.tsv"
    : > "$OUTDIR/final/complete_ITS.fasta"
    : > "$OUTDIR/final/ITS1.fasta"
    : > "$OUTDIR/final/5_8S.fasta"
    : > "$OUTDIR/final/ITS2.fasta"
    : > "$OUTDIR/final/readback_coverage.tsv"
    : > "$OUTDIR/final/path_support.tsv"
    : > "$OUTDIR/final/residual_variants.tsv"
    : > "$OUTDIR/final/ncbi_SSU_hits.tsv"
    : > "$OUTDIR/final/ncbi_ITS_hits.tsv"
    : > "$OUTDIR/final/ncbi_LSU_hits.tsv"
    : > "$OUTDIR/final/ncbi_blast_top_hits.tsv"
    : > "$OUTDIR/final/taxonomy_validation.tsv"
    : > "$OUTDIR/final/rrna_locus_summary.tsv"
    : > "$OUTDIR/final/master_summary.tsv"
    : > "$OUTDIR/master_summary.csv"
fi

for graph in assembly_graph.fastg assembly_graph_with_scaffolds.gfa assembly_graph.gfa; do
    if [[ -s "$OUTDIR/assembly/spades/$graph" ]]; then
        cp -- "$OUTDIR/assembly/spades/$graph" "$OUTDIR/final/$graph"
    fi
done

read -r ALL_CONTIG_COUNT ALL_BP < <(fasta_stats "$ALL_CONTIGS")
read -r CANDIDATE_COUNT CANDIDATE_BP < <(fasta_stats "$CANDIDATE_CONTIGS")
read -r COMPLETE_ITS_COUNT COMPLETE_ITS_BP < <(fasta_stats "$OUTDIR/final/complete_ITS.fasta")
read -r GRAPH_PATH_COUNT GRAPH_PATH_BP < <(fasta_stats "$OUTDIR/final/graph_candidate_paths.fasta")
read -r PARTIAL_LOCUS_COUNT PARTIAL_LOCUS_BP < <(fasta_stats "$OUTDIR/final/partial_locus_contigs.fasta")
VARIANT_COUNT=$(awk 'END { print NR > 0 ? NR - 1 : 0 }' "$OUTDIR/final/residual_variants.tsv")
EXPECTED_PASS_COUNT=$(awk -F '\t' 'NR > 1 && ($3 == "PASS_expected_taxon" || $3 == "PARTIAL_expected_taxon") { n++ } END { print n + 0 }' \
    "$OUTDIR/final/taxonomy_validation.tsv")
NON_TARGET_COUNT=$(awk -F '\t' 'NR > 1 && ($3 == "PASS_non_target_taxon" || $3 == "PARTIAL_non_target_taxon") { n++ } END { print n + 0 }' \
    "$OUTDIR/final/taxonomy_validation.tsv")
COHERENT_TAXONOMY_COUNT=$(awk -F '\t' 'NR > 1 && ($3 == "PASS_taxonomically_coherent" || $3 == "PARTIAL_taxonomically_assigned") { n++ } END { print n + 0 }' \
    "$OUTDIR/final/taxonomy_validation.tsv")
TAXONOMIC_CHIMERA_COUNT=$(awk -F '\t' 'NR > 1 && $3 == "FAIL_taxonomic_chimera" { n++ } END { print n + 0 }' \
    "$OUTDIR/final/taxonomy_validation.tsv")
UNRESOLVED_TAXONOMY_COUNT=$(awk -F '\t' 'NR > 1 && ($3 == "UNRESOLVED" || $3 == "PARTIAL_unresolved") { n++ } END { print n + 0 }' \
    "$OUTDIR/final/taxonomy_validation.tsv")
FINAL_PAIR_COUNT=$(fastq_count "$FINAL_R1")
FINAL_SINGLE_COUNT=$(fastq_count "$FINAL_SINGLE")
ELAPSED=$(format_duration "$SECONDS")

printf '%s\n' \
    "ITSME version: $VERSION" \
    "Raw paired templates: $RAW_PAIR_COUNT" \
    "Raw single reads: $RAW_SINGLE_COUNT" \
    "Seed-accepted templates: $(awk 'NR == 2 { print $5 }' "$METRICS")" \
    "Final accepted template names: $ACCEPTED_COUNT" \
    "Final accepted fraction: $ACCEPTED_FRACTION" \
    "Accepted post-QC paired templates: $FINAL_PAIR_COUNT" \
    "Accepted post-QC single reads: $FINAL_SINGLE_COUNT" \
    "Accepted word-extension rounds: $ROUNDS_ACCEPTED" \
    "Accepted inward rounds: $INWARD_ROUNDS_ACCEPTED" \
    "Outward pass accepted: $OUTWARD_ACCEPTED" \
    "Outward pass result: $OUTWARD_REASON" \
    "Stopping reason: $STOP_REASON" \
    "Maximum accepted templates: $MAX_ACCEPTED_TEMPLATES" \
    "Maximum accepted fraction: $MAX_ACCEPTED_FRACTION" \
    "rRNA database directory: ${DB_DIR:-explicit database files}" \
    "Extension direction: $EXTENSION_DIRECTION" \
    "Inward mode: $INWARD" \
    "Outward mode: $OUTWARD" \
    "Adaptive graph stop: $ADAPTIVE_GRAPH_STOP" \
    "Checkpoint graph k-mer: $GRAPH_CHECK_K" \
    "Checkpoint branch threshold: $GRAPH_STOP_BRANCHES" \
    "Checkpoint cycle threshold: $GRAPH_STOP_CYCLES" \
    "Checkpoint component-node ceiling: $GRAPH_STOP_NODES" \
    "Outward minimum word hits: $OUTWARD_MIN_WORD_HITS" \
    "Outward new-template ceiling: $MAX_OUTWARD_TEMPLATES" \
    "Word size: $WORD_SIZE" \
    "Minimum word hits: $MIN_WORD_HITS" \
    "All assembled contigs: $ALL_CONTIG_COUNT" \
    "All assembled bp: $ALL_BP" \
    "Seed-anchored candidate contigs: $CANDIDATE_COUNT" \
    "Seed-anchored candidate bp: $CANDIDATE_BP" \
    "Graph-derived candidate paths: $GRAPH_PATH_COUNT" \
    "Graph-derived candidate-path bp: $GRAPH_PATH_BP" \
    "Dual-anchor contigs: $DUAL_COUNT" \
    "Validated and oriented dual-anchor contigs: $ORIENTED_COUNT" \
    "Ambiguous dual-anchor contigs: $AMBIGUOUS_COUNT" \
    "Complete ITS regions: $COMPLETE_ITS_COUNT" \
    "Complete ITS bp: $COMPLETE_ITS_BP" \
    "Partial or ambiguous loci retained: $PARTIAL_LOCUS_COUNT" \
    "Partial or ambiguous locus bp: $PARTIAL_LOCUS_BP" \
    "Residual unphased variants: $VARIANT_COUNT" \
    "Expected taxonomy: ${EXPECTED_TAXONOMY:-NA}" \
    "Expected-taxonomy loci: $EXPECTED_PASS_COUNT" \
    "Non-target loci: $NON_TARGET_COUNT" \
    "Taxonomically coherent loci without an expected-clade test: $COHERENT_TAXONOMY_COUNT" \
    "Taxonomic-chimera loci: $TAXONOMIC_CHIMERA_COUNT" \
    "Unresolved-taxonomy loci: $UNRESOLVED_TAXONOMY_COUNT" \
    "NCBI targeted BLAST enabled: $RUN_NCBI_BLAST" \
    "NCBI BLAST database directory: ${NCBI_DB_DIR:-NA}" \
    "NCBI taxdump directory: ${TAXDUMP_DIR:-NA}" \
    "Total analysis time (HH:MM:SS): $ELAPSED" \
    > "$OUTDIR/run_summary.txt"

if [[ "$KEEP_BAM" != true ]]; then
    rm -f -- "$SEED_BAM"
fi

log "Finished. Candidate contigs: $CANDIDATE_CONTIGS"
log "Graph-derived candidate paths: $OUTDIR/final/graph_candidate_paths.fasta"
log "All complete and partial locus candidates: $OUTDIR/final/all_locus_candidates.fasta"
log "Dual-anchor contigs: $DUAL_CONTIGS"
log "Oriented validated contigs: $OUTDIR/final/oriented_dual_anchor_contigs.fasta"
log "Complete ITS regions: $OUTDIR/final/complete_ITS.fasta"
log "Validation report: $OUTDIR/final/contig_validation.tsv"
log "Path support: $OUTDIR/final/path_support.tsv"
log "Residual unphased variants: $OUTDIR/final/residual_variants.tsv"
log "NCBI top hits: $OUTDIR/final/ncbi_blast_top_hits.tsv"
log "Taxonomy validation: $OUTDIR/final/taxonomy_validation.tsv"
log "Integrated locus summary: $OUTDIR/final/rrna_locus_summary.tsv"
log "Master summary: $OUTDIR/final/master_summary.tsv"
log "Master summary CSV: $OUTDIR/master_summary.csv"
log "Recruitment history: $METRICS"
log "Recruitment graph checks: $GRAPH_METRICS"
log "Total analysis time: $ELAPSED (HH:MM:SS)"
