#!/usr/bin/env bash

# Build the database directory expected by ITSME.

set -Eeuo pipefail
IFS=$'\n\t'

VERSION="1.0.0"
DB_DIR="$PWD"
THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')"
FORCE=false

LEFT_NAME="silva-euk-18s-id95.fasta"
RIGHT_NAME="silva-euk-28s-id98.fasta"
LEFT_URL="https://datacache.galaxyproject.org/indexes/rRNA_databases/silva-euk-18s-id95/silva-euk-18s-id95.fasta"
RIGHT_URL="https://datacache.galaxyproject.org/indexes/rRNA_databases/silva-euk-28s-id98/silva-euk-28s-id98.fasta"
TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

usage() {
    cat <<'EOF'
ITSME database setup v1.0.0

Download and arrange the seed, NCBI BLAST, and taxonomy databases used by
itsme.sh.

Usage:
  setup_db.sh [options]

Options:
  -d, --db-dir DIR    Database directory [current working directory]
  -t, --threads INT   Parallel database-download workers [available CPUs]
      --force         Redownload existing database files
  -h, --help          Show this help
      --version       Show the version

Examples:
  mkdir -p itsme_db && cd itsme_db
  /path/to/setup_db.sh

  setup_db.sh --db-dir /home/ark/databases/itsme_db --threads 8
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

log() {
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"
}

need_value() {
    (( $# >= 2 )) || die "Option $1 requires a value."
}

while (( $# > 0 )); do
    case "$1" in
        -d|--db-dir)
            need_value "$@"; DB_DIR="$2"; shift 2 ;;
        -t|--threads)
            need_value "$@"; THREADS="$2"; shift 2 ;;
        --force)
            FORCE=true; shift ;;
        -h|--help)
            usage; exit 0 ;;
        --version)
            printf '%s\n' "$VERSION"; exit 0 ;;
        *)
            die "Unknown option: $1" ;;
    esac
done

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer."
command -v update_blastdb >/dev/null 2>&1 || \
    die "update_blastdb was not found. Activate the ITSME environment first."
command -v tar >/dev/null 2>&1 || die "tar was not found."

if command -v curl >/dev/null 2>&1; then
    DOWNLOADER="curl"
elif command -v wget >/dev/null 2>&1; then
    DOWNLOADER="wget"
else
    die "Either curl or wget is required."
fi

mkdir -p "$DB_DIR"
DB_DIR=$(cd "$DB_DIR" && pwd)
NCBI_DIR="$DB_DIR/NCBI_rRNA_BLAST"
TAX_DIR="$DB_DIR/taxdump"
mkdir -p "$NCBI_DIR" "$TAX_DIR"

download_file() {
    local url="$1"
    local destination="$2"
    local temporary="${destination}.part"

    if [[ -s "$destination" && "$FORCE" == false ]]; then
        log "Keeping existing file: $destination"
        return
    fi

    log "Downloading $(basename "$destination")."
    rm -f "$temporary"
    if [[ "$DOWNLOADER" == "curl" ]]; then
        curl --fail --location --retry 3 --continue-at - \
            --output "$temporary" "$url"
    else
        wget --tries=3 --continue --output-document="$temporary" "$url"
    fi
    [[ -s "$temporary" ]] || die "Download produced an empty file: $url"
    mv -f "$temporary" "$destination"
}

download_file "$LEFT_URL" "$DB_DIR/$LEFT_NAME"
download_file "$RIGHT_URL" "$DB_DIR/$RIGHT_NAME"

log "Downloading NCBI SSU, ITS, LSU, and taxonomy BLAST databases."
UPDATE_ARGS=(--blastdb_version 5 --decompress --num_cores "$THREADS")
[[ "$FORCE" == true ]] && UPDATE_ARGS+=(--force)
(
    cd "$NCBI_DIR"
    update_blastdb "${UPDATE_ARGS[@]}" \
        SSU_eukaryote_rRNA \
        ITS_eukaryote_sequences \
        LSU_eukaryote_rRNA \
        taxdb
)

TAX_ARCHIVE="$TAX_DIR/taxdump.tar.gz"
download_file "$TAXDUMP_URL" "$TAX_ARCHIVE"
log "Extracting NCBI taxonomy ranks."
tar -xzf "$TAX_ARCHIVE" -C "$TAX_DIR" nodes.dmp names.dmp merged.dmp

for required in \
    "$DB_DIR/$LEFT_NAME" \
    "$DB_DIR/$RIGHT_NAME" \
    "$TAX_DIR/nodes.dmp" \
    "$TAX_DIR/names.dmp"; do
    [[ -s "$required" ]] || die "Required file is missing or empty: $required"
done

for database in SSU_eukaryote_rRNA ITS_eukaryote_sequences LSU_eukaryote_rRNA; do
    compgen -G "$NCBI_DIR/${database}.n*" >/dev/null || \
        die "NCBI BLAST database was not created: $database"
done

cat > "$DB_DIR/ITSME_DATABASES.txt" <<EOF
ITSME database setup completed: $(date -Is)
Database root: $DB_DIR
18S seed: $LEFT_NAME
28S seed: $RIGHT_NAME
NCBI BLAST databases: NCBI_rRNA_BLAST/
NCBI taxonomy dump: taxdump/
EOF

log "Database setup complete: $DB_DIR"
log "Run ITSME with: --db-dir $DB_DIR"
