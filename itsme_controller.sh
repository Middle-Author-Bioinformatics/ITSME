#!/usr/bin/env bash

# Batch controller for ITSME v1.
# Discovers paired FASTQs, creates one output directory per library, runs
# itsme.sh sequentially, and combines per-library master summaries.

set -Eeuo pipefail
IFS=$'\n\t'

VERSION="1.3.0"
INPUT_DIR=""
OUTPUT_DIR=""
ITSME="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/itsme.sh"
RESUME=false
STOP_ON_ERROR=false
DRY_RUN=false
ITSME_ARGS=()

usage() {
    cat <<'EOF'
ITSME controller v1.3.0

Run itsme.sh sequentially for every paired FASTQ library in a directory.

Usage:
  itsme_controller.sh -i FASTQ_DIR -o OUTPUT_DIR [controller options] \
      [itsme.sh options]

Required controller arguments:
  -i, --input-dir DIR       Directory containing paired FASTQ files
  -o, --output-dir DIR      Parent directory for per-library results

Controller options:
      --itsme FILE          ITSME executable [itsme.sh beside this script]
      --resume              Skip libraries with completed existing results
      --stop-on-error       Stop after the first failed library
      --dry-run             Print commands without running them
  -h, --help                Show this help
      --version             Show the controller version

ITSME options:
  All options not owned by the controller are passed directly to itsme.sh.
  In particular, provide --db-dir normally; no separator is required.

Input naming:
  The controller discovers files containing _R1 and ending in .fastq.gz.
  It replaces the first _R1 with _R2 to obtain the mate filename. The library
  name is everything before the first _R1.

Example:
  itsme_controller.sh \
      -i /data/fastqs \
      -o /data/itsme_results \
      --resume \
      --db-dir /home/ark/databases/itsme_db

Outputs:
  OUTPUT_DIR/LIBRARY/                Complete ITSME result for one library
  OUTPUT_DIR/logs/LIBRARY.log        Combined stdout and stderr for that run
  OUTPUT_DIR/batch_status.tsv        Status and runtime for every library
  OUTPUT_DIR/batch_master_summary.csv
                                     Combined successful master summaries

Notes:
  * Runs are sequential to avoid oversubscribing CPUs and memory.
  * Do not pass -1, -2, or a per-library -o; the controller supplies them.
  * With --resume, only directories containing both run_summary.txt and
    master_summary.csv are treated as completed.
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

format_duration() {
    local seconds="$1"
    printf '%02d:%02d:%02d' \
        "$((seconds / 3600))" \
        "$(((seconds % 3600) / 60))" \
        "$((seconds % 60))"
}

while (( $# > 0 )); do
    case "$1" in
        -i|--input-dir)
            need_value "$@"; INPUT_DIR="$2"; shift 2 ;;
        -o|--output-dir)
            need_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
        --itsme)
            need_value "$@"; ITSME="$2"; shift 2 ;;
        --resume)
            RESUME=true; shift ;;
        --stop-on-error)
            STOP_ON_ERROR=true; shift ;;
        --dry-run)
            DRY_RUN=true; shift ;;
        -h|--help)
            usage; exit 0 ;;
        --version)
            printf '%s\n' "$VERSION"; exit 0 ;;
        --)
            # Accepted for backward compatibility, but no longer required.
            shift
            ITSME_ARGS+=("$@")
            break ;;
        *)
            ITSME_ARGS+=("$1"); shift ;;
    esac
done

[[ -n "$INPUT_DIR" ]] || die "Missing -i/--input-dir."
[[ -n "$OUTPUT_DIR" ]] || die "Missing -o/--output-dir."
[[ -d "$INPUT_DIR" ]] || die "Input directory not found: $INPUT_DIR"
[[ -x "$ITSME" ]] || die "ITSME executable not found or not executable: $ITSME"

INPUT_DIR=$(cd "$INPUT_DIR" && pwd)
mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/logs"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
ITSME=$(cd "$(dirname "$ITSME")" && pwd)/$(basename "$ITSME")

# Prevent the controller from silently overriding arguments it owns.
for argument in "${ITSME_ARGS[@]}"; do
    case "$argument" in
        -1|--reads1|-2|--reads2|-o|--output-dir)
            die "Do not pass $argument to itsme.sh; the controller supplies input and output paths." ;;
    esac
done

R1_FILES=()
while IFS= read -r -d '' file; do
    R1_FILES+=("$file")
done < <(find "$INPUT_DIR" -maxdepth 1 -type f \
    -name '*_R1*.fastq.gz' -print0 | LC_ALL=C sort -z)

(( ${#R1_FILES[@]} > 0 )) || \
    die "No files matching *_R1*.fastq.gz were found in $INPUT_DIR"

declare -A SEEN_LIBRARIES=()
LIBRARIES=()
R2_FILES=()

# Validate every pair and output name before starting any expensive work.
for r1 in "${R1_FILES[@]}"; do
    filename=${r1##*/}
    library=${filename%%_R1*}
    suffix=${filename#*_R1}
    r2="$INPUT_DIR/${library}_R2${suffix}"

    [[ -n "$library" && "$filename" != "$library" ]] || \
        die "Could not derive a library prefix from: $filename"
    [[ -s "$r2" ]] || die "Missing or empty R2 mate for $filename: ${r2##*/}"
    [[ -z "${SEEN_LIBRARIES[$library]:-}" ]] || \
        die "Multiple R1 files produce the same library prefix: $library"

    SEEN_LIBRARIES[$library]=1
    LIBRARIES+=("$library")
    R2_FILES+=("$r2")
done

STATUS_FILE="$OUTPUT_DIR/batch_status.tsv"
printf 'library\tr1\tr2\toutput_directory\tstatus\texit_code\telapsed_hh_mm_ss\tstarted\tfinished\n' \
    > "$STATUS_FILE"

successes=0
failures=0
skipped=0

log "Discovered ${#LIBRARIES[@]} paired libraries in $INPUT_DIR."
log "Per-library results will be written beneath $OUTPUT_DIR."

for index in "${!LIBRARIES[@]}"; do
    library=${LIBRARIES[$index]}
    r1=${R1_FILES[$index]}
    r2=${R2_FILES[$index]}
    library_output="$OUTPUT_DIR/$library"
    run_log="$OUTPUT_DIR/logs/$library.log"
    started=$(date '+%Y-%m-%d %H:%M:%S')
    start_seconds=$SECONDS

    if [[ -d "$library_output" && \
          -s "$library_output/run_summary.txt" && \
          -s "$library_output/master_summary.csv" && \
          "$RESUME" == true ]]; then
        finished=$(date '+%Y-%m-%d %H:%M:%S')
        printf '%s\t%s\t%s\t%s\tskipped_complete\t0\t00:00:00\t%s\t%s\n' \
            "$library" "${r1##*/}" "${r2##*/}" "$library_output" \
            "$started" "$finished" >> "$STATUS_FILE"
        log "Skipping completed library: $library"
        skipped=$((skipped + 1))
        continue
    fi

    if [[ -e "$library_output" && \
          -n "$(find "$library_output" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
        finished=$(date '+%Y-%m-%d %H:%M:%S')
        printf '%s\t%s\t%s\t%s\tblocked_existing_output\t2\t00:00:00\t%s\t%s\n' \
            "$library" "${r1##*/}" "${r2##*/}" "$library_output" \
            "$started" "$finished" >> "$STATUS_FILE"
        printf 'ERROR: Output for %s already exists but is not resumable: %s\n' \
            "$library" "$library_output" >&2
        failures=$((failures + 1))
        [[ "$STOP_ON_ERROR" == true ]] && break
        continue
    fi

    command=("$ITSME" -1 "$r1" -2 "$r2" -o "$library_output" "${ITSME_ARGS[@]}")
    log "Starting library $library ($((index + 1))/${#LIBRARIES[@]})."

    if [[ "$DRY_RUN" == true ]]; then
        printf 'DRY RUN:'
        printf ' %q' "${command[@]}"
        printf '\n'
        status="dry_run"
        exit_code=0
        skipped=$((skipped + 1))
    elif "${command[@]}" 2>&1 | tee "$run_log"; then
        status="success"
        exit_code=0
        successes=$((successes + 1))
    else
        exit_code=$?
        status="failed"
        failures=$((failures + 1))
    fi

    elapsed=$(format_duration "$((SECONDS - start_seconds))")
    finished=$(date '+%Y-%m-%d %H:%M:%S')
    printf '%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' \
        "$library" "${r1##*/}" "${r2##*/}" "$library_output" \
        "$status" "$exit_code" "$elapsed" "$started" "$finished" \
        >> "$STATUS_FILE"

    log "Library $library finished with status=$status in $elapsed."
    if [[ "$status" == "failed" && "$STOP_ON_ERROR" == true ]]; then
        break
    fi
done

# Combine available per-library summaries. Each row receives a leading library
# column so identical contig names from different samples remain distinguishable.
python3 - "$OUTPUT_DIR" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
destination = root / "batch_master_summary.csv"
summaries = sorted(
    path for path in root.glob("*/master_summary.csv")
    if path.is_file() and path.stat().st_size > 0
)

header = None
with destination.open("w", newline="") as output:
    writer = csv.writer(output, lineterminator="\n")
    for path in summaries:
        with path.open(newline="") as source:
            reader = csv.reader(source)
            try:
                current_header = next(reader)
            except StopIteration:
                continue
            if header is None:
                header = current_header
                writer.writerow(["library", *header])
            elif current_header != header:
                raise SystemExit(
                    f"Summary columns differ for {path}; cannot combine batch CSV."
                )
            for row in reader:
                writer.writerow([path.parent.name, *row])
PY

log "Batch complete: successes=$successes, skipped=$skipped, failures=$failures."
log "Status table: $STATUS_FILE"
log "Combined summary: $OUTPUT_DIR/batch_master_summary.csv"

(( failures == 0 ))
