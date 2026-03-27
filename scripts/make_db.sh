#!/bin/bash

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Default values
input_fasta=""
outdir="./output"
outdb_prefix="rsvDB"

# Get script directory and set scripts path
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPTS_DIR="${SCRIPT_DIR}/../scripts"

# Usage function
usage() {
    cat << EOF
Usage: $0 [options]

Options:
  -i <input_fasta>    Input FASTA file (required)
  -o <outdir>         Output directory (default: ./output)
  -d <outdb_prefix>   Output database prefix (default: rsvDB)
  -h                  Show this help message

Example:
  $0 -i sequences_20250725.fasta -o ./output -d rsvDB-20250725

Requirements:
  - conda activate rsv-analyzer-env before running the script

EOF
    exit 1
}

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Error handling
error_exit() {
    log "ERROR: $1"
    exit "${2:-1}"
}

# Run command
run_cmd() {
    log "Running: $*"
    "$@" || error_exit "Command failed: $*"
}

# Check file exists
check_file() {
    local file="$1"
    local description="${2:-file}"

    [[ -f "$file" ]] || error_exit "$description not found: $file"
    [[ -s "$file" ]] || error_exit "$description is empty: $file"

    log "$description OK: $file"
}

# Skip step if output exists
check_and_run() {
    local output_file="$1"
    local description="$2"
    shift 2

    if [[ -s "$output_file" ]]; then
        log "Skipping $description (exists): $output_file"
        return 0
    fi

    log "Starting $description"
    run_cmd "$@"
    check_file "$output_file" "Output file"
    log "Completed $description"
}

# Parse args
while getopts "i:o:d:h" opt; do
    case $opt in
        i) input_fasta="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        d) outdb_prefix="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate input
[[ -n "$input_fasta" ]] || { echo "ERROR: -i required" >&2; usage; }
check_file "$input_fasta" "Input FASTA"
check_file "${SCRIPTS_DIR}/filterSeqs.py" "filterSeqs.py"

# Create output dir
mkdir -p "$outdir"

# File naming
fasta_prefix="${input_fasta%.fasta}"

# ---------------------------------------------
# Step 1: Filter sequences
# ---------------------------------------------
log "=== Step 1: Filtering FASTA ==="

filtered_fasta="${fasta_prefix}_filter.fasta"

check_and_run "$filtered_fasta" "sequence filtering" \
    python "${SCRIPTS_DIR}/filterSeqs.py" \
    -f "$input_fasta" \
    -o "$filtered_fasta" \
    --minlen 10000 --maxlen 20000 --maxambigs 0

# ---------------------------------------------
# Step 2: CD-HIT clustering
# ---------------------------------------------
log "=== Step 2: CD-HIT clustering ==="

clustered_fasta="${fasta_prefix}.cluster0.99_rep_seq.fasta"
cluster_log="${fasta_prefix}.cluster0.99.log.txt"

if [[ ! -s "$clustered_fasta" ]]; then
    log "Running CD-HIT"
    cd-hit-est \
        -i "$filtered_fasta" \
        -o "$clustered_fasta" \
        -c 0.99 -d 0 -g 1 -M 0 -T 8 -p 1 -sc 1 \
        > "$cluster_log" 2>&1 \
        || error_exit "CD-HIT failed"
else
    log "Skipping CD-HIT (exists)"
fi

# ---------------------------------------------
# Step 3: Mash sketch
# ---------------------------------------------
log "=== Step 3: Mash sketch ==="

mash_output="${outdir}/${outdb_prefix}.msh"

check_and_run "$mash_output" "mash sketch" \
    mash sketch -i "$clustered_fasta" -o "${outdir}/${outdb_prefix}"

# ---------------------------------------------
# Step 4: Organize output
# ---------------------------------------------
log "=== Step 4: Organizing output ==="

final_outdir="${outdir}/${outdb_prefix}"
mkdir -p "$final_outdir"

cp "$clustered_fasta" "${final_outdir}/sequences.fasta"
cp "$mash_output" "${final_outdir}/sequences.msh"

# Summary
log "Generating summary"
cat > "${final_outdir}/processing_summary.txt" << EOF
RSV Database Processing Summary
Generated: $(date)

Input file: $input_fasta
Original sequences: $(grep -c "^>" "$input_fasta" || echo "unknown")
Filtered sequences: $(grep -c "^>" "$filtered_fasta" || echo "unknown")
Clustered sequences: $(grep -c "^>" "$clustered_fasta" || echo "unknown")

Output directory: $final_outdir
EOF

log "Pipeline completed successfully"
log "Output: $final_outdir"
