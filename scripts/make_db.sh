#!/bin/bash

# ---------------------------------------------
# The process to make RSV A, B mash database
# Sequence data, XML, and CSV files were downloaded from:
# https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=human%20respiratory%20syncytial%20virus,%20taxid:11250#
# Usage: sh make_db.sh -i sequences_20250725.fasta -o ./output -d rsvDB-20250725
# ---------------------------------------------

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Default values
input_fasta=""
outdir="./output"
outdb_prefix="rsvDB"
ENV_NAME="virus_env"

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
  - Conda environment '$ENV_NAME' with: mash, seqkit, nextclade, biopython, nextflow, cd-hit
  - Python scripts in ../scripts/ directory

EOF
    exit 1
}

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit "${2:-1}"
}

# Function to check if conda environment exists
check_conda_env() {
    local env_name="$1"
    if conda info --envs | awk -v env="$env_name" '$1 == env {found=1} END {exit !found}'; then
        log "Environment '$env_name' exists"
        return 0
    else
        log "Environment '$env_name' does not exist"
        return 1
    fi
}

# Function to activate conda environment
activate_env() {
    local env_name="$1"
    
    # Find conda installation
    local conda_base
    if command -v conda >/dev/null 2>&1; then
        conda_base="$(conda info --base 2>/dev/null)"
    else
        # Try common conda locations
        for conda_path in "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda" "/usr/local/miniconda3" "/usr/local/anaconda3"; do
            if [[ -d "$conda_path" ]]; then
                conda_base="$conda_path"
                break
            fi
        done
    fi
    
    [[ -n "$conda_base" ]] || error_exit "Could not find conda installation"
    
    # Initialize conda for this session
    log "Initializing conda from: $conda_base"
    # shellcheck source=/dev/null
    source "${conda_base}/etc/profile.d/conda.sh" || error_exit "Could not initialize conda"
    
    # Verify conda is now available
    command -v conda >/dev/null 2>&1 || error_exit "Conda still not available after initialization"
    
    # Activate environment if not already active
    if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
        log "Activating conda environment '$env_name'"
        conda activate "$env_name" || error_exit "Failed to activate conda environment '$env_name'"
        log "Successfully activated environment '$env_name'"
    else
        log "Environment '$env_name' already active"
    fi
}

# Function to run command with logging and error handling
run_cmd() {
    local cmd="$*"
    log "Running: $cmd"
    if ! eval "$cmd"; then
        error_exit "Command failed: $cmd"
    fi
}

# Function to check file exists and is not empty
check_file() {
    local file="$1"
    local description="${2:-file}"
    
    if [[ ! -f "$file" ]]; then
        error_exit "$description not found: $file"
    fi
    
    if [[ ! -s "$file" ]]; then
        error_exit "$description is empty: $file"
    fi
    
    log "$description found and non-empty: $file"
}

# Function to check if output file exists and skip if so
check_and_run() {
    local output_file="$1"
    local cmd="$2"
    local description="${3:-Processing}"
    
    if [[ -f "$output_file" ]] && [[ -s "$output_file" ]]; then
        log "Skipping $description - output already exists: $output_file"
        return 0
    fi
    
    log "Starting $description"
    run_cmd "$cmd"
    check_file "$output_file" "Output file"
    log "Completed $description"
}

# Parse command-line arguments
while getopts "i:o:d:h" opt; do
    case $opt in
        i) input_fasta="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        d) outdb_prefix="$OPTARG" ;;
        h) usage ;;
        *) 
            echo "Invalid option: -$OPTARG" >&2
            usage 
            ;;
    esac
done

# Validate required arguments
if [[ -z "$input_fasta" ]]; then
    echo "ERROR: Input FASTA file (-i) is required." >&2
    usage
fi

# Validate input file
check_file "$input_fasta" "Input FASTA file"

# Check required Python scripts
check_file "${SCRIPTS_DIR}/filterSeqs.py" "Filter sequences script"

# Set up conda environment
log "Setting up conda environment"
if ! check_conda_env "$ENV_NAME"; then
    error_exit "Environment '$ENV_NAME' does not exist. Create it with:
    mamba create -n '$ENV_NAME' mash=2.3 seqkit=2.10.1 nextclade=3.16.0 biopython=1.84 nextflow cd-hit=4.8.1 -y"
fi

activate_env "$ENV_NAME"

# Create output directory
log "Creating output directory: $outdir"
mkdir -p "$outdir" || error_exit "Failed to create output directory: $outdir"

# Get base filename without extension
fasta_base="$(basename "$input_fasta" .fasta)"
fasta_prefix="${input_fasta%.fasta}"

# ---------------------------------------------
# Step 1: Filter and reformat FASTA sequences
log "=== Step 1: Filtering FASTA sequences ==="
filtered_fasta="${fasta_prefix}_filter.fasta"
filter_cmd="python '${SCRIPTS_DIR}/filterSeqs.py' -f '$input_fasta' -o '$filtered_fasta' --minlen 10000 --maxlen 20000 --maxambigs 0"

check_and_run "$filtered_fasta" "$filter_cmd" "sequence filtering"

# ---------------------------------------------
# Step 2: CD-HIT clustering to remove redundancy
log "=== Step 2: CD-HIT clustering ==="
clustered_fasta="${fasta_prefix}.cluster0.99_rep_seq.fasta"
cluster_log="${fasta_prefix}.cluster0.99.log.txt"

cluster_cmd="cd-hit-est -i '$filtered_fasta' -o '$clustered_fasta' -c 0.99 -d 0 -g 1 -M 0 -T 8 -p 1 -sc 1 > '$cluster_log' 2>&1"

check_and_run "$clustered_fasta" "$cluster_cmd" "CD-HIT clustering"

# ---------------------------------------------
# Step 3: Mash sketch creation
log "=== Step 3: Creating Mash sketch ==="
mash_output="${outdir}/${outdb_prefix}.msh"
mash_cmd="mash sketch -i '$clustered_fasta' -o '${outdir}/${outdb_prefix}'"

check_and_run "$mash_output" "$mash_cmd" "Mash sketching"

# ---------------------------------------------
# Step 4: Organize final output
log "=== Step 4: Organizing final output ==="
final_outdir="${outdir}/${outdb_prefix}"
mkdir -p "$final_outdir" || error_exit "Failed to create final output directory: $final_outdir"

# Copy files to final location
log "Copying clustered sequences to final location"
cp "$clustered_fasta" "${final_outdir}/sequences.fasta" || error_exit "Failed to copy sequences"

log "Copying Mash sketch to final location"
cp "$mash_output" "${final_outdir}/sequences.msh" || error_exit "Failed to copy Mash sketch"

# Generate summary report
log "Generating summary report"
cat > "${final_outdir}/processing_summary.txt" << EOF
RSV Database Processing Summary
Generated: $(date)

Input file: $input_fasta
Original sequences: $(grep -c "^>" "$input_fasta" || echo "unknown")
Filtered sequences: $(grep -c "^>" "$filtered_fasta" || echo "unknown")
Clustered sequences: $(grep -c "^>" "$clustered_fasta" || echo "unknown")

Output directory: $final_outdir
Files created:
- sequences.fasta (clustered representative sequences)
- sequences.msh (Mash sketch database)
- processing_summary.txt (this file)

Log files:
- filterSeq.log.txt (if generated)
- $cluster_log
EOF

log "Pipeline completed successfully!"
log "Output directory: $final_outdir"
log "Summary report: ${final_outdir}/processing_summary.txt"