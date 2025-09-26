#!/bin/bash

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Global variables
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly RSV_PROG_BASE="${SCRIPT_DIR}/.."
readonly PROG_BASE="/nfs/APL_Genomics/apps/production"
readonly path_to_qc_pipeline="${PROG_BASE}/qc_pipelines/nf-qc-illumina"
readonly path_to_viralrecon="${PROG_BASE}/viralrecon"
readonly ENV_NAME="virus_env"
readonly RESULTS_DIR="results"

# Command line argument handling
usage() {
    cat << EOF
Usage: $0 <raw_directory>

Arguments:
    raw_directory    Path to directory containing raw nanopore data

Example:
    $0 /path/to/raw_data
EOF
}


# Check command line arguments
if [[ $# -ne 1 ]]; then
    usage
    exit 1
fi

readonly RAW_DIR="$1"

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
        return 1
    fi
}

# Function to activate conda environment
activate_env() {
    local env_name="$1"

    # Ensure conda is initialized (only once)
    if ! type conda > /dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh" || {
            echo "❌ Could not initialize conda." >&2
            return 1
        }
    fi

    # Check current env
    if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
        echo "🔁 Activating '$env_name'..."
        conda activate "$env_name" || {
            echo "❌ Failed to activate '$env_name'" >&2
            return 1
        }
    else
        echo "✅ Environment '$env_name' already active."
    fi
}

# Function to run command with logging
run_cmd() {
    log "Running: $*"
    "$@" || error_exit "Command failed: $*"
}



# Function to check file exists and is not empty
check_file() {
    local file="$1"
    local description="${2:-file}"

    if [[ ! -f "$file" ]]; then
        error_exit "$description not found: $file"
    fi

    if [[ ! -s "$file" ]]; then
        log "WARNING: $description is empty: $file"
        return 1
    fi

    return 0
}

# Validate raw directory
if [[ ! -d "$RAW_DIR" ]]; then
    error_exit "Raw directory does not exist: $RAW_DIR"
fi

# Initialize conda
eval "$(conda shell.bash hook)" || error_exit "Failed to initialize conda"

# Main execution starts here
log "Starting RSV nanopore analysis pipeline"
log "Working directory: $(pwd)"
log "Script directory: $SCRIPT_DIR"
log "Raw data directory: $RAW_DIR"

# Check conda environment
if ! check_conda_env "$ENV_NAME"; then
    error_exit "Environment '$ENV_NAME' does not exist. Create it with:
    mamba create -n '$ENV_NAME' mash=2.3 seqkit=2.10.1 nextclade=3.16.0 biopython=1.84 nextflow cd-hit=4.8.1 -y"
fi

activate_env "$ENV_NAME"

# Create results directory
mkdir -p "$RESULTS_DIR"

# ==========================================
# STEP 1: Generate samplesheet
# ==========================================
log "=== STEP 1: Generating samplesheet ==="

readonly SAMPLESHEET="$RESULTS_DIR/samplesheet_to_qc.csv"
if [[ ! -f "$SAMPLESHEET" ]]; then
    log "Creating samplesheet from raw directory: $RAW_DIR"
    run_cmd python "$SCRIPT_DIR/make_samplesheet.py" "$RAW_DIR" > "$SAMPLESHEET"
else
    log "Samplesheet already exists: $SAMPLESHEET"
fi

check_file "$SAMPLESHEET" "samplesheet"

# ==========================================
# STEP 2: Quality Control
# ==========================================
log "=== STEP 2: Running QC pipeline ==="

run_cmd nextflow ${path_to_qc_pipeline}/main.nf \
-profile singularity \
--input $SAMPLESHEET \
--outdir $RESULTS_DIR/nf-qc-illumina \
-resume

# ==========================================
# STEP 3: RSV A/B Classification
# ==========================================
log "=== STEP 3: RSV A/B classification ==="

run_cmd python $SCRIPT_DIR/screen_rsv_mash.py \
-i $RESULTS_DIR/nf-qc-illumina/qc_reads \
-d ${RSV_PROG_BASE}/resource/db/mash_screen/sequences.msh \
-o $RESULTS_DIR

# ==========================================
# STEP 4: Process RSV-A
# ==========================================
log "=== STEP 4: Processing RSV-A ==="

readonly RSVA_SAMPLESHEET="$RESULTS_DIR/samplesheet_rsvA.csv"
if check_file "$RSVA_SAMPLESHEET" "RSV-A samplesheet"; then
    log "Running RSV-A viralrecon"
    run_cmd nextflow ${path_to_viralrecon}/main.nf \
    -profile singularity \
    -c ${SCRIPT_DIR}/rsvA_illumina.config \
    --input $RSVA_SAMPLESHEET \
    --outdir $RESULTS_DIR/rsvA \
    --platform illumina \
    -resume

    log "Copy RSV-A amplicon genome coverage to the results base dir for better readability"
    run_cmd cp -r $RESULTS_DIR/rsvA/variants/bowtie2/mosdepth $RESULTS_DIR/rsvA/nf-ampgenomecov/

    log "Generating RSV-A consensus statistics and copy consensus files to consensus directory in the base directory"
    run_cmd mkdir -p  $RESULTS_DIR/rsvA/consensus/

    run_cmd cp  $RESULTS_DIR/rsvA/variants/bcftools/consensus/bcftools/*.consensus.fa $RESULTS_DIR/rsvA/consensus/

    run_cmd python $SCRIPT_DIR/consensus_stats.py \
        $RESULTS_DIR/rsvA/consensus/ \
        $RESULTS_DIR/rsvA/consensus/all_consensus_stats.tsv

    log "Running RSV-A nextclade analysis"
    # Check if consensus files exist before running nextclade
    if ls "$RESULTS_DIR/rsvA/consensus"/*.fa >/dev/null 2>&1; then
        run_cmd nextclade run $RESULTS_DIR/rsvA/consensus/*.fa \
            --input-dataset $RSV_PROG_BASE/resource/db/nextclade/db_rsvA_2025-08-22 \
            --output-all $RESULTS_DIR/rsvA/nextclade

    else
        log "WARNING: No RSV-A consensus files found for nextclade analysis"
    fi
else
    log "No RSV-A samples found, skipping RSV-A processing"
fi

# ==========================================
# STEP 5: Process RSV-B
# ==========================================

log "=== STEP 5: Processing RSV-B ==="

readonly RSVB_SAMPLESHEET="$RESULTS_DIR/samplesheet_rsvB.csv"
if check_file "$RSVB_SAMPLESHEET" "RSV-B samplesheet"; then
    log "Running RSV-B viralrecon"
    run_cmd nextflow ${path_to_viralrecon}/main.nf \
    -profile singularity \
    -c  ${SCRIPT_DIR}/rsvB_illumina.config \
    --input $RSVB_SAMPLESHEET \
    --outdir $RESULTS_DIR/rsvB \
    --platform illumina \
    -resume

    log "Copy RSV-B amplicon genome coverage to the results base dir for better readability"
    run_cmd cp -r $RESULTS_DIR/rsvB/variants/bowtie2/mosdepth $RESULTS_DIR/rsvB/nf-ampgenomecov/

    log "Generating RSV-B consensus statistics and copy consensus files to consensus directory in the base directory"
    run_cmd mkdir -p  $RESULTS_DIR/rsvB/consensus/

    run_cmd cp  $RESULTS_DIR/rsvB/variants/bcftools/consensus/bcftools/*.consensus.fa $RESULTS_DIR/rsvB/consensus/

    run_cmd python $SCRIPT_DIR/consensus_stats.py \
        $RESULTS_DIR/rsvB/consensus/ \
        $RESULTS_DIR/rsvB/consensus/all_consensus_stats.tsv

    log "Running RSV-B nextclade analysis"
    # Check if consensus files exist before running nextclade
    if ls "$RESULTS_DIR/rsvB/consensus"/*.fa >/dev/null 2>&1; then
        run_cmd nextclade run $RESULTS_DIR/rsvB/consensus/*.fa \
            --input-dataset $RSV_PROG_BASE/resource/db/nextclade/db_rsvB_2025-08-22 \
            --output-all $RESULTS_DIR/rsvB/nextclade
    else
        log "WARNING: No RSV-B consensus files found for nextclade analysis"
    fi
else
    log "No RSV-B samples found, skipping RSV-B processing"
fi

# ==========================================
# CLEANUP AND SUMMARY
# ==========================================
log "=== PIPELINE COMPLETED SUCCESSFULLY ==="
log "Results directory: $RESULTS_DIR"
log "Finished at: $(date)"

# Ensure we're not in any conda environment when script ends
if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    conda deactivate
fi

exit 0
