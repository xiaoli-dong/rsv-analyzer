#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, and pipe failures

# ==========================================================
# Paths & constants
# ==========================================================
readonly prod_prog_base="/nfs/APL_Genomics/apps/production"
readonly deve_prog_base="/nfs/Genomics_DEV/projects/xdong/deve"

readonly path_to_viralassembly="${prod_prog_base}/viralassembly"
readonly path_to_qc_pipeline="${deve_prog_base}/nf-qcflow"
readonly path_to_covflow="${deve_prog_base}/nf-covflow"

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly RSV_PROG_BASE="${SCRIPT_DIR}/.."
readonly VERSION="$(cat "$SCRIPT_DIR/VERSION" 2>/dev/null || echo "unknown")"
readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
readonly RSV_SLURM_CONFIG="${RSV_PROG_BASE}/conf/slurm.config"
readonly ENV_NAME="virus_env"
#readonly RESULTS_DIR="results"

# Genome references
readonly rsvA_ref="${RSV_PROG_BASE}/resource/primerschemes/rsvA/v3/reference.fasta"
readonly rsvA_bed="${RSV_PROG_BASE}/resource/primerschemes/rsvA/v3/scheme.bed"
readonly rsvB_ref="${RSV_PROG_BASE}/resource/primerschemes/rsvB/v3/reference.fasta"
readonly rsvB_bed="${RSV_PROG_BASE}/resource/primerschemes/rsvB/v3/scheme.bed"
readonly MASH_DB="${RSV_PROG_BASE}/resource/db/mash_screen/sequences.msh"

# Viralassembly config (default + override)
DEFAULT_QCFLOW_CONFIG_FILE="${RSV_PROG_BASE}/conf/qcflow.config"
QCFLOW_CONFIG_FILE="$DEFAULT_QCFLOW_CONFIG_FILE"

# Viralassembly config (default + override)
DEFAULT_VIRALASSEMBLY_CONFIG_FILE="${RSV_PROG_BASE}/conf/viralassembly.config"
VIRALASSEMBLY_CONFIG_FILE="$DEFAULT_VIRALASSEMBLY_CONFIG_FILE"

# ==========================================================
# Help / Version
# ==========================================================
show_version() {
    echo "$SCRIPT_NAME version $VERSION"
}

show_help() {
    cat << EOF
$SCRIPT_NAME - RSV Nanopore Analysis Pipeline

Usage:
   $SCRIPT_NAME [options] <samplesheet.csv> <results_dir>


Options:
  --viralassembly-config FILE  Path to config file (default: ${SCRIPT_DIR}/rsv_nanopore.config)
  -h, --help                   Show this help
  -v, --version                Show version

Example:
  $SCRIPT_NAME samplesheet.csv
  $SCRIPT_NAME samplesheet.csv myresults --viralassembly-config custom.config

EOF
}

# ==========================================================
# Argument parsing
# ==========================================================
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) show_help; exit 0 ;;
        -v|--version) show_version; exit 0 ;;
        --viralassembly-config)
            VIRALASSEMBLY_CONFIG_FILE="${2:-}"
            shift 2
            ;;
        --*) echo "Unknown option: $1" >&2; exit 1 ;;
        *) POSITIONAL_ARGS+=("$1"); shift ;;
    esac
done

set -- "${POSITIONAL_ARGS[@]}"
[[ $# -eq 2 ]] || { show_help; exit 1; }

readonly INPUT_SAMPLESHEET="$1"
readonly RESULTS_DIR="$(mkdir -p "$2" && cd "$2" && pwd)"

# ==========================================================
# Logging & helpers
# ==========================================================
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$SCRIPT_NAME-$VERSION] $*" >&2
}

error_exit() {
    log "ERROR: $1"
    exit "${2:-1}"
}

run_cmd() {
    log "Running: $*"
    "$@" || error_exit "Command failed: $*"
}

check_file() {
    [[ -s "$1" ]] || error_exit "File missing or empty: $1"
}

check_conda_env() {
    conda info --envs | awk -v e="$1" '$1==e{f=1}END{exit !f}'
}

activate_env() {
    eval "$(conda shell.bash hook)"
    [[ "${CONDA_DEFAULT_ENV:-}" == "$1" ]] || conda activate "$1"
}

# ==========================================================
# Validation
# ==========================================================
check_file "$INPUT_SAMPLESHEET"
check_file "$VIRALASSEMBLY_CONFIG_FILE"

log "Input samplesheet      : $INPUT_SAMPLESHEET"
log "Results directory      : $RESULTS_DIR"
log "Viralassembly config file : $VIRALASSEMBLY_CONFIG_FILE"

# ==========================================================
# Conda
# ==========================================================
check_conda_env "$ENV_NAME" || error_exit "Conda env '$ENV_NAME' not found"
activate_env "$ENV_NAME"

# ==========================================================
# STEP 1: Prepare samplesheet
# ==========================================================
log "=== STEP 1: Preparing samplesheet ==="
readonly SAMPLESHEET="$RESULTS_DIR/samplesheet_to_qc.csv"
cp "$INPUT_SAMPLESHEET" "$SAMPLESHEET"
check_file "$SAMPLESHEET"

# ==========================================================
# STEP 2: QC pipeline
# ==========================================================
log "=== STEP 2: Running QC pipeline ==="
run_cmd nextflow run "${path_to_qc_pipeline}/main.nf" \
    -profile singularity,slurm \
    -c "$QCFLOW_CONFIG_FILE" \
    --input "$SAMPLESHEET" \
    --platform nanopore \
    --outdir "$RESULTS_DIR/nf-qcflow" \
    -resume

# ==========================================================
# STEP 3: RSV A/B classification
# ==========================================================
log "=== STEP 3: RSV A/B classification ==="
run_cmd python "$SCRIPT_DIR/screen_rsv_mash.py" \
    -i "$RESULTS_DIR/nf-qcflow/qcreads" \
    -d "$MASH_DB" \
    -o "$RESULTS_DIR"

# ==========================================================
# STEP 4: Process virus function
# ==========================================================
process_virus() {
    local virus="$1"
    local ref="$2"
    local bed="$3"

    local ss="$RESULTS_DIR/samplesheet_${virus}.csv"
    [[ -s "$ss" ]] || { log "No $virus samples found, skipping"; return; }

    log "=== Processing $virus ==="

    # Viral assembly
    run_cmd nextflow run "${path_to_viralassembly}/main.nf" \
        -profile singularity,slurm \
        -c "$RSV_SLURM_CONFIG" \
        -c "$VIRALASSEMBLY_CONFIG_FILE" \
        --input "$ss" \
        --outdir "$RESULTS_DIR/$virus" \
        --scheme "$virus" \
        -resume

    # Prepare for nf-covflow
    mkdir -p "$RESULTS_DIR/$virus/bam_to_covflow"
    cp "$RESULTS_DIR/$virus/bam/"*.primertrimmed.rg.sorted.bam* "$RESULTS_DIR/$virus/bam_to_covflow" || true

    # Build covflow samplesheet
    run_cmd python "$SCRIPT_DIR/build_samplesheet_for_covflow.py" \
        --input_dir "$RESULTS_DIR/$virus/bam_to_covflow" \
        --ref_fasta "$ref" \
        --bed_file "$bed" \
        -o "$RESULTS_DIR/samplesheet_to_covflow_${virus}.csv"

    # Run covflow
    run_cmd nextflow run "${path_to_covflow}/main.nf" \
        -profile singularity,slurm \
        --input "$RESULTS_DIR/samplesheet_to_covflow_${virus}.csv" \
        --outdir "$RESULTS_DIR/$virus/nf-covflow" \
        -resume

    # Consensus stats
    mkdir -p "$RESULTS_DIR/$virus/consensus"
    run_cmd python "$SCRIPT_DIR/consensus_stats.py" \
        "$RESULTS_DIR/$virus/consensus" \
        "$RESULTS_DIR/$virus/consensus/all_consensus_stats.tsv"

    # Nextclade
    last_char="${virus: -1}"
    nextclade_db="./nextstrain/rsv/${last_char,,}"
    run_cmd nextclade dataset get \
      --name "nextstrain/rsv/${last_char,,}" \
      --output-dir "$nextclade_db"

    if ls "$RESULTS_DIR/$virus/consensus/"*.fasta >/dev/null 2>&1; then
        run_cmd nextclade run "$RESULTS_DIR/$virus/consensus/"*.fasta \
            --input-dataset "$nextclade_db" \
            --output-all "$RESULTS_DIR/$virus/nextclade"
    else
        log "WARNING: No $virus consensus FASTA files found for nextclade"
    fi
}

process_virus rsvA "$rsvA_ref" "$rsvA_bed"
process_virus rsvB "$rsvB_ref" "$rsvB_bed"

# ==========================================================
# STEP 5: Final report
# ==========================================================
log "=== Generating final summary report ==="
final_report_dir="${RESULTS_DIR}/summary_report"
mkdir -p "$final_report_dir"

for virus in rsvA rsvB; do
    mkdir -p "${final_report_dir}/${virus}"
    cp "$RESULTS_DIR/samplesheet_${virus}.csv" "${final_report_dir}/"
    cp "$RESULTS_DIR/$virus/consensus/"* "${final_report_dir}/${virus}/" || true
    cp "$RESULTS_DIR/$virus/bam/"*.primertrimmed.rg.sorted.bam* "${final_report_dir}/${virus}/" || true
    cp -r "$RESULTS_DIR/$virus/nextclade" "${final_report_dir}/${virus}/" || true
    cp -r "$RESULTS_DIR/$virus/nf-covflow/report/"* "${final_report_dir}/${virus}/" || true
done

# Copy QC reports
cp "$RESULTS_DIR/nf-qcflow/report/reads_nanopore.qc_report.csv" "$final_report_dir/" || true
cp "$RESULTS_DIR/nf-qcflow/report/reads_nanopore.topmatches.csv" "$final_report_dir/" || true
cp -r "$RESULTS_DIR/mash_screen" "$final_report_dir/" || true

# Generate master summary
conda deactivate || true
cd "$final_report_dir" || error_exit "Failed to cd to $final_report_dir"
run_cmd python "$SCRIPT_DIR/make_summary_report.py" \
    --qc reads_nanopore.qc_report.csv \
    --consensusA rsvA/all_consensus_stats.tsv \
    --depthA rsvA/chromosome_coverage_depth_summary.tsv \
    --nextcladeA rsvA/nextclade/nextclade.tsv \
    --consensusB rsvB/all_consensus_stats.tsv \
    --depthB rsvB/chromosome_coverage_depth_summary.tsv \
    --nextcladeB rsvB/nextclade/nextclade.tsv \
    --samplesheetA ./samplesheet_rsvA.csv \
    --samplesheetB ./samplesheet_rsvB.csv \
    --out rsv_master.tsv

# ==========================================================
# Cleanup
# ==========================================================

log "=== PIPELINE COMPLETED SUCCESSFULLY ==="
log "Pipeline version: $VERSION"
log "Results directory: $RESULTS_DIR"
log "Finished at: $(date)"

exit 0
