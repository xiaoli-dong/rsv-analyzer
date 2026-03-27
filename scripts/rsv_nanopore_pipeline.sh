#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, and pipe failures

# ==========================================================
# Paths & constants
# ==========================================================
readonly prod_prog_base="/nfs/APL_Genomics/apps/production"
readonly deve_prog_base="/nfs/Genomics_DEV/projects/xdong/deve"
#readonly path_to_viroflow="${deve_prog_base}/nf-viroflow"
readonly path_to_viroflow="${prod_prog_base}/viroflow_pipeline/nf-viroflow"
readonly path_to_qc_pipeline="${prod_prog_base}/qcflow_pipeline/nf-qcflow"
#readonly path_to_covflow="${deve_prog_base}/nf-covflow"
readonly path_to_covflow="${prod_prog_base}/covflow_pipeline/nf-covflow"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly RSV_PROG_BASE="${SCRIPT_DIR}/.."
readonly VERSION="$(cat "$SCRIPT_DIR/VERSION" 2>/dev/null || echo "unknown")"
readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
readonly RSV_SLURM_CONFIG="${RSV_PROG_BASE}/conf/slurm.config"

# # virorecon requires 6 columns bed file while align_trim in viralassembly and viroflow requires 7 columns bed file
readonly rsvA_ref="${RSV_PROG_BASE}/resource/primerschemes_plsa/rsvA/v3/reference.fasta"
readonly rsvA_bed="${RSV_PROG_BASE}/resource/primerschemes_plsa/rsvA/v3/primer.bed"
readonly rsvB_ref="${RSV_PROG_BASE}/resource/primerschemes_plsa/rsvB/v3/reference.fasta"
readonly rsvB_bed="${RSV_PROG_BASE}/resource/primerschemes_plsa/rsvB/v3/primer.bed"
readonly MASH_DB="${RSV_PROG_BASE}/resource/db/mash_screen/sequences.msh"

# nf-qcflow config (default + override)
DEFAULT_QCFLOW_CONFIG_FILE="${RSV_PROG_BASE}/conf/qcflow.config"
QCFLOW_CONFIG_FILE="$DEFAULT_QCFLOW_CONFIG_FILE"

# nf-viroflow config (default + override)
DEFAULT_VIROFLOW_CONFIG_FILE="${RSV_PROG_BASE}/conf/viroflow.config"
VIROFLOW_CONFIG_FILE="$DEFAULT_VIROFLOW_CONFIG_FILE"

# ==========================================================
# Help / Version
# ==========================================================
show_version() {
    echo "$SCRIPT_NAME version $VERSION"
}

show_help() {
    cat << EOF
$0 - RSV Nanopore Analysis Pipeline

Usage:
    bash $0 <samplesheet.csv> <results_dir> [options]

REQUIRED:
    samplesheet.csv       Input samplesheet (CSV)
    results_dir           Output directory for pipeline results

Options:
    -h, --help
    -v, --version
    --qcflow-config FILE   Custom qcflow config
    --viroflow-config FILE   Custom viroflow config

EXMPLES:
    bash $0 samplesheet.csv results_dir

    bash $0 samplesheet.csv results_dir \
        --qcflow-config qcflow.config \
        --viroflow-config viroflow.config

CHECK HELP WITHOUT RUNNING:
    bash $0 --help

CHECK VERSION WITHOUT SUBMITTING:
    bash $0 --version


EOF
}


# Allow help without sbatch
case "${1:-}" in
    -h|--help)
        show_help
        exit 0
        ;;
    -v|--version)
        show_version
        exit 0
    ;;
esac

# ==========================================================
# Required arguments
# ==========================================================

if [[ $# -lt 2 ]]; then
    echo "ERROR: Missing required arguments"
    show_help
    exit 1
fi

SAMPLESHEET="$1"
RESULTS_DIR="$2"
shift 2

# Optional named arguments
QCFLOW_CONFIG=""
VIROFLOW_CONFIG=""

# ==========================================================
# Parse named options
# ==========================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --qcflow-config)
            QCFLOW_CONFIG="${2:-}"
            shift 2
            ;;
        --viroflow-config)
            VIROFLOW_CONFIG="${2:-}"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        -v|--version)
            show_version
            exit 0
        ;;
        *)
            echo "ERROR: Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# ==========================================================
# Validation
# ==========================================================

if [[ ! -s "$SAMPLESHEET" ]]; then
    echo "ERROR: Samplesheet missing or empty: $SAMPLESHEET"
    exit 1
fi

validate_config() {
    local cfg="$1"
    local name="$2"

    if [[ -n "$cfg" ]]; then
        if [[ ! -f "$cfg" ]]; then
            echo "ERROR: $name config not found: $cfg"
            exit 1
        fi
        if [[ ! -s "$cfg" ]]; then
            echo "WARNING: $name config is empty: $cfg"
        fi
        echo "$(cd "$(dirname "$cfg")" && pwd)/$(basename "$cfg")"
    fi
}

QCFLOW_CONFIG="$(validate_config "$QCFLOW_CONFIG" "QCflow")"
VIROFLOW_CONFIG="$(validate_config "$VIROFLOW_CONFIG" "Viroflow")"

# After validation
if [[ -n "$VIROFLOW_CONFIG" ]]; then
    VIROFLOW_CONFIG_FILE="$VIROFLOW_CONFIG"
fi

if [[ -n "$QCFLOW_CONFIG" ]]; then
    QCFLOW_CONFIG_FILE="$QCFLOW_CONFIG"
fi


# Resolve absolute paths
SAMPLESHEET="$(cd "$(dirname "$SAMPLESHEET")" && pwd)/$(basename "$SAMPLESHEET")"
RESULTS_DIR="$(mkdir -p "$RESULTS_DIR" && cd "$RESULTS_DIR" && pwd)"

# ==========================================================
# Logging & helpers
# ==========================================================
LOGFILE="$RESULTS_DIR/pipeline.log"

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [$SCRIPT_NAME-$VERSION] $*"
    echo "$msg" >> "$LOGFILE"
}

error_exit() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [$SCRIPT_NAME-$VERSION] ERROR: $1"
    echo "$msg" >&2
    exit "${2:-1}"
}

run_cmd() {
    log "Running: $*"

    {
        echo "----- CMD START: $(date) -----"
        echo "$*"
    } >> "$LOGFILE"

    "$@"
    local status=$?

    {
        echo "----- CMD END: $(date) [exit=$status] -----"
    } >> "$LOGFILE"

    [[ $status -eq 0 ]] || error_exit "Command failed: $*"
}

# ==========================================================
# Validation
# ==========================================================

log "Input samplesheet      : $SAMPLESHEET"
log "Results directory      : $RESULTS_DIR"
log "QCflow config file     : $QCFLOW_CONFIG_FILE"
log "Viroflow config file : $VIROFLOW_CONFIG_FILE"

# ==========================================================
# STEP 1: Prepare samplesheet
# ==========================================================
log "=== STEP 1: Preparing samplesheet ==="

ss="$RESULTS_DIR/samplesheet_to_qc.csv"
cp "$SAMPLESHEET" "$ss"

[[ -s "$ss" ]] || { log "No $virus samples found, skipping"; return; }

# ==========================================================
# STEP 2: QC pipeline
# ==========================================================
log "=== STEP 2: Running QC pipeline ==="
run_cmd nextflow run "${path_to_qc_pipeline}/main.nf" \
    -profile singularity,slurm \
    -c "$QCFLOW_CONFIG_FILE" \
    --input "$ss" \
    --platform nanopore \
    --outdir "$RESULTS_DIR/nf-qcflow" \
    -resume

# ==========================================================
# STEP 3: RSV A/B classification
# produced samplesheet "sample,fastq_1,fastq_2"
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

    # # Declare associative arrays
    # declare -A ref_plsa=(
    #     [rsvA]="$rsvA_ref"
    #     [rsvB]="$rsvB_ref"
    # )

    # declare -A bed_plsa=(
    #     [rsvA]="$rsvA_bed"
    #     [rsvB]="$rsvB_bed"
    # )

    # # Lookup
    # ref="${ref_plsa[$virus]}"
    # bed="${bed_plsa[$virus]}"

    # Validate
    # if [[ -z "$ref" || -z "$bed" ]]; then
    #     echo "ERROR: Invalid virus type or undefined variables: $virus"
    #     exit 1
    # fi

    # Run
    run_cmd python "$SCRIPT_DIR/build_samplesheet_for_viroflow.py" \
        -i "$RESULTS_DIR/samplesheet_${virus}.csv" \
        --ref "$ref" \
        --bed "$bed" \
        -o "$RESULTS_DIR/samplesheet_${virus}_to_viroflow.csv"



    viroflow_outdir="$RESULTS_DIR/$virus/viroflow"

    log "=== Processing $virus ==="

    # Viral assembly
    run_cmd nextflow run "${path_to_viroflow}/main.nf" \
        -profile singularity,slurm \
        -c "$RSV_SLURM_CONFIG" \
        -c "$VIROFLOW_CONFIG_FILE" \
        --input "$RESULTS_DIR/samplesheet_${virus}_to_viroflow.csv" \
        --outdir "${viroflow_outdir}" \
        --input_type amplicon \
        -resume

    run_cmd cat "${viroflow_outdir}/consensus/"*.filtered_consensus.fasta \
        > "${viroflow_outdir}/consensus/all_consensus.${virus}.fasta"

    # Consensus stats
    #mkdir -p "$RESULTS_DIR/$virus/consensus"
    run_cmd  python "$SCRIPT_DIR/fasta_stats.py" \
        "${viroflow_outdir}/consensus/all_consensus.${virus}.fasta" \
        -o "${viroflow_outdir}/consensus/all_consensus.${virus}_stats.tsv"

    # Prepare for nf-covflow
    run_cmd mkdir -p "${viroflow_outdir}/bam_to_covflow"
    run_cmd cp "${viroflow_outdir}/bam/"*.primertrimmed.sorted.bam* \
        "${viroflow_outdir}/bam_to_covflow" || true

    # Build covflow samplesheet
    run_cmd python "$SCRIPT_DIR/build_samplesheet_for_covflow.py" \
        --input_dir "${viroflow_outdir}/bam_to_covflow" \
        --ref_fasta "$ref" \
        --bed_file "$bed" \
        -o "$RESULTS_DIR/samplesheet_to_covflow_${virus}.csv"

    # Run covflow
    run_cmd nextflow run "${path_to_covflow}/main.nf" \
        -profile singularity,slurm \
        --input "$RESULTS_DIR/samplesheet_to_covflow_${virus}.csv" \
        --outdir "$RESULTS_DIR/$virus/nf-covflow" \
        -resume

    # Nextclade
    last_char="${virus: -1}"
    nextclade_db="./nextstrain/rsv/${last_char,,}"
    run_cmd nextclade dataset get \
      --name "nextstrain/rsv/${last_char,,}" \
      --output-dir "$nextclade_db"

    if ls "${viroflow_outdir}/consensus/"*.fasta >/dev/null 2>&1; then
        run_cmd nextclade run "${viroflow_outdir}/consensus/"*.fasta \
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
run_cmd mkdir -p "$final_report_dir"


for virus in rsvA rsvB; do
    run_cmd mkdir -p "${final_report_dir}/${virus}"
    run_cmd mkdir -p "$final_report_dir/$virus/nextclade"

    viroflow_outdir="$RESULTS_DIR/$virus/viroflow"

    #cp "$RESULTS_DIR/samplesheet_${virus}.csv" "${final_report_dir}/"
    run_cmd cp "${viroflow_outdir}/consensus/"all_consensus.* \
        "${final_report_dir}/${virus}/" || true

    run_cmd cp "${viroflow_outdir}/bam/"*.primertrimmed.sorted.bam* \
        "${final_report_dir}/${virus}/" || true

    for ext in csv gff json ndjson tsv nwk tbl; do
        file="$RESULTS_DIR/$virus/nextclade/nextclade.$ext"
        [ -f "$file" ] && cp "$file" "$final_report_dir/$virus/nextclade"
    done

    run_cmd cp -r "$RESULTS_DIR/$virus/nf-covflow/report/"* \
        "${final_report_dir}/${virus}/" || true
done

# Copy QC reports
run_cmd cp "$RESULTS_DIR/nf-qcflow/report/reads_nanopore.qc_report.csv" \
    "$final_report_dir/" || true

run_cmd cp "$RESULTS_DIR/nf-qcflow/report/reads_nanopore.topmatches.csv" \
    "$final_report_dir/" || true

run_cmd cp -r "$RESULTS_DIR/mash_screen" "$final_report_dir/" || true

# ==========================================================
# STEP 6: make_summary_report
# ==========================================================

run_cmd cd "$final_report_dir" || error_exit "Failed to cd to $final_report_dir"

run_cmd python "$SCRIPT_DIR/make_summary_report.py" \
    --qc reads_nanopore.qc_report.csv \
    --consensusA rsvA/all_consensus.rsvA_stats.tsv \
    --depthA rsvA/chromosome_coverage_depth_summary.tsv \
    --nextcladeA rsvA/nextclade.tsv \
    --consensusB rsvB/all_consensus.rsvB_stats.tsv \
    --depthB rsvB/chromosome_coverage_depth_summary.tsv \
    --nextcladeB rsvB/nextclade.tsv \
    --samplesheetA ../samplesheet_rsvA.csv \
    --samplesheetB ../samplesheet_rsvB.csv \
    --out rsv_master.tsv

# ==========================================================
# Cleanup
# ==========================================================

log "=== PIPELINE COMPLETED SUCCESSFULLY ==="
log "Pipeline version: $VERSION"
log "Results directory: $RESULTS_DIR"
log "Finished at: $(date)"

exit 0
