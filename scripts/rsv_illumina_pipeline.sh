#!/bin/bash
set -euo pipefail

# ==========================================================
# Paths & constants
# ==========================================================
readonly prod_prog_base="/nfs/APL_Genomics/apps/production"
readonly deve_prog_base="/nfs/Genomics_DEV/projects/xdong/deve"

readonly path_to_viralrecon="${prod_prog_base}/viralrecon"
readonly path_to_qc_pipeline="${deve_prog_base}/nf-qcflow"
readonly path_to_covflow="${deve_prog_base}/nf-covflow"

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly RSV_PROG_BASE="${SCRIPT_DIR}/.."
readonly VERSION="$(cat "$SCRIPT_DIR/VERSION" 2>/dev/null || echo "unknown")"
readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
readonly RSV_SLURM_CONFIG="${RSV_PROG_BASE}/conf/slurm.config"
readonly ENV_NAME="virus_env"

# Genome references
readonly rsvA_ref="${RSV_PROG_BASE}/resource/primerschemes/rsvA/v3/reference.fasta"
readonly rsvA_bed="${RSV_PROG_BASE}/resource/primerschemes/rsvA/v3/scheme.bed"
readonly rsvB_ref="${RSV_PROG_BASE}/resource/primerschemes/rsvB/v3/reference.fasta"
readonly rsvB_bed="${RSV_PROG_BASE}/resource/primerschemes/rsvB/v3/scheme.bed"
readonly MASH_DB="${RSV_PROG_BASE}/resource/db/mash_screen/sequences.msh"

# Viralrecon config (default + override)
readonly DEFAULT_QCFLOW_CONFIG_FILE="${RSV_PROG_BASE}/conf/qcflow.config"
QCFLOW_CONFIG_FILE="$DEFAULT_QCFLOW_CONFIG_FILE"

# Viralrecon config (default + override)
readonly DEFAULT_VIRALRECON_CONFIG_FILE="${RSV_PROG_BASE}/conf/viralrecon.config"
VIRALRECON_CONFIG_FILE="$DEFAULT_VIRALRECON_CONFIG_FILE"

# ==========================================================
# Help / Version
# ==========================================================
show_version() {
    echo "$SCRIPT_NAME version $VERSION"
}

show_help() {
    cat << EOF
$0 - RSV Illumina analysis pipeline

USAGE:
 bash $0 <samplesheet.csv> <results_dir> [options]

REQUIRED:
  samplesheet.csv       Input samplesheet (CSV)
  results_dir           Output directory for pipeline results

Options:
  -h, --help
  -v, --version
  --qcflow-config FILE   Custom qcflow config
  --viralrecon-config FILE   Custom viralrecon config

EXAMPLES:
  sh $0 samplesheet.csv results_2025_01_16
  sh $0 samplesheet.csv results_2025_01_16 \
      --qcflow-config qcflow.config \
      --viralrecon-config viralrecon.config

CHECK HELP WITHOUT SUBMITTING:
  bash $0 --help

NOTES:
  • This script is intended to be run in a SLURM environment
EOF
}

# Allow help without sbatch
case "${1:-}" in
    -h|--help)
        show_help
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
VIRALRECON_CONFIG=""

# ==========================================================
# Parse named options
# ==========================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --qcflow-config)
            QCFLOW_CONFIG="${2:-}"
            shift 2
            ;;
        --viralrecon-config)
            VIRALRECON_CONFIG="${2:-}"
            shift 2
            ;;
        -h|--help)
            show_help
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
VIRALRECON_CONFIG="$(validate_config "$VIRALRECON_CONFIG" "Viralrecon")"

# Resolve absolute paths
SAMPLESHEET="$(cd "$(dirname "$SAMPLESHEET")" && pwd)/$(basename "$SAMPLESHEET")"
RESULTS_DIR="$(mkdir -p "$RESULTS_DIR" && cd "$RESULTS_DIR" && pwd)"

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

log "Input samplesheet      : $SAMPLESHEET"
log "Results directory      : $RESULTS_DIR"
log "QCflow config file     : $QCFLOW_CONFIG_FILE"
log "Viralrecon config file : $VIRALRECON_CONFIG_FILE"

# ==========================================================
# Conda
# ==========================================================
check_conda_env "$ENV_NAME" || error_exit "Conda env '$ENV_NAME' not found"
activate_env "$ENV_NAME"

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
  --platform illumina \
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
# STEP 4: Process virus
# ==========================================================
process_virus() {
    local virus="$1" ref="$2" bed="$3"
    local ss="$RESULTS_DIR/samplesheet_${virus}.csv"

    [[ -s "$ss" ]] || { log "No $virus samples found, skipping"; return; }

    viralrecon_outdir="$RESULTS_DIR/$virus/viralrecon"

    run_cmd nextflow run "${path_to_viralrecon}/main.nf" \
      -profile singularity,slurm \
      -c "$RSV_SLURM_CONFIG" \
      -c "$VIRALRECON_CONFIG_FILE" \
      --input "$ss" \
      --outdir "${viralrecon_outdir}" \
      --protocol amplicon \
      --platform illumina \
      --primer_bed "$bed" \
      --fasta "$ref" \
      -resume

    run_cmd cat "${viralrecon_outdir}/variants/ivar/consensus/bcftools/"*.consensus.fa | \
        python "$SCRIPT_DIR/reformat_fasta.py" \
        -o "${viralrecon_outdir}/variants/ivar/consensus/bcftools/all_consensus.${virus}.fa"



    run_cmd python "$SCRIPT_DIR/fasta_stats.py" \
      "${viralrecon_outdir}/variants/ivar/consensus/bcftools/all_consensus.${virus}.fa" \
      -o "${viralrecon_outdir}/variants/ivar/consensus/bcftools/all_consensus.${virus}_stats.tsv"

    # run_cmd python "$SCRIPT_DIR/consensus_stats.py" \
    #   "${viralrecon_outdir}/variants/ivar/consensus/bcftools/" \
    #   "${viralrecon_outdir}/variants/ivar/consensus/bcftools/all_consensus_stats.tsv"

    mkdir -p "${viralrecon_outdir}/bam_to_covflow"
    cp "${viralrecon_outdir}/variants/bowtie2/"*.ivar_trim.sorted.bam* \
       "${viralrecon_outdir}/bam_to_covflow"

    run_cmd python "$SCRIPT_DIR/build_samplesheet_for_covflow.py" \
      --input_dir "${viralrecon_outdir}/bam_to_covflow" \
      --ref_fasta "$ref" \
      --bed_file "$bed" \
      -o "$RESULTS_DIR/samplesheet_to_covflow_${virus}.csv"

    run_cmd nextflow run "${path_to_covflow}/main.nf" \
      -profile singularity,slurm \
      --input "$RESULTS_DIR/samplesheet_to_covflow_${virus}.csv" \
      --outdir "$RESULTS_DIR/$virus/nf-covflow" \
      -resume

    last_char="${virus: -1}"
    nextclade_db="./nextstrain/rsv/${last_char,,}"

    run_cmd nextclade dataset get \
      --name "nextstrain/rsv/${last_char,,}" \
      --output-dir "$nextclade_db"

    run_cmd nextclade run \
      "$RESULTS_DIR/$virus/viralrecon/variants/ivar/consensus/bcftools/"*.fa \
      --input-dataset "$nextclade_db" \
      --output-all "$RESULTS_DIR/$virus/nextclade"
}

process_virus rsvA "$rsvA_ref" "$rsvA_bed"
process_virus rsvB "$rsvB_ref" "$rsvB_bed"

# ==========================================================
# STEP 6: Final report
# ==========================================================
log "=== Generating final summary report ==="
final_report_dir="$RESULTS_DIR/summary_report"
mkdir -p "$final_report_dir"

for virus in rsvA rsvB; do
    mkdir -p "$final_report_dir/$virus"
    #cp "$RESULTS_DIR/samplesheet_${virus}.csv" "$final_report_dir/"
    viralrecon_outdir="$RESULTS_DIR/$virus/viralrecon"
    cp "${viralrecon_outdir}/variants/ivar/consensus/bcftools/"all_consensus.* "$final_report_dir/$virus/"
    #cp "${viralrecon_outdir}/variants/ivar/consensus/bcftools/"*.filtered.vcf.gz* "$final_report_dir/$virus/"
    cp "${viralrecon_outdir}/bam_to_covflow/"*.ivar_trim.sorted.bam* "$final_report_dir/$virus/"
    cp "$RESULTS_DIR/$virus/nextclade/"* "$final_report_dir/$virus/"
    cp -r "$RESULTS_DIR/$virus/nf-covflow/report/"* "$final_report_dir/$virus/"
done

cp "$RESULTS_DIR/nf-qcflow/report/reads_illumina.qc_report.csv" "$final_report_dir/"
cp "$RESULTS_DIR/nf-qcflow/report/reads_illumina.topmatches.csv" "$final_report_dir/"
cp -r "$RESULTS_DIR/mash_screen" "$final_report_dir/"

# ==========================================================
# STEP 7: make_summary_report
# ==========================================================
conda deactivate || true
cd "$final_report_dir" || error_exit "Failed to cd to summary report dir"

run_cmd python "$SCRIPT_DIR/make_summary_report.py" \
  --qc reads_illumina.qc_report.csv \
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
