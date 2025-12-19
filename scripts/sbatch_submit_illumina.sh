#!/bin/bash
#SBATCH --job-name=rsv-analyzer-illumina
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=vm-cpu,big-ram,huge-ram
#SBATCH --time=3-00:00:00
#SBATCH --mem=24G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

# ==========================================================
# Help & usage
# ==========================================================

show_help() {
    cat << EOF
RSV Illumina SLURM submission script

USAGE:
 sbatch $0 <samplesheet.csv> <results_dir> [options]

REQUIRED:
  samplesheet.csv       Input samplesheet (CSV)
  results_dir           Output directory for pipeline results

OPTIONS:
  --qcflow-config FILE     Custom nf-qcflow config file
  --viralrecon-config FILE Custom viralrecon config file
  -h, --help               Show this help

EXAMPLES:
  sbatch $0 samplesheet.csv results_2025_01_16
  sbatch $0 samplesheet.csv results_2025_01_16 \
      --qcflow-config qcflow.config \
      --viralrecon-config viralrecon.config

CHECK HELP WITHOUT SUBMITTING:
  bash $0 --help

NOTES:
  • This script must be submitted with sbatch
  • SLURM stdout/stderr will be written to:
      slurm-<jobid>.out
      slurm-<jobid>.err
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
# Logging
# ==========================================================

echo "=================================================="
echo "SLURM job ID        : ${SLURM_JOB_ID:-N/A}"
echo "Job name            : ${SLURM_JOB_NAME:-N/A}"
echo "Node(s)             : ${SLURM_NODELIST:-N/A}"
echo "CPUs per task       : ${SLURM_CPUS_PER_TASK:-N/A}"
echo "Samplesheet         : $SAMPLESHEET"
echo "Results directory   : $RESULTS_DIR"
[[ -n "$QCFLOW_CONFIG" ]] && echo "QCflow config       : $QCFLOW_CONFIG"
[[ -n "$VIRALRECON_CONFIG" ]] && echo "Viralrecon config   : $VIRALRECON_CONFIG"
echo "Start time          : $(date)"
echo "=================================================="

# ==========================================================
# Run pipeline
# ==========================================================

PIPELINE="/nfs/APL_Genomics/virus_rsv/pipelines/rsv-analyzer/scripts/rsv_illumina_pipeline.sh"

if [[ ! -x "$PIPELINE" ]]; then
    echo "ERROR: Pipeline script not executable: $PIPELINE"
    exit 1
fi

CMD=(bash "$PIPELINE" "$SAMPLESHEET" "$RESULTS_DIR")

[[ -n "$QCFLOW_CONFIG" ]] && CMD+=(--qcflow-config "$QCFLOW_CONFIG")
[[ -n "$VIRALRECON_CONFIG" ]] && CMD+=(--viralrecon-config "$VIRALRECON_CONFIG")

echo "Running pipeline command:"
printf '  %q' "${CMD[@]}"
echo

"${CMD[@]}"
exit_code=$?

# ==========================================================
# Finish
# ==========================================================

echo "=================================================="
echo "Pipeline exit code  : $exit_code"
echo "End time            : $(date)"
echo "=================================================="

exit $exit_code
