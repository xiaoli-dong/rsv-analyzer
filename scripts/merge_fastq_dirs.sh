#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <dir1> <dir2> <outdir> <samplesheet>" >&2
    exit 1
fi

DIR1="$1"
DIR2="$2"
OUTDIR="$3"
SAMPLESHEET="$4"

mkdir -p "$OUTDIR"
mkdir -p "$(dirname "$SAMPLESHEET")"

# Create samplesheet header
echo "sample,fastq_1,fastq_2,long_fastq" > "$SAMPLESHEET"

# Loop through files in DIR1
for f in "$DIR1"/*.fastq.gz; do
    base=$(basename "$f")
    sample="${base%%.*}"  # everything before first dot
    f2="$DIR2/$base"

    # Merge if exists in DIR2
    if [[ -f "$f2" ]]; then
        echo "Merging $base"
        zcat "$f" "$f2" | gzip > "$OUTDIR/$base"
    else
        echo "Copying $base (no match in DIR2)"
        cp "$f" "$OUTDIR/$base"
    fi

    # Add row to samplesheet
    echo "${sample},NA,NA,${OUTDIR}/${base}" >> "$SAMPLESHEET"
done
