#!/usr/bin/env python3
import os
import csv
import argparse
from collections import defaultdict
import textwrap

def main():
    example_text = textwrap.dedent("""
    Example:
      python build_samplesheet_for_covflow.py \\
          -i /path/to/bam_dir \\
          -o samplesheet.csv \\
          --ref_fasta reference.fasta \\
          --bed_file scheme.bed

    Sample output (samplesheet.csv):
      sample,bam,ref_fasta,bed_file
      250115_S_I_012-S2,/path/to/250115_S_I_012-S2.ivar_trim.sorted.bam,reference.fasta,scheme.bed
    """)

    parser = argparse.ArgumentParser(
        description="Generate a sample sheet CSV from BAM files in a directory.",
        epilog=example_text,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Directory containing BAM files"
    )

    parser.add_argument(
        "-o", "--output",
        default="samplesheet.csv",
        help="Output CSV file (default: samplesheet.csv)"
    )

    parser.add_argument(
        "--ref_fasta",
        required=True,
        help="Reference FASTA file (will be added to each row)"
    )

    parser.add_argument(
        "--bed_file",
        required=True,
        help="BED file (will be added to each row)"
    )

    parser.add_argument(
        "--delimiter",
        default=".",
        help="Delimiter used to split sample name from BAM filename (default: '.')"
    )

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.input_dir):
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")

    # List BAM files
    bam_files = sorted(f for f in os.listdir(args.input_dir) if f.endswith(".bam"))
    if not bam_files:
        print(f"No BAM files found in {args.input_dir}")
        return

    # Group BAM files by sample
    samples = defaultdict(list)
    for bam in bam_files:
        sample = bam.split(args.delimiter)[0]
        bam_path = os.path.join(args.input_dir, bam)
        samples[sample].append(bam_path)

    # Write sample sheet
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample", "bam", "ref_fasta", "bed_file"])

        for sample, bams in sorted(samples.items()):
            bam_list = ",".join(sorted(bams))
            writer.writerow([sample, bam_list, args.ref_fasta, args.bed_file])

    print(f"Sample sheet written to: {args.output}")


if __name__ == "__main__":
    main()
