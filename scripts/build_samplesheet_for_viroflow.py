#!/usr/bin/env python3
import argparse
import csv
import textwrap
import sys

def main():
    example_text = textwrap.dedent("""
    Example:
      python build_samplesheet_for_viroflow.py \\
          -i input_samplesheet.csv \\
          -o output_samplesheet.csv \\
          --ref reference.fasta \\
          --bed scheme.bed

    Input format:
      sample,fastq_1,fastq_2
      sampleA,/path/to_R1.fastq.gz,/path/to_R2.fastq.gz

    Sample output:
      sample,fastq_1,fastq_2,ref_genome,bed_file
      250116_S_N_012-barcode02,/path/to.fastq.gz,,reference.fasta,scheme.bed
    """)

    parser = argparse.ArgumentParser(
        prog="build_samplesheet_for_viroflow.py",
        description="Build a Viroflow-compatible samplesheet by adding reference genome and BED file columns.",
        epilog=example_text,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV (must contain: sample, fastq_1, fastq_2)"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file"
    )

    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA file (added to all rows)"
    )

    parser.add_argument(
        "--bed",
        required=True,
        help="BED file (added to all rows)"
    )

    args = parser.parse_args()

    fieldnames = ["sample", "fastq_1", "fastq_2", "ref_genome", "bed_file"]

    with open(args.input) as infile, open(args.output, "w", newline="") as outfile:
        reader = csv.DictReader(infile)

        # Validate required columns
        required_cols = {"sample", "fastq_1", "fastq_2"}
        if not reader.fieldnames or not required_cols.issubset(reader.fieldnames):
            print(
                f"Error: Input file must contain columns: {', '.join(required_cols)}",
                file=sys.stderr
            )
            sys.exit(1)

        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            writer.writerow({
                "sample": row.get("sample", ""),
                "fastq_1": row.get("fastq_1", ""),
                "fastq_2": row.get("fastq_2", ""),
                "ref_genome": args.ref,
                "bed_file": args.bed
            })

    print(f"Sample sheet written to: {args.output}")


if __name__ == "__main__":
    main()
