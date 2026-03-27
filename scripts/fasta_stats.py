#!/usr/bin/env python3
import subprocess
import pathlib
import sys
import argparse
import textwrap

# def run_seqkit(fasta_path):
#     """Run seqkit fx2tab on a fasta file and return the parsed rows (without header)."""
#     cmd = [
#         "seqkit", "fx2tab",
#         "--only-id", "--name", "--length",
#         "-C", "ATCG",
#         "-C", "RYSWKMBDHV",
#         "-C", "N",
#         "-H", str(fasta_path)
#     ]
#     result = subprocess.run(
#         cmd,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE,
#         text=True,
#         check=True
#     )
#     return result.stdout.strip().splitlines()
def run_seqkit(fasta_path):
    cmd = [
        "seqkit", "fx2tab",
        "--only-id", "--name", "--length",
        "-C", "ATCG",
        "-C", "RYSWKMBDHV",
        "-C", "N",
        "-H", str(fasta_path)
    ]

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,   # ✅ Python 3.6 compatible
            check=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"seqkit failed: {e.stderr}")

    return result.stdout.strip().splitlines()

def main():
    example_text = textwrap.dedent("""
    Examples:
      Process multiple FASTA files:
        python fasta_stats.py sample1.fasta sample2.fasta sample3.fna

      Use a glob pattern:
        python fasta_stats.py *.fasta -o results.tsv

      Write to a custom output file:
        python fasta_stats.py sample.fasta -o my_stats.tsv

    Output columns:
      file            Input FASTA file path
      #id             Sequence ID
      length          Sequence length
      ATCG            Count of A/T/C/G bases
      RYSWKMBDHV      Count of ambiguous (non-ATCG, non-N) bases
      N               Count of N bases
      coverage        (ATCG + ambiguous) / length
      completeness    ATCG / length
    """)

    parser = argparse.ArgumentParser(
        prog="fasta_stats.py",
        description="Compute FASTA sequence statistics using seqkit (fx2tab).",
        epilog=example_text,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "fastas",
        nargs="+",
        help="Input FASTA files (.fa, .fasta, .fna). Supports multiple files or glob patterns."
    )

    parser.add_argument(
        "-o", "--output",
        default="fasta_stats.tsv",
        help="Output TSV file (default: fasta_stats.tsv)"
    )

    args = parser.parse_args()

    fasta_files = [pathlib.Path(f).resolve() for f in args.fastas]

    # Validate inputs
    for fasta in fasta_files:
        if not fasta.exists():
            print(f"ERROR: File not found: {fasta}", file=sys.stderr)
            sys.exit(1)
        if fasta.suffix.lower() not in [".fa", ".fasta", ".fna"]:
            print(f"WARNING: Skipping non-FASTA file: {fasta}", file=sys.stderr)

    output_file = pathlib.Path(args.output)

    with open(output_file, "w") as out:
        # Write header
        out.write("file\t#id\tlength\tATCG\tRYSWKMBDHV\tN\tcoverage\tcompleteness\n")

        for fasta in fasta_files:
            try:
                rows = run_seqkit(fasta)
                for row in rows[1:]:  # skip seqkit header
                    parts = row.split("\t")
                    if len(parts) < 5:
                        continue

                    seq_id = parts[0]
                    length = int(parts[1])
                    atcg = int(parts[2])
                    rysw = int(parts[3])
                    ncount = int(parts[4])

                    coverage = (atcg + rysw) / length if length > 0 else 0
                    completeness = atcg / length if length > 0 else 0

                    out.write(
                        f"{fasta}\t{seq_id}\t{length}\t{atcg}\t{rysw}\t{ncount}\t"
                        f"{coverage:.4f}\t{completeness:.4f}\n"
                    )
            except subprocess.CalledProcessError as e:
                print(f"ERROR running seqkit on {fasta}: {e.stderr}", file=sys.stderr)

    print(f"Stats written to: {output_file}")

if __name__ == "__main__":
    main()
