#!/usr/bin/env python3
import subprocess
import pathlib
import sys

def run_seqkit(fasta_path):
    """Run seqkit fx2tab on a fasta file and return the parsed rows (without header)."""
    cmd = [
        "seqkit", "fx2tab",
        "--only-id", "--name", "--length",
        "-C", "ATCG",
        "-C", "RYSWKMBDHV",
        "-C", "N",
        "-H", str(fasta_path)
    ]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    return result.stdout.strip().splitlines()

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input_dir> [output_tsv]", file=sys.stderr)
        sys.exit(1)

    input_dir = pathlib.Path(sys.argv[1]).resolve()
    output_file = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else pathlib.Path("fasta_stats.tsv")

    # Match all fasta extensions
    fasta_files = [f for f in input_dir.rglob("*") if f.suffix.lower() in [".fa", ".fasta", ".fna"]]

    if not fasta_files:
        print(f"ERROR: No FASTA files found under {input_dir}", file=sys.stderr)
        sys.exit(1)

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

    print(f"✅ Stats written to {output_file}")

if __name__ == "__main__":
    main()
