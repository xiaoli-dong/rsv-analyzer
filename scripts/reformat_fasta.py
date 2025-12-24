#!/usr/bin/env python3

import argparse
import sys

def reformat_stream(in_stream, out_stream):
    for line in in_stream:
        if line.startswith(">"):
            header = line[1:].strip()
            parts = header.split()

            if len(parts) == 2:
                out_stream.write(f">{parts[0]}|ref_{parts[1]}\n")
            else:
                out_stream.write(line)
        else:
            out_stream.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Reformat FASTA headers from 'ID ACC' to 'ID|ref_ACC'",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""Examples:
  Read from a file and write to a file:
    python reformat_fasta.py -i input.fasta -o output.fasta

  Read from stdin (pipe) and write to stdout:
    cat *.fasta | python reformat_fasta.py > reformatted.fasta

  Read from stdin and write to a file:
    cat *.fasta | python reformat_fasta.py -o output.fasta
"""
    )

    parser.add_argument(
        "-i", "--input",
        help="Input FASTA file (default: stdin)"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output FASTA file (default: stdout)"
    )

    args = parser.parse_args()

    infile = open(args.input, "r") if args.input else sys.stdin
    outfile = open(args.output, "w") if args.output else sys.stdout

    try:
        reformat_stream(infile, outfile)
    finally:
        if args.input:
            infile.close()
        if args.output:
            outfile.close()

if __name__ == "__main__":
    main()
