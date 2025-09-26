#!/usr/bin/env python
import sys
import argparse
import re
import os

__author__ = 'Xiaoli Dong'
__version__ = '0.1.0'
__maintainer__ = 'Xiaoli Dong'
__email__ = 'xiaoli.dong@albertaprecisionlabs.ca'
__status__ = 'Dev'

def count_ambiguous_bases(sequence):
    ambiguous_bases = "RYSMKWBVDHNX"
    return sum(1 for base in sequence.upper() if base in ambiguous_bases)

def remove_whitespace_regex(input_string):
    return re.sub(r'\s+', '', input_string)

def determine_subtype(header_line):
    header_lower = header_line.lower()
    if re.search(r'rsv/?a|virus a|rsva|subgroup a', header_lower):
        return "RSV_A"
    elif re.search(r'rsv/?b|virus b|rsvb|subgroup b', header_lower):
        return "RSV_B"
    else:
        return "UNKNOWN"

def validate_sequences(file_path, output_file, seq_minlen, seq_maxlen, maxambigs):
    seqids_above_cutoff = []
    include_pattern = re.compile(r'rsv/A|rsv/B|virus A|virus B|rsva|rsvb|Subgroup', re.IGNORECASE)

    with open(file_path, 'r') as input_file, open(output_file, 'w') as output:
        sequence = ""
        header = None

        for line in input_file:
            if line.startswith('>'):
                # Process previous record
                if header:
                    sequence = remove_whitespace_regex(sequence)
                    seqlen = len(sequence)
                    ambiguous_count = count_ambiguous_bases(sequence)

                    if (seq_minlen <= seqlen <= seq_maxlen) and (ambiguous_count <= maxambigs):
                        subtype = determine_subtype(header)
                        raw_id = header.split()[0][1:]  # remove '>' and get accession
                        rest = " ".join(header.split()[1:])
                        output.write(f">{raw_id}|{subtype} len={seqlen} {rest}\n")
                        output.write(f"{sequence}\n")
                        seqids_above_cutoff.append(raw_id)
                    else:
                        print(f"validate_sequences:: {header}", file=sys.stderr)
                        print(f"validate_sequences:: sequence length: {seqlen} is outside the acceptable range [{seq_minlen}, {seq_maxlen}] or ambiguous bases: {ambiguous_count} > {maxambigs}", file=sys.stderr)

                header = line.strip()
                if not include_pattern.search(header):
                    print(f"validate_sequences:: Skipping {header} (does not match target pattern)", file=sys.stderr)
                    header = None
                    sequence = ""
                    continue

                sequence = ""

            else:
                sequence += line.strip()

        # Process final record
        if header and include_pattern.search(header):
            sequence = remove_whitespace_regex(sequence)
            seqlen = len(sequence)
            ambiguous_count = count_ambiguous_bases(sequence)
            if (seq_minlen <= seqlen <= seq_maxlen) and (ambiguous_count <= maxambigs):
                subtype = determine_subtype(header)
                raw_id = header.split()[0][1:]
                rest = " ".join(header.split()[1:])
                output.write(f">{raw_id}|{subtype} len={seqlen} {rest}\n")
                output.write(f"{sequence}\n")
                seqids_above_cutoff.append(raw_id)
            else:
                print(f"validate_sequences:: {header}", file=sys.stderr)
                print(f"validate_sequences:: sequence length: {seqlen} is outside the acceptable range [{seq_minlen}, {seq_maxlen}] or ambiguous bases: {ambiguous_count} > {maxambigs}", file=sys.stderr)

    return seqids_above_cutoff

def main():
    script_name = os.path.basename(__file__)
    parser = argparse.ArgumentParser(
        description=f"Filter RSV sequences by length, ambiguity, and match header patterns.\n\nExample:\npython {script_name} -f input.fasta -o filtered.fasta --minlen 15000 --maxlen 16000 --maxambigs 5",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-f", "--fasta", help="Input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Output file for filtered sequences", required=True)
    parser.add_argument("--minlen", type=int, default=10000, help="Minimum sequence length (default: 10000)")
    parser.add_argument("--maxlen", type=int, default=20000, help="Maximum sequence length (default: 20000)")
    parser.add_argument("--maxambigs", type=int, default=0, help="Maximum number of ambiguous bases allowed (default: 0)")

    args = parser.parse_args()
    validate_sequences(args.fasta, args.output, args.minlen, args.maxlen, args.maxambigs)

if __name__ == "__main__":
    main()
