#!/usr/bin/env python
import sys
import argparse
from collections import defaultdict
import re
from csv import DictReader
import os

__author__ = 'Xiaoli Dong'
__version__ = '0.1.0'
__maintainer__ = 'Xiaoli Dong'
__email__ = 'xiaoli.dong@albertaprecisionlabs.ca'
__status__ = 'Dev'

def count_ambiguous_bases(sequence):
    ambiguous_bases = "RYSMKWBVDHNX"  # Set of ambiguous base codes
    count = 0

    # Loop through each base in the sequence and check if it's ambiguous
    for base in sequence.upper():  # Convert to uppercase to handle both cases
        if base in ambiguous_bases:
            count += 1
    return count
def remove_whitespace_regex(input_string):
    return re.sub(r'\s+', '', input_string)  # '\s+' matches all whitespace characters

"""
Filter fasta:
    seq len < seq_minlen
    seq len > seq_maxlen
    The ambiguous number is > maxambigs
    add seq length to the end of header
"""
def validate_sequences(file_path, output_file, seq_minlen, seq_maxlen, maxambigs):
    seqids_above_cutoff = []
    with open(file_path, 'r') as input_file, open(output_file, 'w') as output:
        sequence = ""
        header = None
        for line in input_file:
            if line.startswith('>'):
                if header:
                    sequence = remove_whitespace_regex(sequence)
                    seqlen = len(sequence)
                    ambiguous_count = count_ambiguous_bases(sequence)

                    if (seq_minlen <= seqlen <= seq_maxlen) and (ambiguous_count <= maxambigs):

                        # Write the sequence with the previous header, including its length
                        output.write(f"{header} len={seqlen}\n")
                        output.write(f"{sequence}\n")
                        # Extract the sequence ID without the version number
                        seqid = header.split(' ')[0][1:]
                        seqids_above_cutoff.append(seqid)
                    else:
                        print(f"validate_sequences:: {header}", file=sys.stderr)
                        print(f"validate_sequences:: sequence length:{seqlen} is outside the acceptable range [{seq_minlen}, {seq_maxlen}] or ambiguous bases:{ambiguous_count} > {maxambigs}", file=sys.stderr)
                # Write the current header (without any sequence yet)
                header = line.strip()  # Remove newline characters
                # Use regex to remove the version number (e.g., ".1", ".2", etc.)
                header = re.sub(r'\.\d+', '', header)
                sequence = ""  # Reset sequence
            else:
                sequence += line.strip()  # Add sequence lines, stripping newlines

        # After the loop, write the last sequence with its header and length
        if header:
            sequence = remove_whitespace_regex(sequence)
            seqlen = len(sequence)
            ambiguous_count = count_ambiguous_bases(sequence)
            if seq_minlen <= seqlen <= seq_maxlen and ambiguous_count < maxambigs:
                # Write the sequence with the previous header, including its length
                output.write(f"{header} len={seqlen}\n")
                output.write(f"{sequence}\n")
                # Extract the sequence ID without the version number
                seqid = header.split(' ')[0][1:]
                seqids_above_cutoff.append(seqid)
            else:
                print(f"validate_sequences:: {header}", file=sys.stderr)
                print(f"validate_sequences:: {seqlen} is outside the acceptable range [{seq_minlen}, {seq_maxlen}] or {ambiguous_count} > {maxambigs}", file=sys.stderr)

    return  seqids_above_cutoff



def main():
    script_name = os.path.basename(__file__)
    # Initialize parser
    parser = argparse.ArgumentParser(
        description=f"Example of {script_name} usage.",
        epilog=(f"Example commands:\npython {script_name}\n\t--fasta seq_nt.fasta\n\t--bvbrc BVBRC_genome.csv\n\t> reformat.fasta 2> log.txt\n")
        ,formatter_class=argparse.RawTextHelpFormatter
    )

    # Adding optional argument
    parser.add_argument("-f", "--fasta", help = "fasta", required=True)
    parser.add_argument("-o", "--output", help = "output file", required=True)
    parser.add_argument("--minlen", type=int, default=10000, help = "min sequence length", required=False)
    parser.add_argument("--maxlen", type=int, default=20000, help = "max sequence length", required=False)

    parser.add_argument("--maxambigs", type=int, default=0, help = "max number of ambigs bases", required=False)
   
    
    # Read arguments from command line
    args = parser.parse_args()
    
    ids_above_cutoff = validate_sequences(args.fasta, args.output, args.minlen, args.maxlen, args.maxambigs)
   

if __name__ == "__main__":
    main()
