#!/usr/bin/env python3
import argparse
from collections import defaultdict
import re

def readFasta(fasta):
    id2seq = defaultdict(str)

    with open(fasta, 'r') as in_f:
        seqid = ""
        for line in in_f:
            if line.startswith(">"):
                seqid = line[1:].strip().split()[0]   # keep only first token
                id2seq[seqid] = ""
            else:
                id2seq[seqid] += re.sub(r'\s+', '', line)

    return id2seq


# remove contigs with too many ambiguous bases
def filterSeq(id2seq, maxambigs):
    allowed_bases = set("ATGC")

    filtered = {}
    for seqid, seq in id2seq.items():
        seq = seq.upper()
        seqlen = len(seq)

        if seqlen == 0:
            continue

        valid_count = sum(1 for b in seq if b in allowed_bases)
        ambig_ratio = (seqlen - valid_count) / seqlen

        if ambig_ratio <= maxambigs:
            filtered[seqid] = seq

    return filtered


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by ambiguous base ratio (keep original IDs)."
    )
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA")
    parser.add_argument("--maxambigs", type=float, default=0.25,
                        help="Max ratio of ambiguous bases (default: 0.25)")

    args = parser.parse_args()

    id2seq = readFasta(args.fasta)
    id2seq = filterSeq(id2seq, args.maxambigs)

    # output with original IDs
    for seqid in sorted(id2seq):
        print(f">{seqid}")
        print(id2seq[seqid])


if __name__ == "__main__":
    main()
