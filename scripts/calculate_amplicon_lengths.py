import csv
import re
import sys
from collections import defaultdict

def parse_amplicon_id(primer_name):
    match = re.search(r'_(\d+)_', primer_name)
    return int(match.group(1)) if match else None

def load_primers(bed_file):
    left_primers = {}
    right_primers = {}

    with open(bed_file, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 6:
                continue
            chrom, start, end, name, _, strand = row
            start, end = int(start), int(end)
            amp_id = parse_amplicon_id(name)

            if "_LEFT" in name:
                left_primers[amp_id] = start
            elif "_RIGHT" in name:
                right_primers[amp_id] = end

    return left_primers, right_primers

def compute_amplicon_sizes(left_primers, right_primers):
    sizes = {}
    for amp_id in left_primers:
        if amp_id in right_primers:
            size = right_primers[amp_id] - left_primers[amp_id] + 1
            sizes[amp_id] = size
    return sizes

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python amplicon_size_calc.py <primers.bed>")
        sys.exit(1)

    bed_file = sys.argv[1]
    left, right = load_primers(bed_file)
    sizes = compute_amplicon_sizes(left, right)

    print("AmpliconID\tSize(bp)")
    for amp_id in sorted(sizes):
        print(f"{amp_id}\t\t{sizes[amp_id]}")
