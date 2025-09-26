#!/usr/bin/env python3
import os
import sys
import csv
import re
from collections import defaultdict

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <fastq_directory>", file=sys.stderr)
        sys.exit(1)

    fastq_dir = os.path.abspath(sys.argv[1])

    if not os.path.isdir(fastq_dir):
        print(f"Error: {fastq_dir} is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    # Store reads: {sample_id: {'R1': '', 'R2': ''}}
    samples = defaultdict(lambda: {'R1': '', 'R2': ''})

    for fname in sorted(os.listdir(fastq_dir)):
        if not fname.endswith((".fastq", ".fastq.gz")):
            continue
        fpath = os.path.join(fastq_dir, fname)

        # Paired-end: e.g. 10_S10_L001_R1_001.fastq.gz
        m = re.match(r"(\d+)_S(\d+)_.*_(R1|R2)_\d+\.fastq(?:\.gz)?$", fname)
        if m:
            _, sample_num, read_type = m.groups()
            sample_id = f"S{sample_num}"
            samples[sample_id][read_type] = fpath
            continue

        # Single-end: e.g. barcode49.fastq.gz
        m2 = re.match(r"(barcode\d+).*\.fastq(?:\.gz)?$", fname)
        if m2:
            sample_id = m2.group(1)
            samples[sample_id]['R1'] = fpath
            continue

    # Output to stdout
    writer = csv.writer(sys.stdout)
    writer.writerow(["sample", "fastq_1", "fastq_2"])
    for sample_id in sorted(samples.keys(), key=lambda x: (x.startswith("S"), int(re.sub(r'\D', '', x)) if re.search(r'\d+', x) else x)):
        writer.writerow([sample_id, samples[sample_id]['R1'], samples[sample_id]['R2']])

if __name__ == "__main__":
    main()
