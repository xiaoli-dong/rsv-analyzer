#!/usr/bin/env python3
import os
import subprocess
import csv
import argparse
import re
from collections import defaultdict

# ---------------------------
# Helper Functions
# ---------------------------

def find_sample_pairs(fastq_dir):
    """
    Scan the FASTQ directory and identify sample pairs.
    Supports files ending with .fastq/.fq and .gz
    """
    files = [f for f in os.listdir(fastq_dir) if f.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq'))]
    files.sort()
    sample_dict = defaultdict(dict)

    for f in files:
        path = os.path.join(fastq_dir, f)

        # Extract sample key
        m = re.match(r'^(.*S\d+)_', f)
        if m:
            sample_key = m.group(1)
        else:
            sample_key = f.split(".")[0]

        # Determine read direction
        if re.search(r'(clean_1|_1\.)', f):
            read = 'R1'
        elif re.search(r'(clean_2|_2\.)', f):
            read = 'R2'
        else:
            read = 'R1'

        sample_dict[sample_key][read] = path

    # Build list of sample tuples: (sample_name, [R1, R2])
    sample_pairs = []
    for sample, reads in sample_dict.items():
        r1 = reads.get('R1', '')
        r2 = reads.get('R2', '')
        sample_pairs.append((sample, [r1, r2]))

    # Sort by numeric part of sample name if exists
    def sort_key(x):
        m = re.match(r'S(\d+)', x[0])
        return int(m.group(1)) if m else x[0]
    sample_pairs.sort(key=sort_key)

    return sample_pairs


def run_mash_screen(mash_db, sample_name, fastq_files, outdir):
    """
    Run mash screen for a given sample.
    Returns sets of detected RSV_A and RSV_B.
    Python 3.6 compatible.
    """
    cmd = ["mash", "screen", "-i", "0.95", "-w", "-v", "0.1", mash_db] + [f for f in fastq_files if f]
    print(f"[LOG] Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,  # text=True equivalent
            check=True
        )

        # Save mash output per sample
        mash_outdir = os.path.join(outdir, "mash_screen")
        os.makedirs(mash_outdir, exist_ok=True)
        output_file = os.path.join(mash_outdir, f"{sample_name}.mash.tsv")
        with open(output_file, "w") as f:
            f.write(result.stdout)

        # Parse for RSV_A and RSV_B matches
        matches_A = set()
        matches_B = set()
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            ref_id = parts[4]
            if "|RSV_A" in ref_id:
                matches_A.add("RSV_A")
            if "|RSV_B" in ref_id:
                matches_B.add("RSV_B")

        return matches_A, matches_B

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] mash screen failed for {sample_name}: {e.stderr}")
        return set(), set()


def write_csv(output_file, sample_list):
    """
    Write sample list to CSV: sample, fastq_1, fastq_2
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "fastq_1", "fastq_2"])
        for sample, paths in sample_list:
            r1 = paths[0] if len(paths) > 0 else ""
            r2 = paths[1] if len(paths) > 1 else ""
            writer.writerow([sample, r1, r2])


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(description="Mash screen RSV sequences and split by subtype.")
    parser.add_argument("-i", "--input", required=True, help="Input directory with FASTQ files")
    parser.add_argument("-d", "--db", required=True, help="Mash sketch DB (e.g. rsv.msh)")
    parser.add_argument("-o", "--outdir", default="rsv_screen_output", help="Output directory for mash screen output")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    rsv_a_samples = []
    rsv_b_samples = []

    # Find all sample pairs
    sample_pairs = find_sample_pairs(args.input)
    print("[LOG] Found sample pairs:")
    for s, reads in sample_pairs:
        print(f"  {s}: R1={reads[0]}  R2={reads[1]}")

    # Run mash screen for each sample
    for sample, fastq_paths in sample_pairs:
        print(f"[✓] Processing: {sample}")
        matches_A, matches_B = run_mash_screen(args.db, sample, fastq_paths, args.outdir)

        if matches_A:
            rsv_a_samples.append((sample, fastq_paths))
        if matches_B:
            rsv_b_samples.append((sample, fastq_paths))

    # Write output CSVs
    write_csv(os.path.join(args.outdir, "samplesheet_rsvA.csv"), rsv_a_samples)
    write_csv(os.path.join(args.outdir, "samplesheet_rsvB.csv"), rsv_b_samples)

    print("[✓] Mash screen complete.")
    print(f"[✓] rsvA samplesheet: {os.path.join(args.outdir, 'samplesheet_rsvA.csv')}")
    print(f"[✓] rsvB samplesheet: {os.path.join(args.outdir, 'samplesheet_rsvB.csv')}")
    print(f"[✓] Mash results saved per sample in: {os.path.join(args.outdir, 'mash_screen')}")

if __name__ == "__main__":
    main()
