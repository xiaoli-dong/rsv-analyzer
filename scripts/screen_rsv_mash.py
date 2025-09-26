#!/usr/bin/env python3
import os
import subprocess
import csv
import argparse
import re
from collections import defaultdict

def find_sample_pairs(fastq_dir):
    files = [f for f in os.listdir(fastq_dir) if f.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq'))]
    files.sort()
    sample_dict = defaultdict(dict)

    for f in files:
        path = os.path.join(fastq_dir, f)

       
        m = re.match(r'^(S\d+)_', f)
        if m:
            sample_key = m.group(1)
        else:
           
            #sample_key = os.path.splitext(f)[0]
            sample_key = f.split(".")[0]

       
        if re.search(r'(clean_1|_1\.)', f):
            read = 'R1'
        elif re.search(r'(clean_2|_2\.)', f):
            read = 'R2'
        else:
            read = 'R1'  

        sample_dict[sample_key][read] = path

    sample_pairs = []
    for sample, reads in sample_dict.items():
        r1 = reads.get('R1', '')
        r2 = reads.get('R2', '')
        sample_pairs.append((sample, [r1, r2]))

    sample_pairs.sort(key=lambda x: int(re.match(r'S(\d+)', x[0]).group(1)) if re.match(r'S(\d+)', x[0]) else x[0])
    return sample_pairs

def run_mash_screen(mash_db, sample_name, fastq_files, outdir):
    cmd = ["mash", "screen", "-i", "0.95", "-w", "-v", "0.1", mash_db] + [f for f in fastq_files if f]
    print(f"[LOG] Running: {' '.join(cmd)}")

    try:
        print(cmd)
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        mash_outdir = os.path.join(outdir, "mash_screen")
        os.makedirs(mash_outdir, exist_ok=True)
        output_file = os.path.join(mash_outdir, f"{sample_name}.mash.tsv")
        with open(output_file, "w") as f:
            f.write(result.stdout)

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
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "fastq_1", "fastq_2"])
        for sample, paths in sample_list:
            writer.writerow([sample, paths[0], paths[1] if len(paths) > 1 else ""])

def main():
    parser = argparse.ArgumentParser(description="Mash screen RSV sequences and split by subtype.")
    parser.add_argument("-i", "--input", required=True, help="Input directory with FASTQ files")
    parser.add_argument("-d", "--db", required=True, help="Mash sketch DB (e.g. rsv.msh)")
    parser.add_argument("-o", "--outdir", default="rsv_screen_output", help="Output directory for mash screen output")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    rsv_a_samples = []
    rsv_b_samples = []

    sample_pairs = find_sample_pairs(args.input)
    print("[LOG] Found sample pairs:")
    for s, reads in sample_pairs:
        print(f"  {s}: R1={reads[0]}  R2={reads[1]}")

    for sample, fastq_paths in sample_pairs:
        print(f"[✓] Processing: {sample}")
        matches_A, matches_B = run_mash_screen(args.db, sample, fastq_paths, args.outdir)

        if matches_A:
            rsv_a_samples.append((sample, fastq_paths))
        if matches_B:
            rsv_b_samples.append((sample, fastq_paths))

    write_csv(os.path.join(args.outdir, "samplesheet_rsvA.csv"), rsv_a_samples)
    write_csv(os.path.join(args.outdir, "samplesheet_rsvB.csv"), rsv_b_samples)

    print("[✓] Mash screen complete.")
    print(f"[✓] rsvA samplesheet written to: {os.path.join(args.outdir, 'samplesheet_rsvA.csv')}")
    print(f"[✓] rsvB samplesheet written to: {os.path.join(args.outdir, 'samplesheet_rsvB.csv')}")
    print("[✓] Mash results saved per sample as .tsv files into mash_results directory")

if __name__ == "__main__":
    main()
