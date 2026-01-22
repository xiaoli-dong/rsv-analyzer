#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import re
import os

def process_subgroup(consensus_path, depth_path, nextclade_path, samplesheet_path, subgroup_label):
    """
    Process one RSV subgroup and return a DataFrame of downstream results only:
    - Each sample may have multiple references.
    - QC-only info is not included here (merged later).
    """
    # --- Sample sheet → subgroup
    if not os.path.exists(samplesheet_path):
        print(f"Warning: Samplesheet {samplesheet_path} not found. Returning empty DataFrame.")
        return pd.DataFrame()

    samples = (
        pd.read_csv(samplesheet_path)
        .rename(columns={"sample": "sample_id"})
        .drop_duplicates(subset=["sample_id"])
    )
    if samples.empty:
        return pd.DataFrame()

    samples["subgroup"] = subgroup_label
    lookup = samples[["sample_id", "subgroup"]]

    # --- Consensus
    cons_df = pd.DataFrame()
    if os.path.exists(consensus_path) and os.path.getsize(consensus_path) > 0:
        cons = pd.read_csv(consensus_path, sep="\t")
        if not cons.empty:
            cons_df = (
                cons.rename(columns={"#id": "raw_id",
                                     "coverage": "consensus_coverage/10x",
                                     "completeness": "consensus_completeness/10x"})
                .assign(
                    sample_id=lambda df: df["raw_id"].str.split("|").str[0],
                    ref=lambda df: df["raw_id"].str.split("|").str[1].str.replace("^ref_", "", regex=True)
                )
            )[["sample_id", "ref", "consensus_coverage/10x", "consensus_completeness/10x"]]

    # --- Depth
    depth_df = pd.DataFrame()
    if os.path.exists(depth_path) and os.path.getsize(depth_path) > 0:
        depth = pd.read_csv(depth_path, sep="\t")
        if not depth.empty:
            depth_df = depth.rename(columns={"sample": "sample_id", "chrom": "ref"})[["sample_id", "ref", "mapped_reads", "mean_depth"]]

    # --- Nextclade
    nextclade_df = pd.DataFrame()
    if os.path.exists(nextclade_path) and os.path.getsize(nextclade_path) > 0:
        nc = pd.read_csv(nextclade_path, sep="\t")
        if not nc.empty:
            nextclade_df = (
                nc.rename(columns={"seqName": "raw_id"})
                .assign(
                    sample_id=lambda df: df["raw_id"].str.split().str[0].str.split("-").str[0] + "-" +
                                           df["raw_id"].str.split().str[0].str.split("-").str[-1],
                    ref=lambda df: df["raw_id"].str.split().str[1]
                )[["sample_id", "ref", "clade"]]
            )

    # --- Merge downstream tables on sample_id + ref
    merged = cons_df.merge(depth_df, on=["sample_id", "ref"], how="outer") \
                    .merge(nextclade_df, on=["sample_id", "ref"], how="outer")

    # Attach subgroup from sample sheet
    merged = merged.merge(lookup, on="sample_id", how="left")

    return merged

def extract_trailing_number(sample_id):
    match = re.search(r'(\d+)$', sample_id)
    return int(match.group(1)) if match else 999999

def main():
    parser = argparse.ArgumentParser(description="Merge RSV-A and RSV-B outputs with QC info added later.")
    parser.add_argument("--qc", required=True, help="reads.qc_report.csv (all samples)")
    parser.add_argument("--consensusA", required=True)
    parser.add_argument("--depthA", required=True)
    parser.add_argument("--nextcladeA", required=True)
    parser.add_argument("--consensusB", required=True)
    parser.add_argument("--depthB", required=True)
    parser.add_argument("--nextcladeB", required=True)
    parser.add_argument("--samplesheetA", required=True)
    parser.add_argument("--samplesheetB", required=True)
    parser.add_argument("--out", required=True, help="Output master TSV")
    args = parser.parse_args()

    # --- Load QC table (all samples)
    qc = pd.read_csv(args.qc).rename(columns={"sample": "sample_id"}).drop_duplicates(subset=["sample_id"])
    qc_cols = ["sample_id", "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs"]
    qc = qc[qc_cols]

    # --- Process downstream analysis for A and B
    print("Processing RSV-A...")
    masterA = process_subgroup(args.consensusA, args.depthA, args.nextcladeA, args.samplesheetA, "A")
    print(f"  RSV-A downstream rows: {len(masterA)}")

    print("Processing RSV-B...")
    masterB = process_subgroup(args.consensusB, args.depthB, args.nextcladeB, args.samplesheetB, "B")
    print(f"  RSV-B downstream rows: {len(masterB)}")

    # --- Concatenate A + B downstream results
    dfs = []
    if not masterA.empty:
        dfs.append(masterA)
    if not masterB.empty:
        dfs.append(masterB)
    downstream = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    # --- Merge QC info after concatenation
    master = qc.merge(downstream, on="sample_id", how="left")

    # --- Ensure subgroup column exists (QC-only samples will have NaN)
    #master["subgroup"] = master["subgroup"].fillna("NA")

    # --- Sort by sample_id numeric, then subgroup
    master["_num"] = master["sample_id"].apply(extract_trailing_number)
    master = master.sort_values(by=["_num", "subgroup"], na_position="last").drop(columns=["_num"])

    # --- Reorder columns
    column_order = [
        "sample_id", "subgroup", "clade", "ref",
        "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs",
        "consensus_coverage/10x", "consensus_completeness/10x",
        "mapped_reads", "mean_depth"
    ]
    master = master[column_order]

    # --- Split sample_id into runid and sample
    sample_split = master["sample_id"].str.split("-", n=1, expand=True)
    master.insert(0, "runid", sample_split[0])
    master.insert(1, "sample", sample_split[1])
    master = master.drop(columns=["sample_id"])

    # --- Write output
    master.to_csv(args.out, sep="\t", index=False, na_rep="NaN")
    print(f"✔ Master TSV written: {args.out}")
    print(f"Total rows: {len(master)}")

if __name__ == "__main__":
    main()
