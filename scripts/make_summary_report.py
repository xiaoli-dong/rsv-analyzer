#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import re
import os

def process_subgroup(consensus_path, depth_path, nextclade_path, samplesheet_path, subgroup_label):
    """
    Processes a subgroup (A or B). Returns a clean DataFrame where each row is a unique contig.
    """
    # 1. Load Samplesheet & Create Clean Lookup
    if not os.path.exists(samplesheet_path):
        print(f"--- Subgroup {subgroup_label}: Samplesheet not found. Skipping.")
        return pd.DataFrame()

    samples = pd.read_csv(samplesheet_path).rename(columns={"sample": "sample_id"})
    lookup = samples[["sample_id"]].drop_duplicates()
    lookup["subgroup"] = subgroup_label

    # 2. Process Consensus
    cons_df = pd.DataFrame()
    if os.path.exists(consensus_path) and os.path.getsize(consensus_path) > 0:
        try:
            cons = pd.read_csv(consensus_path, sep="\t")
            if not cons.empty:
                cons_df = cons.rename(columns={
                    "#id": "contigid",
                    "coverage": "consensus_coverage/10x",
                    "completeness": "consensus_completeness/10x"
                })[["contigid", "consensus_coverage/10x", "consensus_completeness/10x"]].drop_duplicates()
        except Exception as e:
            print(f"⚠ Warning: Failed to read consensus file {consensus_path}: {e}")

    # 3. Process Depth
    depth_df = pd.DataFrame()
    if os.path.exists(depth_path) and os.path.getsize(depth_path) > 0:
        try:
            depth = pd.read_csv(depth_path, sep="\t")
            if not depth.empty:
                depth["contigid"] = depth["sample"].astype(str) + "|ref|" + depth["chrom"].astype(str)
                depth_df = depth[["contigid", "mapped_reads", "mean_depth"]].drop_duplicates()
        except Exception as e:
            print(f"⚠ Warning: Failed to read depth file {depth_path}: {e}")

    # 4. Process Nextclade
    nextclade_df = pd.DataFrame()
    if os.path.exists(nextclade_path) and os.path.getsize(nextclade_path) > 0:
        try:
            nc = pd.read_csv(nextclade_path, sep="\t")
            # Accept multiple possible column names
            seq_col_candidates = ["seqName", "sequenceName", "SeqName", "name"]
            found_col = next((c for c in seq_col_candidates if c in nc.columns), None)
            if found_col:
                nextclade_df = nc.rename(columns={found_col: "contigid"})
                if "clade" not in nextclade_df.columns:
                    nextclade_df["clade"] = pd.NA
                nextclade_df = nextclade_df[["contigid", "clade"]].drop_duplicates()
            else:
                print(f"⚠ Warning: No recognized sequence column in {nextclade_path}. Skipping Nextclade merge for {subgroup_label}.")
        except Exception as e:
            print(f"⚠ Warning: Failed to read Nextclade file {nextclade_path}: {e}")

    # 5. Merge Analysis Results (Contig Level)
    merged_analysis = pd.DataFrame()
    if not cons_df.empty or not depth_df.empty or not nextclade_df.empty:
        merged_analysis = cons_df if not cons_df.empty else pd.DataFrame(columns=["contigid"])
        if not depth_df.empty:
            merged_analysis = pd.merge(merged_analysis, depth_df, on="contigid", how="outer") if not merged_analysis.empty else depth_df
        if not nextclade_df.empty:
            merged_analysis = pd.merge(merged_analysis, nextclade_df, on="contigid", how="outer") if not merged_analysis.empty else nextclade_df

        # Extract sample_id from contigid to join with subgroup lookup
        if "contigid" in merged_analysis.columns:
            merged_analysis["sample_id"] = merged_analysis["contigid"].str.split("|").str[0]

        # Join with subgroup label
        final_subgroup_df = merged_analysis.merge(lookup, on="sample_id", how="inner") if "sample_id" in merged_analysis.columns else pd.DataFrame()
        return final_subgroup_df

    return pd.DataFrame()

def extract_trailing_number(sample_id):
    if pd.isna(sample_id): return 999999
    match = re.search(r'(\d+)$', str(sample_id))
    return int(match.group(1)) if match else 999999

def main():
    parser = argparse.ArgumentParser(description="Consolidate RSV A/B Analysis with QC Stats.")
    parser.add_argument("--qc", required=True, help="reads.qc_report.csv")
    parser.add_argument("--consensusA", required=True)
    parser.add_argument("--depthA", required=True)
    parser.add_argument("--nextcladeA", required=True)
    parser.add_argument("--samplesheetA", required=True)
    parser.add_argument("--consensusB", required=True)
    parser.add_argument("--depthB", required=True)
    parser.add_argument("--nextcladeB", required=True)
    parser.add_argument("--samplesheetB", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # --- 1. Load QC Stats
    try:
        qc = pd.read_csv(args.qc).rename(columns={"sample": "sample_id"})
        qc = qc[["sample_id", "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs"]].drop_duplicates(subset=["sample_id"])
    except Exception as e:
        print(f"⚠ Error reading QC file {args.qc}: {e}")
        qc = pd.DataFrame(columns=["sample_id", "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs"])

    # --- 2. Process A and B Subgroups
    df_a = process_subgroup(args.consensusA, args.depthA, args.nextcladeA, args.samplesheetA, "A")
    df_b = process_subgroup(args.consensusB, args.depthB, args.nextcladeB, args.samplesheetB, "B")

    # --- 3. Combine Downstream Analysis
    downstream = pd.concat([df_a, df_b], ignore_index=True) if not df_a.empty or not df_b.empty else pd.DataFrame()

    # --- 4. Final Merge: QC Stats + Analysis
    master = downstream.merge(qc, on="sample_id", how="outer") if not downstream.empty else qc.copy()

    # --- 5. Formatting & Sorting
    master["subgroup"] = master.get("subgroup", pd.Series(["No Hit"]*len(master))).fillna("No Hit")
    master["_num"] = master["sample_id"].apply(extract_trailing_number)
    master = master.sort_values(by=["_num", "subgroup"], na_position="last").drop(columns=["_num"])

    # Define final column order
    column_order = [
        "sample_id", "contigid", "subgroup", "clade",
        "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs",
        "consensus_coverage/10x", "consensus_completeness/10x",
        "mapped_reads", "mean_depth"
    ]
    existing_cols = [c for c in column_order if c in master.columns]
    master = master[existing_cols]

    # Split sample_id for run-level reporting
    if "sample_id" in master.columns:
        sample_split = master["sample_id"].str.split("-", n=1, expand=True)
        master.insert(0, "runid", sample_split[0])
        master.insert(1, "sample", sample_split[1] if sample_split.shape[1] > 1 else "")
        master = master.drop(columns=["sample_id"])

    # --- 6. Save
    master.to_csv(args.out, sep="\t", index=False, na_rep="NaN")
    print(f"✔ Final master report written to: {args.out}")
    print(f"Total Rows: {len(master)} (includes A/B hits and QC-only samples)")

if __name__ == "__main__":
    main()
