#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import re

def process_subgroup(
    qc_df,
    consensus_path,
    depth_path,
    nextclade_path,
    samplesheet_path,
    subgroup_label,
):
    """
    Process one RSV subgroup (A or B) and return a merged DataFrame
    with one row per (sample_id, subgroup).
    """

    # --- Sample sheet → subgroup lookup
    samples = (
        pd.read_csv(samplesheet_path)
        .rename(columns={"sample": "sample_id"})
        .drop_duplicates(subset=["sample_id"])
    )
    samples["subgroup"] = subgroup_label
    lookup = samples[["sample_id", "subgroup"]]

    # --- QC anchor for this subgroup
    qc_sub = qc_df.merge(lookup, on="sample_id", how="inner")

    # --- Consensus stats
    cons = (
        pd.read_csv(consensus_path, sep="\t")
        .rename(
            columns={
                "#id": "sample_id",
                "coverage": "consensus_coverage/10x",
                "completeness": "consensus_completeness/10x",
            }
        )
        .drop_duplicates(subset=["sample_id"])
    )[
        ["sample_id", "consensus_coverage/10x", "consensus_completeness/10x"]
    ]

    # --- Depth summary
    depth = (
        pd.read_csv(depth_path, sep="\t")
        .rename(columns={"sample": "sample_id"})
        .drop_duplicates(subset=["sample_id"])
    )[
        ["sample_id", "mapped_reads", "mean_depth"]
    ]

    # --- Nextclade
    nextclade = (
        pd.read_csv(nextclade_path, sep="\t")
        .rename(columns={"seqName": "sample_id"})
        .drop_duplicates(subset=["sample_id"])
    )
    nextclade["sample_id"] = nextclade["sample_id"].str.split().str[0]
    nextclade = nextclade[["sample_id", "clade"]]

    # --- Merge all subgroup-level tables
    merged = (
        qc_sub
        .merge(cons, on="sample_id", how="left")
        .merge(depth, on="sample_id", how="left")
        .merge(nextclade, on="sample_id", how="left")
    )

    # --- Keep only one row per sample_id within this subgroup
    merged = merged.drop_duplicates(subset=["sample_id"])

    return merged

def extract_trailing_number(sample_id):
    """
    Extract the trailing number at the end of sample_id.
    Returns integer, or a large number if not found (to sort at the end).
    Examples:
      - "250115_S_I_012-S1" -> 1
      - "250115_S_I_012-barcode01" -> 1
      - "sampleX" -> 999999
    """
    match = re.search(r'(\d+)$', sample_id)
    if match:
        return int(match.group(1))
    else:
        return 999999

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge RSV-A and RSV-B outputs into one master TSV "
            "(coinfection-safe, single run)"
        )
    )

    parser.add_argument("--qc", required=True,
                        help="reads.qc_report.csv (all samples)")

    parser.add_argument("--consensusA", required=True)
    parser.add_argument("--depthA", required=True)
    parser.add_argument("--nextcladeA", required=True)

    parser.add_argument("--consensusB", required=True)
    parser.add_argument("--depthB", required=True)
    parser.add_argument("--nextcladeB", required=True)

    parser.add_argument("--samplesheetA", required=True)
    parser.add_argument("--samplesheetB", required=True)

    parser.add_argument("--out", required=True,
                        help="Output master TSV")

    args = parser.parse_args()

    # --- QC anchor table (all samples)
    qc = (
        pd.read_csv(args.qc)
        .rename(columns={"sample": "sample_id"})
        .drop_duplicates(subset=["sample_id"])
    )[
        ["sample_id", "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs"]
    ]

    # --- Process RSV-A and RSV-B independently
    masterA = process_subgroup(
        qc,
        args.consensusA,
        args.depthA,
        args.nextcladeA,
        args.samplesheetA,
        "A",
    )

    masterB = process_subgroup(
        qc,
        args.consensusB,
        args.depthB,
        args.nextcladeB,
        args.samplesheetB,
        "B",
    )

    # --- Combine subgroups
    master = pd.concat([masterA, masterB], ignore_index=True)

    if master.empty:
        sys.exit("ERROR: No samples found after merging")

    # --- Numeric sort by trailing number in sample_id
    master["_num"] = master["sample_id"].apply(extract_trailing_number)
    master["subgroup"] = pd.Categorical(
        master["subgroup"],
        categories=["A", "B"],
        ordered=True,
    )
    master = master.sort_values(by=["_num", "subgroup"])
    master = master.drop(columns=["_num"])

    # --- Reorder columns: sample_id, subgroup, clade, then the rest
    cols = list(master.columns)
    for col in ["subgroup", "clade"]:
        if col in cols:
            cols.remove(col)
    insert_pos = cols.index("sample_id") + 1
    cols.insert(insert_pos, "subgroup")
    cols.insert(insert_pos + 1, "clade")
    master = master[cols]

    # --- Write output
    master.to_csv(args.out, sep="\t", index=False)
    print(f"✔ Master TSV written: {args.out}")

if __name__ == "__main__":
    main()
