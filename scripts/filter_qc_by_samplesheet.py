#!/usr/bin/env python3
import argparse
import pandas as pd

def main():

    parser = argparse.ArgumentParser(
        description="Filter reads.qc_report.csv using sample IDs from samplesheet CSVs."
    )

    parser.add_argument(
        "--samplesheet", "-s",
        required=True,
        nargs="+",
        help="One or more samplesheet CSV files (must contain a 'sample' column)."
    )

    parser.add_argument(
        "--qc",
        required=True,
        help="QC report CSV file (must contain a 'sample' column)."
    )

    parser.add_argument(
        "--out", "-o",
        required=True,
        help="Output filtered CSV file."
    )

    args = parser.parse_args()

    # -------------------------
    # 1. Read sample IDs
    # -------------------------
    sample_ids = set()

    for sheet in args.samplesheet:
        df = pd.read_csv(sheet)
        if "sample" not in df.columns:
            raise ValueError(f"ERROR: File {sheet} does not contain a 'sample' column.")
        sample_ids.update(df["sample"].astype(str))

    print(f"Loaded {len(sample_ids)} unique sample IDs from samplesheets.")

    # -------------------------
    # 2. Read QC report
    # -------------------------
    qc = pd.read_csv(args.qc)

    if "sample" not in qc.columns:
        raise ValueError(f"ERROR: QC file {args.qc} does not contain a 'sample' column.")

    print(f"Loaded {qc.shape[0]} rows from QC report.")

    # -------------------------
    # 3. Filter
    # -------------------------
    qc_filtered = qc[qc["sample"].astype(str).isin(sample_ids)]

    print(f"Filtered down to {qc_filtered.shape[0]} matching rows.")

    # -------------------------
    # 4. Write output
    # -------------------------
    qc_filtered.to_csv(args.out, index=False)
    print(f"Filtered QC report written to: {args.out}")


if __name__ == "__main__":
    main()
