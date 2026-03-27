import pandas as pd
import sys
import re

def natural_sort_key(s):
    """Sorts string-number combinations (S2 before S29) correctly."""
    if pd.isna(s): return s
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]

def pivot_coverage(input_file):
    try:
        # 1. Read input - Assuming header exists: chrom, start, end, region, coverage, sample
        df = pd.read_csv(input_file, sep='\t')

        # 2. Rename for clarity (Adjust these based on your exact TSV header)
        # mapping: region -> Primer, coverage -> Depth, sample -> SampleID
        df = df.rename(columns={
            'region': 'Primer',
            'coverage': 'Depth',
            'sample': 'SampleID'
        })

        # 3. Drop footer/metadata rows (removes that trailing "sample,coverage,,,,")
        df = df[~df['SampleID'].str.contains('sample|coverage|total', case=False, na=False)]

        # 4. Pivot
        pivot_df = df.pivot_table(index='SampleID',
                                  columns='Primer',
                                  values='Depth',
                                  aggfunc='first')

        # 5. Remove the 'Primer' label from the top-left corner
        pivot_df.columns.name = None
        pivot_df.index.name = None

        # 6. Natural Sort Rows (SampleID) and Columns (Primers)
        sorted_indices = sorted(pivot_df.index, key=natural_sort_key)
        pivot_df = pivot_df.reindex(sorted_indices)

        sorted_primers = sorted(pivot_df.columns, key=natural_sort_key)
        pivot_df = pivot_df[sorted_primers]

        # 7. Move SampleID from the index to a proper first column
        pivot_df.reset_index(inplace=True)
        pivot_df.rename(columns={'index': 'SampleID'}, inplace=True)

        # 8. Output to CSV (Comma Delimited)
        # index=False ensures we don't get the extra 0, 1, 2... row numbers
        print(pivot_df.to_csv(index=False, sep=','))

    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        pivot_coverage(sys.argv[1])
    else:
        print("Usage: python script.py your_input_file.tsv")
