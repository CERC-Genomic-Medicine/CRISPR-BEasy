#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Combine CSVs of failed guides and report their count.")
    parser.add_argument('--csv', nargs='+', required=True, metavar='CSV_FILE',
                        help='List of input CSV files (tab-separated)')
    parser.add_argument('--output', required=True, metavar='OUTPUT_PREFIX',
                        help='Output prefix for the combined CSV file')
    args = parser.parse_args()

    # Read and combine CSV files
    csv_frames = [pd.read_csv(f, sep='\t') for f in args.csv]
    csv_combined = pd.concat(csv_frames,ignore_index=True).drop_duplicates()

    # Write combined output
    output_file = f"{args.output}"
    csv_combined.to_csv(output_file, index=False)

    # Print message
    print(f"::message:: {len(csv_combined)} sgRNA were removed by the CFD threshold specified")

if __name__ == '__main__':
    main()

