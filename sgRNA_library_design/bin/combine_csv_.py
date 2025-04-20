#!/usr/bin/env python3
import pandas as pd
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Combine guide CSVs and generate ID mappings.")
    parser.add_argument('--csv', nargs='+', required=True, metavar='CSV_FILE',
                        help='List of input CSV files (tab-separated)')
    parser.add_argument('--output', required=True, metavar='BASE_NAME',
                        help='Base name for output files (without extension)')
    args = parser.parse_args()

    # Combine CSV files
    csv_frames = [pd.read_csv(f) for f in args.csv]
    csv_combined = pd.concat(csv_frames,ignore_index=True).drop_duplicates()

    # Sort and generate ID
    csv_combined = csv_combined.sort_values(['ID'])
    csv_combined['Protein']=[i.split('_')[0] for i in csv_combined['ID']]
    num = csv_combined.groupby(['Protein']).cumcount() + 1
    csv_combined['ID'] = csv_combined['Protein'] + '_' + num.astype(str)
    csv_combined = csv_combined[['ID'] + [col for col in csv_combined.columns if col != 'ID']]
    csv_combined=csv_combined.drop('spacer', axis=1)
    csv_combined=csv_combined.drop('Protospacer', axis=1)

    # Output combined CSV
    csv_combined.to_csv(f"{args.output}.csv", index=False)

    # Output ID dictionary
    dict_df = csv_combined[['ID', 'Protein', 'start', 'strand']]
    dict_df.to_csv(f"{args.output}_id_dic.tsv", sep='\t', index=False)

if __name__ == '__main__':
    main()

