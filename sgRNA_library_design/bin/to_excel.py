#!/usr/bin/env python3
import os
import glob
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='This script produce a lolipop plot (with domain/feature annotations) base on a MaGeCK, a Variant Effect Prediction file (VEP), and bed style file')
parser.add_argument('-c', '--csv', metavar='FILE', dest='csv', required=True, type=str, help='')
parser.add_argument('-v', '--vep', metavar='FILE', dest='vep_files', required=True, type=str, nargs='+', help='tab delimited variant effect prediction file (VEP) with id corresponding to input')
parser.add_argument('-o', '--out', metavar='FILE', dest='output', required=True, type=str, help='')


def save_to_excel(main_df, other_dfs, output_file):
    # Initialize a writer for the Excel file
    with pd.ExcelWriter(output_file) as writer:
        # Write the main DataFrame to the "General" sheet
        main_df.to_excel(writer, sheet_name='Library', engine='openpyxl' ,index=False)
        
        # Loop over the list of other DataFrames
        for df in other_dfs:
            # Extract the unique value in the 'editor' column
            editor_name = df['editor'].iloc[0]
            
            # Write the DataFrame to a sheet named "editor - X"
            sheet_name = f"editor - {editor_name}"
            df.to_excel(writer, sheet_name=sheet_name, engine='openpyxl', index=False)


if __name__ == '__main__':
    args = parser.parse_args()
    csv=pd.read_csv(args.csv)
    veps=[pd.read_csv(i, sep='\t') for i in args.vep_files]
    save_to_excel(csv, veps, args.output)
