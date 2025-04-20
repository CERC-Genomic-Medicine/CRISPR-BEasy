#!/usr/bin/env python3
import os
import glob
import pandas as pd
import argparse
import xlsxwriter

parser = argparse.ArgumentParser(description='This script produce an excel file based on the general annotation of library and the VEP annotation')
parser.add_argument('-c', '--csv', metavar='FILE', dest='csv', required=True, type=str, help='')
parser.add_argument('-v', '--vep', metavar='FILE', dest='vep_files', required=True, type=str, nargs='+', help='tab delimited variant effect prediction file (VEP) with id corresponding to input')
parser.add_argument('-o', '--out', metavar='FILE', dest='output', required=True, type=str, help='')


def is_excel_compatible(df):
    max_rows = 1048576
    max_cols = 16384
    return df.shape[0] <= max_rows and df.shape[1] <= max_cols

def save_to_excel(main_df, other_dfs, output_file):
    # Initialize a writer for the Excel file
    with pd.ExcelWriter(output_file) as writer:
        # Write the main DataFrame to the "General" sheet
        if is_excel_compatible(main_df):
                main_df.to_excel(writer, sheet_name='Library', engine='xlsxwriter' ,index=False)
        else :
                print('::message:: Library is too large for excel look for .csv file in the Auxiliary files')
        
        # Loop over the list of other DataFrames
        for df in other_dfs:
            # Extract the unique value in the 'editor' column
            editor_name = df['editor'].iloc[0]
            
            # Write the DataFrame to a sheet named "editor - X"
            sheet_name = f"editor - {editor_name}"
            if not is_excel_compatible(df):
                print(f'::message:: Annotation {editor_name} of this Library is too large for excel look for vep.tsv file in the Auxiliary files')
            else :
                df.to_excel(writer, sheet_name=sheet_name, engine='xlsxwriter', index=False)


if __name__ == '__main__':
    args = parser.parse_args()
    csv=pd.read_csv(args.csv, low_memory=False)
    veps=[pd.read_csv(i, sep='\t',low_memory=False) for i in args.vep_files]
    save_to_excel(csv, veps, args.output)
