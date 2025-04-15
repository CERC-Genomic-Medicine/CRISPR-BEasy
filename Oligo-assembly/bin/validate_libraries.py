#!/usr/bin/env python3

import os
import pandas as pd
import argparse

# Required columns
REQUIRED_COLUMNS = {'ID','Protospacer','PAM','gRNA_seq_POSstrand','Chromosome','POSstart','strand'}


def read_and_validate_file(library, error_list):
    """
    Reads a CSV file into a DataFrame and checks if it has the required columns.
    Appends errors to the error_list if the file is invalid.
    """
    try:
        if os.path.getsize(library[0]) == 0:
            return None
        df = pd.read_excel(library[0], sheet_name="Library")
        if REQUIRED_COLUMNS.issubset(df.columns):
            df['library_oligomer'] = library[1]
            return df
        else:
            error_list.append(f"File '{filepath}' is missing required columns: {REQUIRED_COLUMNS - set(df.columns)}")
            return None
    except Exception as e:
        error_list.append(f"File '{filepath}' of the wrong format see FAQ for details")
        return None


def concatenate_dataframes(dataframes):
    """Concatenates a list of DataFrames into one."""
    if dataframes:
        return pd.concat(dataframes, ignore_index=True)
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no valid DataFrames are found


def check_duplicates(df):
    """Checks for duplicates based on 'ID' and combination of 'chromosome' and 'position'."""
    if df.empty:
        return pd.DataFrame()  # No duplicates if DataFrame is empty

    # Check for duplicated IDs

    # Check for duplicated chromosome-position combinations
    df['chr_pos'] =df['strand'].astype(str) + '_' + df['Chromosome'].astype(str) + '_' + df['POSstart'].astype(str)
    duplicated_chr_pos = df[df.duplicated(subset=['chr_pos'], keep=False)]
    duplicated_id = df[df.duplicated(subset=['ID'], keep=False)]


    # Combine both types of duplicated rows
    duplicated_rows = pd.concat([duplicated_id, duplicated_chr_pos]).drop_duplicates()
    duplicated_rows=duplicated_rows.loc[:,['ID','Chromosome','POSstart']]

    return duplicated_rows


def save_duplicates(duplicated_rows, output_file):
    """Saves duplicated rows to a CSV file."""
    if not duplicated_rows.empty:
        duplicated_rows.to_csv(f'{output_file}', sep='\t', index=False)


def save_errors(error_list, output_file):
    """Saves the list of errors to a file."""
    with open(f'{output_file}', 'w') as f:
        for error in error_list:
            f.write(f"{error}\n")

def validate_and_check_overlap(files):
    all_sheets = []
    editor_sheets = []
    
    # Process each file
    for file in files:
        # Load Excel file and get sheet names
        if os.path.getsize(file) == 0 : 
            continue
        xlsx = pd.ExcelFile(file,engine='openpyxl')
        sheet_names = xlsx.sheet_names
            
        # Skip file if it only has a "library" sheet or is empty
        if not sheet_names or sheet_names == ["library"]:
            continue
            
            # Collect all sheet names for overlap check
        all_sheets.append(set(sheet_names))
            
            # Collect "editor -" sheets
        editor_sheets.append([name for name in sheet_names if name.startswith("editor - ")])
        
    # Check if we have multiple files to compare
    
    # Check for complete overlap in all sheets
    complete_overlap = all(s == all_sheets[0] for s in all_sheets)
    print("Complete overlap in sheet names." if complete_overlap else "No complete overlap in sheet names.")
    
    # Check for overlap in "editor -" sheets specifically, if multiple lists
    editor_overlap = []
    if len(editor_sheets) > 1:
        editor_overlap = [set(ed) for ed in editor_sheets]
        editor_complete_overlap = all(s == editor_overlap[0] for s in editor_overlap)
        if not editor_complete_overlap :
            with open(f'editor_overlap.warning', 'w') as f:
                f.write(f"The editors do not overlap perfectly are you sure.")


def main(target,positive,negative, output_prefix):
    """Main function to orchestrate the reading, validation, concatenation, and duplicate check."""
    error_list = []  # List to store errors
    files = [[target, 'Target_library'],[positive, 'Positive_control'],[negative, 'Negative_control']]
    dataframes = [read_and_validate_file(file, error_list) for file in files if read_and_validate_file(file, error_list) is not None]

    # If there are no errors, proceed with concatenation and duplicate checking
    if not error_list:
        validate_and_check_overlap([target, positive, negative])
        concatenated_df = concatenate_dataframes(dataframes)
        # Check for duplicates and save them
        duplicated_rows = check_duplicates(concatenated_df)
        output_file = f"{output_prefix}_duplicate.txt"
        save_duplicates(duplicated_rows, output_file)
        if duplicated_rows.empty:
            for i in set(concatenated_df['library_oligomer']):
                concatenated_df.loc[[j==i for j in concatenated_df['library_oligomer']] ,:].to_csv(f"{i}.csv", index=False)
    # If there are errors, save the error log
    else:
        error_file = f"{output_prefix}.err"
        save_errors(error_list, error_file)


if __name__ == "__main__":
    # Setting up argument parser
    parser = argparse.ArgumentParser(description="Concatenate CSV files and check for duplicates.")
    parser.add_argument('-t', metavar = 'filename', dest = 'target', type = str, required = False, default=None, help = 'Editor type (will be examined against reference)')
    parser.add_argument('-p', metavar = 'filename', dest = 'positive', type = str, required = False, default=None, help = 'Editor type (will be examined against reference)')
    parser.add_argument('-n', metavar = 'filename', dest = 'negative', type = str, required = False, default=None, help = 'Editor type (will be examined against reference)')
    parser.add_argument("-o", "--output", required=True, help="Output file prefix for the result files.")

    args = parser.parse_args()

    # Call the main function with the provided directory and output prefix
    main(args.target, args.positive, args.negative, args.output)

