#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import re

# Required columns
REQUIRED_COLUMNS = {'#Uploaded_variation', 'editor', 'Location', 'SYMBOL', 'Consequence'}




def concatenate_dataframes(dataframes):
    """Concatenates a list of DataFrames into one."""
    concat=pd.DataFrame(columns=['#Uploaded_variation','SYMBOL','Consequence','Nchange','editor'])
    if dataframes:
        for sheet_name, df in dataframes.items():
            editor_name=sheet_name.replace('editor - ','')
            df=df.loc[:,['#Uploaded_variation','SYMBOL','Consequence',editor_name+'_Nchange','editor']]
            df.columns=['#Uploaded_variation','SYMBOL','Consequence','Nchange','editor']
            concat=pd.concat([concat,df])
        return concat
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no valid DataFrames are found


def check_duplicates(dataframes):
    """Checks for duplicates based on 'ID' and combination of 'chromosome' and 'position'."""
    if dataframes == {}:
        return pd.DataFrame()  # No duplicates if DataFrame is empty

    # Check for duplicated IDs
    duplicated_rows=pd.DataFrame(columns=['editor','#Uploaded_variation','Location','Allele'])
    for sheet_name, df in dataframes.items():
    # Check for duplicated chromosome-position combinations
        df['ID'] = df['#Uploaded_variation'].astype(str)
        df['chr_pos'] =df['Location'].astype(str)+ '_'+df['Allele']
        duplicated_chr_pos = df[df.duplicated(subset=['chr_pos'], keep=False)]
        duplicated_id = df[df.duplicated(subset=['ID'], keep=False)]
    # Combine both types of duplicated rows
        duplicated_rows = pd.concat([duplicated_rows,duplicated_id, duplicated_chr_pos],axis=0, ignore_index=True).drop_duplicates()
        duplicated_rows=duplicated_rows.loc[:,['editor','#Uploaded_variation','Location','Allele']]
        #['editor','#Uploaded_variation','Location','Allele']

    return duplicated_rows


def save_duplicates(duplicated_rows, output_file):
    """Saves duplicated rows to a CSV file."""
    if not duplicated_rows.empty:
        duplicated_rows.to_csv(output_file,sep='\t', index=False)
        print(f"Duplicated rows have been saved to '{output_file}'.")


def save_errors(error_list, output_file):
    """Saves the list of errors to a file."""
    with open(output_file, 'w') as f:
        for error in error_list:
            f.write(f"{error}\n")
    print(f"Errors have been saved to '{output_file}'.")

def validate_ids(main_df, other_df_dictionnairy, out):
    # Load the main CSV file and extract the IDs
    main_ids = set(main_df.loc[:,'ID'])

    # Prepare a list to collect errors
    errors = []

    # Loop over each TSV file
    for sheet_name, df in other_df_dictionnairy.items():
        # Load the TSV file
        for id in df['#Uploaded_variation']:
            if id not in main_ids:
                errors.append(f"{id}\t{sheet_name}")

    # Write errors or create validated.ok
    if errors:
        with open(f'correspondance_{out}.err', 'w') as error_file:
            error_file.write("ID\tFilename\n")
            error_file.write("\n".join(errors))
        print("Validation failed. Errors written to errors.txt.")

def check_for_pick_flag(dataframes):
    if dataframes =={} :
        return False
    for sheet_name, df in dataframes.items():
        if 'PICK' in df.columns :
                return True
    return False

def validate_sheets(file_path):
    # Load the Excel file and get all sheet names
    xlsx = pd.ExcelFile(file_path)
    sheet_names = xlsx.sheet_names

    # Initialize a list to collect errors
    Errors = []

    # Validate all sheets are either "General" or "editor - {editor_name}"
    for sheet_name in sheet_names:
        if sheet_name != "Library" and not re.match(r"^editor - .+$", sheet_name):
            Errors.append(f"Unexpected sheet name found: '{sheet_name}'")

    if Errors:
        with open(f'Sheets_{args.output}.err', 'w') as error_file:
            for error in Errors:
                error_file.write("\n".join(errors))

    # Extract sheets with names matching "editor - {editor_name}"
    editor_sheets = {
        sheet_name: xlsx.parse(sheet_name)
        for sheet_name in sheet_names
        if sheet_name.startswith("editor - ")
    }
    Errors= []
    for sheet_name, df in editor_sheets.items():
        missing_columns = REQUIRED_COLUMNS - set(df.columns)
        if missing_columns:
            Errors.append(
                f"{args.output}'s Sheet '{sheet_name}' is missing columns: {', '.join(missing_columns)}"
            )
    # Report if there are any unexpected sheet names
    if Errors:
        with open(f'Sheets_{args.output}.err', 'a') as error_file:
            for error in Errors:
                error_file.write("\n".join(Errors))
    
    return editor_sheets

def filter_pick_in_editor_sheets(editor_sheets):
    """
    Filters rows in each DataFrame within the editor_sheets dictionary,
    keeping only rows where the 'consequence' column contains 'PICK'.
    
    Parameters:
        editor_sheets (dict): Dictionary of DataFrames with sheet names as keys.
        
    Returns:
        dict: A new dictionary with filtered DataFrames.
    """
    filtered_sheets = {}
    print("::message:: VEP's PICK flag was detected, Annotations will be filtered for this flag. \n If other annotations are desired, remove the PICK column and deduplicate the annotations according to preference.")
    
    for sheet_name, df in editor_sheets.items():
        if 'PICK' in df.columns:
            filtered_df = df.loc[df['PICK']==1,:]
            filtered_sheets[sheet_name] = filtered_df
        else:
            filtered_sheets[sheet_name] = df
            
    return filtered_sheets

def main(xlsx ,output_prefix):
    """Main function to orchestrate the reading, validation, concatenation, and duplicate check."""
    error_list = []  # List to store errors
    df = pd.read_excel(xlsx, sheet_name="Library")
    dataframes  = validate_sheets(xlsx)
    validate_ids(df, dataframes, output_prefix)

    # If there are no errors, proceed with concatenation and duplicate checking
    if not error_list:
        if check_for_pick_flag(dataframes):
            dataframes = filter_for_pick(dataframes)
        duplicated_rows = check_duplicates(dataframes)
        concatenated_df = concatenate_dataframes(dataframes)
        # Check for duplicates and save them
        if concatenated_df.empty:
            with open(f'{output_prefix}_annotations.csv', 'w') as f:
                pass
        elif duplicated_rows.empty:
            concatenated_df=concatenated_df.rename(columns = {'#Uploaded_variation':'ID'})
            concatenated_df.to_csv(f"{output_prefix}_annotations.csv", index=False)
        else :
            output_file = f"{output_prefix}_duplicate_annotations.err"
            save_duplicates(duplicated_rows, output_file)
    # If there are errors, save the error log
    else:
        error_file = f"{output_prefix}_annotations.err"
        save_errors(error_list, error_file)


if __name__ == "__main__":
    # Setting up argument parser
    parser = argparse.ArgumentParser(description="Concatenate CSV files and check for duplicates.")
    parser.add_argument("--xlsx", required=True, dest='xlsx', help="Directory containing the CSV files.")
    parser.add_argument("-o", "--output", required=True, dest='output', help="Output file prefix for the result files.")

    args = parser.parse_args()

    # Call the main function with the provided directory and output prefix
    main(args.xlsx, args.output)

