#!/usr/bin/env python3
import pandas as pd
import argparse

def map_joined_column(row,col, id_map, sep=","):
    data = row[col].split(sep)
    mapped = [id_map[i] for i in data]
    return ','.join(mapped)

def main():
    parser = argparse.ArgumentParser(description="Combine and remap VEP annotation files using a provided ID mapping.")
    parser.add_argument('--vep', nargs='+', required=True, metavar='VEP_FILE',
                        help='Input VEP TSV files to combine')
    parser.add_argument('--ids', required=True, metavar='ID_MAPPING_FILE',
                        help='Path to the TSV file with Protein, start, strand, and ID columns')
    parser.add_argument('--output', required=True, metavar='OUTPUT_FILE',
                        help='Output filename (e.g. TP53_BE4.tsv)')
    args = parser.parse_args()

    # Read and merge VEP files
    vep_frames = [pd.read_csv(f, sep='\t') for f in args.vep]
    vep_combined = pd.concat(vep_frames,ignore_index=True).drop_duplicates()

    # Load ID mapping
    dict_df = pd.read_csv(args.ids, sep='\t')
    dict_df['key'] = dict_df['Protein'].astype(str) + '_' + dict_df['start'].astype(str) + '_' + dict_df['strand'].astype(str)
    id_map = dict_df.set_index('key')['ID'].to_dict()

    # Remap variation IDs
    editor = vep_combined['editor'][0]
    print(vep_combined[editor])
    vep_combined[editor] = vep_combined.apply(lambda row: map_joined_column(row, editor, id_map, sep=","), axis=1)
    #vep_combined=vep_combined.drop(editor, axis=1)
    vep_combined['#Uploaded_variation'] = vep_combined['#Uploaded_variation'].map(id_map).fillna(vep_combined['#Uploaded_variation'])


    # Write output
    vep_combined.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
