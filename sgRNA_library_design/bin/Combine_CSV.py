#!/usr/bin/env python3

import os
import glob
import pandas as pd
import argparse
argparser = argparse.ArgumentParser(description = 'Combines All CSV in the workspace.')
argparser.add_argument('-e', '--editor', metavar = 'name', dest = 'out_editor', type = str, required = True, help = 'Prefix for output files.')

if __name__ == '__main__':
    args = argparser.parse_args()
    dfs = []
    for file in os.listdir(os.getcwd()):
        if file.endswith(".csv"):
            df = pd.read_csv(file)
            print(file)
            dfs.append(df)
    merged_df = pd.concat(dfs, ignore_index=True)
    remove = merged_df[merged_df.duplicated('protospacer', keep=False)]
    CSV = merged_df.drop_duplicates(subset='protospacer', keep=False)
    output_file = args.out_editor + "_final.csv"
    CSV.to_csv(output_file, index=False)
    remove.to_csv(args.out_editor + 'duplicate_remove.txt')
