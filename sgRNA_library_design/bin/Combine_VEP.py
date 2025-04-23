#!/usr/bin/env python3

import os
import glob
import pandas as pd
import argparse
argparser = argparse.ArgumentParser(description = 'Merges all vep.tsv files within the workspace.')
argparser.add_argument('-e', '--editor', metavar = 'name', dest = 'out_editor', type = str, required = True, help = 'Prefix for output files.')
argparser.add_argument('-r', '--remove', metavar = 'name', dest = 'remove', type = str, required = True, help = 'ID to remove.')

if __name__ == '__main__':
    args = argparser.parse_args()
    remove = pd.read_csv(args.remove)
    dfs = []
    for file in os.listdir(os.getcwd()):
        if file.endswith(".tsv"):
            df = pd.read_csv(file, sep='\t')
            df.loc[:,"library"] = file.split("_")[0] + "_Library"
            dfs.append(df)
    merged_df = pd.concat(dfs, ignore_index=True)
    VEP = merged_df.loc[~merged_df['ID'].isin(remove['ID'])]
    output_file = args.out_editor + "_final_vep.txt"
    VEP.to_csv(output_file, index=False)
