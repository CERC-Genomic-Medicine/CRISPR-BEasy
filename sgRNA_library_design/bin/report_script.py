#!/usr/bin/env python3

import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import argparse
from scipy.stats import rankdata
from scipy.stats import binomtest
import sys
import re
import warnings



parser = argparse.ArgumentParser(description='This script produce a lolipop plot (with domain/feature annotations) base on a MaGeCK, a Variant Effect Prediction file (VEP), and bed style file')
parser.add_argument('-c', '--csv', metavar='FILE', dest='csv', required=True, type=str, help='')
parser.add_argument('-v', '--vep', metavar='FILE', dest='vep_files', required=True, type=str, nargs='+', help='tab delimited variant effect prediction file (VEP) with id corresponding to input')



plt.rcParams.update({'font.size': 18})

variant_consequences_mapping = {
'missense_variant': 'missense',
'intron_variant': 'non_coding',
'downstream_gene_variant': 'non_coding',
'NMD_transcript_variant': 'non_coding',
'upstream_gene_variant': 'non_coding',
'3_prime_UTR_variant': 'non_coding',
'synonymous_variant': 'synonymous',
'non_coding_transcript_exon_variant': 'non_coding',
'splice_region_variant': 'non_coding',
'splice_polypyrimidine_tract_variant': 'non_coding',
'stop_gained': 'nonsense',
'5_prime_UTR_variant': 'non_coding',
'regulatory_region_variant': 'non_coding',
'splice_donor_variant': 'splice',
'splice_acceptor_variant': 'splice',
'non_coding_transcript_variant': 'non_coding',
'splice_donor_region_variant': 'non_coding',
'splice_donor_5th_base_variant': 'non_coding',
'TF_binding_site_variant': 'non_coding',
'start_lost': 'nonsense',
'incomplete_terminal_codon_variant': 'non_coding',
'stop_lost':'nonsense'
}


def get_most_damaging_consequence(consequence_str, variant_consequences_mapping):
    # Split the comma-separated consequence string
    consequences = consequence_str.split(',')
    
    # Map each consequence to the corresponding value using the variant_consequences_mapping
    mapped_consequences = [variant_consequences_mapping.get(c.strip(), 'non_coding') for c in consequences]
    
    # Define the damage hierarchy
    damage_hierarchy = ['splice', 'nonsense', 'missense', 'synonymous', 'non_coding']
    
    # Sort the mapped consequences based on the damage hierarchy and return the most damaging one
    most_damaging = sorted(mapped_consequences, key=lambda x: damage_hierarchy.index(x))[0]
    
    return most_damaging


def createBED(df):
    # Ensure the dataframe has the necessary columns
    required_columns = ['ID', 'Chromosome', 'POSstart', 'PAM','gRNA_seq_POSstrand', 'strand']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in dataframe")
    
    # Process each row and construct the output
    df['protein'] = df['ID'].apply(lambda x: x.split('_')[0])
    df['length']=df['gRNA_seq_POSstrand'].apply(len)
    df['POSend'] = df['POSstart'] + df['length']
    
    # Create the required output columns
    output_df = pd.DataFrame({
        'Chromosome': df['Chromosome'],
        'POSstart': df['POSstart'],
        'POSend': df['POSend'],
        'ID': df['ID'],
        'length': df['length'],
        'strand': df['strand']
    })
    
    return output_df

def Base_summary(df):
    # Ensure the dataframe has the necessary columns including 'library'
    required_columns = ['ID', 'Chromosome', 'POSstart', 'library','PAM','gRNA_seq_POSstrand']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in dataframe")

    # Group by 'protein', count instances, get min and max POSstart/POSen
    df['POSend']=df['PAM']+df['gRNA_seq_POSstrand']
    df['POSend'] = df['POSstart'] + df['length']
    df['protein'] = df['ID'].apply(lambda x: x.split('_')[0])
    summary_df = df.groupby('protein').agg(
        sgRNA=('protein', 'count'),
        chrom=('Chromosome', 'first'),
        min_POSstart=('POSstart', 'min'),
        max_POSend=('POSend', 'max'),
        library=('library', 'first')  # Assuming one unique value per group
    ).reset_index()
     
    # Create the range column
    summary_df['range'] = summary_df['chrom'] + ":" + summary_df['min_POSstart'].astype(str) + "-" + summary_df['max_POSend'].astype(str)
    
    # Rename split_ID to protein
    
    # Reorder columns
    summary_df = summary_df[['protein','range' , 'sgRNA', 'library']]
    return summary_df

def process_vep_summary(editor, summary_df, df_VEP, variant_consequences_mapping):
    # Step 1: Map the 'Consequence' column in df_VEP to the most damaging consequence
    if "PICK" in df_VEP.columns:
        df_VEP=df_VEP[df_VEP['PICK']=='1']
        Filter='VEP PICK'
    elif "CANONICAL" in df_VEP.columns:
        df_VEP=df_VEP[df_VEP['CANONICAL']=='YES']
        Filter='CANONICAL'
    else :
        Filter='ALL annotations'
    df_VEP['most_damaging_consequence'] = df_VEP['Consequence'].apply(
        lambda x: get_most_damaging_consequence(x, variant_consequences_mapping)
    )
    
    # Step 2: Extract 'protein' from the 'ID' column (using split('_')[0])
    df_VEP['protein'] = df_VEP['ID'].apply(lambda x: x.split('_')[0])
    
    # Step 3: Count the number of mapped consequences per type per protein
    consequence_counts = df_VEP.groupby(['protein', 'most_damaging_consequence']).size().unstack(fill_value=0)
    
    # Step 4: Merge the consequence counts with the summary dataframe based on the 'protein' column
    summary_df_library=summary_df.pop('library')
    merged_df = pd.merge(summary_df, consequence_counts, how='left', left_on='protein', right_index=True)
    
    # Fill NaN values resulting from the merge with 0 (for proteins that have no corresponding consequences)
    merged_df.fillna(0, inplace=True)
    
    # Step 5: Rename the consequence columns to include the editor name
    for consequence in consequence_counts.columns:
        merged_df.rename(columns={consequence: f'{editor}_{consequence}'}, inplace=True)
    
    merged_df.loc[:,'filter']=Filter
    merged_df.loc[:,'library']=summary_df_library
    return merged_df

if __name__ == '__main__':
    args = parser.parse_args()
    #### Input Parsing
    csv = pd.read_csv(args.csv)
    bed = createBED(csv)
    Summary = Base_summary(csv)
    for vep_file in args.vep_files:
        df_VEP=pd.read_csv(vep_file)
        Summary = process_vep_summary(vep_file.split("_")[0], Summary, df_VEP, variant_consequences_mapping)
    Summary.columns=[i.replace('.vep', ' ') for i in Summary.columns]
    cols = [col for col in Summary.columns if col != 'filter'] + ['filter']
    Summary = Summary[cols]
    cols = [col for col in Summary.columns if col != 'library'] + ['library']
    Summary = Summary[cols]
    Summary.rename(columns={'protein': 'Gene', 'sgRNA': 'Total sgRNAs'}, inplace=True)
    Summary.to_json('report_summary.json',orient='records')
    bed.to_csv('guides.bed', header=False, index=False, sep='\t')
    # To display the output dataframe
