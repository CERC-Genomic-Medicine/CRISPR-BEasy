#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
import argparse
import re
import pysam
import sys
from collections import defaultdict
argparser = argparse.ArgumentParser(
    description='This Software is part of a tool to design Guiding RNA. \n This script helps design oligomer/concatemere for array Crispr experiment ')



Splice_consequence = ['splice_donor_variant','splice_acceptor_variant']



## Main Library

argparser.add_argument('-L', '--Library', metavar='name',
                       dest='Target', type=str, required=True, help='Target Study')
argparser.add_argument('--annotation', metavar='name',
                       dest='target_VEP', type=str, required=False,default='', help='Annotation')

## Positive Controls

argparser.add_argument('-P', '--Library_pos', metavar='name',
                       dest='Positive', type=str, required=False, help='Target Study')
argparser.add_argument('--annotation_positive', metavar='name',
                       dest='positive_VEP', type=str, required=False,default='', help='Annotation')
argparser.add_argument('--positive_instructions', metavar='name',
                       dest='Positive_instructions', type=str, required=False,default='', help='Annotation')

## Negative Controls

argparser.add_argument('-N', '--Library_neg', metavar='name',
                       dest='Negative', type=str, required=False, help='Target Study')
argparser.add_argument('--annotation_negative', metavar='name',
                       dest='negative_VEP', type=str, required=False,default='', help='Annotation')
argparser.add_argument('--negative_number', metavar='name',
                       dest='negative_number', type=int, required=False, default=0, help='Annotation')

### General

argparser.add_argument('--First_enzyme', metavar='str', dest='First_enzyme',
                       type=str, required=False, default='CGTCTC', help='Restriction Site')
argparser.add_argument('--Second_enzyme', metavar='str', dest='Sec_enzyme',
                       type=str, required=False, default='GAATTC', help='Restriction Site')
argparser.add_argument('--Scaffold', metavar='str', dest='Scaffold',
                       type=str, required=False, default='GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT', help='Scaffold')
argparser.add_argument('--PRIMERS', metavar=str, dest='Primers', type=str,
                       required=True, help='Primer forward,reverse')

def create_oligomer(row, First_enzyme, Sec_enzyme, Scaffold, Primer_for, Primer_rev):
    oligo = Primer_for + First_enzyme + 'G' + row['Protospacer'] + row['PAM'] + Scaffold + Sec_enzyme + Primer_rev
    #oligo = Primer_for + row['Protospace'] + Scaffold + row['ContextSequence'] + Primer_rev
    return oligo

def deduplicate_dataframes(dataframes: list, column_name: str, second_column_name: str):
    """
    Deduplicates a list of dataframes across all dataframes by the specified column, incorporating reverse complements,
    and collects duplicates into a dictionary.

    Parameters:
        dataframes (list): List of pandas DataFrames to process.
        column_name (str): Column name to check for duplicates across dataframes.
        Second_column_name (str): Column name whose values will be grouped by duplicates.

    Returns:
        tuple: (*dataframes, dictionary of {duplicate: [second_column_values]}).
    """
    # Dictionary to store duplicates and their corresponding second column values
    duplicate_dict = defaultdict(list)

    # Precompute Positive_protospace for all dataframes
    #for df in dataframes:
    #    df['Positive_protospace'] = df[column_name]
    #    df.loc[df['strand'] == '-', 'Positive_protospace'] = df.loc[df['strand'] == '-', column_name].map(
    #        lambda x: str(Seq(x).reverse_complement())
    #    )

    # Combine all dataframes into a single dataframe to identify global duplicates
    combined_df = pd.concat(dataframes, ignore_index=True)

    # Identify duplicate sequences using duplicated()
    duplicates = set(combined_df.loc[combined_df.duplicated(subset=[column_name], keep=False), column_name])

    # Collect the second column values for duplicates efficiently
    duplicate_dict = combined_df[combined_df[column_name].isin(duplicates)]\
        .groupby(column_name)[second_column_name]\
        .apply(lambda x: ', '.join(map(str, x)))\
        .to_dict()

    # Deduplicate each dataframe by removing rows with duplicates in Positive_protospace
    deduplicated_dataframes = [df[~df[column_name].isin(duplicates)].reset_index(drop=True) for df in dataframes]

    return (*deduplicated_dataframes, duplicate_dict)


def Positive_library_sel(instructions,Positive_L,pos_Annotation):
    Used = []
    acceptable=[]
    error_list_positive=[]
    Positve=Positive_L.copy()
    if not instructions.empty :
        for index, inst in instructions.iterrows():
            log.append(f" -  {inst.N} guide(s), of predicted {' or '.join(inst['Consequence'].split(','))} consequence(s) with {inst.editor} editor  requested from positive control library")
            accepted =[]
            input_consequence = inst['Consequence'].split(',')
            for i in input_consequence :
                if i.lower()=="splice_altering":
                    accepted.extend(Splice_consequence)
                else :
                    accepted.extend([s.lower() for s in input_consequence])
            if 'none' in accepted:
                acceptable_guides=pos_Annotation.loc[[row.editor == inst.editor],'ID']
            else :
                acceptable_guides = pos_Annotation[
                                        (pos_Annotation['editor'] == inst.editor) &
                                        (pos_Annotation['Consequence'].apply(lambda x: any(i.lower() in accepted for i in x.split(','))))
                                        ]['ID']
            log.append(f"    -   {len(acceptable_guides)} found")
            if inst.N > len(acceptable_guides):
                error_list_positive.append(f'Too many guides were asked in positive controls library instructions {inst["editor"]} \n Reminder each row are processed independantly and a protospace can only be selected once.')
                continue
            else :
                selected=acceptable_guides.sample(inst.N, random_state=11)
                Used.append(Positive_L.loc[selected.array,:])
                Positive_L=Positive_L.drop(selected.array)
                pos_Annotation = pos_Annotation.drop(pos_Annotation[pos_Annotation['ID'].isin(selected.array)].index)
                acceptable.extend(selected.array)
        if error_list_positive :
            return [],[],error_list_positive
        else : 
            Unused=Positive_L.loc[[P in acceptable for P in Positive_L['ID']], :]
            return Used, Unused, error_list_positive
            

def filter_by_enzyme_oligo(df: pd.DataFrame, first_enzyme: str, sec_enzyme: str) -> pd.DataFrame:
    rev_first = str(Seq(first_enzyme).reverse_complement())
    rev_sec = str(Seq(sec_enzyme).reverse_complement())

    cond = (
        (df['Oligo'].str.count(first_enzyme) == 1) &
        (df['Oligo'].str.count(sec_enzyme) == 1)
    )

    if rev_first != first_enzyme:
        cond &= df['Oligo'].str.count(rev_first) == 0
    if rev_sec != sec_enzyme:
        cond &= df['Oligo'].str.count(rev_sec) == 0

    return df[cond]


if __name__ == '__main__':
    error_list=[]
    log=[]
    args = argparser.parse_args()
    Target = pd.read_csv(args.Target)
    primers_reverse = str(Seq(args.Primers.split(',')[1]).reverse_complement())
    primers_forward = args.Primers.split(',')[0]
    Target["Oligo"] = Target.apply(lambda row: create_oligomer(row, args.First_enzyme, args.Sec_enzyme, args.Scaffold, primers_forward, primers_reverse), axis=1)
    Target = filter_by_enzyme_oligo(Target, args.First_enzyme, args.Sec_enzyme)
    Target['editor']='NA'
    Library_list = []
    if args.Positive: 
        Positive = pd.read_csv(args.Positive)
        Positive["Oligo"] = Positive.apply(lambda row: create_oligomer(row, args.First_enzyme, args.Sec_enzyme, args.Scaffold, primers_forward, primers_reverse), axis=1)
        Positive = filter_by_enzyme_oligo(Positive, args.First_enzyme, args.Sec_enzyme)
        instructions=pd.read_csv(args.Positive_instructions, sep=' ', names=['editor' ,'N', 'Consequence'])
        Positive_length= sum([int(i) for i in instructions['N']])
    else :
        Positive = pd.DataFrame(columns = ['ID','Protospacer','Chromosome', 'POSstart', 'strand','Oligo'])
        Positive_length=0
    if args.Negative:     
        Negative = pd.read_csv(args.Negative)
        Negative["Oligo"] = Negative.apply(lambda row: create_oligomer(row, args.First_enzyme, args.Sec_enzyme, args.Scaffold, primers_forward, primers_reverse), axis=1)
        Negative = filter_by_enzyme_oligo(Negative, args.First_enzyme, args.Sec_enzyme)
        N=len(Negative.ID) if args.negative_number == 0 else args.negative_number
    else :  
        Negative = pd.DataFrame(columns = ['ID','Protospacer','Chromosome', 'POSstart', 'strand','Oligo'])
        N=0
    if args.target_VEP:
        t_VEP = pd.read_csv(args.target_VEP)
        t_VEP = t_VEP.loc[t_VEP['ID'].isin(Target['ID']), :]
    if args.positive_VEP:
        p_VEP = pd.read_csv(args.positive_VEP)
        p_VEP = p_VEP.loc[p_VEP['ID'].isin(Positive['ID']), :]
    if args.negative_VEP:
        n_VEP = pd.read_csv(args.negative_VEP)
        n_VEP = n_VEP.loc[n_VEP['ID'].isin(Negative['ID']), :]
    Target, Positive, Negative, protospacer_overlap = deduplicate_dataframes([Target,Positive,Negative], 'Protospacer', 'ID')
    if protospacer_overlap :
        df = pd.DataFrame.from_dict(protospacer_overlap, orient='index').reset_index()
        df.columns = ['Protospace', 'Associated_IDs']
        df.to_csv('overlap.txt', mode='a', sep='\t', header=False, index=False)
    Target.index=Target.ID
    Positive.index=Positive.ID
    Negative.index=Negative.ID
    N=len(Negative.ID) if args.negative_number == 0 else args.negative_number
    Length_libraries = len(Target['ID']) + N + Positive_length
    log.append("<b> Overall statistics </b>")
    log.append(f"Primers (forw,reverse) : {args.Primers}")
    log.append(f"{Length_libraries} sgRNA requested in total")
    if protospacer_overlap :
        log.append(f"{str(sum(len(values) for values in protospacer_overlap.values()))} guides were discarded due to their duplicate protospacers with other guides")
    if args.negative_number>=0 and not Negative.empty :
        log.append(f"\t{N} requested in from negative control library")
        if N > len(Negative.ID) :
            error_list.append(f'Too many guides were asked in negative controls library')
        elif N <= len(Negative['ID']) :
            Library_list.append(Negative.sample(N,random_state=11,axis=0))
        else :
            error_list.append(f" Negative library was too short to overcome the burden of completing concatamer \n Means there would be no Negative controls")
    if  not (Positive.empty  or  p_VEP.empty) :
        log.append(f"{sum([int(i) for i in instructions['N']])} total guides requested from positive control library (pooled)")
        Library_list_Positive, Unused, errors = Positive_library_sel(instructions,Positive,p_VEP)
        if errors :
            error_list.extend(errors)
        else :
            Used=pd.concat(Library_list_Positive, ignore_index=False)
            Library_list.append(Used)
    Library_list.append(Target)
    log.append(f" {Length_libraries} sgRNA requested in total")
    if error_list :
        with open("Oligomer_errors.err", "w") as file:
            for item in error_list:
                file.write(item + "\n")
    else :
        guides = pd.concat(Library_list, ignore_index=False)
        guides = guides.sample(frac=1,random_state=11)
        with open('Concatemere.txt', 'w') as out:
            for index, row in guides.iterrows():
                out.write('\t'.join([row['Oligo'], row['ID']]) + ' \n')
        guides.to_csv('Prepared_libary.csv')
    with open("Oligomer.log", "w") as file:
            for item in log:
                file.write(item + "\n")
