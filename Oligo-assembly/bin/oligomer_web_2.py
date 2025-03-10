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

argparser.add_argument('--nGuidesPerConcatamer', metavar='int', dest='nGperC', type=int,
                       required=True, help='Number of Guides per Concatamer')

argparser.add_argument('--BSMBI', metavar='str', dest='BSMBI',
                       type=str, required=False, default='CGTCTC', help='Restriction Site')
argparser.add_argument('--PRIMERS', metavar=str, dest='Primers', type=str,
                       required=True, help='Primer forward,reverse')
argparser.add_argument('-F', '--fragments', metavar='Annotation', dest='fragments', type=str, required=False, nargs='*', default=[
                       'CGTCTCACACCG', 'GTTTTGAGACGgactgcCGTCTCcCACCG', 'GTTTaGAGACGggactaCGTCTCgCACCG', 'GTTTcGAGACGcttctcCGTCTCtCACCG', 'GTTTgGAGACG'], help='')

def Assert_Fragments_BSMBI(BSMBI, fragments):
    #makes sure fragments meet restriction requierement
    flags = ''
    BSMBI_reverse = str(Seq(BSMBI).reverse_complement())
    if fragments[0][0:len(BSMBI)] != BSMBI:
        flags = flags + \
            'Fragment 1 doesn\'t contain restriction site (BSMBI) at its beginning \n'
    if len(fragments)>2 :
        for i in range(1, len(fragments)-1):
            if fragments[i].find(BSMBI) == -1:
                flags = flags + 'Fragment  ' + \
                    str(i+1) + ' doesn\'t contain forward restriction site (BSMBI) \n'
            if fragments[i].find(BSMBI_reverse) == -1:
                flags = flags + 'Fragment  ' + \
                    str(i+1) + ' doesn\'t contain reverse restriction site (BSMBI) \n'
    last = fragments[len(fragments)-1]
    if last.find(BSMBI_reverse) != (len(last)-len(BSMBI)):
        flags = flags + \
            'Last fragment doesn\'t contain restriction site (BSMBI) at its end \n'
    return flags

def filter_restriction_sgrna (df, library_type) :
    #test if guides induce additional restriction sites
    test_restriction_sgrna=args.fragments[0] + df.Protospacer +args.fragments[1]
    test_answer=[((i.count(args.BSMBI) +  i.count(BSMBI_reverse))==3) for i in test_restriction_sgrna]
    removed=len(test_restriction_sgrna)-sum(test_answer)
    return df.iloc[test_answer], removed

def deduplicate_dataframes(dataframes: list, column_name: str, second_column_name: str):
    """
    Deduplicates a list of dataframes across all dataframes by the specified column, incorporating reverse complements,
    and collects duplicates into a dictionary.

    Parameters:
        dataframes (list): List of pandas DataFrames to process.
        column_name (str): Column name to check for duplicates across dataframes.
        second_column_name (str): Column name whose values will be grouped by duplicates.

    Returns:
        tuple: (*dataframes, dictionary of {duplicate: [second_column_values]}).
    """
    # Dictionary to store duplicates and their corresponding second column values
    duplicate_dict = defaultdict(list)

    # Precompute Positive_protospace for all dataframes
    for df in dataframes:
        df['Positive_protospace'] = df[column_name]
        df.loc[df['strand'] == '-', 'Positive_protospace'] = df.loc[df['strand'] == '-', column_name].map(
            lambda x: str(Seq(x).reverse_complement())
        )

    # Combine all dataframes into a single dataframe to identify global duplicates
    combined_df = pd.concat(dataframes, ignore_index=True)

    # Identify duplicate sequences using duplicated()
    duplicates = set(combined_df.loc[combined_df.duplicated(subset=['Positive_protospace'], keep=False), 'Positive_protospace'])

    # Collect the second column values for duplicates efficiently
    duplicate_dict = combined_df[combined_df['Positive_protospace'].isin(duplicates)]\
        .groupby('Positive_protospace')[second_column_name]\
        .apply(lambda x: ', '.join(map(str, x)))\
        .to_dict()

    # Deduplicate each dataframe by removing rows with duplicates in Positive_protospace
    deduplicated_dataframes = [df[~df['Positive_protospace'].isin(duplicates)].reset_index(drop=True) for df in dataframes]

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
            





if __name__ == '__main__':
    error_list=[]
    log=[]
    args = argparser.parse_args()
    if len(args.fragments) < (args.nGperC + 1):
        raise Exception("The number of Guides per concatament is more than what is possible with the fragments given")
    frags = args.fragments[0:(args.nGperC)]
    frags.append(args.fragments[-1])
    flags = Assert_Fragments_BSMBI(args.BSMBI, frags) # test BSMBI and fragment are appropriate (fragments contains BSMBI)
    if flags:
        raise Exception(flags)
    BSMBI_reverse = str(Seq(args.BSMBI).reverse_complement())
    primers_reverse = str(Seq(args.Primers.split(',')[1]).reverse_complement())
    primers_forward = args.Primers.split(',')[0]
    Library_list = []
    Target = pd.read_csv(args.Target)
    Target, removed=filter_restriction_sgrna(Target, "target library")
    Target['editor']='NA'
    if args.Positive: 
        Positive = pd.read_csv(args.Positive)
        Positive, removed_positive = filter_restriction_sgrna(Positive, "positive control library")
        instructions=pd.read_csv(args.Positive_instructions, sep=' ', names=['editor' ,'N', 'Consequence'])
        Positive_length= sum([int(i) for i in instructions['N']])
    else :
        Positive = pd.DataFrame(columns = ['ID','Protospacer','Chromosome', 'POSstart', 'strand'])
        Positive_length=0
    if args.Negative:     
        Negative = pd.read_csv(args.Negative)
        Negative, removed_negative = filter_restriction_sgrna(Negative, "negative control library")
        N=len(Negative.ID) if args.negative_number == 0 else args.negative_number
    else :  
        Negative = pd.DataFrame(columns = ['ID','Protospacer','Chromosome', 'POSstart', 'strand'])
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
    SguidePerConcat = len(frags)-1
    remainder = SguidePerConcat-(Length_libraries % SguidePerConcat)
    remainder = remainder if remainder != SguidePerConcat else 0

    log.append("<b> Overall statistics </b>")
    log.append("Fagments : " + ",".join(frags))
    log.append(f"Primers (forw,reverse) : {args.Primers}")
    log.append(f"{len(frags)-1} sgRNA per oligomer")
    log.append(f"{Length_libraries} sgRNA requested in total")
    log.append(f"{len(Target['ID'])} guides in target library. {removed} removed du to restriction sites.")
    if protospacer_overlap :
        log.append(f"{str(sum(len(values) for values in protospacer_overlap.values()))} guides were discarded due to their duplicate protospacers with other guides to be tested (including Reverse Complement)")
    if args.negative_number>=0 and not Negative.empty :
        log.append(f"  {len(Negative['ID'])} guides in negative library. {removed_negative} removed du to restriction sites")
        log.append(f"\t{N} requested in from negative control library")
        if N > len(Negative.ID) :
            error_list.append(f'Too many guides were asked in negative controls library')
        elif (remainder + N) <= len(Negative['ID']) :
            if remainder !=0 and not error_list :
                log.append(f"{remainder} Negative library guides were added to complete the concatemer")
            Library_list.append(Negative.sample(N + remainder,random_state=11,axis=0))
            remainder=0
        elif (remainder + N ) > len(Negative['ID']) and not error_list :
            log.append(f"{(SguidePerConcat - remainder)} Negative library guides were removed due to incomplete concatemer \n Normally in this case negative guides are added but too few were provided")
            Library_list.append(Negative.sample(args.negative_number -(SguidePerConcat - remainder),random_state=11,axis=0))
            remainder=0
        else :
            error_list.append(f" Negative library was too short to overcome the burden of completing concatamer \n Means there would be no Negative controls")
    if  not (Positive.empty  or  p_VEP.empty) :
        log.append(f"{len(Positive['ID'])} guides in positive control library (total). {removed_positive} removed du to restriction sites")
        log.append(f"{sum([int(i) for i in instructions['N']])} total guides requested from positive control library (pooled)")
        Library_list_Positive, Unused, errors = Positive_library_sel(instructions,Positive,p_VEP)
        if errors :
            error_list.extend(errors)
        else :
            Used=pd.concat(Library_list_Positive, ignore_index=False)
            if remainder == 0:
                Library_list.append(Used)
            elif len(Unused['ID'])> remainder:
                Library_list.append(Unused.sample(remainder),random_state=11,axis=0)
                Library_list.append(Used)
                if not error_list:
                    log.append(f"{(remainder)} Positive library guides were added at random (amongs specified editors/consqueneces) due to incomplete concatemer \n Normally in this case negative guides are added but too few were provided")
                remainder=0
            elif len(Used['ID']) > (SguidePerConcat - remainder) :
                Used=Used.sample( len(Used['ID'])- (SguidePerConcat - remainder),random_state=11,axis=0)
                Library_list.append(Used)
                remainder=0
                if not error_list :
                    log.append(f"{(SguidePerConcat - remainder)} Positive library guides were removed at random due to incomplete concatemer \n Normally in this case negative guides are added but too few were provided")
            else :
                error_list.append(f" Positive control libraries was too short to overcome the burden of completing concatamer \n Means there would be no positive controls")
    if remainder != 0 and not error_list:
        Target=Target.sample(len(Target['ID'])- (SguidePerConcat - remainder),random_state=11,axis=0)
        Library_list.append(Target)
        log.append(f"{(SguidePerConcat - remainder)} Target library guides were removed at random due to incomplete concatemer \n Normally we would focus on control library guides but not enough were provided")
    elif remainder >= len(Target['ID']):
        error_list.append(f"Library was too short to overcome the burden of completing concatamer \n Means there would be no guides")
    else :
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
            for i in range(0, len(guides.index), SguidePerConcat):
                GuidesNames = []
                concat = [primers_forward]
                for j in range(0, SguidePerConcat):
                    concat.extend([frags[j], guides.iloc[i+j].Protospacer])
                    GuidesNames.append(guides.iloc[i+j].ID)
                Names = ','.join(GuidesNames)
                concat.extend([frags[SguidePerConcat],primers_reverse])
                concatemere = ''.join(concat)
                out.write('\t'.join([concatemere, Names]) + ' \n')
        log.append(f"  {len(guides.ID) / SguidePerConcat} Oligomers produced.")
        guides.to_csv('Prepared_libary.csv')
    with open("Oligomer.log", "w") as file:
            for item in log:
                file.write(item + "\n")
