#!/usr/bin/env python3

''' 

'''


import pandas as pd
import os
import argparse
import itertools
import pickle
import re
import numpy as np
import rs3
import sys
import itertools
from rs3.seq import predict_seq

argparser = argparse.ArgumentParser(description = 'This Software is a tool to design Guiding RNA. \n This script determines the offtarget score of a guide and the RS3 score.')
argparser.add_argument('-I','--input', metavar = 'file', dest = 'Input', type = str, required = True, help = 'CrisprVerse input')
argparser.add_argument('-A','--align', metavar = 'file', dest = 'align', type = str, required = True, help = 'CrisprVerse alignement input')
argparser.add_argument('--CFD', dest = 'CFD', action='store_true', help = 'Whether or not to Calculate_CFD')

argparser.add_argument('-O','--out', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'Ouput file name without extentions')


if __name__ == '__main__':
    args = argparser.parse_args()
    CrisprVerse=pd.read_csv(args.Input,sep='\t')
    CrisprVerseAlign = pd.read_csv(args.align, sep='\t')
    if len(CrisprVerse["spacer"])==0 :
        print(f'::error:: No spacers were availlable for on target evaluation')
        sys.exit(1)
    if len(CrisprVerse["spacer"][0]) == 20 and len(CrisprVerse["PAM"][0]) == 3 :
        context_seqs=CrisprVerse["ExtendedSqueuence"]
        hsu2013 = predict_seq(context_seqs, sequence_tracr='Hsu2013')
        Chen2013 = predict_seq(context_seqs, sequence_tracr='Chen2013')
        CrisprVerse["RS3_hsu2013"]=hsu2013
        CrisprVerse["RS3_Chen2013"]=Chen2013
    if args.CFD and 'CFD_score' in CrisprVerseAlign.columns:
#                CrisprVerse_Spacer_site = CrisprVerse['spacer'] +'_'+ CrisprVerse['pam_site'].astype(str)
#                CrisprVerseAlign_Spacer_site = CrisprVerseAlign['spacer'] +'_'+ CrisprVerseAlign['pam_site'].astype(str)
#                CrisprVerseAlign = CrisprVerseAlign.loc[~CrisprVerseAlign_Spacer_site.isin(CrisprVerse_Spacer_site)]
        print(CrisprVerse)
        score_dict = CrisprVerseAlign.groupby('spacer')['CFD_score'].sum().to_dict()
        print(score_dict)
        CrisprVerse['sgRNA_CFD_score'] = [
                "NA" if sp not in score_dict.keys() or  score_dict[sp] == 0 else 100.0 / score_dict[sp] 
                for sp in CrisprVerse['spacer']
                ]
    CrisprVerse.to_csv(f'{args.Output}.scored.tsv',sep='\t',header=True,index=False)
    CrisprVerseAlign.to_csv(f'{args.Output}_alignements.tsv',sep='\t',header=True,index=False)
