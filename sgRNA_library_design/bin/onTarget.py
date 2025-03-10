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
import itertools
from rs3.seq import predict_seq

argparser = argparse.ArgumentParser(description = 'This Software is a tool to design Guiding RNA. \n This scrip allow for the filtering based on off-target activity prediction.')
argparser.add_argument('-I','--input', metavar = 'file', dest = 'Input', type = str, required = True, help = 'CrisprVerse input')
argparser.add_argument('-A','--align', metavar = 'file', dest = 'align', type = str, required = True, help = 'CrisprVerse alignement input')
argparser.add_argument('-C','--cas', metavar = 'str', dest = 'Cas', type = str, help = 'CAS variant')

argparser.add_argument('-O','--out', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'Ouput file name without extentions')


cas_length_dict = {
    "SpCas9": 20,
    "SpCas9-NG": 20,
    "SpCas9-NRRH": 20,
    "SpCas9-NRTH": 20,
    "SpCas9-NRCH": 20,
    "SpCas9(VQR/VRQR)": 20,
    "xCas9": 20,
    "SpG": 20,
    "SpRY": 20,
    "SpCas9-CP1012": 20,
    "SpCas9-CP1028": 20,
    "SpCas9-CP1041": 20,
    "SpCas9-1249": 20,
    "SpCas9-(HF1/HF2)": 20,
    "eSpCas9": 20,
    "HypaCas9": 20,
    "Sniper-Cas9": 20,
    "evoCas9": 20,
    "ScCas9": 20,
    "ScCas9+": 20,
    "Spymac": 20,
    "iSpymac": 20,
    "SaCas9": 22,
    "SaCas9-KKH": 22,
    "St1Cas9-LMD9": 22,
    "SauriCas9": 22,
    "CjCas9": 22,
    "LbCas12a": 23,
    "LbCas12a-RVRR": 23,
    "AsCas12a": 23,
    "enAsCas12a": 23
}

if __name__ == '__main__':
    args = argparser.parse_args()
    CrisprVerse=pd.read_csv(args.Input,sep='\t')
    CrisprVerseAlign = pd.read_csv(args.align, sep='\t')
    if len(CrisprVerse["ExtendedSqueuence"][0]) == 30 :
        if args.Cas != "SpCas9" :
            print('::notice:: RS3 scoring is based on spCAS9 (Hsu et al, 2013) and (Chen et al. 2013), since your cas has the same length it is Calculable and thus we provide it on an Experimental basis')
        context_seqs=CrisprVerse["ExtendedSqueuence"]
        hsu2013 = predict_seq(context_seqs, sequence_tracr='Hsu2013')
        Chen2013 = predict_seq(context_seqs, sequence_tracr='Chen2013')
        CrisprVerse["RS3_hsu2013"]=hsu2013
        CrisprVerse["RS3_Chen2013"]=Chen2013
        CrisprVerse_Spacer_site = CrisprVerse['spacer'] +'_'+ CrisprVerse['pam_site'].astype(str)
        CrisprVerseAlign_Spacer_site = CrisprVerseAlign['spacer'] +'_'+ CrisprVerseAlign['pam_site'].astype(str)
        CrisprVerseAlign = CrisprVerseAlign.loc[~CrisprVerseAlign_Spacer_site.isin(CrisprVerse_Spacer_site)]
        score_dict = CrisprVerseAlign.groupby('spacer')['CFD_score'].sum().to_dict()
        CrisprVerse['sgRNA_CFD_score'] = [ 100.* 100.0 / score_dict[sp] for sp in CrisprVerse['spacer'] ]
    else :
        print(f'::notice::Off Target for {args.cas} is not currently calculable, in the future we endeavour to implement all Scores')
    CrisprVerse.to_csv(f'{args.Output}.scored.tsv',sep='\t',header=True,index=False)
    CrisprVerseAlign.to_csv(f'{args.Output}_alignements.tsv',sep='\t',header=True,index=False)
