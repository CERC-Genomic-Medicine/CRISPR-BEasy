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
    if len(CrisprVerse["spacer"][0]) == 20 and len(CrisprVerse["PAM"][0]) == 3 :
        if args.Cas != "SpCas9" :
            print('::notice:: Rule set 3 (Hsu et al ; Chen et al. 2013) was developed for SpCas9. We enable scoring for all SpCas9 variants on an experimental basis')
        context_seqs=CrisprVerse["ExtendedSqueuence"]
        hsu2013 = predict_seq(context_seqs, sequence_tracr='Hsu2013')
        Chen2013 = predict_seq(context_seqs, sequence_tracr='Chen2013')
        CrisprVerse["RS3_hsu2013"]=hsu2013
        CrisprVerse["RS3_Chen2013"]=Chen2013
    else :
        print(f'::notice:: Rule set 3 (Hsu et al ; Chen et al. 2013)  is not currently calculated for {args.Cas}  since it is too different from SpCas9 (the experimental basis for CFD scoring). We enable scoring for all SpCas9 variants on an experimental basis')
    if len(CrisprVerse["spacer"][0]) == 20 and 'CFD_score' in CrisprVerseAlign.columns :
#                CrisprVerse_Spacer_site = CrisprVerse['spacer'] +'_'+ CrisprVerse['pam_site'].astype(str)
#                CrisprVerseAlign_Spacer_site = CrisprVerseAlign['spacer'] +'_'+ CrisprVerseAlign['pam_site'].astype(str)
#                CrisprVerseAlign = CrisprVerseAlign.loc[~CrisprVerseAlign_Spacer_site.isin(CrisprVerse_Spacer_site)]
        score_dict = CrisprVerseAlign.groupby('spacer')['CFD_score'].sum().to_dict()
        CrisprVerse['sgRNA_CFD_score'] = [
                "NA" if score_dict[sp] == 0 else 100.0 * 100.0 / score_dict[sp] 
                for sp in CrisprVerse['spacer']
                ]
        if args.Cas != "SpCas9" :
            print('::notice:: CFD scoring was developed for SpCas9 (Doench et al., Nat Biotechnol, 2016). We enable scoring for all SpCas9 variants by adjusting the PAM mismatch penalty matrix.')
    else :
        print(f'::notice:: CFD Off-Target is not currently calculated for {args.Cas} since it is too different from SpCas9 (the experimental basis for CFD scoring). We enable scoring for all SpCas9 variants by adjusting the PAM mismatch penalty matrix.')
    CrisprVerse.to_csv(f'{args.Output}.scored.tsv',sep='\t',header=True,index=False)
    CrisprVerseAlign.to_csv(f'{args.Output}_alignements.tsv',sep='\t',header=True,index=False)
