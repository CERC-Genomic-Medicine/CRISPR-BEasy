#!/usr/bin/env python3

''' 

'''


import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
import itertools
import pickle
import re
import numpy as np



argparser = argparse.ArgumentParser(description = 'This Software is a tool to design Guiding RNA. \n This scrip allow for the filtering based on off-target activity prediction.')
argparser.add_argument('-I','--OffTarget', metavar = 'file', dest = 'OffTarget', type = str, required = True, help = 'Off Target file')
argparser.add_argument('-P','--PAM', metavar = 'str', dest = 'PAM', type = str, required = False, default = 'NGG', help = 'PAM to ajust CFD score if not NGG')
argparser.add_argument('-M','--MM', metavar = 'str', dest = 'Pickle', type = str, required = False, default = '', help = 'Pickle MM matrix file location in pickle format')
argparser.add_argument('-S','--Score', metavar = 'file', dest = 'Score', type = str, required = True, help = 'Score')
argparser.add_argument('-C','--CountTreshold', metavar = 'integer', dest = 'N', type = int, default=5, required = False, help = 'Off Target count threshold')
argparser.add_argument('-T','--Treshold', metavar = 'integer', dest = 'Thres', type = float, default=1.0, required = False, help = 'Threshold CFD to consider')
argparser.add_argument('-O','--out', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'Ouput file name without extentions')

import itertools

def PAM_to_DICT(pattern):
    """
    Given a nucleotide pattern string that may include ambiguous nucleotide codes,
    return a dictionary with keys equal to every possible sequence of A, T, C, G
    of the same length as 'pattern'. The value for a key is 1.0 if that sequence 
    fulfills (matches) the input pattern, and 0.0 otherwise.

    Ambiguous codes defined in this function:
      - A, T, C, G: represent themselves.
      - N: any nucleotide (A, T, C, G)
      - R: A or G (purines)
      - K: G or T
      - M: A or C

    Parameters:
        pattern (str): A string representing the nucleotide pattern (e.g., 'N', 'AR', 'TKM')

    Returns:
        dict: Dictionary where keys are strings (e.g., 'ATC') and values are floats (1.0 or 0.0).
    """
    # Ensure the pattern is uppercase.
    pattern = pattern.upper()
    
    # Define the mapping from ambiguous codes to allowed nucleotides.
    mapping = {
        'A': ['A'],
        'T': ['T'],
        'C': ['C'],
        'G': ['G'],
        'N': ['A', 'T', 'C', 'G'],
        'R': ['A', 'G'],
        'K': ['G', 'T'],
        'M': ['A', 'C']
    }
    
    # Generate the set of concrete sequences that match the pattern.
    try:
        # For each letter in pattern, look up the allowed bases.
        allowed_bases = [mapping[letter] for letter in pattern]
    except KeyError as e:
        raise ValueError(f"Unrecognized nucleotide code in pattern: {e}")
    
    # All sequences that fulfill the input pattern.
    valid_sequences = {''.join(seq) for seq in itertools.product(*allowed_bases)}
    
    # Generate all possible combinations of A, T, C, G for sequences of this length.
    all_sequences = [''.join(seq) for seq in itertools.product("ATCG", repeat=len(pattern))]
    
    # Create the dictionary: assign 1.0 if the sequence is valid, 0.0 otherwise.
    result_dict = {seq: 1.0 if seq in valid_sequences else 0.0 for seq in all_sequences}
    return result_dict


#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

def calc_cfd(wt,sg,pam, mm_scores, mm_PAM):
	m_wt = re.search('[^ATCG]',wt)
	m_off = re.search('[^ATCG]',sg)
	if (m_wt is None) and (m_off is None):
		score = 1.0
		sg = sg.replace('T','U')
		wt = wt.replace('T','U')
		s_list = list(sg)
		wt_list = list(wt)
		for i,sl in enumerate(s_list):
			if wt_list[i] == sl:
				score*=1
			else:
				key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
				score*= mm_scores[key]
		return (score*mm_PAM[pam])
	else :
		return(-1)


if __name__ == '__main__':
    args = argparser.parse_args()
    OffTarget=pd.read_csv(args.OffTarget,sep='\t',usecols=['cfdOfftargetScore','seqId','guideId','guideSeq','offtargetSeq'], dtype={'cfdOfftargetScore':float,'seqId':str,'guideId':str,'guideSeq':str, 'offtargetSeq':str},na_values='None')
    if args.PAM != 'NGG':
        print('::notice:: CFD scoring is normally based on spCAS9 (Doench et al. 2016), currently the only ajustement is to the protospace adjacent motif mismatch penalty matrix. \n The modification makes any acceptable (i.e. fulfiling the motif) PAM 100% and any non acceptable 0%.')
        mm_PAM=PAM_to_DICT(args.PAM)
        mm_scores = pickle.load(open(args.Pickle,'rb'))
        CFD_recalc = [calc_cfd(row['guideSeq'][:20],row['offtargetSeq'][:20],row['offtargetSeq'][-3:],mm_scores,mm_PAM) for index, row in OffTarget.iterrows()]
        OffTarget['cfdOfftargetScore']=CFD_recalc
    df=OffTarget.loc[OffTarget.cfdOfftargetScore>=args.Thres]
    df.index=df.seqId+df.guideId
    impossible=sum(OffTarget.cfdOfftargetScore==-1)
    if impossible > 0 :
        print(f'::notice:: {impossible} alignements were impossible to compare (assumed less than threshold)')
    figB, axB= plt.subplots(1, 1,figsize=(15,15))
    axB.hist(df.index.value_counts(), color=plt.cm.Paired(0))
    figB.savefig(f'{args.Output}_histogram_offtarget_treshold.pdf',format="pdf",bbox_inches="tight")    
    filtered=df.index.value_counts()>=args.N
    Excluded=df.guideSeq[filtered].unique()
    if len(Excluded) > 0 :
        print(f'::notice:: {len(Excluded)} guides were remove du to CFD score threshold and count')
    Score=pd.read_csv(args.Score,sep='\t')
    Filtered = Score[~Score['targetSeq'].isin(Excluded)]
    pd.Series(Excluded).to_csv(f'{args.Output}.txt',sep='\t',header=False,index=False)
    Filtered.to_csv(f'{args.Output}.filtered_score',sep='\t',header=True,index=False)
