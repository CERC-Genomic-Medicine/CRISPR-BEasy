#!/usr/bin/env python3
'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.1
YEAR: 2024

'''
import sys
print (sys. version)

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
import argparse
import re
import gffutils
import pysam
import gzip
import pyranges as pr
from natsort import natsorted, index_natsorted, order_by_index


argparser = argparse.ArgumentParser(description = 'This Software is a part of a pipeline tool to design Guiding RNA. \n This script creates an summary of sgRNA and optionnaly can provide editor specific annotation (VCF file)')
argparser.add_argument('-E','--Editor', metavar = 'file name', dest = 'Editor', type = str, required = False, default=None, help = 'Editor type (will be examined against reference)')
argparser.add_argument('-S','--ScoreGuide', metavar = 'file name', dest = 'scoreGuide', type = str, required = True, help = 'sgRNA scored file')
argparser.add_argument('-O','--Output', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'prefix of the vcf')
argparser.add_argument('-n','--Name', metavar = 'string', dest = 'Name', type = str, required = False, default='Target', help = 'Library name')
argparser.add_argument('-X','--exclude', metavar = 'file', dest = 'excludes', type = str, required = False, default ='', help = 'List of Guides to exclude')
argparser.add_argument('-V','--Per_Variant', dest = 'per_variant', action='store_true', help = 'flag to produce a per variant VCF suitable for VEP')
argparser.add_argument('-R','--Per_sgRNA', dest = 'per_guide', action='store_true', help = 'flag to produce a per guide vcf file suitable for VEP')
argparser.add_argument('-B','--bed', metavar = 'file', dest = 'bed', type = str, required = True, help = 'bedFile protein per region')
argparser.add_argument('--gc', dest = 'gc', action='store_true', required = False, help = 'flag not Consider C in GC as affected')
argparser.add_argument('-G','--Genome', metavar = 'file', dest = 'Genome_file', type = str, required = True, help = 'Genome fasta file')
argparser.add_argument('--head', metavar = 'file', dest = 'header', type = str, required = True, help = 'VCF formats contig file')

def trim_window(row):
    seq1 = row.editing_windowSeq
    seq2 = row.editing_window_mutated
    start_diff = next((j for j in range(len(seq1)) if seq1[j] != seq2[j]), None)
    end_diff = next((j for j in range(len(seq1)-1, -1, -1) if seq1[j] != seq2[j]), None)
    if start_diff is None or end_diff is None or start_diff > end_diff:
        # No mutation
        return pd.Series({
            'editing_windowSTART': row.editing_windowSTART,
            'editing_windowEND': row.editing_windowEND,
            'editing_windowSeq': seq1,
            'editing_window_mutated': seq2,
            'nchange': 0
        })
    new_start = row.editing_windowSTART + start_diff
    new_end = row.editing_windowSTART + end_diff + 1
    return pd.Series({
        'editing_windowSTART': new_start,
        'editing_windowEND': new_end,
        'editing_windowSeq': seq1[start_diff:end_diff + 1],
        'editing_window_mutated': seq2[start_diff:end_diff + 1],
        'nchange': sum(1 for a, b in zip(seq1[start_diff:end_diff + 1], seq2[start_diff:end_diff + 1]) if a != b)
    })


def MutateWindow(rower,edit):
#        print(str(rower.editing_windowSTART-1))
#        print(type(rower.editing_windowSTART - 1))
#        print(rower.editing_windowEND)
#        print(type(rower.editing_windowEND))
        windowSeq=str(Genome_dict[str(rower.Chromosome)][(rower.editing_windowSTART-1):(rower.editing_windowEND)].seq).upper()
        if args.gc and ((rower.strand=='+' and edit['before+'].upper()=='C') or (rower.strand=='-' and edit['before-'].upper()=='G')) :
                returned=""
                for nuc in range(rower.editing_windowSTART-1,rower.editing_windowEND):
                        if rower.strand=='+' and str(Genome_dict[str(rower.Chromosome)][nuc]).upper() ==edit['before+'].upper() and Genome_dict[str(rower.Chromosome)][nuc-1].upper()!='G' :
                                returned=returned+edit['after+']
                        elif rower.strand=='-' and str(Genome_dict[str(rower.Chromosome)][nuc]).upper() ==edit['before-'].upper() and Genome_dict[str(rower.Chromosome)][nuc+1].upper()!='C' :
                                returned=returned+edit['after-']
                        else:
                                returned=returned+Genome_dict[str(rower.Chromosome)][nuc]
        else :
                returned=windowSeq.replace(edit['before+'],edit['after+']) if rower.strand=='+' else windowSeq.replace(edit['before-'],edit['after-'])
        return windowSeq, returned

if __name__ == '__main__':
        args = argparser.parse_args()
        if not args.Editor == None :
                editor=pd.read_csv(args.Editor,index_col='name',sep='\\s+', names=['name','window_start','window_end','before+','after+'])
                editor['after-']=[str(Seq(y).reverse_complement()) for y in editor['after+']]
                editor['before-']=[str(Seq(y).reverse_complement()) for y in editor['before+']]
#                print(editor)
        if args.gc :
                Genome_dict = SeqIO.to_dict(SeqIO.parse(args.Genome_file, "fasta"))
        scoreGuides=pd.read_csv(args.scoreGuide,sep='\t',header=0,index_col=None)
        scoreGuides.drop_duplicates(subset='spacer', inplace=True, keep=False)
        scoreGuides['targetSeq_plusStrand']=[str(Seq(row.protospacer).reverse_complement()) if '-' in  row['strand']  else row['protospacer'] for index, row in scoreGuides.iterrows()]
        ### Find protein corresponding
        scoreGuides['start_seqId']=scoreGuides['start']
#        print(scoreGuides)
        ### Find proiten
        bed=pr.read_bed(args.bed)
        protein=[]
        for index, row in scoreGuides.iterrows():
                df=pd.DataFrame({"Chromosome": [row['Chromosome']], "Start": [row.start],"End":[row.end]})
                position= pr.PyRanges(df)
                names=[i for i in position.join(bed).Name]
                if len(set(names))>1:
                        print('::error:: One protein overlaps with another.')
                        raise Exception('One protein overlaps with anoter.')
                else :
                        protein.extend(list(set(names)))
        scoreGuides['Protein']=protein
        prots=set(protein)
        scoreGuides['ID'] = scoreGuides['Protein'] + '_' + scoreGuides['start'].astype(str) +'_'+scoreGuides['strand']
        scoreGuides.index=scoreGuides['ID']
        with  open(args.Output+'_general.csv', 'wt') as general :
                # Columns to include (preserving original spelling)
                base_columns = [
                'ID', 'Protospacer', 'PAM', 'gRNA_seq_POSstrand',
                'Chromosome', 'POSstart', 'strand', 'library', 'Protein'
                ]
                # Additional non-conflicting columns (exclude those already in base_columns)
                extra_columns = [
                'spacer', 'protospacer', 'start', 'end', 'pam_site', 'cut_site', 'percentGC',
                'polyA', 'polyC', 'polyG', 'polyT', 'startingGGGGG',
                'n0', 'n1', 'n2', 'n3', 'ContextSequence',
                'RS3_hsu2013', 'RS3_Chen2013', 'sgRNA_CFD_score'
                ]
                final_columns = base_columns + [col for col in extra_columns if col not in base_columns and col in scoreGuides.columns]
                general.write(','.join(final_columns) + '\n')
                for index, row in scoreGuides.iterrows():
                        row_dict = {
                        'ID': row.ID,
                        'Protospacer': row['protospacer'],
                        'PAM': row['PAM'],
                        'gRNA_seq_POSstrand': str(Seq(row['protospacer']).reverse_complement()) if row['strand'] == '-' else row['protospacer'],
                        'Chromosome': row['Chromosome'],
                        'POSstart': str(row['start']),
                        'strand': row['strand'],
                        'library': args.Name
                        }
                        # Add extra fields if available
                        for col in final_columns:
                                if col not in row_dict:
                                        row_dict[col] = str(row[col]) if col in row and not pd.isna(row[col]) else ''
                        general.write(','.join(str(row_dict[col]) for col in final_columns) + '\n')
        ### producting Edditor specific files
        if not args.Editor == None :
                for index, i in editor.iterrows() :
                        editor_df=scoreGuides.copy()
                        editor_df['editing_windowSTART']=[editor_df['start'][j]+i.window_start-1 if '+' in editor_df.strand[j] else editor_df['end'][j] -i.window_end+1 for j in editor_df.index]
                        editor_df['editing_windowEND']=editor_df['editing_windowSTART'] + (i.window_end-i.window_start)
                        editor_df[['editing_windowSeq', 'editing_window_mutated']] = editor_df.apply(lambda row: MutateWindow(row,i), axis=1).apply(pd.Series)
                        editor_df['MutationFullWindow']=editor_df['editing_windowSeq'].astype(str) +'/'+ editor_df['editing_window_mutated'].astype(str)
                        editor_df[['editing_windowSTART', 'editing_windowEND', 'editing_windowSeq', 'editing_window_mutated', 'nchange']] = editor_df.apply(trim_window, axis=1)
                        editor_df=editor_df.reindex(index=order_by_index(editor_df.index, index_natsorted(zip(editor_df.Chromosome, editor_df.start))))
                        with open(args.Output+'_'+i.name +'_empties' + '.txt', 'w') as empties:
                                for ind, row in editor_df.iterrows():
                                        if row.editing_windowSeq ==row.editing_window_mutated :
                                                empties.write(str(row.name) +'\n')
                        if args.per_guide :
                                with open(args.header, 'r') as infile, open(args.Output+'_'+i.name+'.per_guides.vcf','wt') as vcf :
                                        vcf.write('##fileformat=VCFv4.1\n')
                                        for line in infile:
                                            vcf.write(line)
                                        vcf.write('##INFO=<ID=Protospacer,Number=.,Type=String,Description="Protospacer"> \n')
                                        vcf.write('##INFO=<ID=STRAND,Number=.,Type=String,Description="Strand"> \n')
                                        vcf.write('##INFO=<ID=PAM,Number=.,Type=String,Description="Protospacer adjacent motif"> \n')
                                        vcf.write('##INFO=<ID=Nchange,Number=.,Type=Integer,Description="Nummber of nucleotide changes"> \n')
                                        vcf.write('##INFO=<ID=MutationFullWindow,Number=.,Type=String,Description="Mutation REF/ALT"> \n')
                                        vcf.write('##INFO=<ID=Library,Number=.,Type=String,Description="Library type"> \n')
                                        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                                        for ind, row in editor_df.iterrows():
                                                vcf.write('{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{END}\t{FILTER}\t{INFO}\n'.format(
                                                        chrom = str(row['Chromosome']),
                                                        pos = str(row['editing_windowSTART']) ,
                                                        vid=str(row.ID),
                                                        ref= str(row['editing_windowSeq']) ,
                                                        alt= str(row['editing_window_mutated']) ,
                                                        END='.',
                                                        FILTER ='.' ,
                                                        INFO=f'Protospacer={row['protospacer']};PAM={row['PAM']};STRAND={row.strand};Nchange= {str(row.nchange)} ;MutationFullWindow={row.MutationFullWindow};Library={args.Name}'))
