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
argparser.add_argument('-L','--length', metavar = 'int', dest = 'length', type = int, required = False, default ='20', help = 'length of the GuideRNA without PAM')
argparser.add_argument('-B','--bed', metavar = 'file', dest = 'bed', type = str, required = True, help = 'bedFile protein per region')
argparser.add_argument('--gc', dest = 'gc', action='store_true', required = False, help = 'flag not Consider C in GC as affected')
argparser.add_argument('-G','--Genome', metavar = 'file', dest = 'Genome_file', type = str, required = True, help = 'Genome fasta file')
argparser.add_argument('--head', metavar = 'file', dest = 'header', type = str, required = True, help = 'VCF formats contig file')

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
                        print("::error:: One protein overlaps with another. This is not allowed due to collision in the interpretation downstream")
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
                        i.window_start=i.window_start
                        i.window_end=i.window_end
                        editor_df['editing_windowSTART']=[editor_df['start'][j]-args.length+i.window_start-1 if '+' in editor_df.strand[j] else editor_df['start'][j] +len(editor_df.protospacer[j])-i.window_end for j in editor_df.index]
                        editor_df['editing_windowEND']=editor_df['editing_windowSTART'] + (i.window_end-i.window_start)
                        editor_df[['editing_windowSeq', 'editing_window_mutated']] = editor_df.apply(lambda row: MutateWindow(row,i), axis=1).apply(pd.Series)
                        editor_df['nchange']=[sum([not rowed.editing_windowSeq[j]==rowed.editing_window_mutated[j] for j in range(0,len(rowed.editing_window_mutated))]) for indexed, rowed in editor_df.iterrows()]
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
                                        vcf.write('##INFO=<ID=Protospacer,Number=.,Type=String,Description="Protospace"> \n')
                                        vcf.write('##INFO=<ID=STRAND,Number=.,Type=String,Description="Strand"> \n')
                                        vcf.write('##INFO=<ID=PAM,Number=.,Type=String,Description="Protospacer adjacent motif"> \n')
                                        vcf.write('##INFO=<ID=Nchange,Number=.,Type=Integer,Description="Nummber of nucleotide changes"> \n')
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
                                                        INFO='Protospacer='+row['protospacer']+ ';PAM='+ row['PAM'] +';STRAND='+ row.strand + ';Nchange=' + str(row.nchange) + f';Library={args.Name}'))
                        if args.per_variant :
                                ranges=editor_df.copy()[['Chromosome','editing_windowSTART','editing_windowEND','strand']]
                                ranges.rename(columns={'Chromosome':"Chromosome","editing_windowSTART":'Start', 'editing_windowEND':'End'}, inplace=True)
                                ranges['Start']=ranges.Start
                                ranges['End']= ranges.End +1      # To include the ending
                                ranges.index=range(0,len(ranges))
                                ranges['Sequence']=[j for j in editor_df.protospacer]
                                ranges['PAM']=[j for j in editor_df.PAM]
                                ranges['strand']=[j for j in editor_df.strand]
                                ranges['ID']=[j for j in editor_df.ID]
                                unmerged=pr.PyRanges(ranges)
                                with open(args.header, 'r') as infile, open(args.Output+'_'+i.name+'.per_variant.vcf','wt') as vcf :
                                        vcf.write('##fileformat=VCFv4.1\n')
                                        for line in infile:
                                            vcf.write(line)
                                        vcf.write('##INFO=<ID=Library,Number=.,Type=String,Description="Library type"> \n')
                                        vcf.write('##INFO=<ID=Protospacer,Number=.,Type=String,Description="Protospace"> \n')
                                        vcf.write('##INFO=<ID=STRAND,Number=.,Type=String,Description="Strand"> \n')
                                        vcf.write('##INFO=<ID=PAM,Number=.,Type=String,Description="Protospacer adjacent motif"> \n')
                                        if args.gc:
                                                vcf.write('##INFO=<ID=GC_FLAG,Number=.,Type=bool,Description="Mutation effectiveness compromised by GC pattern"> \n')
                                        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                                        output=set()
                                        for ind, row in editor_df.iterrows():
                                                for nuc in range(0,len(row.editing_windowSeq)):
                                                        if row.editing_windowSeq[nuc]==i['before+'] and row.strand=='+':
                                                                output.add((row.Chromosome, row.editing_windowSTART+nuc, i['before+'], i['after+'],row.strand))
                                                        if row.editing_windowSeq[nuc]==i['before-'] and row.strand=='-':
                                                                output.add((row.Chromosome, row.editing_windowSTART+nuc, i['before-'], i['after-'],row.strand))
                                        df = pd.DataFrame(list(output),columns =['chrom', 'pos', 'REF','ALT','strand'],index=range(0,len(output)))
                                        output_df=df.reindex(index=order_by_index(df.index, index_natsorted(zip(df.chrom, df.pos))))
                                        for ind, row in output_df.iterrows():
                                                for nuc in range(0,len(row.editing_windowSeq)):
                                                        if row.editing_windowSeq[nuc]==i['before+'] and row.strand=='+':
                                                                output.add((row.Chromosome, row.editing_windowSTART+nuc, i.BE[0], i.BE[1],row.strand))
                                                        if row.editing_windowSeq[nuc]==i.BE_RC[0] and row.strand=='-':
                                                                output.add((row.Chromosome, row.editing_windowSTART+nuc, i.BE_RC[0], i.BE_RC[1],row.strand))
                                        df = pd.DataFrame(list(output),columns =['chrom', 'pos', 'REF','ALT','strand'],index=range(0,len(output)))
                                        output_df=df.reindex(index=order_by_index(df.index, index_natsorted(zip(df.chrom, df.pos))))
                                        for ind, row in output_df.iterrows():
                                                df=pd.DataFrame({"Chromosome": [str(row.chrom)], "Start": [row.pos],"End":[row.pos+1]})
                                                position= pr.PyRanges(df)
                                                overlap=position.join(unmerged).as_df()
                                                overlap=overlap.loc[overlap['strand']==row['strand']]
                                                if args.gc :
                                                        if row['strand']=='+' and str(Genome_dict[row.chrom][row.pos-1]).upper()=='C' and str(Genome_dict[row.chrom][row.pos-2]).upper()=='G':
                                                               gc_flag=';GC_FLAG=True'
                                                        elif  row['strand']=='-' and str(Genome_dict[row.chrom][row.pos-1]).upper()=='G' and str(Genome_dict[row.chrom][row.pos]).upper()=='C':
                                                                gc_flag=';GC_FLAG=True'
                                                        else :
                                                                gc_flag=';GC_FLAG=False'
                                                else :
                                                        gc_flag=''
                                                if not overlap.empty :
                                                        vcf.write('{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{END}\t{FILTER}\t{INFO}\n'.format(
                                                                chrom = str(row['chrom']),
                                                                pos = str(row['pos']) ,
                                                                vid=str(row.chrom)+'_'+str(row['pos'])+"_"+str(row['REF'])+"_"+str(row['ALT']),
                                                                ref= str(row['REF']) ,
                                                                alt= str(row['ALT']) ,
                                                                END='.',
                                                                FILTER ='.' ,
                                                                INFO='STRAND='+ row.strand+';GuideId='+'|'.join(overlap.ID)+';PAM='+ '|'.join(overlap['PAM'])+';Library={args.Name}'+';protospacers='+'|'.join(overlap.Sequence)+gc_flag))
