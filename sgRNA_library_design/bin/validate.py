#!/usr/bin/env python3
'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.1
YEAR: 2024
'''

import sys
print(sys.version)
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
from pathlib import Path
import argparse
import re
import gffutils
import pyranges as pr
import warnings

argparser = argparse.ArgumentParser(
    description='This software produces a bed file corresponding to the regions of interest as defined by a protein (corresponding to a genome) and the desired Feature. This represent the first step (potentially optional) in producing a library design crispr Array')
argparser.add_argument('-T', '--TargetFile', metavar='file', dest='target_file',
                       type=str, required=True, help='List of protein to be examined')
argparser.add_argument('-P', '--PositiveFile', metavar='file', dest='positive_file',
                       type=str, required=True, help='List of protein to be examined as Postivie Controls')
argparser.add_argument('-N', '--NegativeFile', metavar='file', dest='negative_file',
                       type=str, required=True, help='List of protein to be examined as Negative Controls')
argparser.add_argument('-G', '--Genome_database', metavar='file', dest='Genome', type=str, required=True,
                       help='Genome feature database file (either gff3 gtf gff format or gffutil database (if already created))')
argparser.add_argument('-e', '--Genome_bed_encode', metavar='file', dest='encode_bed', type=str, required=False,
                       help='Bed with non N chromosome regions')
argparser.add_argument('-c', '--Genome_bed_chromosome', metavar='file', dest='chromosome_bed', type=str, required=True,
                       help='Bed with chrosome ranges')
argparser.add_argument('--editors', metavar='file', dest='editors', type=str, required=True,
                       help='Editors with their caracteristics')
argparser.add_argument('-o', '--out', metavar='file', dest='Output',
                       type=str, required=False, default='out', help='Output')
argparser.add_argument('-B', '--Border', metavar='bp', dest='border', type=int, required=False,
                       default=30, help='Border to be included around feature as define in feature options')
argparser.add_argument('-F', '--Features', metavar='Annotation', dest='features', type=str,
                       required=False, nargs='*', default=['CDS'], help='Feature Types to be analysed (default exons)')
argparser.add_argument('-L','--limit', dest = 'limit', type=int, required = False, help = 'base pair Limit')
argparser.add_argument('-I','--isoform', dest = 'isoform', type=str, required = False, help = 'Which Isofrom filter shoudl be applied (None, Canonical MANE ')
argparser.add_argument('--protist', dest = 'protist', action='store_true', required = False, help = 'Protist do not have transcript level because there is no alternate isoform')

def validate_editor_file(file_path):
    """
    Validates that the editor file is in the correct format:
    - 5 tab-delimited columns
    - 2nd and 3rd columns are whole numbers (can be negative)
    - 4th and 5th columns are nucleotides (A, C, T, G) and 4th != 5th

    If validation fails, errors are logged to 'editors.err'. If valid, an empty 'editor.ok' file is created.
    If the file is blank, an empty 'editor.blank' file is created.

    :param file_path: Path to the editor file to validate
    """
    errors_edit = []
    editors = []
    notaccepted = ['\\', '_']

    with open(file_path, 'r') as file:
        lines = file.readlines()

        if len(lines) == 0:
            # If the file is blank, create editor.blank
            with open('editor.blank', 'w') as blank_file:
                pass
                return
        for line_num, line in enumerate(lines, start=0):
            columns = str(line.strip()).split(" ")
            # Check if there are exactly 5 columns
            if len(columns) != 5:
                errors_edit.append(f"Line {line_num}: Expected 5 columns, found {len(columns)}")
                errors_edit.append(f"{line}")
                errors_edit.append(f"{columns}")
                continue
                # Check if 2nd and 3rd columns are whole numbers (can be negative)
            if any(ele in columns[0] for ele in notaccepted):
                errors_edit.append(f"Line {line_num}: editor name cannot contain \'_\' or \'\\\' found \'{columns[0]}\'")
            editors.append(columns[0])
            try:
                int(columns[1])
                int(columns[2])
            except:
                errors_edit.append(f"Line {line_num}: Columns 2 and 3 must be whole numbers, found '{columns[1]}' and '{columns[2]}'")
                continue
                # Check if 4th and 5th columns are valid nucleotides (A, C, T, G)
            valid_nucleotides = {'A', 'C', 'T', 'G'}
            if columns[3] not in valid_nucleotides or columns[4] not in valid_nucleotides:
                errors_edit.append(f"Line {line_num}: Columns 4 and 5 must be nucleotides (A, C, T, G), found '{columns[3]}' and '{columns[4]}'")
                continue
                # Check if 4th column is not equal to 5th column
            if columns[3] == columns[4]:
                errors_edit.append(f"Line {line_num}: Columns 4 and 5 must not be equal, found '{columns[3]}' and '{columns[4]}'")
        # If there are errors, write them to editors.err
    if len(editors) != len(set(editors)):
        duplicates = [item for item in editors if strings.count(item) > 1]
        duplicates = set(duplicates)  # Remove repeated duplicates in the result
        errors_edit.append(f"editors listed contains duplicates {duplicates}")
    if errors_edit:
        with open('editors.err', 'w') as error_file:
            error_file.write('\n'.join(errors_edit))
    else:
        # If no errors, create an empty editor.ok file
        with open('editor.ok', 'w') as ok_file:
            pass



def check_bed_overlap(encode_blacklist, chromosome_ranges, chr, start, end, Library_type):
    """
    Function to check if a BED record (chr, start, end) is:
    1. Within the chromosome range (and the chromosome is listed).
    2. Within the ENCODE blacklist.

    :param encode_blacklist: PyRanges object for the ENCODE blacklist.
    :param chromosome_ranges: PyRanges object for the chromosome ranges.
    :param chr: Chromosome of the BED record.
    :param start: Start position of the BED record.
    :param end: End position of the BED record.
    :return: String -> Expected error if region not in assembly
             None -> Everything okay
    """
    returned = None
    # Check if the chromosome exists in the chromosome ranges
    chromosomes_in_bed = chromosome_ranges.chromosomes
    if chr not in chromosomes_in_bed:
        returned= f'custom: {chr}:{start}-{end} ({Library_type}) \t There is no chromosome {chr} in the assembly'

    # Create a PyRanges object for the BED record to check
    bed_record = pr.PyRanges(chromosomes=[chr], starts=[start], ends=[end])

    # Check if the record is within the chromosome range
    within_chrom = chromosome_ranges.intersect(bed_record, how= 'containment').empty

    # Check if the record is within the ENCODE blacklist
    within_blacklist = encode_blacklist.intersect(bed_record,how='containment').empty
    if not within_chrom :
        returned = f'custom: {chr}:{start}-{end} ({Library_type}) \t Region not within assembly'
    elif not within_blacklist:
        returned = f'custom: {chr}:{start}-{end} ({Library_type}) \t Region within poorly defined regions'
    return returned



def transform_func(x):
    # adds some text to the end of transcript IDs
    if 'transcript_id' in x.attributes:
        x.attributes['transcript_id'][0] += '_transcript'
    return x

def Feature_Annotation(Feature, listed):
    returned = []
    if isinstance(Feature, str):
        h = list(
            set([y for x in listed for y in db.children(x, featuretype=Feature)]))
    FeatureType = Feature
    for rec in h:
        Start_edit_freature = rec.start-args.border
        End_edit_freature = rec.end+args.border
        returned.append([rec.seqid, Start_edit_freature,
                         End_edit_freature])
    return returned

def Feature_Annotation_protist(Feature, protein):
    returned = []
    h = list(set([f for f in db.children(db[protein].attributes['ID'][0], featuretype=Feature)]))
    FeatureType = Feature
    for rec in h:
        Start_edit_freature = rec.start-args.border
        End_edit_freature = rec.end+args.border
        returned.append([rec.seqid, Start_edit_freature,
                         End_edit_freature])
    return returned

def fetch_bed(file,encode_blacklist,chromosome_ranges, db, Library_type):
    prot = open(file, 'r')
    Lines = prot.readlines()
    gen_errors=[]
    fetch_error=[]
    returned = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Gene'])
    if os.stat(file).st_size == 0:
        return gen_errors, fetch_error, returned
    try :
        lines = pd.read_csv(file, sep=" ", header=None)
        del lines
    except : 
        gen_errors=errors+[f'file input for {Library_type} is of the wrong format']
        return gen_errors, fetch_error, returned
    for line in Lines:
        protein = str(line.strip())
        if protein.startswith('custom:'):
            print(protein)
            match = re.match(r'custom:(\w+):(\d+)-(\d+)', protein)
            if not match :
                fetch_error = fetch_error + [f'{protein} \t does not match the custom pattern (custom:chromosome:start-end)']
                continue
            chrom, start, end = match.groups()
            chrom, start, end = chrom, int(start), int(end)
            invalidity = check_bed_overlap(encode_blacklist, chromosome_ranges, chrom, start,end,Library_type)
            if invalidity :
                fetch_error = fetch_error + [invalidity]
                continue
            df={'Chromosome':[chrom],'Start':[start],'End' : [end], 'Gene':[protein]}
            returned=pd.concat([returned,pd.DataFrame(df)], axis=0, ignore_index=True)
        else:
            PAM_occurences = []
            if args.protist:
                if args.isoform == 'MANE':
                    fetch_error = fetch_error + ['protist do not have MANE annotations']
                PAM_occurences.extend(Feature_Annotation_protist(Feature,protein))
                pyranges = pr.PyRanges(pd.DataFrame(PAM_occurences, columns=['Chromosome', 'Start', 'End'])).merge()
                df = pyranges.as_df()
                df['Gene'] = protein
                returned=pd.concat([df,returned])
            elif args.isoform == 'MANE':
                if args.genome != 'hg38':
                    fetch_error = fetch_error + ['Only hg38 assembly possess MANE annotations']
                    break
                try : q=[f.attributes['ID'][0] for f in db.children(db[protein].attributes['ID'][0]) if 'tag' in f.attributes.keys() and 'MANE_Select' in f.attributes['tag']]
                except gffutils.exceptions.FeatureNotFoundError :
                    try :
                        # Trying but not keeping
                        g = [f.attributes['ID'][0] for f in db.children(db[protein].attributes['ID'][0])]
                        continue
                    except gffutils.exceptions.FeatureNotFoundError :
                        fetch_error=fetch_error+[f'{protein} ({Library_type}) was not found in database']
                        continue
                    else :
                        fetch_error=fetch_error+[f'{protein} ({Library_type}) did not have a MANE isoform in database']
                        continue
                if not q :
                        fetch_error=fetch_error+[f'{protein} ({Library_type}) did not have a MANE isoform in database']
                        continue
            elif args.isoform == 'Canonical':
                try : q=[f.attributes['ID'][0] for f in db.children(db[protein].attributes['ID'][0]) if 'tag' in f.attributes.keys() and 'Ensembl_canonical' in f.attributes['tag']]
                except gffutils.exceptions.FeatureNotFoundError :
                    try :
                        #Trying not keeping
                        g = [f.attributes['ID'][0] for f in db.children(db[protein].attributes['ID'][0])]
                        continue
                    except gffutils.exceptions.FeatureNotFoundError :
                        fetch_error=fetch_error+[f'{protein} ({Library_type}) was not found in database']
                        continue
                    else :
                        fetch_error=fetch_error+[f'{protein} ({Library_type}) did not have a Canonical isoform in database']
                        continue
                if not q :
                    fetch_error=fetch_error+[f'{protein} ({Library_type}) did not have a Canonical isoform in database']
                    continue
            else:
                try : q = [f.attributes['ID'][0] for f in db.children(db[protein].attributes['ID'][0])]
                except gffutils.exceptions.FeatureNotFoundError :
                    fetch_error=fetch_error+[f'{protein} ({Library_type}) was not found in database']
                    continue
            if not args.protist:
                for Feature in args.features:
                    PAM_occurences.extend(Feature_Annotation(Feature, q))
                    pyranges = pr.PyRanges(pd.DataFrame(PAM_occurences, columns=[
                                   'Chromosome', 'Start', 'End'])).merge()
                    df = pyranges.as_df()
                    df['Gene'] = protein
                    if returned.empty:
                        returned=df
                    else:
                        returned=pd.concat([df,returned])
    return gen_errors, fetch_error, returned


def merge_overlapping_genes(df):
    """
    Takes a BED-like DataFrame with 'Chromosome', 'Start', 'End', and 'Gene' columns.
    Merges overlapping and continuous regions, concatenates gene names, and removes empty gene segments.
    """
    # Remove empty gene entries
    df = df.dropna(subset=["Gene"])
    df = df[df["Gene"] != ""]
    
    # Sort by chromosome and start position
    df = df.sort_values(by=["Chromosome", "Start"]).reset_index(drop=True)

    merged_intervals = []
    current_start, current_end, current_genes = df.iloc[0]["Start"], df.iloc[0]["End"], set(df.iloc[0]["Gene"].split("-"))

    # Iterate through regions
    for i in range(1, len(df)):
        row = df.iloc[i]
        row_genes = set(row["Gene"].split("-"))

        # If the current region overlaps or is continuous with the previous one, merge it
        if row["Start"] <= current_end:
            current_end = max(current_end, row["End"])
            current_genes.update(row_genes)
        else:
            # Append the last merged region
            merged_intervals.append([df.iloc[0]["Chromosome"], current_start, current_end, ":".join(sorted(current_genes))])
            # Start a new merged region
            current_start, current_end, current_genes = row["Start"], row["End"], row_genes

    # Append the last merged region
    merged_intervals.append([df.iloc[0]["Chromosome"], current_start, current_end, ":".join(sorted(current_genes))])

    # Convert back to DataFrame
    merged_df = pd.DataFrame(merged_intervals, columns=["Chromosome", "Start", "End", "Gene"])

    return merged_df

def find_duplicates(file1, file2, file3):
    files = [file1, file2, file3]
    seen = set()
    duplicates = set()
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line in seen:
                    duplicates.add(line)
                else:
                    seen.add(line)
    return duplicates

if __name__ == '__main__':
    args = argparser.parse_args()
    validate_editor_file(args.editors)
    Genome_file=Path(args.Genome).stem
    chrom_df = pd.read_csv(args.chromosome_bed, sep='\t', header=None, names=['Chromosome', 'End'])
    chrom_df['Start'] = 0
    chromosome_ranges = pr.PyRanges(chrom_df[['Chromosome', 'Start', 'End']])
    if args.encode_bed :
        encode_blacklist=pr.read_bed(args.encode_bed)
    else :
        empty_df = pd.DataFrame(columns=["Chromosome", "Start", "End"])
        encode_blacklist = pr.PyRanges(empty_df)
    errors=[]
    fetch_errors=[]
    if find_duplicates(args.target_file, args.positive_file, args.negative_file):
        fetch_errors=fetch_errors+[f'Duplicate gene entries were found']
    if not os.path.exists(args.Genome):
        raise ValueError(f'{args.Genome} does not exist')
    try :
        db = gffutils.FeatureDB(args.Genome)
    except:
        try :
            warnings.warn(f"Creating the database {Genome_file}.db this may take a few minutes. \n In the future this database can be used in the -G argument directly Avoiding re-computation")
            db = gffutils.create_db(args.Genome, Genome_file + '.db', id_spec={'gene': 'gene_name', 'transcript': "transcript_id"}, merge_strategy="create_unique", transform=transform_func, keep_order=True)
        except:
            raise ValueError(f'Genome freature database {args.Genome} is not in the correct format')
    target_error, fetch_target_error, target_df = fetch_bed(args.target_file ,encode_blacklist,chromosome_ranges, db, 'Target Library')
    positive_error , fetch_positive_error, positive_df = fetch_bed(args.positive_file,encode_blacklist,chromosome_ranges, db, 'Positive Library')
    negative_error , fetch_negative_error, negative_df = fetch_bed(args.negative_file,encode_blacklist,chromosome_ranges, db, 'Negative Library')
    errors = errors+ target_error + positive_error + negative_error
    fetch_errors=fetch_errors+fetch_target_error+fetch_negative_error+fetch_positive_error
    if fetch_errors :
        with open('fetch.err', 'w') as file:
            file.write('\n'.join(fetch_errors))
    elif errors :
        with open('errors.err', 'w') as file:
            file.write('\n'.join(errors))
    else :
        target_pr = pr.PyRanges(target_df[['Chromosome', 'Start', 'End']])
        positive_pr = pr.PyRanges(positive_df[['Chromosome', 'Start', 'End']])
        negative_pr = pr.PyRanges(negative_df[['Chromosome', 'Start', 'End']])
        length = target_pr.lengths().sum() + negative_pr.lengths().sum() + positive_pr.lengths().sum()
        if length > args.limit:
            print(f'::error:: libraries queries regions total is too high ({length} > {args.limit}')
            sys.exit("Error exit")
        if target_pr.join(positive_pr) :
            print(f'::error:: The Target library and Positive library overlaps {target_pr.join(positive_pr)}')
            sys.exit("Error exit")
        if target_pr.join(negative_pr) :
            print(f'::error:: The Target library and Positive library overlaps {target_pr.join(negative_pr)}')
            sys.exit("Error exit")
        if negative_pr.join(positive_pr) :
            print(f'::error:: The Target library and Positive library overlaps {negative_pr.join(positive_pr)}')
            sys.exit("Error exit")
        target_df.to_csv('Study_Target_library.bed',sep='\t', index=False, header=False)
        if not positive_df.empty:
            positive_df.to_csv('Positive_Controls_library.bed',sep='\t', index=False, header=False)
        if not negative_df.empty:
            negative_df.to_csv('Negative_Controls_library.bed',sep='\t', index=False, header=False)
