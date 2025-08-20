#!/usr/bin/env python3

import gffutils
import os
import argparse
import re
import gzip
import json
from collections import Counter

argparser = argparse.ArgumentParser(
    description='install gff3')

argparser.add_argument('-G', '--Genome_database', metavar='file', dest='Genome', type=str, required=True,
                       help='Genome feature database file')
argparser.add_argument('-O', metavar='string', dest='output', type=str, required=True,
                       help='output')

def transform_func(x):
    # adds some text to the end of transcript IDs
    if 'Name' in x.attributes and 'gene_id' in x.attributes and x.attributes['Name'][0] in duplicate_gene:
        x.attributes['Unique_name'] = x.attributes['gene_id']
    if 'Parent' in x.attributes and 'transcript_id' in x.attributes:
        if ',' in x.attributes['Parent'][0] :
            parents= x.attributes['Parent'][0].split(',')
            x.attributes['Parent'][0]=",".join([dicti[parenti] for parenti in parents])
        else :
            x.attributes['Parent'][0] = dicti[x.attributes['Parent'][0]]
    return x



def test_mane_select(db):
    found_mane = False
    for feature in db.all_features():
        tags = feature.attributes.get("tag", [])
        if "MANE_Select" in tags:
            found_mane = True
            print(f"Found MANE_Select on {feature.featuretype} (ID={feature.id}).")
            break
    return found_mane


if __name__ == '__main__':
    args = argparser.parse_args()
    genomes=os.path.splitext(args.Genome)[0].split('.')
    genome=genomes[0]
    dicti ={}
    list_gene=[]
    db = gffutils.create_db(args.Genome, ":memory:",
    disable_infer_transcripts=True,
    disable_infer_genes=True,
    id_spec={'gene': ["Name","gene_id"], 'chromosome':'ID', 'mRNA':'ID', 'ncRNA':'ID', 'pseudogenic_transcript':'ID', 'rRNA':'ID', 'snRNA':'ID', 'snoRNA':'ID', 'tRNA':'ID', 'transposable_element':'ID'}, merge_strategy='create_unique')
    for x in db.all_features():
        if 'Name' in x.attributes and 'gene_id' in x.attributes:
            list_gene.append(x.attributes['Name'][0])
    duplicate_gene = {x for x, c in Counter(list_gene).items() if c > 1}
    for feature in db.all_features():
        if 'Name' in feature.attributes and feature.attributes['Name'][0] in duplicate_gene and "gene_id" in feature.attributes  :
            ided=feature.attributes['gene_id'][0]
        else :
            ided=feature.id
        if 'ID' in list(feature.attributes) and "gene_id" in feature.attributes:
            dicti[feature['ID'][0]]=ided
    del db 
    db1 = gffutils.create_db(args.Genome, args.output + '.db',
    disable_infer_transcripts=True,
    disable_infer_genes=True,
    id_spec={'gene': ["Unique_name","Name","gene_id"], 'chromosome':'ID', 'mRNA':'ID', 'ncRNA':'ID', 'pseudogenic_transcript':'ID', 'rRNA':'ID', 'snRNA':'ID', 'snoRNA':'ID', 'tRNA':'ID', 'transposable_element':'ID'}, merge_strategy='error', transform=transform_func)
    #Test
    total=[i for i in db1.features_of_type("gene")]
    firsts=total[:10]
    firsts_id=[i.id for i in firsts]
    overall_mane=test_mane_select(db1)
    output_data = {
    "mane_select": overall_mane,
    "first_genes": firsts_id
    }
    with open(args.output+"_gff.json", "w") as out_f:
        json.dump(output_data, out_f, indent=2)
