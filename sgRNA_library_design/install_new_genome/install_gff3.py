import gffutils
import os
import argparse
import re
import gzip

argparser = argparse.ArgumentParser(
    description='install gff3')

argparser.add_argument('-G', '--Genome_database', metavar='file', dest='Genome', type=str, required=True,
                       help='Genome feature database file')

def transform_func(x):
    # adds some text to the end of transcript IDs
    if 'transcript_id' in x.attributes:
        x.attributes['transcript_id'][0] += '_transcript'
    return x

def check_gene_column(file_path):
    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as file:  # 'rt' mode ensures text reading
        for line in file:
            columns = line.strip().split()
            if len(columns) >= 9:
                match_Name = re.search(r'\bName=', columns[8])
                match_gene = re.search(r'\bgene=', columns[8])
                match_gene_name = re.search(r'\bgene_name=', columns[8])
                match_gene_id = re.search(r'\bgene_id=', columns[8])
                if match_Name:
                    return 'Name'
                if match_gene:
                    return 'gene'
                elif match_gene_name:
                    return 'gene_name'
                elif match_gene_id:
                    return 'gene_id'
            continue
    return None  # Return None if neither pattern is found in any line

if __name__ == '__main__':
    args = argparser.parse_args()
    genomes=os.path.splitext(args.Genome)[0].split('.')
    genome=genomes[0]
    gene_desc=check_gene_column(args.Genome)
    print(gene_desc)
    db = gffutils.create_db(args.Genome, genome + '.db', id_spec={'gene': gene_desc, 'transcript': "transcript_id"}, merge_strategy="create_unique", transform=transform_func, keep_order=True, force=True)
    #Test
    total=[i for i in db.all_features()]
    firsts=total[:20]
    firsts_id=[i.id for i in firsts]
    print(firsts_id)

