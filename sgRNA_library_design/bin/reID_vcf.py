#!/usr/bin/env python3

import argparse
import gzip
import csv

def load_id_map(dict_path):
    id_map = {}
    with open(dict_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = f"{row['Protein']}_{row['start']}_{row['strand']}"
            id_map[key] = row['ID']
    return id_map

def process_vcf_line(line, id_map):
    fields = line.strip().split('\t')
    if len(fields) < 3:
        return line  # Skip malformed lines
    old_id = fields[2]
    new_id = id_map.get(old_id, old_id)
    fields[2] = new_id
    return '\t'.join(fields)

def main():
    parser = argparse.ArgumentParser(description="Replace VCF ID (3rd column) using mapping from dictionary.")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--dict", required=True, help="Dictionary TSV with columns: ID, Protein, start, strand")
    parser.add_argument("--output", required=True, help="Output .vcf.gz file")
    args = parser.parse_args()

    id_map = load_id_map(args.dict)

    with gzip.open(args.vcf, 'rt') as vcf_in, gzip.open(args.output, 'wt') as vcf_out:
        for line in vcf_in:
            if line.startswith('#'):
                vcf_out.write(line)
            else:
                processed_line = process_vcf_line(line, id_map)
                vcf_out.write(processed_line + '\n')

if __name__ == "__main__":
    main()

