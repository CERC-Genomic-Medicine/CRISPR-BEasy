#!/usr/bin/python3

import argparse
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--template', required=True, help='Path to Jinja2 template file')
    parser.add_argument('--genome', required=True)
    parser.add_argument('--gff3_url', required=True)
    parser.add_argument('--fasta_url', required=True)
    parser.add_argument('--fasta_index_url', required=True)
    parser.add_argument('--fasta_gz_gzi', required=True)
    parser.add_argument('--output', required=True, help='Path to save the rendered output')
    args = parser.parse_args()

    template_path = Path(args.template).resolve()
    env = Environment(loader=FileSystemLoader(template_path.parent))
    template = env.get_template(template_path.name)
    rendered = template.render(
                genome = args.genome,
                fasta_url = args.fasta_url,
                gff3_url = args.gff3_url,
                fasta_index_url = args.fasta_index_url,
                fasta_gz_gzi = args.fasta_gz_gzi
     )


    with open(args.output, 'w') as f:
        f.write(rendered)

if __name__ == '__main__':
    main()

