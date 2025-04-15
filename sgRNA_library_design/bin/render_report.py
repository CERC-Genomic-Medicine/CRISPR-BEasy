#!/usr/bin/python3

import argparse
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--template', required=True, help='Path to Jinja2 template file')
    parser.add_argument('--genome', required=True)
    parser.add_argument('--genome_vep', required=True)
    parser.add_argument('--species_vep', required=True)
    parser.add_argument('--release', required=True)
    parser.add_argument('--baseurl', required=True)
    parser.add_argument('--collection', required=False)
    parser.add_argument('--output', required=True, help='Path to save the rendered output')
    args = parser.parse_args()

    template_path = Path(args.template).resolve()
    env = Environment(loader=FileSystemLoader(template_path.parent))
    template = env.get_template(template_path.name)
    if args.collection :
        rendered = template.render(
                genome=args.genome,
                genome_vep=args.genome_vep,
                species_vep=args.species_vep,
                release=args.release,
                baseurl=args.baseurl,
                collection=args.collection
         )
    else :
        rendered = template.render(
                genome=args.genome,
                genome_vep=args.genome_vep,
                species_vep=args.species_vep,
                release=args.release,
                baseurl=args.baseurl,
         )


    with open(args.output, 'w') as f:
        f.write(rendered)

if __name__ == '__main__':
    main()

