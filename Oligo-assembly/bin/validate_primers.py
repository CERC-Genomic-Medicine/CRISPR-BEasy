#!/usr/bin/env python3

'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.1
YEAR: 2024

Expension possible: Testing annealing temperature, testing selfing etc.
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
import warnings

argparser = argparse.ArgumentParser(
    description='This script validates the primers.')
argparser.add_argument('--Primer_Customs', metavar='string', dest='Primer_custom',
                       type=str, required=True, help='Custom primer forw,rev')
argparser.add_argument('-s', action='store_true', help='Sensor assembly flag')

def is_invalid_sequence(seq):
    return re.fullmatch(r'[atcgATCG]+', seq) is None

def validate_primers(Primer_custom):
    errors = []
    messages = []
    if (Primer_custom == ",") :
        errors.append(f"No primers were selected")
    else :
        list_primer = Primer_custom.split(",")
        if len(list_primer)!=2:
            errors.append(f"Primers are limited to nuclotides (i.e. ACGT)")
        else :
            if (is_invalid_sequence(list_primer[0]) or is_invalid_sequence(list_primer[1])):
                errors.append(f"Primers are limited to nuclotides (i.e. ACGT)")
    if errors:
        with open('primer.err', 'w') as error_file:
            error_file.write('\n'.join(errors))
    else :
        with open('primer.ok', 'w') as ok_file:
            ok_file.write('Primers were accepted')


if __name__ == '__main__':
    args = argparser.parse_args()
    validate_primers(args.Primer_custom)
