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
import warnings

argparser = argparse.ArgumentParser(
    description='This software produces a bed file corresponding to the regions of interest as defined by a protein (corresponding to a genome) and the desired Feature. This represent the first step (potentially optional) in producing a library design crispr Array')
argparser.add_argument('-i', '--instructions', metavar='file', dest='instructions',
                       type=str, required=True, help='set of instruction to validate')
argparser.add_argument('-x', '--xlsx', metavar='file', dest='xlsx',
                       type=str, required=True, help='xlsx corresponding to instruction')


def is_invalid_sequence(seq):
    return re.fullmatch(r'[atcgATCG]+', seq) is None

def validate_instruction_file(file_path, editors_correct):
    """
    Validates that the editor file is in the correct format:
    - 3 space-delimited columns
    - 2nd and 3rd columns are whole numbers (can be negative)

    If validation fails, errors are logged to 'editors.err'. If valid, an empty 'editor.ok' file is created.
    If the file is blank, an empty 'editor.blank' file is created.

    :param file_path: Path to the editor file to validate
    """
    errors = []
    messages = []

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
            if len(columns) ==0:
                continue
            if len(columns) != 3:
                errors.append(f"Line {line_num}: Expected 3 columns, found {len(columns)}")
                continue
            if columns[0] not in editors_correct :
                error.append(f"Line {line_num}: editor {columns[0]} is not found in the list of editors found in positive library's sheet")
            try:
                int(columns[1])
            except:
                errors.append(f"Line {line_num}: Columns 2 must be whole numbers, found '{columns[1]}'")
                continue
            if int(columns[1])<0 :
                errors.append(f"Line {line_num}: Columns 2 cannot be negative numbers, found '{columns[1]}'")
                continue
            elif int(columns[1])==0:
                messages.append(f"Line {line_num}: Columns 2 was 0, All sgRNA meeting criteria will be outputed.")
            valid_consequences = {'none','stop_gained', 'splice_altering', 'start_loss'}
            consequences=columns[2].split(",")
            consequences=[o.lower() for o in consequences]
            if len(consequences)>1 and 'none' in consequences :
                errors.append(f"Line {line_num}: Columns 3, Consequences cannot be combined with 'None'")
            elif not set(consequences).issubset(valid_consequences):
                errors.append(f"Line {line_num}: Columns 3, invalid consequences specified")
    if errors:
        with open('Instructions.err', 'w') as error_file:
            error_file.write('\n'.join(errors))
    else:
        # If no errors, create an empty editor.ok file
        with open('Instructions.ok', 'w') as ok_file:
            ok_file.write('\n'.join(messages))
            pass




if __name__ == '__main__':
    args = argparser.parse_args()
    editor_sheets = []
    xlsx = pd.ExcelFile(args.xlsx)
    sheet_names = xlsx.sheet_names
    editor_sheets.extend([name for name in sheet_names if name.startswith("editor - ")])
    editors = [s.replace("editor - ", "") for s in editor_sheets]
    validate_instruction_file(args.instructions, editors)
