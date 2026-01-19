import argparse
import os.path

import json, time, os, sys, glob

import gzip
import numpy as np
import itertools
import random
from collections import defaultdict

import sys

def parse_args( argv ):
    argv_tmp = sys.argv
    sys.argv = argv
    description = 'do protein sequence design using the MPNN model ...'
    parser = argparse.ArgumentParser( description = description )
    parser.add_argument('-pdbs', type=str, nargs='*', help='name of the input pdb file')
    parser.add_argument('-pdb_list', type=str, help='a list file of all pdb files')
    args = parser.parse_args()
    sys.argv = argv_tmp

    return args

args = parse_args( sys.argv )

if args.pdbs == None:
    assert (args.pdb_list != None)
    with open(args.pdb_list) as f:
        all_pdbs = [line.strip() for line in f]
else:
    all_pdbs = args.pdbs

for pdb in all_pdbs:

    assert(pdb.endswith('pdb.gz'))

    with gzip.open(pdb, 'rt') as f:
        pdblines = f.readlines()

    new_lines = []
    for line in pdblines:
        if 'BLA' in line:
            new_line = line.replace('BLA', 'FDI')
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    with gzip.open(pdb, 'wt') as f:
        for line in new_lines:
            f.write(line)
