#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

import os
import sys
import subprocess


def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    output, error = process.communicate()
    return output, error


COLS = ['CHR', 'SNP', 'BP', 'NMISS', 'BETA', 'SE', 'T', 'P']
MAF_FILTER=0.05
FILE = sys.argv[1]


command = f"cat {FILE} | cut -f 1,3,6,9,10,11,12 -d $'\\t' > {os.path.splitext(FILE)[0]}_.tsv"

print(command)
run(command)
print('Shortened')

data = pd.read_csv(f"{os.path.splitext(FILE)[0]}_.tsv", sep='\t')
print('Readed')

data.minor_AF = data.minor_AF.apply(lambda x: min(x, 1-x))
data = data[data.minor_AF>=MAF_FILTER]
print('Filtered MAF')

data['CHR'] = data['variant'].apply(lambda x: x.split(':')[0])
data['BP'] = data['variant'].apply(lambda x: x.split(':')[1])
data['SNP'] = data['CHR']+':'+data['BP']
data = data.rename(columns={'n_complete_samples':'NMISS', 'beta':'BETA', 'se':'SE', 'tstat':'T', 'pval':'P'})

print('Saving')
data[COLS].to_csv(os.path.splitext(FILE)[0]+'.qassoc', sep='\t', index=False)
print('Done')

command = f"rm {os.path.splitext(FILE)[0]}_.tsv"
print(command)
run(command)
print('Deleted extra file')