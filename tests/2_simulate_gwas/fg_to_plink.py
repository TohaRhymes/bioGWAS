#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

import os
import sys
import subprocess

import math


CHR_ORDER = {'1': 0,
 '2': 1,
 '3': 2,
 '4': 3,
 '5': 4,
 '6': 5,
 '7': 6,
 '8': 7,
 '9': 8,
 '10': 9,
 '11': 10,
 '12': 11,
 '13': 12,
 '14': 13,
 '15': 14,
 '16': 15,
 '17': 16,
 '18': 17,
 '19': 18,
 '20': 19,
 '21': 20,
 '22': 21,
 'X': 22,
 'Y': 23}

N_COLS = ['n_hom_cases',
 'n_hom_ref_cases',
 'n_het_cases',
 'n_hom_controls',
 'n_hom_ref_controls',
 'n_het_controls']


def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    output, error = process.communicate()
    return output, error


COLS = ['CHR', 'SNP', 'BP', 'NMISS', 'BETA', 'SE', 'T', 'P']
MAF_FILTER=0.05
FILE = sys.argv[1]


command = f"cat {FILE} | cut -f 7,9,10,11,14,15,16,17,18,19,20,21 -d $'\\t' > {os.path.splitext(FILE)[0]}_.tsv"

print(command)
run(command)
print('Shortened')

data = pd.read_csv(f"{os.path.splitext(FILE)[0]}_.tsv", sep='\t')
print('Readed')

data.af_alt = data.af_alt.apply(lambda x: min(x, 1-x))
data = data[data.af_alt>=MAF_FILTER]
print('Filtered MAF')

data['CHR'] = data['new_chr'].astype(str).apply(lambda x: x.replace('chr', ''))
data = data[data['CHR'].isin(list(CHR_ORDER.keys()))]
data['BP'] = data['new_coord'].astype(int)
data['chr_order'] = data.CHR.apply(lambda x: CHR_ORDER[x])
data = data.sort_values(by=['chr_order', 'BP'])

data['NMISS'] = data[N_COLS].sum(axis=1).apply(lambda x: math.floor(x))

data['SNP'] = data['CHR']+':'+data['BP'].astype(str)

data['T'] = data.beta/data.sebeta

data = data.rename(columns={'beta':'BETA', 'sebeta':'SE', 'pval':'P'})

print('Saving')
data[COLS].to_csv(os.path.splitext(FILE)[0]+'.qassoc', sep='\t', index=False)
print('Done')

command = f"rm {os.path.splitext(FILE)[0]}_.tsv"

print(command)
run(command)

print('Deleted extra file')