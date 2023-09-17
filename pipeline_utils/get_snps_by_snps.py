#!/usr/bin/env python

import subprocess
import pandas as pd
import re
import os
import sys


def extract_substrings(input_string):
    matches = re.findall(r'\"\"(.+?)\"\"', input_string)
    return '_'.join(matches)

def flatten(l):
    return [item for sublist in l for item in sublist]

# 

# params
## todo possibility to chose snps, not genes
# ========================
# data dir
data_dir = sys.argv[1]
# sim pattern
pattern = sys.argv[2]
# pheno pattern
pheno_pattern = sys.argv[3]
# simulated bed file
vcf_file = sys.argv[4]
# ------------------------
# pathway 2 genes file
snps_file = sys.argv[5]
K = int(float(sys.argv[6]))
# ========================

# annotated snps from bed
annotated_snps_file = os.path.join(data_dir, f"{pattern}_{pheno_pattern}_annotated_snps.tsv")
# chosen snps
selected_snps_file = os.path.join(data_dir, f"{pattern}_{pheno_pattern}_chosen_snps.tsv")


# chose this specific genes in file
bash_command = f"value=$(tail -n 1 {vcf_file} | awk '{{print NF}}') ; col=$(echo $value-3 | bc) ; grep -w -f {snps_file} {vcf_file} | cut -f 1,2,4,$col > {annotated_snps_file}"
process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
output, error = process.communicate()

snps = pd.read_csv(annotated_snps_file, sep='\t', names=['chr', 's', 'e', 'gene'], header=None)
snps.gene = snps.gene.apply(lambda x: extract_substrings(x))
snps = snps.groupby('s').agg({'chr':'first', 'e' : 'first', 'gene' : 'first'}).reset_index()[['chr', 's', 'e', 'gene']].sort_values(by=['chr', 's'])

snps = snps.sample(min(snps.shape[0], K), replace=False).reset_index(drop=True)

snps.e = snps.s+snps.e.apply(lambda x: len(x)-1)
print(snps.shape)
snps = snps.sort_values(['chr', 's'])
snps.to_csv(selected_snps_file, header=False, index=False, sep='\t')
print(f'Succesfully saved to {selected_snps_file}')
