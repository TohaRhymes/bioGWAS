#!/usr/bin/env python

import subprocess
import pandas as pd
import re
import numpy as np
from collections import defaultdict
import random
import os
import sys


def extract_substrings(input_string):
    matches = re.findall(r'\"\"(.+?)\"\"', input_string)
    return '_'.join(matches)


def extract_second(input_string):
    return input_string.split('_')[1]

def choose_random_row(group_df):
    random_idx = np.random.choice(group_df.index)
    return group_df.loc[random_idx]


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
# simulated bed file
bed_file = sys.argv[3]
# ------------------------
# pathway 2 genes file
gmt_file = sys.argv[4]
# file with selected pathways
pathways_file = sys.argv[5]
K = int(sys.argv[6])
k = int(sys.argv[7])
# ========================

#included in bed file genes
included_genes_file = os.path.join(data_dir, f'{pattern}_included_genes_snps.gtf')
# where to save gene causal geneset
causal_genes_file = os.path.join(data_dir, f'{pattern}_causal_geneset_snps.txt')
# annotated snps from bed
annotated_snps_file = os.path.join(data_dir, f"{pattern}_annotated_snps.tsv")
# chosen snps
selected_snps_file = os.path.join(data_dir, f"{pattern}_chosen_snps.tsv")


with open(pathways_file) as f:
    pathways = flatten([line.strip().split(',') for line in f])



bash_command = f"value=$(tail -n 1 {bed_file} | awk '{{print NF}}') ; col=$(echo $value-4 | bc) ; cut -f $col {bed_file} | sort | uniq > {included_genes_file}"
process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
output, error = process.communicate()


included_genes = pd.read_csv(included_genes_file, sep='\t', names=['gene'], header=None)
included_genes.gene = included_genes.gene.apply(lambda x: extract_second(extract_substrings(x)))
included_genes = set(included_genes.gene.unique())


# Read file with pathway-<genes> per line
lines = []
with open(gmt_file, 'r') as gmt:
    lines = list(map(lambda x: x.split('\t'), gmt.readlines()))

# Make dict key = pathway, value = set of genes 
path2genes = defaultdict(str)
for l in lines:
    path2genes[l[0]] = set(map(lambda x: x.replace('\n', ''), l[2:]))
path2genes = dict(path2genes)

# For target pathways make set of target genes and set of other genes 
target_genes = set(flatten([path2genes[target] for target in pathways])) & included_genes
other_genes = set(flatten(path2genes.values())) & included_genes
other_genes = list(other_genes-target_genes)
target_genes = list(target_genes)

target_genes = random.choices(target_genes, k=k)
other_genes = random.choices(other_genes, k=K-k)
print(f'Chosen target genes: {target_genes}')
print(f'Chosen other genes: {other_genes}')
gene_set = target_genes+other_genes

with open(causal_genes_file, 'w') as f:
    for x in gene_set[:-1]:
        f.write(f"{x}\n")
    f.write(f"{gene_set[-1]}")

    
    
# chose this specific genes in file
bash_command = f"value=$(tail -n 1 {bed_file} | awk '{{print NF}}') ; col=$(echo $value-4 | bc) ; grep -w -f <(cut -f 1 {causal_genes_file}) {bed_file} | cut -f 1,2,3,$col > {annotated_snps_file}"
process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
output, error = process.communicate()



snps = pd.read_csv(annotated_snps_file, sep='\t', names=['chr', 's', 'e', 'gene'], header=None)
snps.gene = snps.gene.apply(lambda x: extract_substrings(x))
snps = snps.groupby('gene').apply(lambda x: choose_random_row(x)).reset_index(drop=True)
snps.s = snps.e
snps.to_csv(selected_snps_file, header=False, index=False, sep='\t')
print(f'Succesfully saved to {selected_snps_file}')
