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
    try:
        return input_string.split('_')[1]
    except IndexError:
        return pd.NA

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
# pheno pattern
pheno_pattern = sys.argv[3]
# simulated bed file
vcf_file = sys.argv[4]
# ------------------------
# pathway 2 genes file
gmt_file = sys.argv[5]
# file with selected pathways
pathways_file = sys.argv[6]
K = int(float(sys.argv[7]))
k = int(float(sys.argv[8]))
# ========================

#included in bed file genes
included_genes_file = os.path.join(data_dir, f'{pattern}_{pheno_pattern}_included_genes_snps.txt')
# where to save gene causal geneset
causal_genes_file = os.path.join(data_dir, f'{pattern}_{pheno_pattern}_causal_geneset_snps.txt')
# annotated snps from bed
annotated_snps_file = os.path.join(data_dir, f"{pattern}_{pheno_pattern}_annotated_snps.tsv")
# chosen snps
selected_snps_file = os.path.join(data_dir, f"{pattern}_{pheno_pattern}_chosen_snps.tsv")


with open(pathways_file) as f:
    pathways = flatten([line.strip().split(',') for line in f])


included_genes = pd.read_csv(included_genes_file, sep='\t', names=['gene'], header=None)
included_genes.gene = included_genes.gene.apply(lambda x: extract_second(extract_substrings(x)))
included_genes = included_genes[~included_genes.gene.isna()]
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

target_genes = random.sample(target_genes, k=k)
other_genes = random.sample(other_genes, k=K-k)
gene_set = target_genes+other_genes
print('Genes were chose!')

while len(gene_set) != K:
    target_genes = random.sample(target_genes, k=k)
    other_genes = random.sample(other_genes, k=K-k)
    gene_set = target_genes+other_genes

    
print(f'Chosen target genes: {target_genes}')
print(f'Chosen other genes: {other_genes}')

with open(causal_genes_file, 'w') as f:
    for x in gene_set[:-1]:
        f.write(f"\"{x}\"\n")
    f.write(f"{gene_set[-1]}")

    
print("Chosing these specific genes in file...")
# chose this specific genes in file
bash_command = f"value=$(tail -n 1 {vcf_file} | awk '{{print NF}}') ; col=$(echo $value-3 | bc) ; grep -w -f <(cut -f 1 {causal_genes_file}) {vcf_file} | cut -f 1,2,4,$col > {annotated_snps_file}"
print(bash_command)
process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
output, error = process.communicate()

print("Finishing")
snps = pd.read_csv(annotated_snps_file, sep='\t', names=['chr', 's', 'e', 'gene'], header=None)
snps.gene = snps.gene.apply(lambda x: extract_substrings(x))
snps = snps.groupby('s').agg({'chr':'first', 'e' : 'first', 'gene' : 'first'}).reset_index()[['chr', 's', 'e', 'gene']].sort_values(by=['chr', 's'])

snps = snps.groupby('gene').apply(lambda x: choose_random_row(x)).reset_index(drop=True)
snps.e = snps.s+snps.e.apply(lambda x: len(x)-1)
snps = snps.sort_values(['chr', 's'])
snps.to_csv(selected_snps_file, header=False, index=False, sep='\t')
print(f'Succesfully saved to {selected_snps_file}')
