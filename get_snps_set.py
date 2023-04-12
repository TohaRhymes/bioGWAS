#!/usr/bin/env python

import subprocess
import pandas as pd
import re
import numpy as np
from collections import defaultdict
import random


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

#params
gmt_file = './h.all.v2023.1.Hs.symbols.gmt'
genes_file = './causal_geneset.txt'
included_gene_file = './unique_chr_EUR_sim.gtf'
annotated_snps_file="chr_EUR_sim.sorted.annotated_chosen.tsv"
selected_snps_file="snps_chr_EUR.tsv"
pathways = ['HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_HYPOXIA']
K = 20
k = 10




bash_command = "value=$(tail -n 1 chr_EUR_sim.sorted.annotated.bed | awk '{print NF}') ; col=$(echo $value-4 | bc) ; cut -f $col chr_EUR_sim.sorted.annotated.bed | sort | uniq > unique_chr_EUR_sim.gtf"

process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


included_genes = pd.read_csv(included_gene_file, sep='\t', names=['gene'], header=None)
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
path2genes = dict(path_genes)

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

with open(genes_file, 'w') as f:
    for x in gene_set[:-1]:
        f.write(f"{x}\n")
    f.write(f"{gene_set[-1]}")

    
    



bash_command = "value=$(tail -n 1 chr_EUR_sim.sorted.annotated.bed | awk '{print NF}') ; col=$(echo $value-4 | bc) ; grep -w -f <(cut -f 1 causal_geneset.txt) chr_EUR_sim.sorted.annotated.bed | cut -f 1,2,3,$col > chr_EUR_sim.sorted.annotated_chosen.tsv"

process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
output, error = process.communicate()



snps = pd.read_csv(annotated_snps_file, sep='\t', names=['chr', 's', 'e', 'gene'], header=None)
snps.gene = snps.gene.apply(lambda x: extract_substrings(x))
snps = snps.groupby('gene').apply(lambda x: choose_random_row(x)).reset_index(drop=True)
snps.s = snps.e
snps.to_csv(selected_snps_file, header=False, index=False, sep='\t')
