#!/usr/bin/env python
# coding: utf-8

from itertools import product
from pprint import pprint
import pandas as pd
import os
import subprocess
import re
from collections import defaultdict

from tqdm import tqdm

import glob



from utils_2 import VALIDATION_COMPARISON, IDs, CAUSAL_SNP_FILE, GWAS_FILE, FILE_RESULTS, BFILE

def extract_digit_after_seed(input_string, seed='seed'):
    pattern = f"{re.escape(seed)}(\d+)_gwas"
    match = re.search(pattern, input_string)
    
    if match:
        digit = match.group(1)
        return digit
    else:
        return None


def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    output, error = process.communicate()
    return output, error

def flatten(l):
    return [item for sublist in l for item in sublist]

def get_snps(line):
    if line=='NONE':
        return []
    else:
        return list(map(lambda x: re.sub('\([\d]+\)', '', x), line.split(',')))


print(f"Summary will be saved at: {FILE_RESULTS}")

snps = pd.read_csv(VALIDATION_COMPARISON, sep='\t').drop(['clumps_not_causal', 'causal_not_found'], axis=1)
snps['F1']=2/(1/snps.clumps_causal+1/snps.causal_found)
snps = snps[snps.F1 > 0.9]
params = snps.T.to_dict()

print(f"Amount of parameters sets: {len(params)}")


# ## Make table of results

result_set = []

## SNPs set for every k
for param_key in tqdm(params):
    param = params[param_key]
    K = int(float(param["K"]))
    k = K // 2
    
    m_beta = param["m_beta"]
    sd_beta = param["sd_beta"]
    gen_var = param["gen_var"]
    alpha = param["alpha"]
    theta = param["theta"]
    pIndep = param["pIndep"]
    
    param['seed']="*"
    causal_snp_file = CAUSAL_SNP_FILE.format(**{'K':K})
    gwas_file_pattern = GWAS_FILE.format(**{'K':K}, **param)
    gwas_files = {extract_digit_after_seed(f):f for f in glob.glob(gwas_file_pattern)}
    
    # SET OF CAUSAL SNPS
    causal_snp = pd.read_csv(causal_snp_file, sep='\t', header=None, names=['chr', 'pos', 'pos_e', 'gene'])
    causal_set = set(causal_snp.chr.astype(str)+':'+causal_snp.pos.astype(str))
    print(causal_set)
        
    for seed in gwas_files:
        param['seed']=seed
        gwas_file = GWAS_FILE.format(K, **param)
        # count number of snps
        with open(gwas_file, "rb") as f:
            N = sum(1 for _ in f) - 1

        # make clump file
        FILE_IN = gwas_file.replace('.tsv', '')
        FILE_OUT = FILE_IN + '_clump'
        PVAL_CUTOFF = 0.05/N
        KB = 1000

        try:
            clumps_data = pd.read_fwf(f"{FILE_OUT}.clumped", sep='\t')
        except FileNotFoundError:
            command_format = f"""sed '1{{ s/chr/CHR/; s/rsid/SNP/; s/pos/BP/; s/n/NMISS/; s/beta/BETA/; s/se/SE/; s/r2/R2/; s/t/T/; s/pval/P/;}}' {FILE_IN}.tsv | awk '{{$1=$1}};1'| tr -s ' ' '\t' > {FILE_IN}.qassoc"""
            command_clump = f"""plink         --bfile {BFILE}         --allow-no-sex         --clump {FILE_IN}.qassoc         --clump-p1 {'{:.30f}'.format(PVAL_CUTOFF)}         --clump-p2 {'{:.30f}'.format(PVAL_CUTOFF)}         --clump-r2 0.2         --clump-kb {KB}         --out {FILE_OUT}"""
            print(command_format)
            run(command_format)
            print(command_clump)
            run(command_clump)
            print('done clump')

            # read clump file and make list of snps for every clump
            try:
                clumps_data = pd.read_fwf(f"{FILE_OUT}.clumped", sep='\t')
            except FileNotFoundError:
                results = {"K": K, 
                **param,
                'clumps_total': 0,
                'clumps_causal': 0,
                'clumps_not_causal': 0,
                'causal_found': 0,
                'causal_not_found': 1.0,
              }
                result_set.append(results)
                continue

        clumps_data['SP2_list'] = clumps_data.SP2.apply(get_snps)
        clumps_data.sort_values(by=["CHR", "BP"],
                       inplace = True,)

        snp2list = defaultdict(list)
        last_pos = -1
        last_chrom = -1
        last_SNP = -1
        for key, value in clumps_data.iterrows():
            if last_chrom != value.CHR or (last_chrom == value.CHR and abs(last_pos - value.BP) > KB * 1000):
                if last_SNP != -1:
                    snp2list[last_SNP] = set(flatten(snp2list[last_SNP]))
                last_SNP = value.SNP
            last_chrom = value.CHR
            last_pos = value.BP
            snp2list[last_SNP].append(value.SP2_list+[value.SNP])
        snp2list[last_SNP] = set(flatten(snp2list[last_SNP]))
        snp2list = dict(snp2list)

        # make dict with results for this paramset
        clumps_causal = []
        causal_found = []
        for snp in causal_set:
            for key in snp2list:
                if snp in snp2list[key]:
                    clumps_causal.append(key)
                    causal_found.append(snp)
                    break
        clumps_causal = set(clumps_causal)
        causal_found = set(causal_found)
        clumps_not_causal = set(snp2list.keys()) - set(clumps_causal)
        causal_not_found = causal_set - set(causal_found)


        results = {**param,
                    'clumps_total': len(snp2list),
                    'clumps_causal': round(len(clumps_causal)/len(snp2list), 5),
                    'clumps_not_causal': round(len(clumps_not_causal)/len(snp2list), 5),
                    'causal_found': round(len(causal_found)/K, 5),
                   'causal_not_found': round(len(causal_not_found)/K, 5),
                  }

        result_set.append(results)




result_table = pd.DataFrame(result_set)
result_table.to_csv(FILE_RESULTS, sep='\t', index=False)
