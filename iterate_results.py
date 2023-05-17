#!/usr/bin/env python
# coding: utf-8

# In[165]:


from itertools import product
from pprint import pprint
import pandas as pd
import os
import subprocess
import re
from collections import defaultdict

from tqdm import tqdm

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


# In[24]:


CAUSAL_SNP_FILE = "data/chr_ph_sperm_K{}_chosen_snps.tsv"
GWAS_FILE = "data/chr_ph_sperm_K{}_m{}_sd{}_gv{}_h2s{}_gwas.tsv"


# ## Params for one K

# In[4]:


Ks = [10, 30, 50]

# try these m_betas
m_betas = (
    0.01, 
    0.05, 
    0.1, 
    0.2, 
    0.3, 
    0.5)
# use these sd for m_betas above respectively
sd_betas = (
    (0.001, 0.005),
    (0.001, 0.01),
    (0.001, 0.01, 0.05),
    (0.01, 0.05, 0.1),
    (0.01, 0.05, 0.1),
    (0.01, 0.1, 0.3),
)

assert len(m_betas) == len(sd_betas)
# make products
m_sd_comb = []
for mb, sds in zip(m_betas, sd_betas):
    m_sd_comb += list(product([mb], sds))
# for every of combinations above, use these genvar and sd
gen_vars = (0.05, 0.5)
h2ss = (0.1, 0.5)
gv_h2s_comb = list(product(gen_vars, h2ss))

# print("These combinations of m and sd:")
# pprint(m_sd_comb)
# print("These combinations of gv and h2s:")
# pprint(gv_h2s_comb)

params = [
    {"m_beta": m, "sd_beta": sd, "gen_var": gv, "h2s": h2s}
    for (m, sd), (gv, h2s) in list(product(m_sd_comb, gv_h2s_comb))
]
# pprint(params)


# ## All params set

# In[13]:


print(f'Set of K: {Ks}')
print(f"For every K there are {len(params)} param sets (in variable `params`)")
print("\nSet selection process: ")
print(f"* for every `m_beta` from this set: {m_betas}")
print(f"* `sd_beta` was selected from this set respectively (2-3 `sd_beta` for every `m_beta`):")
pprint(sd_betas)
print(f"* for these {len(m_sd_comb)} combinations of `sd_beta` and `m_beta`, it was used 4 combinations of `gen_var` and `h2s`:")
pprint(gv_h2s_comb)


# ## Make table of results

# In[166]:


# result_set = []

# ## SNPs set for every k
# for param in tqdm(params):
#     for K in Ks:
#         m = param["m_beta"]
#         sd = param["sd_beta"]
#         gv = param["gen_var"]
#         h2s = param["h2s"]

#         causal_snp_file = CAUSAL_SNP_FILE.format(K)
#         gwas_file = GWAS_FILE.format(K, m, sd, gv, h2s)
#         # count number of snps
#         with open(gwas_file, "rb") as f:
#             N = sum(1 for _ in f) - 1
#         # SET OF CAUSAL SNPS
#         causal_snp = pd.read_csv(causal_snp_file, sep='\t', header=None, names=['chr', 'pos', 'pos_e', 'gene'])
#         causal_set = set(causal_snp.chr.astype(str)+':'+causal_snp.pos.astype(str))
#         print(causal_set)

#         # make clump file
#         FILE = gwas_file.replace('.tsv', '')
#         BFILE = "data/chr_filt_sim"
#         PVAL_CUTOFF = 0.05/N
#         KB = 1000

#         command_format = f"""sed '1{{ s/chr/CHR/; s/rsid/SNP/; s/pos/BP/; s/n/NMISS/; s/beta/BETA/; s/se/SE/; s/r2/R2/; s/t/T/; s/pval/P/;}}' {FILE}.tsv | awk '{{$1=$1}};1'| tr -s ' ' '\t' > {FILE}.qassoc"""
#         command_clump = f"""plink         --bfile {BFILE}         --allow-no-sex         --clump {FILE}.qassoc         --clump-p1 {'{:.30f}'.format(PVAL_CUTOFF)}         --clump-p2 {'{:.30f}'.format(PVAL_CUTOFF)}         --clump-r2 0.2         --clump-kb {KB}         --out {FILE}"""
#         print(command_format)
#         run(command_format)
#         print(command_clump)
#         run(command_clump)
#         print('done clump')

#         # read clump file and make list of snps for every clump
#         try:
#             clumps_data = pd.read_fwf(f"{FILE}.clumped", sep='\t')
#         except FileNotFoundError:
#             results = {"K": K, 
#            "m_beta" : m, 
#            "sd_beta" : sd, 
#            "gen_var" : gv, 
#            "h2s" : h2s,
#             'clumps_total': 0,
#             'clumps_causal': 0,
#             'clumps_not_causal': 0,
#             'causal_found': 0,
#            'causal_not_found': 1,
#           }
#             result_set.append(results)
#             continue
            
#         clumps_data['SP2_list'] = clumps_data.SP2.apply(get_snps)
#         clumps_data.sort_values(by=["CHR", "BP"],
#                        inplace = True,)
#         clumps_data.head(20)
#         snp2list = defaultdict(list)
#         last_pos = -1
#         last_chrom = -1
#         last_SNP = -1
#         for key, value in clumps_data.iterrows():
#             if last_chrom != value.CHR or (last_chrom == value.CHR and abs(last_pos - value.BP) > KB * 1000):
#                 if last_SNP != -1:
#                     snp2list[last_SNP] = set(flatten(snp2list[last_SNP]))
#                 last_SNP = value.SNP
#             last_chrom = value.CHR
#             last_pos = value.BP
#             snp2list[last_SNP].append(value.SP2_list+[value.SNP])
#         snp2list[last_SNP] = set(flatten(snp2list[last_SNP]))
#         snp2list = dict(snp2list)

#         # make dict with results for this paramset
#         clumps_causal = []
#         causal_found = []
#         for snp in causal_set:
#             for key in snp2list:
#                 if snp in snp2list[key]:
#                     clumps_causal.append(key)
#                     causal_found.append(snp)
#                     break
#         clumps_not_causal = set(snp2list.keys()) - set(clumps_causal)
#         causal_not_found = causal_set - set(causal_found)


#         results = {"K": K, 
#                    "m_beta" : m, 
#                    "sd_beta" : sd, 
#                    "gen_var" : gv, 
#                    "h2s" : h2s,
#                     'clumps_total': len(snp2list),
#                     'clumps_causal': round(len(clumps_causal)/len(snp2list), 5),
#                     'clumps_not_causal': round(len(clumps_not_causal)/len(snp2list), 5),
#                     'causal_found': round(len(causal_found)/K, 5),
#                    'causal_not_found': round(len(causal_not_found)/K, 5),
#                   }

#         result_set.append(results)


# # In[167]:


# result_table = pd.DataFrame(result_set)
# result_table.to_csv("data/chr_ph_sperm_compare_results.tsv", sep='\t', index=False)
