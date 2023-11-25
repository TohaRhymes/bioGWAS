import os


from utils_1 import N, VALIDATION_COMPARISON
from utils_1 import params as params_launch1

import pandas as pd
import numpy as np

from collections import defaultdict

from copy import deepcopy



def top_n_with_ties(df, column, n=3):
    top_n_values = df[column].nlargest(n).unique() 
    return df[column].isin(top_n_values)


param_names = list(params_launch1[0].keys())

# select params we're going to iterate over the second time
snps = pd.read_csv(VALIDATION_COMPARISON, sep=',')
cols_for_flags = ['F1_avg', 
                  'F1 (K=10)', 
                  'F1 (K=20)', 
                  'F1 (K=30)']
Ks_for_cols = [[10,20,30], 
               [10,], 
               [20,], 
               [30,]]
flags = defaultdict(list)
TOP_K=2

for Ks, col in zip(Ks_for_cols, cols_for_flags):
    for cur_K in Ks:
        flags[cur_K].append(top_n_with_ties(snps, col, TOP_K))

Ks_snps = []
for cur_K in flags:
    cur_snps = deepcopy(snps[np.any(flags[cur_K], axis=0)][param_names])
    cur_snps['K']=cur_K
    Ks_snps.append(cur_snps)
snps = pd.concat(Ks_snps).reset_index()

params = snps.T.to_dict()




# =============
# filenames
# =============
IDs = {'pattern': "test10000",
        'casual_id': "best_hallmark_K{K}",
       'sim_id':"m{m_beta}_sd{sd_beta}_gv{gen_var}_alpha{alpha}_theta{theta}_pIndep{pIndep}_seed{seed}"}



# =============
# filenames
# =============
# aux
DATA_DIR = "data"
PATTERN = IDs['pattern']
CAUSAL_ID = IDs['casual_id']
SIM_ID = IDs['sim_id']
res_CAUSAL_ID = CAUSAL_ID.replace("K{K}", "")
# to import
CAUSAL_SNP_FILE = os.path.join(DATA_DIR, f"{PATTERN}_{CAUSAL_ID}_chosen_snps.tsv")
GWAS_FILE = os.path.join(DATA_DIR, f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas.tsv")
FILE_RESULTS=os.path.join(DATA_DIR, f"{PATTERN}_{res_CAUSAL_ID}_compare_results.tsv")
BFILE=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim")