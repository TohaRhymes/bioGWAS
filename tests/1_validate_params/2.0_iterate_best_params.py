#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint
import random
from tqdm import tqdm

import pandas as pd
import numpy as np


from utils_1 import Ks, N, VALIDATION_COMPARISON
from utils_1 import params as params_launch1
from utils_2 import IDs

def top_n_with_ties(df, column, n=3):
    top_n_values = df[column].nlargest(n).unique() 
    return df[column].isin(top_n_values)


param_names = list(params_launch1[0].keys())

# select params we're going to iterate over the second time
snps = pd.read_csv(VALIDATION_COMPARISON, sep=',')
cols_for_flags = ['F1_avg', 'F1 (K=10)', 'F1 (K=20)', 'F1 (K=30)']
flags = []
TOP_K=3

for col in cols_for_flags:
    flags.append(top_n_with_ties(snps, col, TOP_K))

snps = snps[np.any(flags, axis=0)][param_names]
params = snps.T.to_dict()

pprint(params)
print(f"Amount of parameters sets: {len(params)}")


for K in Ks:
    k = K // 2
    print(f"Started for K={K}")
    for param_key in tqdm(params):
        param = params[param_key]
        
        m_beta = param["m_beta"]
        sd_beta = param["sd_beta"]
        gen_var = param["gen_var"]
        alpha = param["alpha"]
        theta = param["theta"]
        pIndep = param["pIndep"]
        
        
        pat = IDs['pattern']
        cas_id = IDs['casual_id'].format(**{'K':K})

        seeds = random.sample(range(1, 9999999), 20)
        for seed in seeds:
            sim_id = IDs['sim_id'].format(**{'K':K, 'seed':seed}, **param)
            command = f"docker run \
            -v /media/DATA/gwasim/round2/bioGWAS/tests:/wd \
            biogwas \
            /bioGWAS/biogwas.py \
            --dependencies /dependencies.yaml \
            --threads 8 \
            --input_dir /wd/data \
            --data_dir /wd/1_validate_params/data \
            --img_dir /wd/1_validate_params/images \
            --vcf_in_flag \
            --input_list  /wd/data/chr.list \
            --ids_file  /wd/data/EUR_SAMPLES_ID.txt \
            --anno_file /wd/data/gencode.v37.annotation.gtf \
            --gmt_file /wd/data/h.all.v2023.1.Hs.symbols.gmt \
            --use_causal_snps \
            --causal_snps /wd/data/pathways_{K}.csv \
            --maf_filter 0.05 \
            --N {N} \
            --m_beta {m_beta} \
            --sd_beta {sd_beta} \
            --gen_var {gen_var} \
            --alpha {alpha} \
            --theta {theta} \
            --p_independent_genetic {pIndep} \
            --K {K} \
            --k {k} \
            --no-draw_flag \
            --pattern {pat} \
            --causal_id {cas_id}  \
            --sim_id {sim_id} \
            --seed {seed}"
            print(command)
            process = subprocess.Popen(
                 command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
            )
            output, error = process.communicate()
            print("\n=============================\nOUTPUT: \n", output)
            print("\n=============================\nERRORS: \n", error)
        
print('Done with params iterations')