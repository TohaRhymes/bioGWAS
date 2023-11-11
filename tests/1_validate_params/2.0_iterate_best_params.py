#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint
import random
from tqdm import tqdm

import pandas as pd


from utils_2 import VALIDATION_COMPARISON, IDs


snps = pd.read_csv(VALIDATION_COMPARISON, sep='\t').drop(['clumps_not_causal', 'causal_not_found'], axis=1)
snps['F1']=2/(1/snps.clumps_causal+1/snps.causal_found)
snps = snps[snps.F1 > 0.9]
params = snps.T.to_dict()

pprint(params)
print(f"Amount of parameters sets: {len(params)}")


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
    
    causal_id = IDs['casual_id'].format({'K':K})

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
        --gmt_file /wd/data/c2.cp.kegg.v2023.1.Hs.symbols.gmt \
        --causal_pathways /wd/data/pathways.csv \
        --maf_filter 0.05 \
        --N 1000 \
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
