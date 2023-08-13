#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint
import random
from tqdm import tqdm

import pandas as pd


snps = pd.read_csv('data3/PAT_ph_compare_results.tsv', sep='\t').drop(['clumps_not_causal', 'causal_not_found'], axis=1)
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
    h2s = param["h2s"]
    theta = param["theta"]
    pIndep = param["pIndep"]
    phi = param["phi"]
    alpha = param["alpha"]
    causal_id = f"ph_best_K{K}"

    seeds = random.sample(range(1, 9999999), 20)
    for seed in seeds:
        sim_id = f"m{m_beta}_sd{sd_beta}_gv{gen_var}_h2s{h2s}_theta{theta}_pIndep{pIndep}_phi{phi}_alpha{alpha}_seed{seed}"
        command = f"snakemake \
        --cores 8 \
        --config \
        m_beta={m_beta} \
        sd_beta={sd_beta} \
        gen_var={gen_var} \
        h2s={h2s} \
        theta={theta} \
        pIndep={pIndep} \
        phi={phi} \
        alpha={alpha} \
        K={K} \
        k={k} \
        draw_flag={False} \
        sim_id={sim_id} \
        causal_id={causal_id} \
        seed={seed}"
        print(command)
        process = subprocess.Popen(
             command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print("\n=============================\nOUTPUT: \n", output)
        print("\n=============================\nERRORS: \n", error)
        
print('Done with params iterations')
