#!/usr/bin/env python

from itertools import product
from collections import defaultdict
from copy import deepcopy
import subprocess
from pprint import pprint
import random
from tqdm import tqdm
import os
import glob

import pandas as pd
import numpy as np


from utils_1 import N
from utils_2 import params, IDs

DATA_DIR='data'

pprint(params)
print(f"Amount of parameters sets: {len(params)}")

for param_key in tqdm(params):
    param = params[param_key]

    m_beta = param["m_beta"]
    sd_beta = param["sd_beta"]
    gen_var = param["gen_var"]
    alpha = param["alpha"]
    theta = param["theta"]
    pIndep = param["pIndep"]
    param["K"] = int(param["K"])
    K = int(param["K"])
    k = K // 2


    pat = IDs['pattern']
    cas_id = IDs['casual_id'].format(**{'K':K})
    
    sim_id_check = IDs['sim_id'].format(**{'seed':'*'}, **param)
    if len(glob.glob(os.path.join(DATA_DIR, f'*{cas_id}_{sim_id_check}_gwas.tsv'))) == 0:
        seeds = random.sample(range(1, 9999999), 20)
        for seed in seeds:
            sim_id = IDs['sim_id'].format(**{'seed':seed}, **param)
            command = f"docker run \
            -v /media/DATA/gwasim/round2/bioGWAS/tests:/wd \
            biogwas \
            /bioGWAS/biogwas.py \
            --dependencies /dependencies.yaml \
            --threads 8 \
            --data_dir /wd/1_validate_params/data \
            --img_dir /wd/1_validate_params/images \
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