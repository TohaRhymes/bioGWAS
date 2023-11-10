#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint

from utils_1 import params, Ks


# Params for one K -- imported from utils_1.py
pprint(params)
print(f"Amount of parameters sets: {len(params)}")
# Ks itself
pprint(Ks)
print(f"Amount of Ks: {len(Ks)}")

for K in Ks:
    k = K // 2
    print(f"Started for K={K}")
    for param in params:
        m_beta = param["m_beta"]
        sd_beta = param["sd_beta"]
        gen_var = param["gen_var"]
        alpha = param["alpha"]
        theta = param["theta"]
        pIndep = param["pIndep"]
        
        pat = "test1"
        cas_id = f"K{K}"
        sim_id = f"m{m_beta}_sd{sd_beta}_gv{gen_var}_alpha{alpha}_theta{theta}_pIndep{pIndep}"
        
        command = f"""docker run \
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
        """
        print(command)
        process = subprocess.Popen(
             command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print("\n=============================\nOUTPUT: \n", output)
        print("\n=============================\nERRORS: \n", error)
        
print('Done with params iterations')
