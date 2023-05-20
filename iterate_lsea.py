#!/usr/bin/env python

import subprocess
import os
from tqdm import tqdm


# variables
for i in tqdm(range(0, 50)):
    for k in [0, 15]:
        # best constants
        K=30
        m_beta=0.5
        sd_beta=0.1
        gen_var=0.5
        h2s=0.5

        phenos_id = f'ph_{k}_{i}'
        sim_id=f"K{K}_m{m_beta}_sd{sd_beta}_gv{gen_var}_h2s{h2s}"
        LSEA = "../LSEA/LSEA_2.3.py" 
        snake_gwas = f"data/chr_{phenos_id}_{sim_id}_gwas.tsv"
        universe = "../LSEA/data/universe.json"
        out_dir = f"data/lsea_{phenos_id}_{sim_id}"
        plink_dir = "/home/achangalidi/tools/plink/"
        bfile = "data/chr_filt_sim"
        pval = 0.00000000729730
        result_file = os.path.join(out_dir, f"universe_result_{pval}.tsv")


        command_sim = f"snakemake \
                --cores 8 \
                --config \
                m_beta={m_beta} \
                sd_beta={sd_beta} \
                gen_var={gen_var} \
                h2s={h2s} \
                K={K} \
                k={k} \
                draw_flag={False} \
                sim_id={sim_id} \
                phenos_id={phenos_id}"
        command_lsea = f"python3 \
                          {LSEA} \
                         -input {snake_gwas} \
                         -universe {universe} \
                         -out {out_dir} \
                         -plink_dir {plink_dir} \
                         -bfile {bfile} \
                         -column_names chr pos rsid pval \
                         -p {pval}"
        
        
        print('-------------\nROUND 1\n-------------')
        print(command_sim)
        process = subprocess.Popen(
            command_sim, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print(output)
        print(error)
        print('-------------\nROUND 2\n-------------')
        print(command_lsea)
        process = subprocess.Popen(
            command_lsea, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print(output)
        print(error)
        
print('Done with lsea iterations')