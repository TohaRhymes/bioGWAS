#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint
import os
from tqdm import tqdm


# to create path with all genes: echo -e "ALL_GENES\thttps://github.com/TohaRhymes/GWAS_simulator\t$(awk '{print $1}' gencode.v37.annotation.gtf.loc | sort | uniq | paste -s -d '\t')" > all.gmt

files = [
    "/wd/data/path_small.txt",
    "/wd/data/path_medium.txt",
    "/wd/data/path_big.txt",
    "/wd/data/path_random.txt",
]
gmts = [
    "/wd/data/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
    "/wd/data/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
    "/wd/data/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
    '/wd/data/all.gmt',
]
# total amnt of causal SNPs
K=30
# amnt of causal SNPs from pathway
ks = [
    15, 
    15, 
    15, 
    30
]
# amnt of iterations
ITER = 30


for file, gmt, k in zip(files, gmts, ks):
    file_just_name = os.path.splitext(os.path.basename(file))[0]
    print(f"Started for {file}")
    pat = "test10000"
    for i in tqdm(range(ITER)):
        cas_id = f"{file_just_name}_{i}"
        command = f"""docker run \
        -v /media/DATA/gwasim/round2/bioGWAS/tests:/wd \
        biogwas \
        /bioGWAS/biogwas.py \
        --dependencies /dependencies.yaml \
        --threads 8 \
        --input_dir /wd/data \
        --data_dir /wd/3_pathways/in_data \
        --img_dir /wd/3_pathways/images \
        --vcf_in_flag \
        --input_list  /wd/data/chr.list \
        --ids_file  /wd/data/EUR_SAMPLES_ID.txt \
        --anno_file /wd/data/gencode.v37.annotation.gtf \
        --gmt_file {gmt} \
        --causal_pathways {file} \
        --maf_filter 0.05 \
        --N 10000 \
        --m_beta 0.5 \
        --sd_beta 0.01 \
        --gen_var 0.5 \
        --alpha 0.5 \
        --theta 1 \
        --p_independent_genetic 0 \
        --K {K} \
        --k {k} \
        --no-draw_flag \
        --pattern {pat} \
        --causal_id {cas_id}  \
        --sim_id {cas_id} \
        """  

        print(command)
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print("\n=============================\nOUTPUT: \n", output)
        print("\n=============================\nERRORS: \n", error)

print("Done with params iterations")
