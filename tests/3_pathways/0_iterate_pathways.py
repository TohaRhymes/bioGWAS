#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint
import os


# to create path with all genes: echo -e "ALL_GENES\thttps://github.com/TohaRhymes/GWAS_simulator\t$(awk '{print $1}' gencode.v37.annotation.gtf.loc | sort | uniq | paste -s -d '\t')" > all.gmt

files = [
    "in_data/path_small.txt",
    "in_data/path_medium.txt",
    "in_data/path_big.txt",
             'in_data/path_random.txt',
]
gmts = [
    "../../1000genomes/data2/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
    "../../1000genomes/data2/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
    "../../1000genomes/data2/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
            '../../1000genomes/data2/all.gmt',
]

ks = [
    15, 
    15, 
    15, 
    30
]
ITER = 30


for i in range(ITER):
    for file, gmt, k in zip(files, gmts, ks):
        file_just_name = os.path.splitext(os.path.basename(file))[0]
        print(f"Started for {file}")

        cas_id = f"{file_just_name}_{i}"
        command = f"../../1000genomes/./biogwas.py \
        -d ../../1000genomes/./dependencies.yaml \
        --input_dir ../../1000genomes/data2  \
        --data_dir out_data  \
        --img_dir img  \
        --vcf_in_flag \
        -N 1000 \
        --input_list  ../../1000genomes/data2/chr.list \
        --ids_file  ../../1000genomes/data2/EUR_SAMPLES_ID.txt \
        --anno_file ../../1000genomes/data2/gencode.v37.annotation.gtf \
        --gmt_file {gmt}  \
        -cp {file} \
        -K 30 \
        -k {k} \
        --pattern SIM  \
        --causal_id {cas_id}  \
        --sim_id {cas_id} \
        -T 8"

        print(command)
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print("\n=============================\nOUTPUT: \n", output)
        print("\n=============================\nERRORS: \n", error)

print("Done with params iterations")