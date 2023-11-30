#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint
import os

from tqdm import tqdm

import sys


def launch_command(command):
    print(command)
    process = subprocess.Popen(
         command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
    )
    output, error = process.communicate()
    print("\n=============================\nOUTPUT: \n", output)
    print("\n=============================\nERRORS: \n", error)


GWAS_OUT_DATA = "./in_data"
PATH_OUT_DATA = "./data_enrich"
MAGMA_SCRIPTS_DIR="./magma_scripts"

GENE_LOC = "../data/gencode.v37.annotation.gtf"
GENE_SETS = "../data/c2.cp.kegg.v2023.1.Hs.symbols.gmt"
ANNO_FILE = "./data_enrich/magma_anno"


BFILE = "./in_data/test10000_filt_sim"

N = 10000
ITER = 30

models = {
"linreg":"linreg",
"mean":"snp-wise=mean",
"top":"snp-wise=top",
}

total_command=""


for pathname in [
    "path_small", 
    "path_medium", 
    "path_big", 
    "path_random"
]:
    for i in tqdm(range(ITER)):
        cas_id = f"{pathname}_{i}"
        SIM_ID = f"test10000_{cas_id}_{cas_id}"

        PHENO = f"{GWAS_OUT_DATA}/{SIM_ID}_phenos.tsv"
        GWAS = f"{GWAS_OUT_DATA}/{SIM_ID}_gwas.tsv"

        command_bfile = f"plink \
                --bfile {BFILE} \
                --pheno {PHENO} \
                --allow-no-sex \
                --make-bed \
                --out {BFILE}_{SIM_ID}" # bed bim fam с фенотипом
        
        launch_command(command_bfile)
        
        for model in models:

            GENES_OUT = f"{PATH_OUT_DATA}/{SIM_ID}_magma_genes_{model}"
            GENE_SETS_OUT = f"{PATH_OUT_DATA}/{SIM_ID}_magma_sets_{model}"
            SCRIPT_FILE=f"{MAGMA_SCRIPTS_DIR}/{SIM_ID}_{model}.sh"

            # Method 1: using pathfiles

            command1 = f"magma \
                        --bfile {BFILE} \
                        --gene-annot {ANNO_FILE}.genes.annot \
                        --pval {GWAS} \
                        use=rsid,pval \
                        N={N} \
                        --gene-model {models[model]} \
                        --out {GENES_OUT}" # Output: /home/achangalidi/ukb_finngen/gwassim_check/path_sim/pathway_analysis/*_magma_genes.genes.raw

            command2 = f"magma \
                            --gene-results {GENES_OUT}.genes.raw \
                            --set-annot {GENE_SETS}.ssv \
                            --out {GENE_SETS_OUT}" # Output: /home/achangalidi/ukb_finngen/gwassim_check/path_sim/pathway_analysis/*_magma_sets.gsa.out

            command1_bfile = f"magma \
                        --bfile {BFILE}_{SIM_ID} \
                        --gene-annot {ANNO_FILE}.genes.annot \
                        --gene-model {models[model]} \
                        --out {GENES_OUT}"

            command2_bfile = f"magma \
                            --gene-results {GENES_OUT}.genes.raw \
                            --set-annot {GENE_SETS}.ssv \
                            --out {GENE_SETS_OUT}"
            command_rm_bfile = f"rm -rf {BFILE}_{SIM_ID}.*"
            if model =="linreg":
                commands = [
                    command1_bfile, 
                    command2_bfile,
                    command_rm_bfile
                ]
            else:
                commands = [
                    command1, 
                    command2]
            total_command = '\n\n'.join(commands)+'\n'

            with open(SCRIPT_FILE, 'w') as f:
                f.write(total_command)



print('Finished')