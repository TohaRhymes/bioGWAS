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


GWAS_OUT_DATA = "/home/achangalidi/ukb_finngen/gwassim_check/path_sim/out_data"
PATH_OUT_DATA = "/home/achangalidi/ukb_finngen/gwassim_check/path_sim/pathway_analysis"

GENE_LOC = "/home/achangalidi/ukb_finngen/1000genomes/data2/gencode.v37.annotation.gtf"
GENE_SETS = (
    "/home/achangalidi/ukb_finngen/1000genomes/data2/c2.cp.kegg.v2023.1.Hs.symbols.gmt"
)
ANNO_FILE = (
    "/home/achangalidi/ukb_finngen/gwassim_check/path_sim/pathway_analysis/magma_anno"
)


BFILE = f"{GWAS_OUT_DATA}/SIM_filt_sim"

N = 1000


# =====================================================
# Preparation (just one time) -- all commands in bash
# =====================================================

# GWAS_OUT_DATA=/home/achangalidi/ukb_finngen/gwassim_check/path_sim/out_data
# GENE_LOC=/home/achangalidi/ukb_finngen/1000genomes/data2/gencode.v37.annotation.gtf
# GENE_SETS=/home/achangalidi/ukb_finngen/1000genomes/data2/c2.cp.kegg.v2023.1.Hs.symbols.gmt
# ANNO_FILE=/home/achangalidi/ukb_finngen/gwassim_check/path_sim/pathway_analysis/magma_anno
# BFILE=${GWAS_OUT_DATA}/SIM_filt_sim

# ------------
# Step 0
# ------------
# CHANGE \t to spaces in gmt files

# sed 's/\t/ /g' ${GENE_SETS} | awk '{for (i=1; i<=NF; i++) if (i!=2) printf $i (i==NF? RS : FS)}' > ${GENE_SETS}.ssv

# ------------
# Step 1
# ------------
# Make custom .loc file from gtf

# cat ${GENE_LOC} \
# | awk -F'\t' -v OFS='\t' \
# '$3=="gene"{split($9,a,";"); split(a[1],b," "); split(a[3],c," "); $NF="\""c[2]"\""; print}' \
# | awk -F'\t' \
# '!seen[$9]++ {gsub("\"", "", $9); gsub("chr", "", $1); print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
# > ${GENE_LOC}.loc

# ------------
# Step 2
# ------------
# Make magma annotations

# magma \
# --annotate window=1,0.5 \
# --snp-loc ${BFILE}.bim \
# --gene-loc  ${GENE_LOC}.loc \
# --out ${ANNO_FILE}


# pathname = sys.argv[1]  # path_small path_medium path_big path_random
ITER = 30

models = {
"linreg":"linreg",
"all":"multi=all",
"mean":"snp-wise=mean",
"top":"snp-wise=top",
"top15":"snp-wise=top,15",
}

total_command=""

for pathname in [
#     "path_small", 
#     "path_medium", 
#     "path_big", 
    "path_random",
]:
    for i in tqdm(range(ITER)):
        cas_id = f"{pathname}_{i}"
        SIM_ID = f"SIM_{cas_id}_{cas_id}"

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
            SCRIPT_FILE=f"magma_scripts/{SIM_ID}_{model}.sh"

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
            if model in ["linreg", "all"]:
                commands = [
                    command1_bfile, 
                    command2_bfile]
            else:
                commands = [
                    command1, 
                    command2]
            total_command = '\n'.join(commands)+'\n'

            with open(SCRIPT_FILE, 'w') as f:
                f.write(total_command)

print('Finished')