import yaml
import os
from pathlib import Path
import pandas as pd



# ===================================
# READ FILES
# ===================================

module_gwas_final_outputs = [
    os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas.tsv"),
]

rule module_sim_gwas_all:
    input:
        inputd = module_gwas_final_outputs
        
        

# ===================================
# MAKE GWAS
# ===================================



# todo make custom gwas        
rule gwas:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
        pheno=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_phenos.tsv"),
    output:
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas.tsv"),
        qassoc=temp(os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas.qassoc"))
    params:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim"),
        pheno=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_phenos"),
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas")
    message: 
        """
        Description: Performing GWAS on plink binary {params.bed} and {input.pheno}.
        I/O info:    Summary statistics are in {output.gwas}
        Errors:      That can be a plink's gwas error. Check logs and in/out files, their formats and try again.
        """
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/gwas_analysis_binary.sh' if BINARY_PHENO else './pipeline_utils/gwas_analysis.sh')} \
        {{params.bed}} \
        {{params.pheno}} \
        {{params.gwas}} \
        {PLINK_PATH}
        """   
           
