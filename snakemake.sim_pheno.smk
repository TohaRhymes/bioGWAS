import yaml
import os
from pathlib import Path
import pandas as pd



# ===================================
# READ FILES
# ===================================

module_sim_pheno_final_outputs = [
        os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_chosen_snps.tsv"),
        os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_phenos.tsv"),
]

rule module_sim_pheno_all:
    input:
        inputc = module_sim_pheno_final_outputs
        
        

# ===================================
# SIM PHENO 
# ===================================


rule recode_merged:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
    output:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf")
    params:
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim"),
        CAUSAL_MAF_MIN=CAUSAL_MAF_MIN,
        CAUSAL_MAF_MAX=CAUSAL_MAF_MAX 
    shell:
        f"""
        {PLINK2_PATH} --bfile {{params.data}} \
        --recode vcf \
        --maf {{params.CAUSAL_MAF_MIN}} \
        --max-maf {{params.CAUSAL_MAF_MAX}} \
        --out {{params.data}}
        """
    message: "Recoding merged plinks's binary into VCF format {params.data}.vcf."
    onsuccess: "Recoding into VCF format {params.data}.vcf completed."
    onerror: "Error in recoding into VCF format {params.data}.vcf. That is probably plink2 error. Check logs and in/out files, their formats and try again. Logs are in: {log}."
        
        
        
rule transform_gtf:
    input:
        gtf=GTF_IN
    output:
        filt_gtf=os.path.join(DATA_DIR,f"{get_file_name(GTF_IN)}_filt_sort.gtf")
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/script_make_gtf.sh')} {{input.gtf}} {{output.filt_gtf}}
        """
    message: "Transforming+filtering {input.gtf} gtf-file."
    onsuccess: "GTF file {input.gtf} transformed successfully into {output.filt_gtf}."
    onerror: "Error transforming GTF file {input.gtf} into {output.filt_gtf}. Check logs and in/out files, their formats and try again. Logs are in: {log}."
        
        
        
rule annotate_vcf:
    input:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf"),
        gtf=os.path.join(DATA_DIR,f"{get_file_name(GTF_IN)}_filt_sort.gtf")
    output:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno.vcf"),
    shell:
        f"""
        {BEDTOOLS_PATH} closest -a  {{input.vcf}} -b {{input.gtf}} > {{output.vcf}}
        """     
    message: "Annotating VCF file {input.vcf} using annotations in {input.gtf} and bedtools."
    onsuccess: "VCF file {input.vcf} annotated successfully. Annotated file is in: {output.vcf}."
    onerror: "Error annotating VCF file {input.vcf}. That can be problem in bedtools. Check logs and in/out files, their formats and try again. Logs are in: {log}."
        
        
rule find_included:
    input:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno.vcf")
    output:
        included_txt = os.path.join(DATA_DIR, f'{PATTERN}_{CAUSAL_ID}_included_genes_snps.txt'),
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/find_included.sh')} {{input.vcf}} {{output.included_txt}}
        """      
    message: "Identifying included genes and SNPs for {input.vcf}."
    onsuccess: "Successfully identified included genes and SNPs for {input.vcf}."
    onerror: "Error in identifying included genes and SNPs for {input.vcf}. Check logs and in/out files, their formats and try again. Logs are in: {log}."
        
            
# Understand which to execute: select snps per pathway, or use preselected snps
causal_genes_output = []
if not SNPS_PROVIDED:
    causal_genes_output = os.path.join(DATA_DIR, f'{PATTERN}_{CAUSAL_ID}_causal_geneset_snps.txt')
    
gmt_file = []
path_file = []
included_txt_file = []
snps_file = []
if SNPS_PROVIDED:
    snps_file=SNPS_FILE
else:
    gmt_file = GMT_IN
    path_file = PATH_FILE
    included_txt_file = os.path.join(DATA_DIR, f'{PATTERN}_{CAUSAL_ID}_included_genes_snps.txt')
    


rule get_snps_list:
    input:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno.vcf"),
        gmt=gmt_file,
        path=path_file,
        included_txt = included_txt_file,
        snps_file=snps_file
    output:
        causal_genes = temp(causal_genes_output),
        annotated_snps = temp(os.path.join(DATA_DIR, f"{PATTERN}_{CAUSAL_ID}_annotated_snps.tsv")),
        tsv=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_chosen_snps.tsv")
    params:
        data_dir=DATA_DIR,
        pattern=PATTERN,
        pheno_pattern=CAUSAL_ID,
        K=K,
        k=k
    run:
        if SNPS_PROVIDED:
            shell(f"""
            {os.path.join(BIOGWAS_PATH, './pipeline_utils/get_snps_by_snps.py')} \
            {{params.data_dir}} \
            {{params.pattern}} \
            {{params.pheno_pattern}} \
            {{input.vcf}} \
            {{input.snps_file}} \
            {{params.K}}
            """) 
        else:
            shell(f"""
            {os.path.join(BIOGWAS_PATH, './pipeline_utils/get_snps_set.py')} \
            {{params.data_dir}} \
            {{params.pattern}} \
            {{params.pheno_pattern}} \
            {{input.vcf}} \
            {{input.gmt}} \
            {{input.path}} \
            {{params.K}} \
            {{params.k}}
            """)     
    message: "Generating causal SNPs list."
    onsuccess: "Causal SNPs list generated successfully in {output.tsv}."
    onerror: "Error in generating causal SNPs list. Check logs and in/out files, their formats and try again. Logs are in: {log}."




rule extract_snps:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
        tsv=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_chosen_snps.tsv"),
    output:
        bed=temp(os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps.bed")),
        bim=temp(os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps.bim")),
        fam=temp(os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps.fam"))
    params:
        filt_sim_data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim"),
        out_data=os.path.join(DATA_DIR, f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps"),
    shell:
        f"""
        {PLINK2_PATH} \
        --bfile {{params.filt_sim_data}} \
        --extract range {{input.tsv}} \
        --make-bed --out {{params.out_data}}
        """   
    message: "Extracting selected causal SNPs from plinks binary {params.filt_sim_data}."
    onsuccess: "Selected causal SNPs extracted successfully from plinks binary {params.filt_sim_data}. Corresponding plinks binary is: {params.out_data}"
    onerror: "Error in extracting selected causal SNPs from plinks binary {params.filt_sim_data}. That can be a problem in plink2. Check logs and in/out files, their formats and try again. Logs are in: {log}."


rule pheno_sim:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps.fam"),
        K_file=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_chosen_snps.tsv")
    output:
        tsv=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_phenos.tsv")
    params:
        data_dir=DATA_DIR,
        N=N,
        data=os.path.join(DATA_DIR, f"{PATTERN}_{CAUSAL_ID}_filt_sim_snps"),
        M_BETA=M_BETA,
        SD_BETA=SD_BETA,
        GEN_VAR=GEN_VAR,
        H2S=H2S,
        P_INDEPENDENT_GENETIC=P_INDEPENDENT_GENETIC,
        THETA=THETA,
        PHI=PHI,
        ALPHA=ALPHA,
        SEED=SEED
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/pheno_sim.R')} \
        ./ \
        {{params.data}} \
        {{output.tsv}} \
        {{params.N}} \
        {{input.K_file}} \
        {{M_BETA}} \
        {{SD_BETA}} \
        {{GEN_VAR}} \
        {{H2S}} \
        {{P_INDEPENDENT_GENETIC}} \
        {{THETA}} \
        {{PHI}} \
        {{ALPHA}} \
        {{SEED}}
        """   
    message: "Simulating phenotypes from selected causal SNPs ({params.data})."
    onsuccess: "Phenotype simulation from selected causal SNPs ({params.data}) completed}. Simulated phenotype is in {output.tsv}."
    onerror: "Error in simulating phenotypes from selected causal SNPs ({params.data}). That can be an R error, or Phenotype simulator inner error. Check logs and in/out files, their formats and try again. Logs are in: {log}."
