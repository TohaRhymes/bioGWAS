import yaml
import os
from pathlib import Path
import pandas as pd

# ===================================
# UTIL FUNCTIONS
# ===================================

def get_chrom(wildcards):
    return file2chrom[os.path.basename(wildcards.file)]

def get_file_name(full_path: str):
    return os.path.splitext(os.path.basename(full_path))[0]
    
def remove_ext(full_path: str):
    return os.path.splitext(full_path)[0]


# ===================================
# READING CONFIGURATION
# ===================================


BIOGWAS_PATH = config['biogwas_path']

PLINK_PATH = config['plink']
PLINK2_PATH = config['plink2']
BEDTOOLS_PATH = config['bedtools']
HAPGEN2_PATH = config['hapgen2']

VCF_IN_FLAG=config['vcf_in_flag']
VCF_IN_DIR = config['vcf_in_dir']
DATA_DIR = config['data_dir']
IMAGES_DIR = config['images_dir']

# VCFS_LIST -- file with structure: `path_to_vcf,chromosome` per row
VCFS_LIST = config['vcfs_list']
IDS_FILE = config['ids_file']
GTF_IN = config['anno_file']
GMT_IN = config['gmt_file']

# if use_causal_snps (default = True), causal SNPs are taken from `causal_snps` file
# otherwise -- randomly selected from `causal_snps` file.
SNPS_PROVIDED = config['use_causal_snps'] 

PATH_FILE = config['causal_pathways']
SNPS_FILE = config['causal_snps']

# General filter
MAF_FILTER = config['maf_filter']
# Filter for causal SNPs only: defaults are [0.0,1.0]
CAUSAL_MAF_MIN = config['causal_maf_min'] 
CAUSAL_MAF_MAX = config['causal_maf_max']

N = config['N']

K = config['K']
k = config['k']


M_BETA = config['m_beta']
SD_BETA = config['sd_beta']

GEN_VAR = config['gen_var']

H2S = 1.0
P_INDEPENDENT_GENETIC = config['p_independent_genetic']
THETA = config['theta']

PHI = 1.0
ALPHA = config['alpha']

# pattern of all files (starting from genotypes)
PATTERN = config['pattern']
# ID of chosen causal SNPS set for the simulations.
CAUSAL_ID = config['causal_id']
# ID of phenotypes simulation for selected causal SNPs
SIM_ID = config['sim_id']
DRAW_FLAG = config['draw_flag']


SEED = config['seed']



# ===================================
# READ FILES
# ===================================

with open(VCFS_LIST) as f:
    file_chrom = [line.strip().split(',') for line in f]
    file2chrom = {os.path.basename(file):chrom for file, chrom in file_chrom}
    files = list(map(lambda x: get_file_name(x), list(file2chrom.keys())))
    
    
VCFS_LIST_PATTERN = get_file_name(VCFS_LIST)
    
with open(os.path.join(DATA_DIR, f"{VCFS_LIST_PATTERN}_filt_sim.list"), 'w') as f:
    f.write('\n'.join(map(lambda x: f"{os.path.join(DATA_DIR, x)}_filt_sim", files)))
    

# ===================================
# SET OUTPUT FILES
# ===================================

draw_output = list()

if DRAW_FLAG:
    output_pca=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_pca.pdf")
    output_pca_indep=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_indep_pca.pdf")
    mh=os.path.join(IMAGES_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas_mh.pdf")
    qq=os.path.join(IMAGES_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas_qq.pdf")
    draw_output.append(output_pca)
    draw_output.append(output_pca_indep)
    draw_output.append(mh)
    draw_output.append(qq)

rule all:
    priority: 1000
    input:
        filt_bed=expand(os.path.join(DATA_DIR, "{file}_filt.bed"), file=files),
        filt_bim=expand(os.path.join(DATA_DIR, "{file}_filt.bim"), file=files),
        filt_fam=expand(os.path.join(DATA_DIR, "{file}_filt.fam"), file=files),
        filt_sim_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        filt_sim_bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        filt_sim_fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
        filt_anno_vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno.vcf"),
        chosen_snps=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_chosen_snps.tsv"),
        included_txt = os.path.join(DATA_DIR, f'{PATTERN}_{CAUSAL_ID}_included_genes_snps.txt'),
        phenos_tsv=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_phenos.tsv"),
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas.tsv"),
        draw_output = draw_output
        

# ===================================
# RUN SIMULATION
# ===================================



rule init_vcf_bfile:
    input:
        vcf=lambda wildcards: os.path.join(VCF_IN_DIR, os.path.basename(f"{wildcards.file}.vcf")) if VCF_IN_FLAG else []
    output:
        bed="{file}.bed" if VCF_IN_FLAG else [],
        bim="{file}.bim" if VCF_IN_FLAG else [],
        fam="{file}.fam" if VCF_IN_FLAG else []
    shell:
        f'''{PLINK2_PATH} --vcf {{input.vcf}} \
        --max-alleles 2 \
        --set-all-var-ids @:#\$r:\$a --new-id-max-allele-len 5000 missing \
        --make-bed  \
        --out {{wildcards.file}}'''
        
        
rule init_bfile_bfile:
    input:
        vcf=lambda wildcards: os.path.join(VCF_IN_DIR, os.path.basename(f"{wildcards.file}.bed")) if not VCF_IN_FLAG else []
    output:
        bed="{file}.bed" if not VCF_IN_FLAG else [],
        bim="{file}.bim" if not VCF_IN_FLAG else [],
        fam="{file}.fam" if not VCF_IN_FLAG else []
    params:
        bed=os.path.join(VCF_IN_DIR, os.path.basename("{file}"))
    shell:
        f'''{PLINK2_PATH}  --bfile {{params.bed}} \
        --max-alleles 2 \
        --set-all-var-ids @:#\$r,\$a --new-id-max-allele-len 5000 missing \
        --make-bed \
        --out {{wildcards.file}}
        '''
        
        
print(VCF_IN_FLAG)
if VCF_IN_FLAG:
    ruleorder: init_vcf_bfile > init_bfile_bfile
else:
    ruleorder: init_bfile_bfile > init_vcf_bfile
        
        



rule filter_bfile:
    input:
        bed="{file}.bed",
        bim="{file}.bim",
        fam="{file}.fam"
    output:
        filt_bed="{file}_filt.bed",
        filt_bim="{file}_filt.bim",
        filt_fam="{file}_filt.fam"
    params:
        data="{file}"
    shell:
        f'''{PLINK2_PATH} --bfile {{params.data}} \
        --max-alleles 2 \
        --maf {MAF_FILTER} \
        --keep {IDS_FILE}  \
        --make-bed \
        --out {{params.data}}_filt'''
        

rule haps_legend_map_bfile:
    input:
        filt_bed="{file}_filt.bed",
        filt_bim="{file}_filt.bim",
        filt_fam="{file}_filt.fam"
    output:
        filt_haps=temp("{file}_filt.haps"),
        filt_leg=temp("{file}_filt.legend"),
        filt_map=temp("{file}_filt.map"),
        filt_ped=temp("{file}_filt.ped")
    params:
        data="{file}_filt"
    shell:
        f"""{PLINK2_PATH} --bfile {{params.data}} --export ped  --out {{params.data}}
        {PLINK2_PATH}  --bfile {{params.data}}  --export hapslegend  --out {{params.data}}"""
        
        

rule hapgen2:
    input:
        filt_haps="{file}_filt.haps",
        filt_legend="{file}_filt.legend",
        filt_map="{file}_filt.map",
        filt_bim="{file}_filt.bim"
    output:
        controls_gen=temp("{file}_filt_sim.controls.gen"),
        controls_haps=temp("{file}_filt_sim.controls.haps"),
        controls_sample=temp("{file}_filt_sim.controls.sample"),
        filt_sim_legend=temp("{file}_filt_sim.legend"),
    params:
        data="{file}_filt"
    shell:
        f"""
        {HAPGEN2_PATH} \
        -h {{input.filt_haps}} \
        -l {{input.filt_legend}} \
        -m {{input.filt_map}} \
        -o {{params.data}}_sim \
        -dl $(head -1 {{input.filt_bim}} | awk '{{print $4}}') 1 1.5 2.25 \
        -int 0 500000000 \
        -n {N} 0 \
        -Ne 11418 \
        -theta 1
        """
        


rule postprocess_hapgen2:
    input:
        controls_gen="{file}_filt_sim.controls.gen",
        controls_haps="{file}_filt_sim.controls.haps",
        controls_sample="{file}_filt_sim.controls.sample",
    output:
        filt_sim_gen=temp("{file}_filt_sim.gen"),
        filt_sim_sample=temp("{file}_filt_sim.sample"),
        filt_sim_haps=temp("{file}_filt_sim.haps"),
    params:
        data="{file}_filt_sim"
    shell:
        f"""
        rm {{params.data}}.cases.*
        rename 's,sim\.controls,sim,' {{params.data}}.controls.*
        """        


rule change_snp_hapgen:
    input:
        filt_sim_gen="{file}_filt_sim.gen",
        filt_sim_sample="{file}_filt_sim.sample"
    output:
        filt_sim_gen_=temp("{file}_filt_sim_.gen"),
        filt_sim_sample_=temp("{file}_filt_sim_.sample")
    params:
        get_chrom = get_chrom
    shell:
        f"""
        sed s/"snp_"/"{{params.get_chrom}} snp{{params.get_chrom}}_"/g {{input.filt_sim_gen}} > {{output.filt_sim_gen_}}
        mv {{input.filt_sim_sample}} {{output.filt_sim_sample_}}
        """


rule make_bed_from_hapgen2:
    input:
        filt_sim_gen_="{file}_filt_sim_.gen",
        filt_sim_sample_="{file}_filt_sim_.sample"
    output:
        filt_sim_bed=temp("{file}_filt_sim.bed"),
        filt_sim_bim=temp("{file}_filt_sim.bim"),
        filt_sim_fam=temp("{file}_filt_sim.fam")
    params:
        data_="{file}_filt_sim_",
        data="{file}_filt_sim"
    shell:
        f"""
         {PLINK2_PATH} \
         --data {{params.data_}} ref-first \
         --make-bed \
         --out {{params.data}}
        """
        
        
rule merge_chroms:
    input:
        chrom_list=os.path.join(DATA_DIR, f"{VCFS_LIST_PATTERN}_filt_sim.list"),
        input_bed=expand(os.path.join(DATA_DIR, "{file}_filt_sim.bed"), file=files),
        input_bim=expand(os.path.join(DATA_DIR, "{file}_filt_sim.bim"), file=files),
        input_fam=expand(os.path.join(DATA_DIR, "{file}_filt_sim.fam"), file=files)
    output:
        filt_sim_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        filt_sim_bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        filt_sim_fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam")
    params:
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim")
    shell:
        f"""
        {PLINK2_PATH} \
        --pmerge-list {{input.chrom_list}} bfile \
        --set-all-var-ids @:#\$r:\$a --new-id-max-allele-len 5000 missing \
        --make-bed \
        --out {{params.data}}
        """    
        
        
rule recode_merged:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
    output:
        vcf=temp(os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf"))
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
        
        
        
rule transform_gtf:
    input:
        gtf=GTF_IN
    output:
        filt_gtf=temp(os.path.join(DATA_DIR,f"{get_file_name(GTF_IN)}_filt_sort.gtf"))
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/script_make_gtf.sh')} {{input.gtf}} {{output.filt_gtf}}
        """
        
        
        
rule annotate_vcf:
    input:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf"),
        gtf=os.path.join(DATA_DIR,f"{get_file_name(GTF_IN)}_filt_sort.gtf")
    output:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno.vcf"),
        included_txt = os.path.join(DATA_DIR, f'{PATTERN}_{CAUSAL_ID}_included_genes_snps.txt'),
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/bedtools_closest.sh')} {{input.vcf}} {{input.gtf}} {{output.vcf}} {{output.included_txt}} {BEDTOOLS_PATH}
        """        
        
            

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
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/gwas_analysis.sh')} \
        {{params.bed}} \
        {{params.pheno}} \
        {{params.gwas}} \
        {PLINK_PATH}
        """   
             
             
rule pca:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
    output:
        pca_val=temp(os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.eigenval")),
        pca_val_indep=temp(os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep.eigenval")),
        pca_vec=temp(os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.eigenvec")),
        pca_vec_indep=temp(os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep.eigenvec"))
    params:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim"),
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/calc_indep_snps_and_pca.sh')} \
        {{params.bed}} \
        {PLINK2_PATH}
        """           
     
        
             
rule draw_gwas:
    input:
        qassoc=os.path.join(DATA_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas.qassoc"),
        bed_file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed")
    output:
        mh=os.path.join(IMAGES_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas_mh.pdf"),
        qq=os.path.join(IMAGES_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas_qq.pdf")
    params:
        name=f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas",
        bfile=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim")
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/draw_pvals.R')} \
        {{params.name}} \
        {{input.qassoc}} \
        TRUE \
        {{params.bfile}} \
        {PLINK_PATH}
   
        mv QQplot.pval_{{params.name}}.pdf {{output.qq}}
        mv Rectangular-Manhattan.pval_{{params.name}}.pdf {{output.mh}}
        """       
             
rule draw_pca:
    input:
        pca_val=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.eigenval"),
        pca_val_indep=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep.eigenval"),
        pca_vec=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.eigenvec"),
        pca_vec_indep=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep.eigenvec")
    output:
        output_pca=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_pca.pdf"),
        output_pca_indep=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_indep_pca.pdf")
    params:
        file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim"),
        file_indep=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep"),
        pdf=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim"),
        pdf_indep=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_indep")
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/draw_pca.py')} \
        {{params.file}} \
        {{params.pdf}}
        
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/draw_pca.py')} \
        {{params.file_indep}} \
        {{params.pdf_indep}}
        """   
        