import yaml
import os
from pathlib import Path
from distutils.util import strtobool



configfile: "config.yaml"
    
DATA_DIR = config['data_dir']
IMAGES_DIR = config['images_dir']
VCFS_LIST = config['vcfs_list']
IDS_FILE = os.path.join(DATA_DIR, config['ids_file'])
GTF_IN = config['anno_file'].replace('.gtf', '')
GMT_IN = config['gmt_data'].replace('.gmt', '')
PATH_FILE = config['pathways']

MAF_FILTER = config['maf_filter']
N = config['N']

PHENOS_ID = config['phenos_id']
SIM_ID = config['sim_id']

K = config['K']
k = config['k']
GEN_VAR = config['gen_var']
H2S = config['h2s']
SHARED = config['shared']
SD_BETA = config['sd_beta']
M_BETA = config['m_beta']
DRAW_FLAG = bool(strtobool(str(config['draw_flag'])))

print(DRAW_FLAG, GEN_VAR, SIM_ID, config['sim_id'], config['draw_flag'], DRAW_FLAG)


with open(os.path.join(DATA_DIR, VCFS_LIST)) as f:
    file_chrom = [line.strip().split(',') for line in f]
    file2chrom = {k:v for k, v in file_chrom}
    files = list(file2chrom.keys())
    
def get_chrom(wildcards):
    return file2chrom[os.path.basename(wildcards.file)]
    
PATTERN = f"{VCFS_LIST}".replace('.list', '')
    
with open(os.path.join(DATA_DIR, f"{PATTERN}_filt_sim.list"), 'w') as f:
    f.write('\n'.join(map(lambda x: f"{os.path.join(DATA_DIR, x)}_filt_sim", files)))
    
    
draw_output = list()

if DRAW_FLAG:
    output_pca=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_pca.pdf")
    output_pca_indep=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_indep_pca.pdf")
    mh=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas_mh.pdf")
    qq=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas_qq.pdf")
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
        included_txt = os.path.join(DATA_DIR, f'{PATTERN}_{PHENOS_ID}_included_genes_snps.txt'),
        phenos_tsv=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_phenos.tsv"),
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas.tsv"),
        draw_output = draw_output
        


rule filter_vcf:
    input:
        vcf="{file}.vcf"
    output:
        vcf=temp("{file}_filt.vcf")
    priority: 1
    shell:
        f'''plink2 --vcf {{input.vcf}} \
        --max-alleles 2 \
        --maf {MAF_FILTER} \
        --keep {IDS_FILE}  \
        --recode vcf \
        --out {{wildcards.file}}_filt'''
        

rule haps_legend_map_bfile:
    input:
        filt_vcf="{file}_filt.vcf"
    output:
        filt_haps=temp("{file}_filt.haps"),
        filt_leg=temp("{file}_filt.legend"),
        filt_map=temp("{file}_filt.map"),
        filt_ped=temp("{file}_filt.ped"),
        filt_bed="{file}_filt.bed",
        filt_bim="{file}_filt.bim",
        filt_fam="{file}_filt.fam"
    priority: 3
    params:
        data="{file}_filt"
    shell:
        f"""plink2 --vcf {{input.filt_vcf}} --export ped  --out {{params.data}}
        plink2  --vcf {{input.filt_vcf}}  --export hapslegend  --out {{params.data}}
        plink2  --vcf {{input.filt_vcf}}  --set-all-var-ids @:#  --make-bed  --out {{params.data}}"""
        
        

rule hapgen2:
    input:
        filt_haps="{file}_filt.haps",
        filt_legend="{file}_filt.legend",
        filt_map="{file}_filt.map",
        filt_vcf="{file}_filt.vcf",
    output:
        controls_gen="{file}_filt_sim.controls.gen"
    params:
        data="{file}_filt"
    priority: 4
    shell:
        f"""
        hapgen2 \
        -h {{input.filt_haps}} \
        -l {{input.filt_legend}} \
        -m {{input.filt_map}} \
        -o {{params.data}}_sim \
        -dl $(bcftools query -f "%POS\n" {{input.filt_vcf}} | head -n 1) 1 1.5 2.25 \
        -int 0 500000000 \
        -n {N} 0 \
        -Ne 11418 \
        -theta 1"""
        


rule postprocess_hapgen2:
    input:
        controls_gen="{file}_filt_sim.controls.gen"
    output:
        filt_sim_gen=temp("{file}_filt_sim.gen"),
        filt_sim_sample=temp("{file}_filt_sim.sample"),
        filt_sim_haps=temp("{file}_filt_sim.haps"),
        filt_sim_legend=temp("{file}_filt_sim.legend"),
    priority: 5
    params:
        data="{file}_filt_sim"
    shell:
        f"""
        rm {{params.data}}.cases.*
        rename 's,sim\.controls,sim,' {{params.data}}.controls.*
        """        


rule change_snp_hapgen:
    input:
        filt_sim_gen="{file}_filt_sim.gen"
    output:
        filt_sim_gen_="{file}_filt_sim_.gen"
    priority: 6
    params:
        get_chrom = get_chrom
    run:         
        shell(f"""sed s/"snp_"/"{{params.get_chrom}} snp{{params.get_chrom}}_"/g {{input.filt_sim_gen}} > {{output.filt_sim_gen_}}""")


rule make_bed_from_hapgen2:
    input:
        filt_sim_gen_="{file}_filt_sim_.gen",
        filt_sim_gen="{file}_filt_sim.gen",
        filt_sim_sample="{file}_filt_sim.sample"
    output:
        filt_sim_bed="{file}_filt_sim.bed"
    params:
        data="{file}_filt_sim"
    priority: 7
    shell:
        f"""
        rm {{input.filt_sim_gen}}
        mv {{input.filt_sim_gen_}} {{input.filt_sim_gen}}
         plink2 \
         --data {{params.data}} ref-first \
         --make-bed \
         --out {{params.data}}
        """
        
        
rule merge_chroms:
    input:
        input=expand(os.path.join(DATA_DIR, "{file}_filt_sim.bed"), file=files)
    output:
        filt_sim_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        filt_sim_bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        filt_sim_fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam")
    priority: 8
    params:
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim")
    shell:
        f"""
        plink2 \
        --pmerge-list {{params.data}}.list bfile \
        --set-all-var-ids @:# \
        --make-bed \
        --out {{params.data}}
        """    
        
        
rule recode_merged:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
    output:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf")
    priority: 9
    params:
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim")
    shell:
        f"""
        plink2 --bfile {{params.data}} \
        --recode vcf \
        --out {{params.data}}
        """
        
        
        
rule transform_gtf:
    input:
        gtf=os.path.join(DATA_DIR,f"{GTF_IN}.gtf")
    output:
        filt_gtf=temp(os.path.join(DATA_DIR,f"{GTF_IN}_filt_sort.gtf"))
    priority: 10
    params:
        data_dir=DATA_DIR,
        gtf_in=GTF_IN,
    shell:
        f"""
        ./3_2_script_make_gtf.sh {{params.data_dir}} {{params.gtf_in}}
        """
        
        
        
rule annotate_vcf:
    input:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf"),
        gtf=os.path.join(DATA_DIR,f"{GTF_IN}_filt_sort.gtf")
    output:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno.vcf"),
        included_txt = os.path.join(DATA_DIR, f'{PATTERN}_{PHENOS_ID}_included_genes_snps.txt'),
    priority: 11
    shell:
        f"""
        ./3_3_bedtools_closest.sh {{input.vcf}} {{input.gtf}} {{output.vcf}} {{output.included_txt}}
        """        
        
            
        
rule get_snps_list:
    input:
        vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno.vcf"),
        gmt=os.path.join(DATA_DIR,f"{GMT_IN}.gmt"),
        path=os.path.join(DATA_DIR,f"{PATH_FILE}"),
        included_txt = os.path.join(DATA_DIR, f'{PATTERN}_{PHENOS_ID}_included_genes_snps.txt'),
    output:
        causal_genes = temp(os.path.join(DATA_DIR, f'{PATTERN}_{PHENOS_ID}_causal_geneset_snps.txt')),
        annotated_snps = temp(os.path.join(DATA_DIR, f"{PATTERN}_{PHENOS_ID}_annotated_snps.tsv")),
        tsv=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_chosen_snps.tsv")
    params:
        data_dir=DATA_DIR,
        pattern=PATTERN,
        pheno_pattern=PHENOS_ID,
        K=K,
        k=k
    priority: 12
    shell:
        f"""
        ./4_1_get_snps_set.py {{params.data_dir}} {{params.pattern}} {{params.pheno_pattern}} {{input.vcf}} {{input.gmt}} {{input.path}} {{params.K}} {{params.k}}
        """     



rule extract_snps:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
        tsv=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_chosen_snps.tsv"),
    output:
        bed=temp(os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_filt_sim_snps.bed")),
        bim=temp(os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_filt_sim_snps.bim")),
        fam=temp(os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_filt_sim_snps.fam"))
    priority: 13
    params:
        filt_sim_data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim"),
        out_data=os.path.join(DATA_DIR, f"{PATTERN}_{PHENOS_ID}_filt_sim_snps"),
    shell:
        f"""
        plink2 \
        --bfile {{params.filt_sim_data}} \
        --extract range {{input.tsv}} \
        --make-bed --out {{params.out_data}}
        """   
        


rule pheno_sim:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_filt_sim_snps.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_filt_sim_snps.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_filt_sim_snps.fam")
    output:
        tsv=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_phenos.tsv")
    params:
        data_dir=DATA_DIR,
        K=K,
        N=N,
        data=os.path.join(DATA_DIR, f"{PATTERN}_{PHENOS_ID}_filt_sim_snps"),
        GEN_VAR=GEN_VAR,
        H2S=H2S,
        SHARED=SHARED,
        M_BETA=M_BETA,
        SD_BETA=SD_BETA,
    priority: 14
    shell:
        f"""
        ./6_1_pheno_sim.R \
        ./ \
        {{params.data}} \
        {{output.tsv}} \
        {{params.N}} \
        {{params.K}} \
        {{GEN_VAR}} \
        {{H2S}} \
        {{SHARED}}\
        {{M_BETA}} \
        {{SD_BETA}}
        """   
        
# todo make custom gwas        
rule gwas:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
        pheno=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_phenos.tsv"),
    output:
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas.tsv")
    priority: 15
    params:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim"),
        pheno=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_phenos"),
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas")
    shell:
        f"""
        ./7_1_gwas_analysis.sh \
        {{params.bed}} \
        {{params.pheno}} \
        {{params.gwas}}
        """   
             
             
rule pca:
    input:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam"),
    output:
        pca_val=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.eigenval"),
        pca_val_indep=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep.eigenval"),
        pca_vec=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.eigenvec"),
        pca_vec_indep=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep.eigenvec")
    priority: 16
    params:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim"),
    shell:
        f"""
        ./7_2_calc_indep_snps_and_pca.sh \
        {{params.bed}}
        """           
     
        
             
rule draw_gwas:
    input:
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas.tsv")
    output:
        mh=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas_mh.pdf"),
        qq=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas_qq.pdf")
    priority: 17
    params:
        name=f"{PATTERN}_{PHENOS_ID}_{SIM_ID}_gwas"
    shell:
        f"""
        ./8_1_draw_pvals.R \
        {{params.name}} \
        {{input.gwas}} 
   
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
    priority: 18
    params:
        file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim"),
        file_indep=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_indep"),
        pdf=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim"),
        pdf_indep=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_indep")
    shell:
        f"""
        ./8_2_draw_pca.py \
        {{params.file}} \
        {{params.pdf}}
        
        ./8_2_draw_pca.py \
        {{params.file_indep}} \
        {{params.pdf_indep}}
        """   
        
        
        
        
        
        
        
        
        
        