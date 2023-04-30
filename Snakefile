import yaml
import os
from pathlib import Path

with open('config.yaml') as cf_file:
    config = yaml.safe_load( cf_file.read() )
    
DATA_DIR = config['data_dir']
IMAGES_DIR = config['images_dir']
VCFS_LIST = config['vcfs_list']
MAF_FILTER = config['maf_filter']
N = config['N']
IDS_FILE = os.path.join(DATA_DIR, config['ids_file'])

GTF_IN = config['anno_file'].replace('.gtf', '')
GMT_IN = config['gmt_data'].replace('.gmt', '')
PATH_FILE = config['pathways']
PHENOS_ID = config['phenos_id']
K = config['K']
k = config['k']
GEN_VAR = config['gen_var']
H2S = config['h2s']
SHARED = config['shared']


splitted_pattern = os.path.join(DATA_DIR, "{file}_filt_sim.bed")
splitted_init_pattern = os.path.join(DATA_DIR, "{file}_filt.bed")

with open(os.path.join(DATA_DIR, VCFS_LIST)) as f:
    file_chrom = [line.strip().split(',') for line in f]
    file2chrom = {k:v for k, v in file_chrom}
    files = list(file2chrom.keys())
    
def get_chrom(wildcards):
    return file2chrom[os.path.basename(wildcards.file)]
    
PATTERN = f"{VCFS_LIST}".replace('.list', '')
    
with open(os.path.join(DATA_DIR, f"{PATTERN}_filt_sim.list"), 'w') as f:
    f.write('\n'.join(map(lambda x: f"{os.path.join(DATA_DIR, x)}_filt_sim", files)))

rule all:
    priority: 1000
    input:
        output_pca=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_pca.pdf"),
        output_pca_indep=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_indep_pca.pdf"),
        mh=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_gwas_mh.pdf"),
        qq=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_gwas_qq.pdf")


rule filter_vcf:
    input:
        input_file="{file}.vcf"
    output:
        output_file="{file}_filt.vcf"
    priority: 1
    shell:
        f'''plink2 --vcf {{input.input_file}} \
        --max-alleles 2 \
        --maf {MAF_FILTER} \
        --keep {IDS_FILE}  \
        --recode vcf \
        --out {{wildcards.file}}_filt'''
        

rule ped:
    input:
        input_file="{file}_filt.vcf"
    output:
        output_file_haps="{file}_filt.haps",
        output_file_leg="{file}_filt.legend",
        output_file_map="{file}_filt.map"
    priority: 3
    params:
        data="{file}_filt"
    shell:
        f"""plink2 --vcf {{input.input_file}} --export ped  --out {{params.data}}
        plink2  --vcf {{input.input_file}}  --export hapslegend  --out {{params.data}}"""
        
        

rule hapgen2:
    input:
        input_file="{file}_filt.haps"
    output:
        output_file="{file}_filt_sim.controls.gen"
    params:
        data="{file}_filt"
    priority: 4
    shell:
        f"""
        hapgen2 \
        -h {{params.data}}.haps \
        -l {{params.data}}.legend \
        -m {{params.data}}.map \
        -o {{params.data}}_sim \
        -dl $(bcftools query -f "%POS\n" {{params.data}}.vcf | head -n 1) 1 1.5 2.25 \
        -int 0 500000000 \
        -n {N} 0 \
        -Ne 11418 \
        -theta 1"""
        


rule postprocess_hapgen2:
    input:
        input_file="{file}_filt_sim.controls.gen"
    output:
        output_file="{file}_filt_sim.gen"
    priority: 5
    shell:
        f"""
        rm {{wildcards.file}}_filt_sim.cases.*
        rm {{wildcards.file}}_filt.*
        rename 's,sim\.controls,sim,' {{wildcards.file}}_filt_sim.controls.*
        """        


rule change_snp_hapgen:
    input:
        input_file="{file}_filt_sim.gen"
    output:
        output_file="{file}_filt_sim_.gen"
    priority: 6
    params:
        myparam = get_chrom
    run:         
        shell(f"""sed s/"snp_"/"{{params.myparam}} snp{{params.myparam}}_"/g {{input.input_file}} > {{output.output_file}}""")


rule make_bed_from_hapgen2:
    input:
        input_file_="{file}_filt_sim_.gen",
        input_file="{file}_filt_sim.gen"
    output:
        output_file="{file}_filt_sim.bed"
    params:
        data="{file}_filt_sim"
    priority: 7
    shell:
        f"""
        rm {{input.input_file}}
        mv {{input.input_file_}} {{input.input_file}}
         plink2 \
         --data {{params.data}} ref-first \
         --make-bed \
         --out {{params.data}}
        """
        
        
rule merge_chroms:
    input:
        input_file=expand(splitted_pattern, file=files)
    output:
        output_file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed")
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
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed")
    output:
        output_vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf"),
        output_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_ucsc.bed"),
    priority: 9
    params:
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim"),
        f_data=os.path.join(DATA_DIR, f"filtered_{PATTERN}_filt_sim_ucsc")
    shell:
        f"""
        plink2 --bfile {{params.data}} \
        --recode vcf \
        --out {{params.data}}
        
        # vcf to ucsc's bed (vcf2bed sorts it by default)
        vcf2bed < {{output.output_vcf}} > {{output.output_bed}}
        
        # some of rows sont include all samples - filter 'em
        awk -v N={N+10} 'NF==N' {{output.output_bed}} > {{params.f_data}}.bed 
        
        # rename files
        rm {{output.output_bed}}
        mv {{params.f_data}}.bed {{output.output_bed}}
        """
        
        
rule transform_gtf:
    input:
        input_file=os.path.join(DATA_DIR,f"{GTF_IN}.gtf")
    output:
        output_file=os.path.join(DATA_DIR,f"{GTF_IN}_filt_sort.gtf")
    priority: 10
    shell:
        f"""
        ./3_2_script_make_gtf.sh {DATA_DIR} {GTF_IN}
        """
        
        
        
rule annotate_bed:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_ucsc.bed"),
        input_gtf=os.path.join(DATA_DIR,f"{GTF_IN}_filt_sort.gtf")
    output:
        output_file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno_ucsc.bed")
    priority: 11
    shell:
        f"""
        bedtools closest -d -a {{input.input_bed}} -b {{input.input_gtf}}  > {{output.output_file}}
        """        
        
            
        
rule get_snps_list:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_anno_ucsc.bed"),
        input_gmt=os.path.join(DATA_DIR,f"{GMT_IN}.gmt"),
        input_path=os.path.join(DATA_DIR,f"{PATH_FILE}")
    output:
        output_file=os.path.join(DATA_DIR,f"{PATTERN}_chosen_snps.tsv")
    params:
        data_dir=DATA_DIR,
        pattern=PATTERN,
        K=K,
        k=k
    priority: 12
    shell:
        f"""
        ./4_1_get_snps_set.py {{params.data_dir}} {{params.pattern}} {{input.input_bed}} {{input.input_gmt}} {{input.input_path}} {{params.K}} {{params.k}}
        """     


rule extract_snps:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        input_tsv=os.path.join(DATA_DIR,f"{PATTERN}_chosen_snps.tsv"),
    output:
        output_file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_snps.bed")
    priority: 13
    params:
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim_snps")
    shell:
        f"""
        plink2 \
        --bfile {os.path.join(DATA_DIR,f"{PATTERN}_filt_sim")} \
        --extract range {{input.input_tsv}} \
        --make-bed --out {{params.data}}
        """   
        


rule pheno_sim:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_snps.bed")
    output:
        output_file=os.path.join(DATA_DIR,f"phenos_{PHENOS_ID}.tsv")
    params:
        data_dir=DATA_DIR,
        K=K,
        N=N,
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim_snps")
    priority: 14
    shell:
        f"""
        ./6_1_pheno_sim.R \
        ./ \
        {{params.data}} \
        {{output.output_file}} \
        {{params.N}} \
        {{params.K}} \
        {{GEN_VAR}} \
        {{H2S}} \
        {{SHARED}}
        """   
        
# todo make custom gwas        
rule gwas:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        input_pheno=os.path.join(DATA_DIR,f"phenos_{PHENOS_ID}.tsv"),
    output:
        output_gwas=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_gwas.tsv")
    priority: 15
    params:
        bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim"),
        pheno=os.path.join(DATA_DIR,f"phenos_{PHENOS_ID}"),
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_gwas")
    shell:
        f"""
        ./7_1_gwas_analysis.sh \
        {{params.bed}} \
        {{params.pheno}} \
        {{params.gwas}}
        """   
             
             
rule pca:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed")
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
        gwas=os.path.join(DATA_DIR,f"{PATTERN}_{PHENOS_ID}_gwas.tsv")
    output:
        mh=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_gwas_mh.pdf"),
        qq=os.path.join(IMAGES_DIR,f"{PATTERN}_{PHENOS_ID}_gwas_qq.pdf")
    priority: 17
    params:
        name=f"{PATTERN}_{PHENOS_ID}_gwas"
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
        
        
        
        
        
        
        
        
        
        