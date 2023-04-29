import yaml
import os
from pathlib import Path

with open('config.yaml') as cf_file:
    config = yaml.safe_load( cf_file.read() )
    
DATA_DIR = config['data_dir']
VCFS_LIST = config['vcfs_list']
MAF_FILTER = config['maf_filter']
N = config['N']
IDS_FILE = os.path.join(DATA_DIR, config['ids_file'])

GTF_IN = config['anno_file'].replace('.gtf', '')
GMT_IN = config['gmt_data'].replace('.gmt', '')
PATH_FILE = config['pathways']
K = config['K']
k = config['k']

splitted_pattern = os.path.join(DATA_DIR, "{file}_filt_sim.bed")

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
        input_file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_snps.bed")


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
        output_file="{file}_filt.haps"
    priority: 3
    shell:
        f"""plink2 --vcf {{input.input_file}} --export ped  --out {{wildcards.file}}_filt
        plink2  --vcf {{input.input_file}}  --export hapslegend  --out {{wildcards.file}}_filt"""


rule hapgen2:
    input:
        input_file="{file}_filt.haps"
    output:
        output_file="{file}_filt_sim.controls.gen"
    priority: 4
    shell:
        f"""
        hapgen2 \
        -h {{wildcards.file}}_filt.haps \
        -l {{wildcards.file}}_filt.legend \
        -m {{wildcards.file}}_filt.map \
        -o {{wildcards.file}}_filt_sim \
        -dl $(bcftools query -f "%POS\n" {{wildcards.file}}_filt.vcf | head -n 1) 1 1.5 2.25 \
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
    priority: 7
    shell:
        f"""
        rm {{input.input_file}}
        mv {{input.input_file_}} {{input.input_file}}
         plink2 \
         --data {{wildcards.file}}_filt_sim ref-first \
         --make-bed \
         --out {{wildcards.file}}_filt_sim
        """
        
        
rule merge_chroms:
    input:
        input_file=expand(splitted_pattern, file=files)
    output:
        output_file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed")
    priority: 8
    shell:
        f"""
        plink2 \
        --pmerge-list {os.path.join(DATA_DIR, PATTERN)}_filt_sim.list bfile \
        --set-all-var-ids @:# \
        --make-bed \
        --out {os.path.join(DATA_DIR, PATTERN)}_filt_sim
        """        
        
rule recode_merged:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed")
    output:
        output_vcf=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.vcf"),
        output_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_ucsc.bed"),
    priority: 9
    shell:
        f"""
        plink2 --bfile {os.path.join(DATA_DIR, PATTERN)}_filt_sim \
        --recode vcf \
        --out {os.path.join(DATA_DIR, PATTERN)}_filt_sim
        
        # vcf to ucsc's bed (vcf2bed sorts it by default)
        vcf2bed < {{output.output_vcf}} > {{output.output_bed}}
        
        # some of rows sont include all samples - filter 'em
        awk -v N={N+10} 'NF==N' {{output.output_bed}} > {os.path.join(DATA_DIR, "")}filtered_{PATTERN}_filt_sim_ucsc.bed 
        
        # rename files
        rm {{output.output_bed}}
        mv {os.path.join(DATA_DIR, "")}filtered_{PATTERN}_filt_sim_ucsc.bed {{output.output_bed}}
        """
        
        
rule transform_gtf:
    input:
        input_file=os.path.join(DATA_DIR,f"{GTF_IN}.gtf")
    output:
        output_file=os.path.join(DATA_DIR,f"{GTF_IN}_filt_sort.gtf")
    priority: 10
    shell:
        f"""
        ./5_script_make_gtf.sh {DATA_DIR} {GTF_IN}
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
        ./get_snps_set.py {{params.data_dir}} {{params.pattern}} {{input.input_bed}} {{input.input_gmt}} {{input.input_path}} {{params.K}} {{params.k}}
        """     


rule extract_snps:
    input:
        input_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        input_tsv=os.path.join(DATA_DIR,f"{PATTERN}_chosen_snps.tsv"),
    output:
        output_file=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_snps.bed")
    priority: 13
    shell:
        f"""
        plink2 \
        --bfile {os.path.join(DATA_DIR,f"{PATTERN}_filt_sim")} \
        --extract range {{input.input_tsv}} \
        --make-bed --out {os.path.join(DATA_DIR,f"{PATTERN}_filt_sim_snps")}
        """   
        
        
        
        
        
        
        
        
        
        
        
        
        