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

output_pattern = f"{DATA_DIR}/{{file}}_filt_sim.gen"
input_pattern = f"{DATA_DIR}/{{file}}.vcf"

with open(f"{DATA_DIR}/{VCFS_LIST}") as f:
    files = [line.strip() for line in f]
    

rule all:
    priority: 1000
    input:
        expand(output_pattern, file=files)

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
        rm {{wildcards.file}}_filt_sim.controls.*
        rename 's,sim\.controls,sim,' {{wildcards.file}}_filt_sim.controls.*
        """

        
        