import yaml
import os
from pathlib import Path
import pandas as pd



# ===================================
# READ FILES
# ===================================

files_basenames = [get_basename(f) for f in files]
module_sim_geno_final_outputs = [
    os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
    os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
    os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam")
]

rule module_sim_geno_all:
    input:
        inputa = module_sim_geno_final_outputs
        
        

# ===================================
# RUN SIMULATION 
# ===================================

rule prefilter:
    input:
        vcf=lambda wildcards: next(f"{f}.vcf" for f in files if get_basename(f) == wildcards.file)
    output:
        vcf=temp("{file}_prefilt.vcf"),    
        vmiss=temp("{file}_prefilt.vmiss"),    
        smiss=temp("{file}_prefilt.smiss")    
    message:
        """
        Description: Starting preparation of {input.vcf}: filter samples. It also calculates variant and sample rate for filtering. 
        I/O info:    Result will be available in {output.vcf}.
        Errors:      If the error occurs, possible reasons are: problems with MAF or sample filters. Check PLINK requirements if this happened. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        {PLINK2_PATH} --vcf {{input.vcf}} \
        --keep-allele-order \
        --max-alleles 2 \
        --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 5000 missing \
        {f"--keep {IDS_FILE} " if IDS_FILE is not None else ""}\
        --recode vcf  \
        --missing \
        --out {{wildcards.file}}_prefilt
        """       
        

rule get_filtered_variants_samples:
    input:
        vmiss="{file}_prefilt.vmiss",   
        smiss="{file}_prefilt.smiss" 
    output:  
        vmiss=temp("{file}_filt.vmiss"),    
        smiss=temp("{file}_filt.smiss")    
    message:
        """
        Description: Get samples and variants to filter due to call rate. 
        I/O info:    Result will be available in {output.vmiss} and  {output.smiss}.
        Errors: Check PLINK2 requirements if this happened. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        awk 'NR>1 && $4 > 0' \
        {{input.smiss}} \
        > {{output.smiss}}
        
        awk 'NR>1 && $5 > 0' \
        {{input.vmiss}} \
        > {{output.vmiss}}
        """        
     

rule postfilter:
    input:
        vcf="{file}_prefilt.vcf",
        vmiss="{file}_filt.vmiss",    
        smiss="{file}_filt.smiss"  
    output:
        vcf=temp("{file}_filt.vcf") 
    params:
        bfile="{file}"   
    message:
        """
        Description: Starting preparation of {params.bfile}: filter by call rate. 
        I/O info:    Result will be available in {output.vcf}.
        Errors:      If the error occurs, possible reasons are: problems with MAF or sample filters. Check PLINK requirements if this happened. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        {PLINK2_PATH} --vcf {{input.vcf}} \
        --remove {{input.smiss}} \
        --exclude {{input.vmiss}} \
        --maf {MAF_FILTER} \
        --recode vcf  \
        --out {{wildcards.file}}_filt
        """
        
        
        
rule haps_legend_map_bfile:
    input:
        filt_bed="{file}_filt.vcf"
    output:
        filt_haps=temp("{file}_filt.haps"),
        filt_leg=temp("{file}_filt.legend"),
        filt_map=temp("{file}_filt.map"),
        filt_ped=temp("{file}_filt.ped"),
        filt_bed=temp("{file}_filt.bed"),
        filt_bim=temp("{file}_filt.bim"),
        filt_fam=temp("{file}_filt.fam")
    params:
        data="{file}_filt"
    message: 
        """
        Description: Generating haps, legend, map, ped files (for hapgen2) & bed, bim, fam (for user to check) for {wildcards.file}. 
        I/O info:    Input file: {input.filt_bed}, output files start with: {params.data}.*
        Errors:      If the error occurs, this can be plink2 error. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        {PLINK2_PATH} --vcf {{input.filt_bed}} --export ped  --out {{params.data}}
        {PLINK2_PATH}  --vcf {{input.filt_bed}}  --export hapslegend  --out {{params.data}}
        {PLINK2_PATH}  --vcf {{input.filt_bed}}  --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 5000 missing  --make-bed --out {{params.data}}
        """


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
    message: 
        """
        Description: Executing hapgen2 for {params.data}.* to simulate new dataset.
        Errors:      If the error occurs, this can be an inner error of Hapgen2. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        {HAPGEN2_PATH} \
        -h {{input.filt_haps}} \
        -l {{input.filt_legend}} \
        -m {{input.filt_map}} \
        -o {{params.data}}_sim \
        -dl $(head -1 {{input.filt_bim}} | cut -f4) 1 1.5 2.25 \
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
    message: 
        """
        Description: Post-processing hapgen2 results for {params.data}.
        Errors:      If the error occurs check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        rm -f {{params.data}}.cases.*
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
    message: 
        """
        Description: Changing SNP format of hapgen output for {input.filt_sim_gen}.
        Errors:      If the error occurs, check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
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
    message: 
        """
        Description: Creating plinks's binary files from hapgen2 simulation of genotypes for {params.data}.
        Errors:      If the error occurs, that is probably plink2 error. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
         {PLINK2_PATH} \
         --data {{params.data_}} ref-first \
         --make-bed \
         --out {{params.data}}
        """
        
        
rule merge_chroms:
    input:
        chrom_list=VCFS_LIST_SIMULATED,
        input_bed=expand("{file}_filt_sim.bed", file=new_files_pat),
        input_bim=expand("{file}_filt_sim.bim", file=new_files_pat),
        input_fam=expand("{file}_filt_sim.fam", file=new_files_pat)
    output:
        filt_sim_bed=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
        filt_sim_bim=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
        filt_sim_fam=os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam")
    params:
        data=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim")
    message: 
        """
        Merging plinks's binary files by chromosomes to plinks's binary file {params.data}.
        If the error occurs, that is probably plink2 error. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """ 
    shell:
        f"""
        if [ "{TO_MERGE}" = "True" ]; then
            {PLINK2_PATH} \
            --pmerge-list {{input.chrom_list}} bfile \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 5000 missing \
            --make-bed \
            --out {{params.data}}
        else
            cp {{input.input_bed}} {{output.filt_sim_bed}}
            cp {{input.input_bim}} {{output.filt_sim_bim}}
            cp {{input.input_fam}} {{output.filt_sim_fam}}
        fi
        """
        
        
        
        
        
        
        
        
        
        
        
        