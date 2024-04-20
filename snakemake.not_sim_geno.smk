import yaml
import os
from pathlib import Path
import pandas as pd



# ===================================
# READ FILES
# ===================================

files_basenames = [remove_ext(f) for f in files]
module_not_sim_geno_final_outputs = [
    os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bed"),
    os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.bim"),
    os.path.join(DATA_DIR,f"{PATTERN}_filt_sim.fam")
]

rule module_not_sim_geno_all:
    input:
        inputb = module_not_sim_geno_final_outputs
        

# ===================================
# SKIP SIMULATION -- but prepare genos 
# ===================================

rule prefilter:
    input:
        vcf=lambda wildcards: next(
        f+'.vcf' if not BFILE_IN_FLAG else f+'.bed' 
        for f in files if (
        (get_basename(f) == wildcards.file and not BFILE_IN_FLAG) 
        or 
        (remove_ext(get_basename(f)) == wildcards.file and BFILE_IN_FLAG)
        ))
    output:
        bed=temp("{file}_prefilt.bed"),
        bim=temp("{file}_prefilt.bim"),
        fam=temp("{file}_prefilt.fam"),
        vmiss=temp("{file}_prefilt.vmiss"),    
        smiss=temp("{file}_prefilt.smiss")   
    params:
        vcf_bfile = lambda wildcards: next(
            f + ".vcf" if not BFILE_IN_FLAG else f
            for f in files
            if (
                (get_basename(f) == wildcards.file and not BFILE_IN_FLAG)
                or (remove_ext(get_basename(f)) == wildcards.file and BFILE_IN_FLAG)
            )
        )  # Adjust for plink input
    message:
        """
        Description: Starting preparation of {input.vcf}: filter alleles. It also calculates variant and sample rate for filtering. 
        I/O info:    Result will be available in {output} PLINK's bfile.
        Errors:      If the error occurs, possible reasons are: problems with MAF or sample filters. Check PLINK requirements if this happened. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        {PLINK2_PATH} {'--vcf'  if not BFILE_IN_FLAG else '--bfile'} {{params.vcf_bfile}} \
        --keep-allele-order \
        --max-alleles 2 \
        --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 5000 missing \
        {f"--keep {IDS_FILE} " if IDS_FILE is not None else ""}\
        --make-bed  \
        --missing \
        --out {{wildcards.file}}_prefilt
        """
        
rule get_filtered_variants_samples:
    input:
        vmiss="{file}_prefilt.vmiss",   
        smiss="{file}_prefilt.smiss" 
    output:  
        vmiss=temp("{file}_filt_sim.vmiss"),    
        smiss=temp("{file}_filt_sim.smiss")      
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
        bed="{file}_prefilt.bed",
        bim="{file}_prefilt.bim",
        fam="{file}_prefilt.fam",
        vmiss="{file}_filt_sim.vmiss",    
        smiss="{file}_filt_sim.smiss"    
    output:
        bed=temp("{file}_filt_sim.bed"),
        bim=temp("{file}_filt_sim.bim"),
        fam=temp("{file}_filt_sim.fam")  
    params:
        bfile="{file}"
    message:
        """
        Description: Starting preparation of {params.bfile}: filter by call rate. 
        I/O info:    Result will be available in {output} PLINK's bfile.
        Errors:      If the error occurs, possible reasons are: problems with MAF or sample filters. Check PLINK requirements if this happened. Check logs and in/out files, their formats and try again. Logs are in: {log}.
        """
    shell:
        f"""
        {PLINK2_PATH} --bfile {{params.bfile}} \
        --remove {{input.smiss}} \
        --exclude {{input.vmiss}} \
        --maf {MAF_FILTER} \
        --make-bed  \
        --out {{wildcards.file}}_filt_sim
        """        
        
        
        


rule merge_chroms:
    input:
        chrom_list=VCFS_LIST_SIMULATED,
        input_bed=expand("{file}_filt_sim.bed", file=files_basenames),
        input_bim=expand("{file}_filt_sim.bim", file=files_basenames),
        input_fam=expand("{file}_filt_sim.fam", file=files_basenames)
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
        
        
        
        
        
        
        
        
        
        
        
        