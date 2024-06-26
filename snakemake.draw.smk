import yaml
import os
from pathlib import Path
import pandas as pd



# ===================================
# FILES
# ===================================

module_draw_final_outputs = list()

if DRAW_FLAG:
    output_pca=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_pca.pdf")
    output_pca_indep=os.path.join(IMAGES_DIR,f"{PATTERN}_filt_sim_indep_pca.pdf")
    mh=os.path.join(IMAGES_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas_mh.pdf")
    qq=os.path.join(IMAGES_DIR,f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas_qq.pdf")
    module_draw_final_outputs.append(output_pca)
    module_draw_final_outputs.append(output_pca_indep)
    module_draw_final_outputs.append(mh)
    module_draw_final_outputs.append(qq)



# ===================================
# RULES
# ===================================


rule module_draw_all:
    input:
        inpute = module_draw_final_outputs


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
        bed=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim")
    message: 
        """
        Description: Executing PCA for {params.bed}. 
        I/O info:    Result will be available in *.eigenval and *.eigenvec files.
        Errors:      If the error occurs, it can be Python problem. Check logs and in/out files, their formats and try again.
        """
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
    message: 
        """
        Description: Drawing Manhattan and Q-Q plots for {params.gwas}.
        I/O info:    Plots are in: {output.mh} & {output.qq}.
        oErrors:      Error in drawing Manhattan and Q-Q plots. That can be R problem. Check logs and in/out files, their formats and try again.
        """
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/draw_pvals.R')} \
        {{params.name}} \
        {{input.qassoc}} \
        TRUE \
        {{params.bfile}} \
        {PLINK_PATH} \
        TRUE \
        {QQ_WIDTH} \
        {QQ_HEIGHT} \
        {QQ_DPI} \
        {MH_WIDTH} \
        {MH_HEIGHT} \
        {MH_DPI}
   
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
    message: 
        """
        Description: Creating PCA plots using previously computed PC.
        I/O info:    PCA plots created successfully in pdf formats.
        Errors:      Error in creating PCA plots. That can be Python problem. Check logs and in/out files, their formats and try again.
        """
    shell:
        f"""
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/draw_pca.py')} \
        {{params.file}} \
        {{params.pdf}} \
        {PCA_WIDTH} \
        {PCA_HEIGHT} \
        {PCA_DPI}
        
        {os.path.join(BIOGWAS_PATH, './pipeline_utils/draw_pca.py')} \
        {{params.file_indep}} \
        {{params.pdf_indep}} \
        {PCA_WIDTH} \
        {PCA_HEIGHT} \
        {PCA_DPI}
        """   