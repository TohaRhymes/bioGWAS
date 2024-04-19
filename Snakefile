# Include module Snakefiles
# and compile final outputs
final_outputs = []

include: "snakemake.common.smk"

# STEP 1: sim genotypes/skip it
if not SKIP_SIMULATION:
    include: "snakemake.sim_geno.smk"
    final_outputs+=module_sim_geno_final_outputs
else:
    include: "snakemake.not_sim_geno.smk"
    final_outputs+=module_not_sim_geno_final_outputs
    
# STEP 2: sim phenos
include: "snakemake.sim_pheno.smk"
final_outputs+=module_sim_pheno_final_outputs

# STEP 3: sim gwas
include: "snakemake.gwas.smk"
final_outputs+=module_gwas_final_outputs

# STEP 4: draw
include: "snakemake.draw.smk"
final_outputs+=module_draw_final_outputs

# Rule all to define the final outputs
rule all:
    input:
        final_outputs
        
onsuccess:
    print("bioGWAS finished, no errors.")

onerror:
    print("An error occurred. Logs: {log}")
    shell("Create an issue on https://github.com/TohaRhymes/bioGWAS/")