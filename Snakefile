# Include module Snakefiles
# and compile final outputs
final_outputs = []

include: "snakemake.common.smk"

if not SKIP_SIMULATION:
    include: "snakemake.sim_geno.smk"
    final_outputs+=module_sim_geno_final_outputs
else:
    include: "snakemake.not_sim_geno.smk"
    final_outputs+=module_not_sim_geno_final_outputs
# include: "moduleC.snakefile"
# include: "moduleD.snakefile"


# Rule all to define the final outputs
rule all:
    input:
        final_outputs