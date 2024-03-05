# Include module Snakefiles
# and compile final outputs
final_outputs = []

include: "snakemake.common.smk"

if not SKIP_SIMULATION:
    include: "snakemake.a.smk"
    final_outputs+=module_a_final_outputs
else:
    include: "snakemake.b.smk"
    final_outputs+=module_b_final_outputs
# include: "moduleC.snakefile"
# include: "moduleD.snakefile"


# Rule all to define the final outputs
rule all:
    input:
        final_outputs