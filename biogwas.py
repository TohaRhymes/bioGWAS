#!/usr/bin/env python

import argparse
import os
import subprocess
import sys


SNAKE_YAML_TEMPLATE = """
biogwas_path: {biogwas_path}
vcf_in_dir: {vcf_in_dir}
data_dir: {data_dir}
images_dir: {images_dir}

vcf_in_flag: {vcf_in_flag}
vcfs_list: {vcfs_list}
ids_file: {ids_file}
anno_file: {anno_file}
gmt_file: {gmt_file}

use_causal_snps: {use_causal_snps}
causal_pathways: {causal_pathways}
causal_snps: {causal_snps}

maf_filter: {maf_filter}
causal_maf_min: {causal_maf_min}
causal_maf_max: {causal_maf_max}

N: {N}

K: {K}
k: {k}
m_beta: {m_beta}
sd_beta: {sd_beta}
gen_var: {gen_var}
p_independent_genetic: {p_independent_genetic}
theta: {theta}
alpha: {alpha}

pattern: {pattern}
causal_id: {causal_id}
sim_id: {sim_id}
draw_flag: {draw_flag}

seed: {seed}"""


LOGO = """
     _     _        ______        ___    ____  
    | |__ (_) ___  / ___\ \      / / \  / ___| 
    | '_ \| |/ _ \| |  _ \ \ /\ / / _ \ \___ \ 
    | |_) | | (_) | |_| | \ V  V / ___ \ ___) |
    |_.__/|_|\___/ \____|  \_/\_/_/   \_\____/ 

(C) Changalidis et al., 2023
"""


def main(path, args):
    with open(args.dependencies) as f:
        dep_yaml = f.read()

    snake_args = {
        "biogwas_path": os.path.abspath(path),
        "vcf_in_dir": os.path.abspath(args.input_dir),
        "data_dir": os.path.abspath(args.data_dir),
        "images_dir": os.path.abspath(args.img_dir),
        "vcf_in_flag": args.vcf_in_flag,
        "vcfs_list": os.path.abspath(args.input_list),
        "ids_file": os.path.abspath(args.ids_file),
        "anno_file": os.path.abspath(args.anno_file),
        "gmt_file": os.path.abspath(args.gmt_file),
        "use_causal_snps": args.use_causal_snps,
        "causal_pathways": os.path.abspath(args.causal_pathways) if args.causal_pathways is not None else "",
        "causal_snps": os.path.abspath(args.causal_snps) if args.causal_snps is not None else "",
        "maf_filter": args.maf_filter,
        "causal_maf_min": args.causal_maf_min,
        "causal_maf_max": args.causal_maf_max,
        "N": args.N,
        "K": args.K,
        "k": args.k,
        "m_beta": args.m_beta,
        "sd_beta": args.sd_beta,
        "gen_var": args.gen_var,
        "p_independent_genetic": args.p_independent_genetic,
        "theta": args.theta,
        "alpha": args.alpha,
        "pattern": args.pattern,
        "causal_id": args.causal_id,
        "sim_id": args.sim_id,
        "draw_flag": args.draw_flag,
        "seed": args.seed,
    }
    snake_yaml = SNAKE_YAML_TEMPLATE.format(**snake_args)
    config = os.path.abspath(args.config)
    print(config)
    with open(config, "w") as w:
        w.write(dep_yaml + "\n" + snake_yaml)

    command = f'snakemake \
     --nolock \
     -s {os.path.join(path, "Snakefile")} \
     --configfile {config} \
     --cores {args.threads} \
     --directory {os.path.abspath(args.data_dir)}'
    print(command)
    subprocess.call(command, shell=True)

    print(LOGO)


if __name__ == "__main__":
    print(LOGO)
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # SCRIPT SETTINGS
    script_settings = parser.add_argument_group("script settings")
    script_settings.add_argument(
        "-d",
        "--dependencies",
        type=str,
        default="./dependencies.yaml",
        required=False,
        help="Path to the dependencies file with all constant paths for all internal programs.",
    )
    script_settings.add_argument(
        "-cfg",
        "--config",
        type=str,
        default="snakefile_config.yaml",
        required=False,
        help="Path to the future config file.",
    )
    script_settings.add_argument(
        "-T", "--threads", type=int, default=1, required=False, help="N cores."
    )
    # DIRS
    dir_settings = parser.add_argument_group("directories settings")
    dir_settings.add_argument(
        "-id",
        "--input_dir",
        required=True,
        type=str,
        help="Path to the directory which contains input genotypes.",
    )
    dir_settings.add_argument(
        "-dd",
        "--data_dir",
        required=True,
        type=str,
        help="Path to the directory which will be used as a storage for output data.",
    )
    dir_settings.add_argument(
        "-imd",
        "--img_dir",
        required=True,
        type=str,
        help="Path to the directory which will be used as a storage for the images.",
    )
    # INPUT FILES
    input_settings = parser.add_argument_group("input settings")
    input_settings.add_argument(
        "-vcf",
        "--vcf_in_flag",
        action=argparse.BooleanOptionalAction,
        default=False,
        type=bool,
        help="If flag is set, algorithm expect input files be in `vcf` format; othervise - in `plink binary` format (requires `bed`+`bim`+`fam` files).",
    )
    input_settings.add_argument(
        "-il",
        "--input_list",
        required=True,
        type=str,
        help="CSV-file that describes one input file per line with structure: <file>,<chromosome_number>. For plink bfiles, provide just name without extension.",
    )
    input_settings.add_argument(
        "-if",
        "--ids_file",
        required=True,
        type=str,
        help="File with samples ids' (from input files) to use for sampling (one id per line).",
    )
    input_settings.add_argument(
        "-af",
        "--anno_file",
        required=True,
        type=str,
        help="File with chromosome annotations in gtf-format.",
    )
    input_settings.add_argument(
        "-gf",
        "--gmt_file",
        required=False,
        type=str,
        help="File with gene sets (pathways) in gmt-format. Parameter should be used only if use_causal_snps is set to False.",
    )
    # CAUSAL SNPS SETTING
    causal_settings = parser.add_argument_group("causal SNPs settings")
    causal_settings.add_argument(
        "-ucs",
        "--use_causal_snps",
        action=argparse.BooleanOptionalAction,
        default=False,
        type=bool,
        help="Use specific set of SNPs (from the file in `--causal_snps` parameter), or random SNPs from the specific pathways (from the file in `--causal_pathways` parameter and that were annotated using gmt-file from `--gmt_file` parameter).",
    )
    causal_settings.add_argument(
        "-cp",
        "--causal_pathways",
        required=False,
        type=str,
        help="Path to the file with causal pathways/gene sets (causal SNPs will be randomly selected from genes from these pathways). Several sets' names should be separated by comma/new line (both ways are possible). Required when `--use_causal_snps` is set to False.",
    )
    causal_settings.add_argument(
        "-cs",
        "--causal_snps",
        required=False,
        type=str,
        help="Path to the file with preselected SNPs (one SNP per line). Required when `--use_causal_snps` is set to True. The format is: <chrom>:<serial number of basepair>",
    )
    # MAF FILTERS
    maf_settings = parser.add_argument_group("maf settings")
    maf_settings.add_argument(
        "-maf",
        "--maf_filter",
        default=0.05,
        required=False,
        type=float,
        help="Filters out all SNPs with minor allele frequency below the provided threshold. Only SNPs with MAF greater than or equal to provided are used for the subsequent analysis.",
    )
    maf_settings.add_argument(
        "-maf_min",
        "--causal_maf_min",
        default=0.05,
        required=False,
        type=float,
        help="Only SNPs with MAF greater or equal to provided are used as causal SNPs (they influence phenotype). The maximal theshold is set ufing `--causal_maf_max`.",
    )
    maf_settings.add_argument(
        "-maf_max",
        "--causal_maf_max",
        default=0.9999999999,
        required=False,
        type=float,
        help="Only SNPs with MAF less or equal to provided are used as causal SNPs (they influence phenotype). The minimal theshold is set ufing `--causal_maf_max`.",
    )
    # SIMULATED SAMPLES
    samples_settings = parser.add_argument_group("samples settings")
    samples_settings.add_argument(
        "-N", "--N", default=100, type=int, help="Amount of samples to simulate."
    )
    # SIMULATION PHENOTYPES SETTINGS
    sim_pheno_settings = parser.add_argument_group("phenotypes simulation settings")
    sim_pheno_settings.add_argument(
        "-K",
        "--K",
        required=False,
        default=10,
        type=int,
        help="Amount of causal SNPs. In case when SNPs were set (instead of sets), K should be less or equal than amount of given SNPs).",
    )

    sim_pheno_settings.add_argument(
        "-k",
        "--k",
        required=False,
        default=5,
        type=int,
        help="Required when `--use_causal_snps` is set to False. Amount of causal SNPs selected from causal gene sets' genes (other causal SNPs are chosen randomly).",
    )

    sim_pheno_settings.add_argument(
        "-mb",
        "--m_beta",
        required=False,
        default=0.05,
        type=float,
        help="Phenotype simulator parameter: SNP's effect size.",
    )

    sim_pheno_settings.add_argument(
        "-sb",
        "--sd_beta",
        required=False,
        default=0.001,
        type=float,
        help="Phenotype simulator parameter: standard deviation of `--m_beta`",
    )

    sim_pheno_settings.add_argument(
        "-gv",
        "--gen_var",
        required=False,
        default=0.5,
        type=float,
        help="Phenotype simulator parameter: the proportion of variation explained by genetics (the rest is noise).",
    )

    sim_pheno_settings.add_argument(
        "-pig",
        "--p_independent_genetic",
        required=False,
        default=1.0,
        type=float,
        help="Phenotype simulator parameter: proportion of genetic variant effects to have a trait-independent fixed effect",
    )

    sim_pheno_settings.add_argument(
        "-t",
        "--theta",
        required=False,
        default=0.0,
        type=float,
        help="Phenotype simulator parameter: proportion of variance of shared genetic variant effects",
    )

    sim_pheno_settings.add_argument(
        "-a",
        "--alpha",
        required=False,
        default=0.5,
        type=float,
        help="Phenotype simulator parameter: variance of shared observational noise effect.",
    )
    # Filenames patterns
    pattern_settings = parser.add_argument_group("filename patterns settings")
    pattern_settings.add_argument(
        "-p",
        "--pattern",
        required=True,
        type=str,
        help="String, that will be use as a prefix for every output file of this simulation.",
    )
    pattern_settings.add_argument(
        "-cid",
        "--causal_id",
        required=True,
        type=str,
        help="ID which represents each file of this specific SNPs' set (foes after `pattern`).",
    )
    pattern_settings.add_argument(
        "-sid",
        "--sim_id",
        required=True,
        type=str,
        help="ID which represents each output file of this specific phenotype and its association with SNPs (goes after `pattern` and `causal_id`).",
    )
    # Draw
    draw_settings = parser.add_argument_group("draw settings")
    draw_settings.add_argument(
        "-df",
        "--draw_flag",
        action=argparse.BooleanOptionalAction,
        default=False,
        type=bool,
        help="If flag is set, make graphical representation of simulated genotypes and associations (PCA, MH and QQ).",
    )
    # Seed
    seed_settings = parser.add_argument_group("seed settings")
    seed_settings.add_argument(
        "-S",
        "--seed",
        required=False,
        default=566,
        type=int,
        help="Random seed of this phenotype (!!!) simulation.",
    )

    args = parser.parse_args()
    path = os.path.dirname(sys.argv[0])
    main(path, args)