import yaml
import os
from pathlib import Path
import pandas as pd

# ===================================
# UTIL FUNCTIONS
# ===================================

def get_chrom(wildcards):
    return file2chrom[get_basename(wildcards.file)]

def get_basename(full_path: str):
    return os.path.splitext(os.path.basename(full_path))[0]
    
def remove_ext(full_path: str):
    return os.path.splitext(full_path)[0]
    
def count_rows_in_file(filename):
    try:
        with open(filename, 'r') as file:
            count = sum(1 for line in file)
        return count
    except FileNotFoundError:
        return f"File {filename} not found."


# ===================================
# READING CONFIGURATION
# ===================================

# Tools paths
BIOGWAS_PATH = config['biogwas_path']

PLINK_PATH = config['plink']
PLINK2_PATH = config['plink2']
BEDTOOLS_PATH = config['bedtools']
HAPGEN2_PATH = config['hapgen2']

# I/O config
BFILE_IN_FLAG=config['bfile_in_flag']
DATA_DIR = config['data_dir']
IMAGES_DIR = config['images_dir']

# VCFS_LIST -- file with structure: `path_to_vcf,chromosome` per row
VCFS_LIST = config['vcfs_list']
GENO_LIST = config['geno_list']
CHR_LIST = config['chr_list']
IDS_FILE = config['ids_file']
GTF_IN = config['anno_file']
GMT_IN = config['gmt_file']

# if use_causal_snps (default = True), causal SNPs are taken from `causal_snps` file
# otherwise -- randomly selected from `causal_snps` file.
SNPS_PROVIDED = config['use_causal_snps'] 

PATH_FILE = config['causal_pathways']
SNPS_FILE = config['causal_snps']

# General filter
SAMPLE_CALL_RATE = config["sample_call_rate"]
VARIANT_CALL_RATE = config["variant_call_rate"]
MAF_FILTER = config['maf_filter']
# Filter for causal SNPs only: defaults are [0.0,1.0]
CAUSAL_MAF_MIN = config['causal_maf_min'] 
CAUSAL_MAF_MAX = config['causal_maf_max']

# Simulate geno?
SKIP_SIMULATION = config['skip_simulation']
# Constants
N = config['N']

K = config['K']
k = config['k']

M_BETA = config['m_beta']
SD_BETA = config['sd_beta']
GEN_VAR = config['gen_var']
H2S = 1.0
P_INDEPENDENT_GENETIC = config['p_independent_genetic']
THETA = config['theta']
PHI = 1.0
ALPHA = config['alpha']


#IDs
# pattern of all files (starting from genotypes)
PATTERN = config['pattern']
# ID of chosen causal SNPS set for the simulations.
CAUSAL_ID = config['causal_id']
# ID of phenotypes simulation for selected causal SNPs
SIM_ID = config['sim_id']

# draw settings
DRAW_FLAG = config['draw_flag']

PCA_WIDTH = config['pca_width']
PCA_HEIGHT = config['pca_height']
PCA_DPI = config['pca_dpi']
QQ_WIDTH = config['qq_width']
QQ_HEIGHT = config['qq_height']
QQ_DPI = config['qq_dpi']
MH_WIDTH = config['mh_width']
MH_HEIGHT = config['mh_height']
MH_DPI = config['mh_dpi']

SEED = config['seed']

# genotype files and chromosomes -- will be assigned later
files = None
file2chrom = None
VCFS_LIST_SIMULATED= os.path.join(DATA_DIR, f"simulated_filt_sim.list")

# ==========================================
# Make file and chrom lists and dicts 
# ==========================================

if not GENO_LIST and VCFS_LIST:
    with open(VCFS_LIST) as f:
        file_chrom = [line.strip().split(',') for line in f]
        file2chrom = {get_basename(file):chrom for file, chrom in file_chrom}
        files = [remove_ext(file) for file, _ in file_chrom]
        assert len(files)==len(file2chrom), "Genotype base file names should be unique."
elif GENO_LIST and not VCFS_LIST:
    files = GENO_LIST.split(',')
    if CHR_LIST != '':
        chroms = str(CHR_LIST).split(',')
    else:
        chroms = [i+1 for i in range(len(GENO_LIST))]
    assert len(files)==len(chroms), "Length of genotype files and chroms are not equal."
    file2chrom = {get_basename(file):chrom for file, chrom in zip(files, chroms)}
    files = [remove_ext(file) for file in files]
    assert len(files)==len(file2chrom), "Genotype base file names should be unique."
else:
    raise ValueError('You have to provide either `--geno_list`, or `--input_list`.')
    
new_files = [os.path.join(DATA_DIR,f"{get_basename(f)}_filt_sim") for f in files]
new_files_pat = [os.path.join(DATA_DIR,f"{get_basename(f)}") for f in files]

with open(VCFS_LIST_SIMULATED, 'w') as f:
    f.write('\n'.join(new_files))
    
# if more than 1 file, we need to merge at the end
TO_MERGE = len(new_files)>1


# ==========================================
# Check vcf/bfile
# ==========================================

# VCF in -- always ok, 
# BFILE -- only if skip simulation
    
if not SKIP_SIMULATION and BFILE_IN_FLAG:
    raise ValueError('It is not possible to use bfiles when simulating genotypes: you have to input `.vcf` files, OR skip genotypes\' simulation step (--skip-simulation).')

    
    
    
    
    