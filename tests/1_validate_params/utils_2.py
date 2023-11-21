import os

N=10000

# =============
# filenames
# =============
IDs = {'pattern': "test10000",
        'casual_id': "best_hallmark_K{K}",
       'sim_id':"m{m_beta}_sd{sd_beta}_gv{gen_var}_alpha{alpha}_theta{theta}_pIndep{pIndep}_seed{seed}"}



# =============
# filenames
# =============
# aux
DATA_DIR = "data"
PATTERN = IDs['pattern']
CAUSAL_ID = IDs['casual_id']
SIM_ID = IDs['sim_id']
res_CAUSAL_ID = CAUSAL_ID.replace("K{K}", "")
# to import
CAUSAL_SNP_FILE = os.path.join(DATA_DIR, f"{PATTERN}_{CAUSAL_ID}_chosen_snps.tsv")
GWAS_FILE = os.path.join(DATA_DIR, f"{PATTERN}_{CAUSAL_ID}_{SIM_ID}_gwas.tsv")
FILE_RESULTS=os.path.join(DATA_DIR, f"{PATTERN}_{res_CAUSAL_ID}_compare_results.tsv")
BFILE=os.path.join(DATA_DIR, f"{PATTERN}_filt_sim")