from itertools import product
from pprint import pprint

import os


# try these m_betas
m_betas = (
    0.05,
    0.5)
# use these sd for m_betas above respectively
sd_betas = (
    (0.001, 0.01),
    (0.01, 0.1),
)

assert len(m_betas) == len(sd_betas)
# make products
m_sd_comb = []
for mb, sds in zip(m_betas, sd_betas):
    m_sd_comb += list(product([mb], sds))
# for every of combinations above, use these genvar and sd
gen_vars = (0, 0.1, 0.5, 0.9, 1)
alphas = (0, 0.5)
gv_alpha_comb = list(product(gen_vars, alphas))

theta_pIndep_comb = [(0,1), (0.5, 0.5), (1, 0)]


print("These combinations of m and sd:")
pprint(m_sd_comb)
print("These combinations of gv and alpha:")
pprint(gv_alpha_comb)
print("These combinations of theta and pIndep:")
pprint(theta_pIndep_comb)

# =============
# sim params
# =============
params = [
    {"m_beta": m, 
     "sd_beta": sd, 
     "gen_var": gv,
     "alpha":alpha,
     "theta": theta, 
     "pIndep": pIndep}
    for (m, sd), 
    (gv, alpha), 
    (theta, pIndep) 
    in list(product(m_sd_comb,
                    gv_alpha_comb,
                    theta_pIndep_comb))
]



Ks = [
    10,
    20, 
    30
]

N=10000

# =============
# bioGWAS ids
# =============
IDs = {'pattern': "test10000",
      'casual_id': "hallmark_K{K}",
      'sim_id': "m{m_beta}_sd{sd_beta}_gv{gen_var}_alpha{alpha}_theta{theta}_pIndep{pIndep}"}

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



VALIDATION_COMPARISON = 'data/aggregated_sim1.csv'


