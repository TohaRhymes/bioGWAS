#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint

# try these m_betas
m_betas = (
    0.01, 
    0.05, 
    0.1, 
    0.2, 
    0.3, 
    0.5)
# use these sd for m_betas above respectively
sd_betas = (
    (0.001, 0.005),
    (0.001, 0.01),
    (0.001, 0.01, 0.05),
    (0.01, 0.05, 0.1),
    (0.01, 0.05, 0.1),
    (0.01, 0.1, 0.3),
)

# try these m_betas
m_betas = (
    0.2, 
    0.3)
# use these sd for m_betas above respectively
sd_betas = (
    (0.1,),
    (0.1,)
)

Ks = [
     10, 
     30, 
    50
]

assert len(m_betas) == len(sd_betas)
# make products
m_sd_comb = []
for mb, sds in zip(m_betas, sd_betas):
    m_sd_comb += list(product([mb], sds))
# for every of combinations above, use these genvar and sd
gen_vars = (0.05, 0.5)
h2ss = (0.1, 0.5)
gv_h2s_comb = list(product(gen_vars, h2ss))

print("These combinations of m and sd:")
pprint(m_sd_comb)
print("These combinations of gv and h2s:")
pprint(gv_h2s_comb)

params = [
    {"m_beta": m, "sd_beta": sd, "gen_var": gv, "h2s": h2s}
    for (m, sd), (gv, h2s) in list(product(m_sd_comb, gv_h2s_comb))
]
pprint(params)

for K in Ks:
    k = K // 2
    for param in params:
        m_beta = param["m_beta"]
        sd_beta = param["sd_beta"]
        gen_var = param["gen_var"]
        h2s = param["h2s"]
        sim_id = f"m{m_beta}_sd{sd_beta}_gv{gen_var}_h2s{h2s}"
        phenos_id = f"ph_sperm_K{K}"
        command = f"snakemake \
        --cores 8 \
        --config \
        m_beta={m_beta} \
        sd_beta={sd_beta} \
        gen_var={gen_var} \
        h2s={h2s} \
        K={K} \
        k={k} \
        draw_flag={True} \
        sim_id={sim_id} \
        phenos_id={phenos_id}"
        print(command)
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print("\n=============================\nOUTPUT: \n", output)
        print("\n=============================\nERRORS: \n", error)
        
print('Done with params iterations')