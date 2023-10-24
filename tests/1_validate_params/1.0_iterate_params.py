#!/usr/bin/env python

from itertools import product
import subprocess
from pprint import pprint

# try these m_betas
m_betas = (
    0.05,
    0.5)
# use these sd for m_betas above respectively
sd_betas = (
    (0.001, 0.01),
    (0.05, 0.2),
)

Ks = [
     10,
    30
]

assert len(m_betas) == len(sd_betas)
# make products
m_sd_comb = []
for mb, sds in zip(m_betas, sd_betas):
    m_sd_comb += list(product([mb], sds))
# for every of combinations above, use these genvar and sd
gen_vars = (0, 0.1, 0.5, 0.9, 1)
h2ss = (1.0,)
gv_h2s_comb = list(product(gen_vars, h2ss))


theta_pIndep_comb = [(0,1), (0.5, 0.5), (1, 0)]

phi_alpha_comb = [(1, 0), (1, 0.5)]

print("These combinations of m and sd:")
pprint(m_sd_comb)
print("These combinations of gv and h2s:")
pprint(gv_h2s_comb)
print("These combinations of theta and pIndep:")
pprint(theta_pIndep_comb)
print("These combinations of phi and alpha:")
pprint(phi_alpha_comb)

params = [
    {"m_beta": m, 
     "sd_beta": sd, 
     "gen_var": gv, 
     "h2s": h2s, 
     "theta": theta, 
     "pIndep": pIndep, 
     "phi":phi, 
     "alpha":alpha}
    for (m, sd), 
    (gv, h2s), 
    (theta, pIndep), 
    (phi, alpha) in list(product(m_sd_comb, 
                                           gv_h2s_comb, 
                                           theta_pIndep_comb,
                                           phi_alpha_comb))
]
pprint(params)
print(f"Amount of parameters sets: {len(params)}")

for K in Ks:
    k = K // 2
    print(f"Started for K={K}")
    for param in params:
        m_beta = param["m_beta"]
        sd_beta = param["sd_beta"]
        gen_var = param["gen_var"]
        h2s = param["h2s"]
        theta = param["theta"]
        pIndep = param["pIndep"]
        phi = param["phi"]
        alpha = param["alpha"]
        causal_id = f"ph_K{K}"
        sim_id = f"m{m_beta}_sd{sd_beta}_gv{gen_var}_h2s{h2s}_theta{theta}_pIndep{pIndep}_phi{phi}_alpha{alpha}"
        command = f"snakemake \
        --cores 8 \
        --config \
        m_beta={m_beta} \
        sd_beta={sd_beta} \
        gen_var={gen_var} \
        h2s={h2s} \
        theta={theta} \
        pIndep={pIndep} \
        phi={phi} \
        alpha={alpha} \
        K={K} \
        k={k} \
        draw_flag={False} \
        sim_id={sim_id} \
        causal_id={causal_id}"
        print(command)
        process = subprocess.Popen(
             command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        output, error = process.communicate()
        print("\n=============================\nOUTPUT: \n", output)
        print("\n=============================\nERRORS: \n", error)
        
print('Done with params iterations')
