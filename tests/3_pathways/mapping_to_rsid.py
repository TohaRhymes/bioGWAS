#!/usr/bin/env python
# coding: utf-8

# In[57]:


import pandas as pd
import tabix
import numpy as np
import os
import glob

from time import time
from tqdm import tqdm
## step 0: gunzip -c -k -v /media/MIRROR/ukb_finngen/imputed_v3/variants.tsv.bgz > /home/achangalidi/ukb_finngen/gwassim_check/path_sim/in_data/variants.tsv
# step 1: 
DATA_MAP = "/home/achangalidi/ukb_finngen/gwassim_check/path_sim/in_data/variants.tsv"
file_pattern = "path_sim/_out_data/*random*_gwas.tsv"


def map_rsid(x):
    try:
        return _pos2rsid[x]
    except KeyError:
        return np.nan


# In[17]:


# Convert the records to a DataFrame
data = pd.read_csv(DATA_MAP, sep='\t')
data['variant_short'] = data.chr.astype(str)+":"+data.pos.astype(str)
_pos2rsid = data[['variant_short', 'rsid']].set_index('variant_short').to_dict()['rsid']


# In[66]:


# Use glob to get all the filenames that match the pattern
file_list = glob.glob(file_pattern)

# Iterate over the files
for sample_name in tqdm(file_list):
    name_to_save =  os.path.splitext(sample_name)[0]+'_rsid.tsv'
    gwas = pd.read_csv(sample_name, sep='\t')
    gwas = gwas.rename(columns={'rsid':'variant'})
    gwas['variant_short'] = gwas.chr.astype(str)+":"+gwas.pos.astype(str)
    gwas['rsid'] = gwas.variant_short.apply(map_rsid)
    gwas = gwas[~gwas.rsid.isna()]
    gwas[['chr', 'rsid', 'pos', 'n', 'beta', 'se', 'r2', 't', 'pval',]].to_csv(name_to_save, sep='\t', index=False)
