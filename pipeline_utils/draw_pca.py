#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import sys

import warnings
warnings.filterwarnings('ignore')


def draw_PCA(FILE_PREFIX, IMAGES_PREFIX):
    VEC = f"{FILE_PREFIX}.eigenvec"
    VAL = f"{FILE_PREFIX}.eigenval"
    TO = f"{IMAGES_PREFIX}_pca.pdf"

    data = pd.read_csv(VEC, sep='\t')
    eigenval = pd.read_csv(VAL, sep='\t', header=None)
    ev = np.array(eigenval[0])
    ev = ev/ev.sum()*100
    print(f"Dispersion: {ev}")

    plt.figure(figsize=(18,5))

    fig, axs = plt.subplots(nrows=1, ncols=3)
    for i in range(1,4,1):
        ax = sns.scatterplot(data, x=f'PC{i}', y=f'PC{i+1}', ax=axs[i-1])
        ax.set(xlabel=f'PC{i}, {round(ev[i-1], 2)}%', 
               ylabel=f'PC{i+1}, {round(ev[i], 2)}%')
    plt.savefig(TO)
    print(f"Draw to {TO}.")
    plt.show()
    

FILE_PREFIX = sys.argv[1]
IMAGES_PREFIX = sys.argv[2]
draw_PCA(FILE_PREFIX, IMAGES_PREFIX)