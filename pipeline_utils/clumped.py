#!/usr/bin/env python
# coding: utf-8

import subprocess
import sys
import os

def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    output, error = process.communicate()
    return output, error


BFILE=sys.argv[1]
PFILE=sys.argv[2]
PVAL_CUTOFF=float(sys.argv[3])
R2=float(sys.argv[4])
KB=int(sys.argv[5])
OFILE=os.path.splitext(BFILE)[0]+'_clumps'


command_clump = f"""plink \
--bfile {BFILE} \
--allow-no-sex \
--clump {PFILE} \
--clump-p1 {'{:.30f}'.format(PVAL_CUTOFF)} \
--clump-p2 {'{:.30f}'.format(PVAL_CUTOFF)} \
--clump-r2 {R2} \
--clump-kb {KB} \
--out {OFILE}"""
print(command_clump)
run(command_clump)

command_get = f"""awk '{{ print $3 }}' {OFILE}.clumped > {OFILE}.txt"""
print(command_get)
run(command_get)

# command_rm = f"""rm {OFILE}.clumped"""
# print(command_rm)
# run(command_rm)


