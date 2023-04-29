#!/usr/bin/env bash

BFILE_IN=$1
PHENO_IN=$2
SUMSTAT_OUT=$3


# Make a GWA
plink \
--bfile ${BFILE_IN} \
--allow-no-sex \
--assoc \
--pheno ${PHENO_IN}.tsv \
--out ${SUMSTAT_OUT}


# Convert header a bit
sed '1{ s/CHR/chr/; s/SNP/rsid/; s/BP/pos/; s/NMISS/n/; s/BETA/beta/; s/SE/se/; s/R2/r2/; s/T/t/; s/P/pval/;}' ${SUMSTAT_OUT}.qassoc | awk '{$1=$1};1'| tr -s ' ' '\t' > ${SUMSTAT_OUT}.tsv

rm ${SUMSTAT_OUT}.qassoc
