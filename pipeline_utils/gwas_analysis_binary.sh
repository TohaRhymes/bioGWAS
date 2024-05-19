#!/usr/bin/env bash

BFILE_IN=$1
PHENO_IN=$2
SUMSTAT_OUT=$3
PLINK_PATH=$4


# Make a GWA
${PLINK_PATH} \
--bfile ${BFILE_IN} \
--allow-no-sex \
--assoc \
--1 \
--pheno ${PHENO_IN}.tsv \
--out ${SUMSTAT_OUT}


# Convert header a bit
sed '1{ s/CHR/chr/; s/SNP/rsid/; s/BP/pos/; s/NMISS/n/; s/BETA/beta/; s/SE/se/; s/R2/r2/; s/T/t/; s/P/pval/;}' ${SUMSTAT_OUT}.assoc | awk '{$1=$1};1'| tr -s ' ' '\t' > ${SUMSTAT_OUT}.tsv

# Convert qassoc to tsv-file

awk '{$1=$1};1' ${SUMSTAT_OUT}.assoc | tr -s ' ' '\t' > ${SUMSTAT_OUT}_.assoc
rm -f ${SUMSTAT_OUT}.assoc
mv ${SUMSTAT_OUT}_.assoc ${SUMSTAT_OUT}.qassoc



