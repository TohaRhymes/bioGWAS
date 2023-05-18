#!/usr/bin/env bash

BFILE_IN=$1


# Make independent set of SNPs
plink2 \
--bfile ${BFILE_IN} \
--indep-pairwise 1000kb 1 0.8 \
--out ${BFILE_IN}_indep

plink2 \
--bfile ${BFILE_IN} \
--extract ${BFILE_IN}_indep.prune.in  \
--make-bed \
--out ${BFILE_IN}_indep

rm ${BFILE_IN}_indep.prune.in
rm ${BFILE_IN}_indep.prune.out

# Make PCA on independent and full sets of SNPS
plink2 --bfile ${BFILE_IN}_indep --pca 10 approx --out ${BFILE_IN}_indep

plink2 --bfile ${BFILE_IN} --pca 10 approx --out ${BFILE_IN}


