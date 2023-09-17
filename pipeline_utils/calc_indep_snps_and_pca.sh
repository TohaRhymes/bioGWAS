#!/usr/bin/env bash

BFILE_IN=$1
PLINK2_PATH=$2


# Make independent set of SNPs
${PLINK2_PATH} \
--bfile ${BFILE_IN} \
--indep-pairwise 1000kb 1 0.8 \
--out ${BFILE_IN}_indep

${PLINK2_PATH} \
--bfile ${BFILE_IN} \
--extract ${BFILE_IN}_indep.prune.in  \
--make-bed \
--out ${BFILE_IN}_indep

rm ${BFILE_IN}_indep.prune.in
rm ${BFILE_IN}_indep.prune.out

# Make PCA on independent and full sets of SNPS
${PLINK2_PATH} --bfile ${BFILE_IN}_indep --pca 10 approx --out ${BFILE_IN}_indep

${PLINK2_PATH} --bfile ${BFILE_IN} --pca 10 approx --out ${BFILE_IN}


