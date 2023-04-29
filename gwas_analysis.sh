#!/usr/bin/env bash

# Make independent set of SNPs
plink2 \
--bfile chr_EUR_sim \
--indep-pairwise 1000kb 1 0.8 \
--out chr_EUR_sim_indep

plink2 \
--bfile chr_EUR_sim \ 
--extract chr_EUR_sim_indep.prune.in  \
--make-bed \
--out chr_EUR_sim_indep

# Make PCA on independent and full sets of SNPS
plink2 --bfile chr_EUR_sim_indep --pca 10 approx --out chr_EUR_sim_indep

plink2 --bfile chr_EUR_sim --pca 10 approx --out chr_EUR_sim



# Make a GWAS
plink \
--bfile chr_EUR_sim \
--allow-no-sex \
--assoc \
--pheno phenos.tsv \
--out chr_EUR_sim_p1


# Convert header a bit
sed '1{ s/CHR/chr/; s/SNP/rsid/; s/BP/pos/; s/NMISS/n/; s/BETA/beta/; s/SE/se/; s/R2/r2/; s/T/t/; s/P/pval/;}' chr_EUR_sim_p1.qassoc | awk '{$1=$1};1'| tr -s ' ' '\t' > chr_EUR_sim_p1.tsv

