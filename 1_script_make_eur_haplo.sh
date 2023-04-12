#!/bin/bash

DATA_DIR=$1
PATTERN=$2
MAF_FILTER=$3
IDS_FILE=$4

now=$(date)
echo "Started filtering at: $now"

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do    
    echo "Started with $CHR chromosome!" 
    
    NAME_TO_PREF=${PATTERN}${CHR}
    NAME_TO=${DATA_DIR}/${NAME_TO_PREF}.vcf
    CHR_FILE_EUR=${DATA_DIR}/${NAME_TO_PREF}_FILT
    
    # Filter VCF-файл (we want only 2-alleles SNPs and only Europeas + MAF > 0.05)
    plink2 --vcf ${NAME_TO_PREF}.vcf \
    --max-alleles 2 \
    --maf $MAF_FILTER \
    --keep $IDS_FILE  \
    --recode vcf \
    --out ${CHR_FILE_EUR}
    # Recode vcf to hapgen-needed files
    plink2 --vcf ${CHR_FILE_EUR}.vcf  --export ped  --out ${CHR_FILE_EUR}
    plink2 --vcf ${CHR_FILE_EUR}.vcf  --export hapslegend --out  ${CHR_FILE_EUR}
    
    
    now=$(date)
    echo "Done with $CHR chromosome!" 
    echo "Time now: $now"
done

now=$(date)
echo "Finished at: $now"
