#!/bin/bash

DATA_DIR=$1
PATTERN=$2
N_10=$3

cd $DATA_DIR

now=$(date)
echo "[LOG] Started at: $now"

> ${PATTERN}_FILT_sim.list

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do 
	plink2 \
	--data ${PATTERN}${CHR}_FILT_sim ref-first \
	--make-bed \
	--out ${PATTERN}${CHR}_FILT_sim

	echo "${PATTERN}${CHR}_FILT_sim" >> ${PATTERN}_FILT_sim.list
	now=$(date)
	echo "[LOG] Finished with chr $CHR at: $now"

done

now=$(date)
echo "[LOG] Finished converting at: $now"

echo "[LOG] Started merging at: $now"

plink2 \
--pmerge-list ${PATTERN}_FILT_sim.list bfile \
--set-all-var-ids @:#[b37] \
--make-bed \
--out ${PATTERN}_FILT_sim


now=$(date)
echo "[LOG] Started converting to vcf at: $now"

plink2 --bfile ${PATTERN}_FILT_sim \
--recode vcf \
--out ${PATTERN}_FILT_sim

now=$(date)
echo "[LOG] Started converting to ucsc's bed at: $now"

vcf2bed < ${PATTERN}_FILT_sim.vcf > ${PATTERN}_FILT_sim_ucsc.bed

now=$(date)
echo "[LOG] Started filtering ucsc's bed at: $now"

awk -v N=${N_10} 'NF==N' ${PATTERN}_FILT_sim_ucsc.bed > filtered_${PATTERN}_FILT_sim_ucsc.bed

now=$(date)
echo "[LOG] Finished at: $now"

