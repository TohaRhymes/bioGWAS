#!/usr/bin/env bash


INPUT_VCF=$1
INPUT_GTF=$2
OUTPUT_VCF=$3
OUTPUT_TXT=$4

bedtools closest -a ${INPUT_VCF} -b ${INPUT_GTF}  > ${OUTPUT_VCF}

value=$(tail -n 1 ${OUTPUT_VCF} | awk '{print NF}')
col=$(echo $value-3 | bc)
cut -f $col ${OUTPUT_VCF} | uniq > ${OUTPUT_TXT}
