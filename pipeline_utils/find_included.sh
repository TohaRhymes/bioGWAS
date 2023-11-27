#!/usr/bin/env bash

OUTPUT_VCF=$1
OUTPUT_TXT=$2

value=$(tail -n 1 ${OUTPUT_VCF} | awk '{print NF}')
col=$(echo $value-3 | bc)
cut -f $col ${OUTPUT_VCF} | uniq > ${OUTPUT_TXT}
