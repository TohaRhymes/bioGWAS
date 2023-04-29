#!/bin/bash


DATA_DIR=$1
PATTERN=$2
GTF_IN=$3

cd ${DATA_DIR}

bedtools closest -d -a filtered_${PATTERN}_FILT_sim_ucsc.bed -b gen_$2.sorted.gtf  > ${PATTERN}_FILT_sim.sorted.annotated.bed