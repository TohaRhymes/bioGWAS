#!/bin/bash

DATA_DIR=$1
GTF_IN=$2

cd $DATA_DIR

#extract only "gene" features from gtf
> gen_${GTF_IN}
awk -F'\t' 'NR>=6 && $3=="gene"' ${GTF_IN}.gtf >>  gen_${GTF_IN}.gtf
sed -i '1,$s/^chr//' gen_${GTF_IN}.gtf

awk -F'\t' -v OFS='\t' '$3=="gene"{split($9,a,";"); split(a[1],b," "); split(a[3],c," "); $NF="gene_id \"" b[2] "\"; gene_name \"" c[2] "\";"; print}' gen_${GTF_IN}.gtf > gen_${GTF_IN}.gtf


# sort by positions
sort -k1,1 -k4,4n gen_${GTF_IN}.gtf > gen_${GTF_IN}.sorted.gtf

