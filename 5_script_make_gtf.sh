#!/bin/bash

DATA_DIR=$1
GTF_IN=$2

cd $DATA_DIR

#extract only "gene" features from gtf (and delete comments from the beginning)
> gen_${GTF_IN}.gtf
sed -n '/^[^#]/,$p' ${GTF_IN}.gtf | awk -F'\t' '$3=="gene"' >>  gen_${GTF_IN}.gtf
sed -i '1,$s/^chr//' gen_${GTF_IN}.gtf

awk -F'\t' -v OFS='\t' '$3=="gene"{split($9,a,";"); split(a[1],b," "); split(a[3],c," "); $NF="gene_id \"" b[2] "\"; gene_name \"" c[2] "\";"; print}' gen_${GTF_IN}.gtf > selected_gen_${GTF_IN}.gtf


# sort by positions
sort -k1,1 -k4,4n selected_gen_${GTF_IN}.gtf > ${GTF_IN}_filt_sort.gtf

rm gen_${GTF_IN}.gtf
rm selected_gen_${GTF_IN}.gtf

cd ..

