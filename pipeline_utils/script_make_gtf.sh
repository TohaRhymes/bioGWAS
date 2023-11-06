#!/usr/bin/env bash

GTF_IN=$1
GTF_OUT=$2

GTF_MID=${GTF_OUT}_gen.gtf
GTF_MID_SEL=${GTF_OUT}_gen_sel.gtf

#extract only "gene" features from gtf (and delete comments from the beginning)
> ${GTF_MID}
sed -n '/^[^#]/,$p' ${GTF_IN} | awk -F'\t' '$3=="gene"'  | grep -E "^chr[0-9]+" >>  ${GTF_MID}
sed -i '1,$s/^chr//' ${GTF_MID}

awk -F'\t' -v OFS='\t' '$3=="gene"{split($9,a,";"); split(a[1],b," "); split(a[3],c," "); $NF="gene_id \"" b[2] "\"; gene_name \"" c[2] "\";"; print}' ${GTF_MID} > ${GTF_MID_SEL}


# sort by positions
sort -k1,1n -k4,4n ${GTF_MID_SEL}  > ${GTF_OUT}

rm -f ${GTF_MID}
rm -f ${GTF_MID_SEL}
