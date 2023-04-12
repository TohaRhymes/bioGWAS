#extract only "gene" features from gtf
> gen_gencode.v37.annotation.gtf
awk -F'\t' 'NR>=6 && $3=="gene"' gencode.v37.annotation.gtf >>  gen_gencode.v37.annotation.gtf
sed -i '1,$s/^chr//' gen_gencode.v37.annotation.gtf

awk -F'\t' -v OFS='\t' '$3=="gene"{split($9,a,";"); split(a[1],b," "); split(a[3],c," "); $NF="gene_id \"" b[2] "\"; gene_name \"" c[2] "\";"; print}' gen_gencode.v37.annotation.gtf > gen_gencode.v37.annotation.gene.gtf


# sort by positions
sort -k1,1 -k4,4n gen_gencode.v37.annotation.gene.gtf > gen_gencode.v37.annotation.gene.sorted.gtf

