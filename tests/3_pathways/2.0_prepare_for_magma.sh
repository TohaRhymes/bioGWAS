# =====================================================
# Preparation (just one time) -- all commands in bash
# =====================================================

GWAS_OUT_DATA=./data

GENE_LOC=../data/gencode.v37.annotation.gtf
GENE_SETS=../data/c2.cp.kegg.v2023.1.Hs.symbols.gmt
GENE_SETS_ALL=../data/all.gmt
ANNO_FILE=./data_enrich/magma_anno

BFILE=./data/test10000_filt_sim

# ------------
# Step 0
# ------------
# CHANGE \t to spaces in gmt files

sed 's/\t/ /g' ${GENE_SETS} | awk '{for (i=1; i<=NF; i++) if (i!=2) printf $i (i==NF? RS : FS)}' > ${GENE_SETS}.ssv
sed 's/\t/ /g' ${GENE_SETS_ALL} | awk '{for (i=1; i<=NF; i++) if (i!=2) printf $i (i==NF? RS : FS)}' > ${GENE_SETS_ALL}.ssv

# ------------
# Step 1
# ------------
# Make custom .loc file from gtf

cat ${GENE_LOC} \
| awk -F'\t' -v OFS='\t' \
'$3=="gene"{split($9,a,";"); split(a[1],b," "); split(a[3],c," "); $NF="\""c[2]"\""; print}' \
| awk -F'\t' \
'!seen[$9]++ {gsub("\"", "", $9); gsub("chr", "", $1); print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
> ${GENE_LOC}.loc

# ------------
# Step 2
# ------------
# Make magma annotations

magma \
--annotate window=1,0.5 \
--snp-loc ${BFILE}.bim \
--gene-loc  ${GENE_LOC}.loc \
--out ${ANNO_FILE}
