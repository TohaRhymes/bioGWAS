
now=$(date)
echo "[LOG] Started at: $now"

> chr_EUR_sim.list

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do 
	plink2 \
	--data chr${CHR}_EUR_sim ref-first \
	--make-bed \
	--out chr${CHR}_EUR_sim

	echo "chr${CHR}_EUR_sim" >> chr_EUR_sim.list
	now=$(date)
	echo "[LOG] Finished with chr $CHR at: $now"

done

now=$(date)
echo "[LOG] Finished converting at: $now"

echo "[LOG] Started merging at: $now"

plink2 \
--pmerge-list chr_EUR_sim.list bfile \
--set-all-var-ids @:#[b37] \
--make-bed \
--out chr_EUR_sim


now=$(date)
echo "[LOG] Started converting to vcf at: $now"

plink2 --bfile chr_EUR_sim \
--recode vcf \
--out chr_EUR_sim

now=$(date)
echo "[LOG] Started converting to ucsc's bed at: $now"

vcf2bed < chr_EUR_sim.vcf > chr_EUR_sim_ucsc.bed

now=$(date)
echo "[LOG] Started filtering ucsc's bed at: $now"

awk 'NF==5010' chr_EUR_sim_ucsc.bed > filtered_chr_EUR_sim_ucsc.bed

now=$(date)
echo "[LOG] Finished at: $now"

