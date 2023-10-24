../../1000genomes/pipeline_utils/calc_indep_snps_and_pca.sh \
out_data/SIM_filt_sim \
/home/achangalidi/tools/plink2

../../1000genomes/pipeline_utils/draw_pca.py \
out_data/SIM_filt_sim_indep \
images/SIM_filt_sim_indep
../../1000genomes/pipeline_utils/draw_pca.py \
out_data/SIM_filt_sim \
images/SIM_filt_sim




cd ../../1000genomes
ls data3/chr*.bed | sed 's/\.bed$//' | sort -V > data3/chr_filt.list


plink2 \
--pmerge-list data3/chr_filt.list bfile \
--set-all-var-ids @:#\$r,\$a --new-id-max-allele-len 5000 missing \
--max-alleles 2 \
--multiallelics-already-joined \
--make-bed \
--out data3/chr_filt

pipeline_utils/calc_indep_snps_and_pca.sh \
data3/chr_filt \
/home/achangalidi/tools/plink2

pipeline_utils/draw_pca.py \
data3/chr_filt_indep \
images/chr_filt_indep
pipeline_utils/draw_pca.py \
data3/chr_filt \
images/chr_filt



