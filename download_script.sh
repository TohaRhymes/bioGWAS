#!/bin/bash

DATENOW=$( date )
echo "Started at: $DATENOW"


# Declare an array of string with type
declare -a files=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
 
# Iterate the string array using for loop
for file in ${files[@]}; do
   echo ${file}
   while [ $( ps -f -u $USER | grep 'wget ftp' | wc -l ) -ge 7 ]; do sleep 1; done
   echo ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${file}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
   wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${file}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz &
   echo ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${file}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
   wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${file}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi &
done
wait

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz.tbi &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz.tbi &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_vcf_info_annotation.20141104 &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_chrY_calls_20141104 &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/20140625_related_individuals.txt &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220 &
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped &
wait




DATENOW=$( date )
echo "Finished at: $DATENOW"

