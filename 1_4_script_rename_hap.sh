#!/bin/bash

DATA_DIR=$1
PATTERN=$2

now=$(date)
echo "Started at: $now"

cd $DATA_DIR

rm *sim.cases*
echo "Removed *sim.cases*"

rename 's,sim\.controls,sim,' *
echo "Renamed *sim.controls* to *sim*"


for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do 
	sed -i s/"snp_"/"${CHR} snp${CHR}_"/g ${PATTERN}${CHR}_FILT_sim.gen
	now=$(date)
	echo "Finished with chr $CHR at: $now"

done

now=$(date)
echo "Finished at: $now"

