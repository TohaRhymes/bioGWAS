#!/bin/bash

for i in ./in_data/*_rsid.tsv
do
	while [ $( ps -AF | grep 'Pascal' | wc -l ) -ge 5 ] ; do sleep 1 ; done
	TAG=$( basename $i )
	cut -f 2,9 $i | tail -n +2  | grep 'rs' > ${TAG%%.tsv}.txt
	./Pascal --runpathway=on --pval=${TAG%%.tsv}.txt &
done
