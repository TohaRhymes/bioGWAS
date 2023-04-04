
now=$(date)
echo "Started at: $now"

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do 
	now=$(date)
	echo "Started with chr $CHR at: $now"

	CHR_FILE=chr${CHR}
	CHR_FILE_EUR=${CHR_FILE}_EUR
	SNP_CASE=$(tail -1 ${CHR_FILE_EUR}.vcf |  awk '{ print $2 }')

	hapgen2 \
	-h ${CHR_FILE_EUR}.haps \
	-l ${CHR_FILE_EUR}.legend \
	-m ${CHR_FILE_EUR}.map \
	-o ${CHR_FILE_EUR}_sim \
	-dl ${SNP_CASE} 1 1.5 2.25 \
	-int 0 500000000 \
	-n 5000 0 \
	-Ne 11418 \
	-theta 1


	now=$(date)
	echo "Finished with chr $CHR at: $now"

done

now=$(date)
echo "Finished at: $now"

