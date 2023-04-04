now=$(date)
echo "Started at: $now"

# in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do    
    echo "Started with $CHR chromosome!" 
    
    NAME=ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf
    NAME_TO_PREF=chr${CHR}
    NAME_TO=${NAME_TO_PREF}.vcf
    CHR_FILE_EUR=${NAME_TO_PREF}_EUR
    
    if test -f ${NAME}.gz; then
        gzip -d ${NAME}.gz
    fi

    first_pos=`grep -n -e '##contig=<ID=22,' ${NAME}  | cut -f1 -d:`
    mt_pos=`grep -n -e '##contig=<ID=MT,' ${NAME}  | cut -f1 -d:`
    x_pos=`grep -n -e '##contig=<ID=X,' ${NAME}  | cut -f1 -d:`
    y_pos=`grep -n -e '##contig=<ID=Y,' ${NAME}  | cut -f1 -d:`
    end_pos=`grep -n -e '##contig=<ID=hs37d5,' ${NAME}  | cut -f1 -d:`
    ((end_pos++))

    head -${first_pos} ${NAME} > ${NAME_TO}
    head -${mt_pos} ${NAME} | tail -1 >> ${NAME_TO}
    head -${x_pos} ${NAME} | tail -1 >> ${NAME_TO}
    head -${y_pos} ${NAME} | tail -1 >> ${NAME_TO}

    tail -n +${end_pos} ${NAME} >> ${NAME_TO}

    # Gzip origial 1khg file
    gzip ${NAME}
    
    
    # Filter VCF-файл (we want only 2-alleles SNPs and only Europeas + MAF > 0.05)
    plink2 --vcf ${NAME_TO_PREF}.vcf --max-alleles 2 --maf 0.05 --keep /media/MIRROR/ukb_finngen/1000genomes/EUR_SAMPLES_ID.txt  --recode vcf --out ${CHR_FILE_EUR}
    # Recode vcf to hapgen-needed files
    plink2 --vcf ${CHR_FILE_EUR}.vcf  --export ped  --out ${CHR_FILE_EUR}
    plink2 --vcf ${CHR_FILE_EUR}.vcf  --export hapslegend --out  ${CHR_FILE_EUR}
    
    
    now=$(date)
    echo "Done with $CHR chromosome!" 
    echo "Time now: $now"
done

now=$(date)
echo "Finished at: $now"
