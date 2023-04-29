for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do    
    echo $CHR
    echo started!
    NAME=ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf
    NAME_TO_PREF=chr${CHR}
    NAME_TO=${NAME_TO_PREF}.vcf
    
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

    plink --allow-no-sex --vcf ${NAME_TO} --make-bed --out ${NAME_TO_PREF}
    echo $CHR
    echo done! 
done

CHR=X
echo $CHR
echo started!
NAME=ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf
NAME_TO_PREF=chr${CHR}
NAME_TO=${NAME_TO_PREF}.vcf

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

plink --allow-no-sex --vcf ${NAME_TO} --make-bed --out ${NAME_TO_PREF}
echo $CHR
echo done!


CHR=Y
echo $CHR
echo started!
NAME=ALL.chr${CHR}.phase3_integrated_v2b.20130502.genotypes.vcf
NAME_TO_PREF=chr${CHR}
NAME_TO=${NAME_TO_PREF}.vcf

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

plink --allow-no-sex --vcf ${NAME_TO} --make-bed --out ${NAME_TO_PREF}
echo $CHR
echo done!
