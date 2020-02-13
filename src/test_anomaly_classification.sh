#!/bin/sh

# # =====================================================
# # Generate test data
# # =====================================================
# =====================================================
## Select barcode with Abnormal(target_deletion) > 80%
# =====================================================
true > .tmp_/barcodelist
for barcode in $(cat .tmp_/anomaly_classification_revised.txt | cut -f 1 | sort -u); do
    cat .tmp_/anomaly_classification_revised.txt |
    awk -v id=${barcode} '$1 == id' |
    cut -f 3 |
    sort |
    uniq -c |
    awk '{if($2 ~ "target_deletion") td=$1;
        sum+=$1}
        END{print int(td/sum*100)}' |
    xargs -I @ [ @ -gt 80 ] &&
    printf "${barcode}\n" \
    >> .tmp_/barcodelist
done
# ========================================
# ========================================
output_csv="abnormal_class_info.csv"
output_dir="bam/abnormal_class"
#
printf "barcodeID,cluster,number of reads\n" \
> ${output_dir}/${output_csv}
mkdir -p ${output_dir}
#
for barcode in $(cat .tmp_/barcodelist); do
    cat .tmp_/anomaly_classification_revised.txt |
    awk -v id=${barcode} '$1 == id' |
    grep target_deletion |
    cut -f 2 |
    awk '{gsub("-.*","",$1); print}' |
    sort > .tmp_/tmp_id
    #
    # Extract ID, cuttind sites and cutting size
    #
    samtools view bam/${barcode}.bam |
    awk '{gsub("-.*","",$1); print}' |
    sort -k 1,1 |
    join .tmp_/tmp_id - 2>/dev/null |
    awk '$2==0 || $2==16' |
    sed "s/ /\t/g" |
    awk '{print $1,$4,$(NF-1)}' |
    sed -e "s/cs:Z:://g" \
    -e "s/~[a-z]*/ /g" |
    awk '{gsub("[a-z].*","",$NF); print}' |
    sed -e "s/:/ /g" -e "s/\*[a-z][a-z]/+1/g" \
        -e "s/+/ /g" -e "s/-/ -/g" |
    awk '{for(i=3;i<=NF-1;i++) if($i~"-") $i="-"length($i)-1;
        else if($i ~ /[a-z]/) $i=length($i);
        print}' |
    awk '{sum=1; for(i=3;i<=NF-1;i++) sum+=$i
        print $1,$2+sum,$2+sum+$NF, $NF}' \
    > .tmp_/tmp_anomaly_classification
    # =============================
    python ./DAJIN/src/test_anomaly_classification.py
    # =============================
    for cl in $(cat .tmp_/tmp_anomaly_label | awk '{print $NF}' | sort | uniq); do
        #
        cat .tmp_/tmp_anomaly_label |
        awk -v cl=${cl} '$NF==cl' |
        cut -f 1 |
        sort \
        > .tmp_/tmp_
        #
        # if [ ${cl} -eq 0 ]; then
        #     cl="unclassified"
        # else
        #     cl="cluster"${cl}
        # fi
        cl="cluster"${cl}
#        [ ${cl} -eq 0 ] && cl="unclassified"
#        [ ${cl} -ne 0 ] && cl="cluster"${cl}
        #
        samtools view -h bam/${barcode}.bam |
        grep "^@" >.tmp_/tmp_header
        #
        samtools view bam/${barcode}.bam |
        awk '$2==0 || $2==16' |
        awk '{gsub("-.*","",$1); print}' |
        sort |
        join .tmp_/tmp_ - 2>/dev/null |
        head -n 100 |
        sed "s/ /\t/g" \
        >> .tmp_/tmp_header
        #
        samtools sort .tmp_/tmp_header > ${output_dir}/${barcode}_${cl}.bam
        samtools index ${output_dir}/${barcode}_${cl}.bam
        #
        cat .tmp_/tmp_ |
        wc -l |
        sed "s/^/${barcode},${cl},/g" \
        >> ${output_dir}/${output_csv}
    done
done
rm .tmp_/tmp_*
exit 0

