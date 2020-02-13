#!/bin/sh

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

# =====================================================
## Define output files
# =====================================================
# -----------------------------------------------------
output_csv="abnormal_class_info.csv"
output_dir="bam/abnormal_class"
output_dir_seqlogo="bam/abnormal_class/seqlogo"
# -----------------------------------------------------
#
printf "barcodeID,cluster,number of reads\n" \
> ${output_dir}/${output_csv}
mkdir -p ${output_dir} ${output_dir_seqlogo}
#
for barcode in $(cat .tmp_/barcodelist); do
    # =====================================================
    ## Format input data
    # =====================================================
    cat .tmp_/anomaly_classification_revised.txt |
    awk -v id=${barcode} '$1 == id' |
    grep target_deletion |
    cut -f 2 |
    awk '{gsub("-.*","",$1); print}' |
    sort > .tmp_/tmp_id
    # -----------------------------------------------------
    # Extract ID, cuttind sites and cutting size
    # -----------------------------------------------------
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
    awk '{for(i=3;i<=NF-1;i++) if($i~"-") $i=length($i)-1;
        else if($i ~ /[a-z]/) $i=0 # $i="-"length($i)
        print}' |
    awk '{sum=-1; for(i=3;i<=NF-1;i++) sum+=$i
        print $1,$2+sum,$2+sum+$NF, $NF,sum}' \
    > .tmp_/tmp_anomaly_classification
    # =====================================================
    ## Clustering
    # =====================================================
    python ./DAJIN/src/test_anomaly_classification.py
    # =====================================================
    ## Output bam files and read statistics
    # =====================================================
    for cl in $(cat .tmp_/tmp_anomaly_label | awk '{print $NF}' | sort | uniq); do
        #
        cat .tmp_/tmp_anomaly_label |
        awk -v cl=${cl} '$NF==cl' |
        cut -f 1 |
        sort \
        > .tmp_/tmp_cl
        #
        cat .tmp_/tmp_cl |
        wc -l |
        sed "s/^/${barcode},${cl},/g" \
        >> ${output_dir}/${output_csv}
        #
        cl="cluster"${cl}
        #
        # -----------------------------------------------------
        # Generate BAM file
        # -----------------------------------------------------
        samtools view -h bam/${barcode}.bam |
        grep "^@" >.tmp_/tmp_header
        #
        samtools view bam/${barcode}.bam |
        awk '$2==0 || $2==16' |
        awk '{gsub("-.*","",$1); print}' |
        sort |
        join .tmp_/tmp_cl - 2>/dev/null |
        head -n 100 |
        sed "s/ /\t/g" \
        >> .tmp_/tmp_header
        #
        samtools sort .tmp_/tmp_header > ${output_dir}/${barcode}_${cl}.bam
        samtools index ${output_dir}/${barcode}_${cl}.bam
        #
        # =====================================================
        ## Sequence logo
        # =====================================================
        # -----------------------------------------------------
        # Extract mode of cutting site and cutting length
        # -----------------------------------------------------
        cat .tmp_/tmp_anomaly_label |
        awk -v cl=${cl} '$NF==cl' |
        cut -f 2 |
        sort |
        uniq -c |
        sort -nr |
        head -n 1 |
        awk '{print $2}' |
        > .tmp_/tmp_mode_start
        #
        cat .tmp_/tmp_anomaly_label |
        awk -v cl=${cl} '$NF==cl' |
        cut -f 4 |
        sort |
        uniq -c |
        sort -nr |
        head -n 1 |
        awk '{print $2}' |
        > .tmp_/tmp_mode_length

    done
    
    #
    samtools view bam/${barcode}.bam |
    awk '{gsub("-.*","",$1); print}' |
    sort -k 1,1 |
    join .tmp_/tmp_id - 2>/dev/null |
    awk '$2==0 || $2==16' |
    cut -d " " -f1,10 |
    sort |
    # join - .tmp_/tmp_anomaly_classification 2>/dev/null |
    join - .tmp_/tmp_anomaly_label |
    awk '{print $1, substr($2,$(NF-1)-15,30),$NF}' \
    > .tmp_/tmp_anomaly_seqlogo
    #
    for cl in $(cat .tmp_/tmp_anomaly_seqlogo | awk '{print $NF}' | sort | uniq); do
        cat .tmp_/tmp_anomaly_seqlogo |
        awk -v cl=${cl} '$NF==cl' |
        awk '{print ">"$1"\n"$2}' |
        clustalo -t DNA --auto -i - \
        > .tmp_/clustalo.fa 2>/dev/null
    # -----------------------------------------------------------------
    # REMOVE GAP
    # -----------------------------------------------------------------
    output_rmgap=$(echo .tmp_/clustalo.fa | sed -e "s/.fa/_rmgap.fa/g")
    # Extract gap-enriched nucreotide location
    true > .tmp_/remove_gaprow
    seqnum=$(cat .tmp_/clustalo.fa | awk -F "" '{if(NR==2) print length($0)}')
    for i in $(awk -v num=${seqnum} 'BEGIN{for(i=1;i<=num;i++) print i}'); do
        # echo "$i ==============="
        cat .tmp_/clustalo.fa |
        awk -F "" -v i=${i} '{if(NR%2==0) print $i}' |
        sort |
        uniq -c |
        awk -v i=${i} '{sum+=$1; if(max<$1) {max=$1; nuc=$2}}
        END{print i,nuc,max/sum*100}' |
        #Extract nucleotide position with gap "-" > 20%
        awk '$2 == "-" && $3>20' |
        cut -d " " -f 1 >> .tmp_/remove_gaprow
    done
    # Remove gap-enriched nucreotide location
    cat .tmp_/remove_gaprow |
    sed -e "s/^/\$/g" \
    -e 's/$/="";@/g' |
    tr -d "\n" |
    sed -e "s/@/ /g" \
    -e "s/^/{if(NR%2==0){/g" \
    -e "s/$/ print} else print}/g" \
    > .tmp_/remove_gap.awk
    #
    cat .tmp_/clustalo.fa |
    awk -F "" -f .tmp_/remove_gap.awk |
    sed "s/ //g" \
    > ${output_rmgap}
    #
    # Output sequence logo
    ## PNG
    { weblogo --title "${barcode}: Joint sequence" --scale-width no -n 50 --errorbars no -c classic --format png_print \
    < ${output_rmgap} > ${output_dir_seqlogo}/${barcode}_${cl}.png & } 1>/dev/null 2>/dev/null
    ## SVG
    { weblogo --title "${barcode}: Joint sequence" --scale-width no -n 50 --errorbars no -c classic --format svg \
    < ${output_rmgap} > ${output_dir_seqlogo}/${barcode}_${cl}.svg & } 1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null
    done
    


    cat .tmp_/tmp_anomaly_classification
done
rm .tmp_/tmp_*


output_dir="bam/abnormal_class/seqlogo"
mkdir -p ${output_dir}


exit 0

