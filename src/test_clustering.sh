#!/bin/sh
# cat data_for_ml/sequence_MIDS.txt |
# grep "wt" |
# cut -f 1,3 > test_labels.txt

# cat data_for_ml/sequence_MIDS.txt |
# grep "wt" |
# cut -f 2 |
# awk -F "" '{
#     for(i=1; i<=NF; i++){
#         if($i=="I" && $(i+1)!="I") $(i+1)="@"
#         }
#     print $0}' |
# sed -e "s/I//g" -e "s/ //g" \
# > tmp_seq

# paste test_labels.txt tmp_seq \
# > test_tmp_seq


zcat data_for_ml/sequence_MIDS.txt.gz |
grep "barcode02" |
join -1 3 -2 2 - .tmp_/sorted_prediction_result |
cut -d " " -f 1,2 |
sed "s/ /\t/g" > test_labels.txt

zcat data_for_ml/sequence_MIDS.txt.gz |
grep "barcode02" |
join -1 3 -2 2 - .tmp_/sorted_prediction_result |
cut -d " " -f 3 |
awk -F "" '{
    for(i=1; i<=NF; i++){
        if($i=="I" && $(i+1)!="I") $(i+1)="@"
        }
    print $0}' |
sed -e "s/I//g" -e "s/ //g" \
> tmp_seq

input="tmp_seq"
output="test.csv"
#
seqnum=$(cat ${input} |
    grep -v "^$" |
    awk -F "" 'BEGIN{min="inf"} \
    {if(min>length($0)) min=length($0)} \
    END{print min}')
#
cat ${input} |
awk -F "" -v seqnum=${seqnum} \
    '{for(i=1;i<=seqnum;i++) {
    seq[i]=seq[i]$i
}}
END{for(key in seq) print seq[key]}' |
awk -F "" 'BEGIN{OFS=","}{
    numGap=gsub("=","=",$0)
    numM=gsub("M","M",$0)
    numI=gsub("@","@",$0)
    numD=gsub("D","D",$0)
    numS=gsub("S","S",$0)
    num_result=""
    for(i=1; i<=NF; i++){
        if($i=="=") $i=0
        else if($i=="M") $i=0
        else if($i=="@") $i=numI
        else if($i=="D") $i=numD
        else if($i=="S") $i=numD
    }
    # print $0,numGap,numM,numI,numD,numS
    print $0}' \
> ${output}


# HDBSCAN ------------------

Rscript clustering.R
# cat fasta_ont/* > fasta_ont/all_merged.fa
# cat fasta_ont/wt* > fasta_ont/wt_merged.fa
# ./DAJIN/src/igvjs.sh ${genome:-mm10} ${threads:-1}

# Separate BAM files
barcode=barcode02
cat test_result.txt | cut -f 2 | sort | uniq -c
#
rm ${barcode}*.bam* 2>/dev/null
for i in $(cat test_result.txt | cut -f 2 | sort -u);do
    cat test_result.txt | grep ${i}$ | cut -f 1 | sort > tmp_id
    samtools view -h bam/${barcode}.bam | grep "^@" > tmp_header 
    #
    samtools view bam/${barcode}.bam |
    sort |
    join - tmp_id |
    sed "s/ /\t/g" |
    head -n 100 \
    >> tmp_header
    #
    samtools sort tmp_header > ${barcode}_cl${i}.bam
    samtools index ${barcode}_cl${i}.bam
done

samtools view -h bam/${barcode}.bam | grep "^@" > tmp_header 
#
samtools view bam/${barcode}.bam |
sort |
join - tmp_id |
sed "s/ /\t/g" |
head -n 100 \
>> tmp_header
#
samtools sort tmp_header > ${barcode}_all.bam
samtools index ${barcode}_all.bam

rm tmp_*
