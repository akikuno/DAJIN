#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# genome=mm10
genome=${1}
threads=${2:-1}

# ======================================
# Identify Target Exon
# ======================================

minimap2 -t ${threads:-1} -ax splice --cs \
.tmp_/ref.fa .tmp_/target.fa 2>/dev/null |
grep -v "^@" |
awk '{print $3, $(NF-1)}' |
sed -e "s/cs:Z:://g" |
awk -F ":" '{print $2,$3}' |
sed -e "s/-/ /g" -e "s/~[a-z]*/ /g" -e "s/[a-z]*$//g" |
awk 'BEGIN{OFS="\t"}{print $1,$2+$4,$2+$4+$5}' \
> .tmp_/deletion_location.bed
#
wget -qO - ftp://hgdownload.soe.ucsc.edu/goldenPath/${genome}/database/refFlat.txt.gz |
gzip -dc |
awk 'BEGIN{OFS="\t"}{print $3,$5,$6,$0}' |
bedtools intersect -a - -b .tmp_/gggenome_location -wb |
cut -f 13-14 > .tmp_/tmp_1
#
cat .tmp_/tmp_1 | cut -f 1 | sed "s/,/\n/g" | grep -v "^$" > .tmp_/tmp_2
cat .tmp_/tmp_1 | cut -f 2 | sed "s/,/\n/g" | grep -v "^$" > .tmp_/tmp_3
#
chr=$(cat .tmp_/deletion_location.bed | cut -f 1 | sort | uniq)
#
paste .tmp_/tmp_2 .tmp_/tmp_3 | grep -v "^$" | sed "s/^/${chr}\t/g" | sort -k 1,1 -k 2,2n | uniq |
bedtools intersect -a - -b .tmp_/deletion_location.bed \
> .tmp_/target_exon.bed
#
start=$(cat .tmp_/target_exon.bed | cut -f 2)
end=$(cat .tmp_/target_exon.bed | cut -f 3)
#
paste .tmp_/tmp_2 .tmp_/tmp_3 | grep -v "^$" | sed "s/^/${chr}\t/g" | sort -k 1,1 -k 2,2n | uniq |
bedtools intersect -a - -b .tmp_/gggenome_location |
awk -v start=${start} -v end=${end} \
    '{if($2==start && $3==end) print $0"\ttarget"; else print $0"\tothers"}' \
> .tmp_/all_exons.bed
rm .tmp_/tmp*

# ======================================
# Output "Problematic anomaly"
# (Target exon was not cutted, or Other exons were deleted)
# ======================================

cat .tmp_/anomaly_classification.txt |
sort -k 2,2  \
> .tmp_/anomaly_classification_sorted.txt

output_anomaly=".tmp_/anomaly_nonproblematic.txt"
true > ${output_anomaly}
barcode=barcode27
for barcode in $(cut -f 1 .tmp_/abnormal_sequenceids.txt | sort | uniq); do
    echo ${barcode}
    samtools view bam/${barcode}.bam |
    # grep  525565d | # !================================
    sort |
    join -1 1 -2 2 - .tmp_/anomaly_classification_sorted.txt |
    grep Abnormal |
    sed "s/ /\t/g" |
    awk '$2==0||$2==16' \
    > .tmp_/abnormal.sam
    #
    cat .tmp_/abnormal.sam  |
    awk '{for(i=1;i<=NF;i++) if($i ~ "cs:Z::") print $1, $i}' |
    sed -e "s/cs:Z:://g" -e "s/~[a-z]*\([0-9]*\)/ \1\t/" |
    cut -f 1 |
    sed -e "s/:/ /g" -e "s/\*[a-z][a-z]/+1/g" |
    # Remove inversion reads
    awk '{seq=$(NF-1);
        if (seq !~ /\+[a-z]+$/ ) print }' | 
    #
    awk '{id=$1; $1="ID";
        gsub("-"," - ",$0); gsub("+"," + ",$0);
        $1=id; print}' |
    awk '{for(i=2;i<=NF-1;i++) if($i ~ /[a-z]/) $i=length($i); print}' |
    awk '{id=$1; $1="ID";
        gsub("- "," -",$0); gsub("+ "," ",$0);
        $1=id; print}' |
    awk '{sum=1; for(i=2;i<=NF-1;i++) sum+=$i; print $1, sum, $NF}' |
    sort \
    > .tmp_/abnormal_deletionsite
    #
    cat .tmp_/abnormal.sam |
    cut -f 1,3,4 |
    join - .tmp_/abnormal_deletionsite |
    awk 'BEGIN{OFS="\t"}{print $2, $3+$4, $3+$4+$5, $1}' |
    sort -k 1,1 -k 2,2n |
    bedtools intersect -a - -b .tmp_/all_exons.bed -wa -wb |
    awk -v start=${start} -v end=${end} \
    '$2<=start && $3>=end' |
    cut -f 4,8 |
    awk 'BEGIN{OFS="\t"}{print $1,1,100,$2}' |
    sort -k 1,1 -k 2,2n |
    bedtools merge -i - -c 4 -o distinct |
    grep -v other |
    cut -f 1 |
    sed -e "s/^/${barcode}\t/g" \
        -e "s/$/\tAbnormal(target_deletion)/g" \
    >> ${output_anomaly}
done

cat ${output_anomaly} |
sort -k 2,2 |
join -a 1 -1 2 -2 2 .tmp_/anomaly_classification_sorted.txt - |
# grep 01a16fea-5549-41e3-aa4d-160c68280d5d |
awk 'BEGIN{OFS="\t"}{if(NF==5){$3=$5}; print $2,$1,$3}' |
sort \
> .tmp_/anomaly_classification_revised.txt
## confirmation

# cat .tmp_/anomaly_classification.txt |
# grep ${barcode} | cut -f 2 | sort > tmp

# echo ${barcode}
# samtools view -h bam/${barcode}.bam | grep "^@" > .tmp_/header
# samtools view bam/${barcode}.bam |
# sort |
# join - tmp |
# awk '$2=="0" || $2=="16"' |
# sed "s/ /\t/g" >> .tmp_/header
# samtools sort .tmp_/header > .tmp_/bam_abnormal/${barcode}_nonproblem.bam
# samtools index .tmp_/bam_abnormal/${barcode}_nonproblem.bam
# wc -l .tmp_/header
# rm tmp .tmp_/header
