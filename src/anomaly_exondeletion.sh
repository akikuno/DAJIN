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

output_anomaly=".tmp_/anomaly_classification.txt"
true > ${output_anomaly}
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


chr=$(cat .tmp_/deletion_location.bed | cut -f 1 | sort | uniq)

wget -qO - https://hgdownload.soe.ucsc.edu/goldenPath/${genome}/database/refFlat.txt.gz |
gzip -dc |
awk 'BEGIN{OFS="\t"}{print $3,$5,$6,$0}' |
bedtools intersect -a - -b .tmp_/gggenome_location -wb |
cut -f 13-14 > .tmp_/tmp_1
#
cat .tmp_/tmp_1 | cut -f 1 | sed "s/,/\n/g" | grep -v "^$" > .tmp_/tmp_2
cat .tmp_/tmp_1 | cut -f 2 | sed "s/,/\n/g" | grep -v "^$" > .tmp_/tmp_3
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

cat .tmp_/abnormal_sequenceids.txt |
cut -f 1,2 |
sort -k 2,2  \
> .tmp_/abnormal_seqid_sorted.txt

barcode=barcode17
for barcode in $(cut -f 1 .tmp_/abnormal_sequenceids.txt | sort | uniq); do
    echo ${barcode}
    samtools view bam/${barcode}.bam |
    sort |
    join -1 1 -2 2 - .tmp_/abnormal_seqid_sorted.txt |
    sed "s/ /\t/g" |
    awk '$2==0||$2==16' > .tmp_/abnormal.sam
    #
    cat .tmp_/abnormal.sam  | # grep 91606edc | # ! ==================
    awk '{print $(NF-2)}' |
    sed -e "s/cs:Z:://g" -e "s/~[a-z]*\([0-9]*\)/ \1\t/" |
    cut -f 1 |
    sed -e "s/:/ /g" -e "s/\*[a-z][a-z]/+1/g" -e "s/\+/ /g" | #  |
    # awk '{for(i=1;i<=NF;i++) {gsub("+[a-z]*","",$i)}; print}' |
    sed -e "s/\-/ - /g" |
    awk '{for(i=1;i<=NF;i++) if($i ~ /[a-z]/) $i=length($i); print}' |
    sed "s/\- /\-/g" |
    awk '{sum=1; for(i=1;i<=NF-1;i++) sum+=$i; print sum, $NF}' \
    > .tmp_/abnormal_deletionsite
    #
    cat .tmp_/abnormal.sam | # grep 91606edc | # ! ==================
    cut -f 1,3,4 |
    paste - .tmp_/abnormal_deletionsite |
    awk 'BEGIN{OFS="\t"}{print $2, $3+$4, $3+$4+$5, $1}' |
    sort -k 1,1 -k 2,2n |
    bedtools intersect -a - -b .tmp_/all_exons.bed -wb |
    awk -v start=${start} -v end=${end} \
    '$2==start && $3==end' |
    cut -f 4,8 |
    awk 'BEGIN{OFS="\t"}{print $1,1,100,$2}' |
    sort -k 1,1 -k 2,2n |
    bedtools merge -i - -c 4 -o distinct |
    grep -v other |
    cut -f 1 |
    sort |
    sed -e "s/^/${barcode}\t/g" \
        -e "s/$/\tNon_problematic_anomaly/g" \
    >> ${output_anomaly}
done


#確認

cat .tmp_/anomaly_classification.txt |
grep ${barcode} | cut -f 2 | sort > tmp

echo ${barcode}
samtools view -h bam/${barcode}.bam | grep "^@" > .tmp_/header
samtools view bam/${barcode}.bam |
sort |
join - tmp |
sed "s/ /\t/g" >> .tmp_/header
samtools sort .tmp_/header > .tmp_/bam_abnormal/${barcode}_nonproblem.bam
samtools index .tmp_/bam_abnormal/${barcode}_nonproblem.bam
rm tmp .tmp_/header

# ======================================
# Extract abnormal reads
# ======================================




mkdir -p .tmp_/bam_abnormal
cat .tmp_/abnormal_sequenceids.txt |
cut -f 1,2 |
sort -k 2,2  \
> .tmp_/abnormal_seqid_sorted.txt
for barcode in $(cut -f 1 .tmp_/abnormal_sequenceids.txt | sort | uniq); do
    echo ${barcode}
    samtools view -h bam/${barcode}.bam | grep "^@" > .tmp_/header
    samtools view bam/${barcode}.bam |
    sort |
    join - .tmp_/abnormal_seqid_sorted.txt |
    sed "s/ /\t/g" >> .tmp_/header
    samtools sort .tmp_/header > .tmp_/bam_abnormal/${barcode}.bam
    samtools index .tmp_/bam_abnormal/${barcode}.bam
    rm .tmp_/header
done

    # ----------------------------------------
    # None-exon cutting
    # ----------------------------------------
    cat tmp_$$ |
    cut -f 4 |
    sort > tmp2_$$
    #
    cat .tmp_/abnormal.sam | grep cbd69deb -A100 |
    cut -f 1 |
    sort |
    join -v 1 - tmp2_$$ |
    sed -e "s/^/${barcode}\t/g" \
        -e "s/$/\tProblematic_abnormal/g" \
    >> ${output_anomaly}
    # ----------------------------------------
    # Cutting Other Exons, retain Target Exon
    # ----------------------------------------
    cat tmp_$$ |
    bedtools intersect -v -a - -b .tmp_/target_exon.bed |  
    cut -f 4 |
    sort > tmp3_$$
    cat .tmp_/abnormal.sam | grep cbd69deb -A100 |
    cut -f 1 |
    sort |
    join - tmp3_$$ |
    sed -e "s/^/${barcode}\t/g" \
        -e "s/$/\tProblematic_abnormal/g" \
    >> ${output_anomaly}

    # ----------------------------------------
    # Cutting Target Exons AND Other Exons
    # ----------------------------------------


    cat .tmp_/abnormal.sam | grep cbd69deb -A100 |
    cut -f 1 |
    sort |
    join -v 1 - tmp |
    sed -e "s/^/${barcode}\t/g" \
        -e "s/$/\tProblematic_abnormal/g" 

    cat tmp |
    awk -v num=${target_exon_num} '$2!=num' |

    sort |
    join head
