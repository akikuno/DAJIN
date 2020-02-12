#!/bin/sh

# =====================================================
# Generate test data
# =====================================================
printf "
02bbb775\t0
042c068f\t0
0324be19\t0
07bececa\t1
0df976dc\t1
05f9d0e3\t1
0f838c72\t2
0d29b0af\t2
07fa3b6f\t2
" |
sort |
grep -v "^$" \
> tmp_testdata

# ========================================
cat .tmp_/anomaly_classification_revised.txt |
grep ${barcode} | grep target_deletion |
cut -f 2 |
awk '{gsub("\-.*","",$1); print}' |
sort > tmpid

samtools view bam/${barcode}.bam |
awk '{gsub("\-.*","",$1); print}' |
sort |
join tmpid - > tmp

cat tmp | #grep dd3a08ee |
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
> tmp_testdata2

# --------------------
barcode=barcode12
samtools view -h bam/${barcode}.bam |
grep "^@" > tmp_header

samtools view bam/${barcode}.bam |
awk '{gsub("\-.*","",$1); print}' |
sort |
join tmp_testdata - |
sed "s/ /\t/g" |
cut -f 1,3- \
>> tmp_header
samtools sort tmp_header > tmp.bam
samtools index tmp.bam

rm tmp*
for cl in $(cat hoge.txt | awk '{print $NF}' | sort | uniq); do
    echo $cl
    cat hoge.txt |
    awk -v cl=${cl} '$NF==cl' |
    cut -f 1 |
    sort \
    > tmp_${cl}
done

barcode=barcode12
samtools view -h bam/${barcode}.bam |
grep "^@" > tmp_header

samtools view bam/${barcode}.bam |
awk '{gsub("\-.*","",$1); print}' |
sort |
join tmp_-1 - |
head -n 100 |
sed "s/ /\t/g" \
>> tmp_header

samtools sort tmp_header > tmp.bam
samtools index tmp.bam