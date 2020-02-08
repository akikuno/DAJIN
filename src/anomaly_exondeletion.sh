#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

genome=mm10
genome=${1}
wget -qO - https://hgdownload.soe.ucsc.edu/goldenPath/${genome}/database/refFlat.txt.gz |
gzip -dc > tmp.txt
cat tmp.txt | awk 'BEGIN{OFS="\t"}{print $3,$5,$6,$0}' |
bedtools intersect -a - -b .tmp_/gggenome_location -wb > tmp2
cat tmp2 | cut -f 13-14 > tmp3
chr=$(cat tmp2 | cut -f 1)
cat tmp3 | cut -f 1 | sed "s/,/\n/g" | grep -v "^$" > tmp4
cat tmp3 | cut -f 2 | sed "s/,/\n/g" | grep -v "^$" > tmp5
paste tmp4 tmp5 | grep -v "^$" | sed "s/^/${chr}\t/g" | sort -k 1,1 -k 2,2n |
bedtools intersect -a - -b .tmp_/gggenome_location

genePredToBed tmp.txt out.bed

wget -qO - http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/knownGene.txt.gz |
gzip -dc > tmp2.txt
cat tmp2.txt | awk -F"\t" 'BEGIN{OFS="\t";} {$11 = 0; print;}' > tmp3.txt
genePredToBed tmp3.txt tmp_hoge.bed
# genePredToGtf file stdin hg38.knownGene.gtf
