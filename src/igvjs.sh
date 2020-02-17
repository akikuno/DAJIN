#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================

set -u
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ======================================
# Define the functions for printing usage and error message
# ======================================

error_exit() {
    ${2+:} false && echo "${0##*/}: $2" 1>&2
    exit $1
}
warning() {
    ${1+:} false && echo "${0##*/}: $1" 1>&2
}

# ---------------------------------------------------------------------

genome=${1}
threads=${2}

true > .tmp_/gggenome.txt
# left flank
cat .tmp_/wt.fa | sed 1d |
awk -v genome=${genome} '{seq=substr($0, 1, 50);
print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
sh -e | grep chr | cut -f 1-4 >> .tmp_/gggenome.txt

[ $(cat .tmp_/gggenome.txt | wc -l) -ne 1 ] &&
error_exit 1 'No matched sequence found in reference genome.
Check and correct FASTA sequence and reference genome.'

# right flank
cat .tmp_/wt.fa | sed 1d |
awk -v genome=${genome} '{seq=substr($0, length($0)-50, length($0));
print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
sh -e | grep chr | cut -f 1-4 >> .tmp_/gggenome.txt 
[ $(cat .tmp_/gggenome.txt | wc -l) -ne 2 ] &&
error_exit 1 '
No matched sequence found in reference genome: 
Check FASTA sequence and reference genome.'

chromosome=$(cat .tmp_/gggenome.txt | sort -k 3,3n | sed -n 1p | cut -f 1)
start=$(cat .tmp_/gggenome.txt | sort -k 3,3n | sed -n 1p | cut -f 3)
end=$(cat .tmp_/gggenome.txt | sort -k 3,3nr | sed -n 1p | cut -f 4)
strand=$(cat .tmp_/gggenome.txt | sort -k 3,3nr | sed -n 1p | cut -f 2)

printf "${chromosome}\t${start}\t${end}\t${strand}\n" \
> .tmp_/gggenome_location # this file will be used at "knockin search"

wget -q -O - http://togows.org/api/ucsc/${genome}/${chromosome}:${start}-${end}.fasta > .tmp_/ref.fa
[ $(cat .tmp_/ref.fa | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'

chrom_len=$(
    wget -q -O - http://hgdownload.cse.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.chrom.sizes |
    awk -v chrom=${chromosome} '$1 == chrom' | cut -f 2)
[ $(echo ${chrom_len} | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'

reference=.tmp_/ref.fa

#* export alingments in "bam" dir
for input in fasta_ont/*; do
    output=$(echo ${input} |
        sed -e "s#\..*#.bam#g" \
        -e "s#.*/#bam/#g")
    echo "${output} is now generating..."
    ####
    minimap2 -t ${threads:-1} -ax splice --cs=long ${reference} ${input} 2>/dev/null |
    # header: replace chromosome number and length
    awk -v chrom=${chromosome} -v chrom_len=${chrom_len} '{
    if($1=="@SQ") print "@SQ\tSN:"chrom"\tLN:"chrom_len
    else print 
    }' |
    # main: replace chromosome number and start site
    awk -v chrom=${chromosome} -v start=${start} 'BEGIN{OFS="\t"}{
    if($1!~/^@/) {$3=chrom; $4=start+$4-1; print}
    else print
    }' |
    samtools sort -@ ${threads:-1} - 2>/dev/null > ${output}
    samtools index -@ ${threads:-1} ${output}
done

#* export 100 reads in "igvjs/bam_100reads" dir
read_num=100

mkdir -p bam/bam_${read_num}reads results/igvjs/

for input in $(ls bam/*bam) ; do
    output=$(echo ${input} |
    sed -e "s#.*/#bam/bam_${read_num}reads/#g")
    # echo "${output} is now generating..."
    ####
    samtools view -h ${input} | awk '$1 ~ /^@/ || $2 == 0 || $2 == 16'| head -n $((${read_num}+5)) |
    samtools sort -@ ${threads:-1} - 2>/dev/null > ${output}
    samtools index -@ ${threads:-1} ${output}
done
cp -rf bam/bam_${read_num}reads results/igvjs/
# 
ls results/igvjs/bam_${read_num}reads/* | grep -e bam$ | sed -e "s#^.*/bam#bam#g" -e 's#_#\\_#g' > .tmp_/tmp1
ls results/igvjs/bam_${read_num}reads/* | grep -e bam$ | sed -e "s#^.*/bam#bam#g" -e "s#.*/##g" -e "s#.bam##g" -e 's#_#\\_#g' > .tmp_/tmp2
paste .tmp_/tmp1 .tmp_/tmp2 > .tmp_/igvjs_template.txt
rm .tmp_/tmp1 .tmp_/tmp2

./DAJIN/src/mojihame-l -l LABEL DAJIN/src/igvjs_template.html .tmp_/igvjs_template.txt |
sed -e "s/genome_info/${genome}/g" \
-e "s/locus_info/${chromosome}:${start}-${end}/g" \
> results/igvjs/igvjs.html


exit 0
