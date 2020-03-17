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

# ============================================================================
# Arguments
# ============================================================================
genome=${1}
threads=${2}

parent_dir=".DAJIN_temp"
# ============================================================================
# GGGENOME
output_gggenome_temp="${parent_dir}"/data/gggenome.txt
output_gggenome_location="${parent_dir}"/data/gggenome_location 
output_reference="${parent_dir}"/data/ref.fa
# ============================================================================

true > "${output_gggenome_temp}"
# left flank
cat "${parent_dir}"/fasta_conv/wt.fa | sed 1d |
awk -v genome=${genome} '{seq=substr($0, 1, 50);
print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
sh -e | grep chr | cut -f 1-4 >> "${output_gggenome_temp}"

[ $(cat "${output_gggenome_temp}" | wc -l) -ne 1 ] &&
error_exit 1 '
No matched sequence found in reference genome:
Check and correct FASTA sequence and reference genome.'

# right flank
cat "${parent_dir}"/fasta_conv/wt.fa | sed 1d |
awk -v genome=${genome} '{seq=substr($0, length($0)-50, length($0));
print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
sh -e | grep chr | cut -f 1-4 >> "${output_gggenome_temp}"
####
[ $(cat "${output_gggenome_temp}" | wc -l) -ne 2 ] &&
error_exit 1 '
No matched sequence found in reference genome: 
Check FASTA sequence and reference genome.'

chromosome=$(cat "${output_gggenome_temp}" | head -n 1 | cut -f 1)
start=$(cat "${output_gggenome_temp}" | sort -k 3,3n | head -n 1 | cut -f 3)
end=$(cat "${output_gggenome_temp}" | sort -k 3,3nr | head -n 1 | cut -f 4)
strand=$(cat "${output_gggenome_temp}" | head -n 1 | cut -f 2)

printf "${chromosome}\t${start}\t${end}\t${strand}\n" \
> "${output_gggenome_location}" # this file will be used at "knockin search"

wget -q -O - http://togows.org/api/ucsc/${genome}/${chromosome}:${start}-${end}.fasta \
> "${output_reference}"
[ $(cat "${output_reference}" | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'

chrom_len=$(
    wget -q -O - http://hgdownload.cse.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.chrom.sizes |
    awk -v chrom=${chromosome} '$1 == chrom' | cut -f 2)
[ $(echo ${chrom_len} | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'

reference="${output_reference}"

rm ${output_gggenome_temp}
# ============================================================================
# Minimap2
output_bam_all=DAJIN_Report/bam
mkdir -p "${output_bam_all}"
# ============================================================================
for input in "${parent_dir}"/fasta_ont/*; do
    output=$(echo "${input}" | sed -e "s#.*/#${output_bam_all}/#g" -e "s/\.f.*$/.bam/g")
    echo "${output} is now generating..."
    ####
    minimap2 -t ${threads:-1} -ax map-ont --cs=long ${reference} ${input} 2>/dev/null |
    awk -v chrom="${chromosome}" -v chrom_len="${chrom_len}" -v start="${start}" \
    'BEGIN{OFS="\t"}
    $1=="@SQ" {$0="@SQ\tSN:"chrom"\tLN:"chrom_len; print}
    $1!~/^@/ {$3=chrom; $4=start+$4-1; print}' |
    samtools sort -@ ${threads:-1} - 2>/dev/null > "${output}"
    samtools index -@ ${threads:-1} "${output}"
done

# ============================================================================
# IGV.JS
output_igvjs=DAJIN_Report/igvjs
output_bam_100="${output_igvjs}"/bam_100reads
mkdir -p "${output_bam_100}"
# ============================================================================
read_num=100
for input in "${output_bam_all}"/*bam ; do
    output=$(echo ${input} | sed "s#.*/#${output_bam_100}/#g")
    # echo "${output} is now generating..."
    ####
    samtools view -h ${input} 2>/dev/null |
    awk '$1 ~ /^@/ || $2 != 4' |
    head -n $((${read_num}+5)) |
    samtools sort -@ ${threads:-1} - 2>/dev/null > "${output}"
    samtools index -@ ${threads:-1} "${output}"
done
# 
find "${output_bam_100}" | grep -e bam$ | sort | sed -e "s#^.*/bam#bam#g" -e 's#_#\\_#g' > "${parent_dir}"/data/tmp1_$$
find "${output_bam_100}" | grep -e bam$ | sort | sed -e "s#^.*/bam#bam#g" -e "s#.*/##g" -e "s#.bam##g" -e 's#_#\\_#g' > "${parent_dir}"/data/tmp2_$$
paste "${parent_dir}"/data/tmp1_$$ "${parent_dir}"/data/tmp2_$$ > "${parent_dir}"/data/igvjs_template.txt
rm "${parent_dir}"/data/tmp1_$$ "${parent_dir}"/data/tmp2_$$

./DAJIN/src/mojihame-l -l LABEL \
    DAJIN/src/igvjs_template.html \
    "${parent_dir}"/data/igvjs_template.txt |
sed -e "s/genome_info/${genome}/g" \
-e "s/locus_info/${chromosome}:${start}-${end}/g" \
> "${output_igvjs}"/igvjs.html

exit 0
