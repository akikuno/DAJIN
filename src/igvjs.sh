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

# ============================================================================
# Arguments
# ============================================================================
genome=${1}
threads=${2}

# ============================================================================
# ターゲットのゲノム座標を入手する
# ============================================================================
tmp_gggenome=.DAJIN_temp/data/gggenome.txt
output_gggenome_location=.DAJIN_temp/data/gggenome_location 
# ----------------------------------------------------------------

# left flank
cat .DAJIN_temp/fasta_conv/wt.fa |
    sed 1d |
    awk -v genome=${genome} '{seq=substr($0, 1, 50);
    print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
    sh -e |
    grep chr |
    cut -f 1-4 |
cat - > "${tmp_gggenome}"

[ $(cat "${tmp_gggenome}" | wc -l) -ne 1 ] &&
error_exit 1 '
No matched sequence found in reference genome:
Check and correct FASTA sequence and reference genome.'

# right flank
cat .DAJIN_temp/fasta_conv/wt.fa |
    sed 1d |
    awk -v genome=${genome} \
        '{seq=substr($0, length($0)-50, length($0));
        print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
    sh -e |
    grep chr |
    cut -f 1-4 |
cat - >> "${tmp_gggenome}"

[ $(cat "${tmp_gggenome}" | wc -l) -ne 2 ] &&
error_exit 1 '
No matched sequence found in reference genome: 
Check FASTA sequence and reference genome.'

chromosome=$(cat "${tmp_gggenome}" | head -n 1 | cut -f 1)
start=$(cat "${tmp_gggenome}" | sort -k 3,3n | head -n 1 | cut -f 3)
end=$(cat "${tmp_gggenome}" | sort -k 3,3nr | head -n 1 | cut -f 4)
strand=$(cat "${tmp_gggenome}" | head -n 1 | cut -f 2)

rm ${tmp_gggenome}


printf "${chromosome}\t${start}\t${end}\t${strand}\n" \
> "${output_gggenome_location}" # this file will be used at "knockin search"

# ============================================================================
# ターゲットのゲノム座標を入手する
# ============================================================================
output_reference=.DAJIN_temp/data/ref.fa
# ----------------------------------------------------------------

# ----------------------------------------------------------------------------
# Targetが一塩基変異の場合: 
# Cas9の切断部に対してgRNA部自体の欠損およびgRNA長分の塩基挿入したものを異常アレルとして作成する
# ----------------------------------------------------------------------------

# Rerefence FASTA file
url_ucsc_usa="http://genome.ucsc.edu/cgi-bin/das/${genome}/dna?segment=${chromosome}:${start},${end}"
url_ucsc_asia="http://genome-asia.ucsc.edu/cgi-bin/das/${genome}/dna?segment=${chromosome}:${start},${end}"
url_ucsc_euro="http://genome-euro.ucsc.edu/cgi-bin/das/${genome}/dna?segment=${chromosome}:${start},${end}"

if [ $(wget -nv --spider --timeout 5 -t 1 "${url_ucsc_usa}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_ucsc="${url_ucsc_usa}"
elif [ $(wget -nv --spider --timeout 5 -t 1 "${url_ucsc_asia}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_ucsc="${url_ucsc_asia}"
elif [ $(wget -nv --spider --timeout 5 -t 1 "${url_ucsc_euro}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_ucsc="${url_ucsc_euro}"
else
    error_exit 1 'Reference genome can not be obtained due to UCSC server error'
fi

wget -qO - "${url_ucsc}" |
    grep -v "^<" |
    awk '{print toupper($0)}' |
    sed -e "1i >${genome}:${chromosome}:${start}-${end}" |
cat - >> "${output_reference}"

[ $(cat "${output_reference}" | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'

# Rerefence Chromosome length  ---------------------------------------------------
chrom_len=$(
    wget -q -O - http://hgdownload.cse.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.chrom.sizes |
    awk -v chrom=${chromosome} '$1 == chrom' |
    cut -f 2)

[ $(echo ${chrom_len} | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'

# ============================================================================
# Minimap2
output_bam_all=.DAJIN_temp/bam/
mkdir -p "${output_bam_all}"
# ============================================================================
reference="${output_reference}"
for input in .DAJIN_temp/fasta_ont/*; do
    output=$(echo "${input}" |
        sed -e "s#.*/#${output_bam_all}/#g" \
            -e "s/\.f.*$/.bam/g")
    # echo "${output} is now generating..."
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
output_igvjs=.DAJIN_temp/bam/igvjs
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
find "${output_bam_100}" | grep -e bam$ | sort | sed -e "s#^.*/bam#bam#g" -e 's#_#\\_#g' > .DAJIN_temp/data/tmp1_$$
find "${output_bam_100}" | grep -e bam$ | sort | sed -e "s#^.*/bam#bam#g" -e "s#.*/##g" -e "s#.bam##g" -e 's#_#\\_#g' > .DAJIN_temp/data/tmp2_$$
paste .DAJIN_temp/data/tmp1_$$ .DAJIN_temp/data/tmp2_$$ > .DAJIN_temp/data/igvjs_template.txt
rm .DAJIN_temp/data/tmp1_$$ .DAJIN_temp/data/tmp2_$$

./DAJIN/src/mojihame-l -l LABEL \
    DAJIN/src/igvjs_template.html \
    .DAJIN_temp/data/igvjs_template.txt |
    sed -e "s/genome_info/${genome}/g" \
        -e "s/locus_info/${chromosome}:${start}-${end}/g" |
cat - > "${output_igvjs}"/igvjs.html

exit 0
