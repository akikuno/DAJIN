#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
export PATH="$(command -p getconf PATH 2>/dev/null)${PATH+:}${PATH-}"
case $PATH in :*) PATH=${PATH#?};; esac
export UNIX_STD=2003  # to make HP-UX conform to POSIX


################################################################################
#! Define the functions for printing usage and error message
################################################################################

error_exit() {
    ${2+:} false && echo "${0##*/}: $2" 1>&2
    exit $1
}

################################################################################
#! I/O naming
################################################################################

#===============================================================================
#? TEST Aurguments
#===============================================================================

#===============================================================================
#? Aurguments
#===============================================================================
genome=${1}
threads=${2}

#===============================================================================
#? Input
#===============================================================================

#===============================================================================
#? Output
#===============================================================================
gggenome_location=.DAJIN_temp/data/gggenome_location
ref_fa=.DAJIN_temp/data/ref.fa

bam_all=.DAJIN_temp/bam/
mkdir -p "${bam_all}"

#===============================================================================
#? Temporal
#===============================================================================
tmp_gggenome=.DAJIN_temp/data/tmp_gggenome

################################################################################
#! Server response
################################################################################

#===============================================================================
#? GGGenome
#===============================================================================

[ $(wget -nv --spider --timeout 5 -t 1 https://gggenome.dbcls.jp/ja/ 2>&1 | grep -c '200 OK') -eq 0 ] &&
error_exit 1 'No connection to GGGenome server'


#===============================================================================
#? UCSC Genome Browser
#===============================================================================
url_ucsc_usa="http://genome.ucsc.edu/cgi-bin/das/"
url_ucsc_asia="http://genome-asia.ucsc.edu/cgi-bin/das/"
url_ucsc_euro="http://genome-euro.ucsc.edu/cgi-bin/das/"

if [ $(wget -nv --spider --timeout 5 -t 1 "${url_ucsc_usa}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_ucsc="${url_ucsc_usa}"
elif [ $(wget -nv --spider --timeout 5 -t 1 "${url_ucsc_asia}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_ucsc="${url_ucsc_asia}"
elif [ $(wget -nv --spider --timeout 5 -t 1 "${url_ucsc_euro}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_ucsc="${url_ucsc_euro}"
else
    error_exit 1 'Reference genome can not be obtained due to UCSC server error'
fi

#===============================================================================
#? GoldenPath of UCSC Genome Browser
#===============================================================================
url_golden_usa="http://hgdownload.cse.ucsc.edu/goldenPath"
url_golden_asia="http://hgdownload-asia.soe.ucsc.edu/goldenPath"
url_golden_euro="http://hgdownload-euro.soe.ucsc.edu/goldenPath"

if [ $(wget -nv --spider --timeout 5 -t 1 "${url_golden_usa}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_golden="${url_golden_usa}"
elif [ $(wget -nv --spider --timeout 5 -t 1 "${url_golden_asia}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_golden="${url_ucsc_aurl_golden_asiasia}"
elif [ $(wget -nv --spider --timeout 5 -t 1 "${url_golden_euro}" 2>&1 | grep -c '200 OK') -eq 1 ]; then
    url_golden="${url_golden_euro}"
else
    error_exit 1 'Reference genome can not be obtained due to UCSC server error'
fi

################################################################################
#! Obtain Genome coodinate from GGGnome
################################################################################

#===============================================================================
#? Left flank
#===============================================================================

cat .DAJIN_temp/fasta_conv/wt.fa |
    sed 1d |
    awk -v genome=${genome} '{seq=substr($0, 1, 50);
    print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
    sh -e |
    grep chr |
    cut -f 1-4 |
cat > "${tmp_gggenome}"

[ $(cat "${tmp_gggenome}" | wc -l) -ne 1 ] &&
error_exit 1 '
No matched sequence found in reference genome:
Check and correct FASTA sequence and reference genome.'

#===============================================================================
#? Right flank
#===============================================================================

cat .DAJIN_temp/fasta_conv/wt.fa |
    sed 1d |
    awk -v genome=${genome} \
        '{seq=substr($0, length($0)-50, length($0));
        print "wget -q -O - https://gggenome.dbcls.jp/ja/"genome"/"seq".txt"}' |
    sh -e |
    grep chr |
    cut -f 1-4 |
cat >> "${tmp_gggenome}"

[ $(cat "${tmp_gggenome}" | wc -l) -ne 2 ] &&
error_exit 1 '
No matched sequence found in reference genome: 
Check FASTA sequence and reference genome.'

chromosome=$(cat "${tmp_gggenome}" | head -n 1 | cut -f 1)
start=$(cat "${tmp_gggenome}" | sort -k 3,3n | head -n 1 | cut -f 3)
end=$(cat "${tmp_gggenome}" | sort -k 3,3nr | head -n 1 | cut -f 4)
strand=$(cat "${tmp_gggenome}" | head -n 1 | cut -f 2)

printf "${chromosome}\t${start}\t${end}\t${strand}\n" |
cat > "${gggenome_location}" # this file will be used at "knockin search"

################################################################################
#! Obtain Reference fasta file from UCSC Genome browser
################################################################################

url_ucsc_genome="${url_ucsc}/${genome}/dna?segment=${chromosome}:${start},${end}"

wget -qO - "${url_ucsc_genome}" |
    grep -v "^<" |
    awk '{print toupper($0)}' |
    sed -e "1i >${genome}:${chromosome}:${start}-${end}" |
cat > "${ref_fa}"

[ $(cat "${ref_fa}" | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'


################################################################################
#! Mapping for IGV
################################################################################

#===============================================================================
#? Rerefence Chromosome length
#===============================================================================
chrom_len=$(
    wget -q -O - "${url_golden}/${genome}/bigZips/${genome}.chrom.sizes" |
    awk -v chrom=${chromosome} '$1 == chrom' |
    cut -f 2)

[ $(echo ${chrom_len} | wc -l) -eq 0 ] &&
error_exit 1 'Invalid reference genome.'

#===============================================================================
#? Mapping all reads
#===============================================================================

reference="${ref_fa}"
for input in .DAJIN_temp/fasta_ont/*; do
    output=$(echo "${input}" |
        sed -e "s#.*/#${bam_all}/#g" \
            -e "s/\.f.*$/.bam/g")
    # echo "${output} is now generating..."
    ####
    minimap2 -t ${threads:-1} -ax map-ont --cs=long ${reference} ${input} 2>/dev/null |
        awk -v chrom="${chromosome}" -v chrom_len="${chrom_len}" -v start="${start}" \
        'BEGIN{OFS="\t"}
        $1~/@SQ/ {$0="@SQ\tSN:"chrom"\tLN:"chrom_len; print}
        $1!~/^@/ {$3=chrom; $4=start+$4-1; print}' |
        samtools sort -@ ${threads:-1} - 2>/dev/null |
    cat > "${output}"
    samtools index -@ ${threads:-1} "${output}"
done

# #===============================================================================
# #? 100 reads
# #===============================================================================
# for input in .DAJIN_temp/bam/*.bam; do
#     output=$(echo "${input}" |
#         sed "s#.DAJIN_temp/bam#.DAJIN_temp/bam/reads100#g")
#     # echo "${output} is now generating..."

#     header_num=$(samtools view -H -@ ${threads} ${input} | wc -l)
#     bam_num=$((100+${header_num}))

#     samtools view -h -@ ${threads} ${input} |
#         head -n "${bam_num}" |
#     samtools sort -@ ${threads} > "${output}"
#     samtools index -@ "${threads}" "${output}"
# done

# # ============================================================================
# # IGV.JS
# output_igvjs=.DAJIN_temp/bam/igvjs
# output_bam_100="${output_igvjs}"/bam_100reads
# mkdir -p "${output_bam_100}"
# # ============================================================================
# read_num=100
# for input in "${bam_all}"/*bam ; do
#     output=$(echo ${input} | sed "s#.*/#${output_bam_100}/#g")
#     # echo "${output} is now generating..."
#     ####
#     samtools view -h ${input} 2>/dev/null |
#     awk '$1 ~ /^@/ || $2 != 4' |
#     head -n $((${read_num}+5)) |
#     samtools sort -@ ${threads:-1} - 2>/dev/null > "${output}"
#     samtools index -@ ${threads:-1} "${output}"
# done
# # 
# find "${output_bam_100}" | grep -e bam$ | sort | sed -e "s#^.*/bam#bam#g" -e 's#_#\\_#g' > .DAJIN_temp/data/tmp1_$$
# find "${output_bam_100}" | grep -e bam$ | sort | sed -e "s#^.*/bam#bam#g" -e "s#.*/##g" -e "s#.bam##g" -e 's#_#\\_#g' > .DAJIN_temp/data/tmp2_$$
# paste .DAJIN_temp/data/tmp1_$$ .DAJIN_temp/data/tmp2_$$ > .DAJIN_temp/data/igvjs_template.txt
# rm .DAJIN_temp/data/tmp1_$$ .DAJIN_temp/data/tmp2_$$

# ./DAJIN/src/mojihame-l -l LABEL \
#     DAJIN/src/igvjs_template.html \
#     .DAJIN_temp/data/igvjs_template.txt |
#     sed -e "s/genome_info/${genome}/g" \
#         -e "s/locus_info/${chromosome}:${start}-${end}/g" |
# cat > "${output_igvjs}"/igvjs.html

rm ${tmp_gggenome}

exit 0
