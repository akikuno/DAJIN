#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
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
ref_fa=.DAJIN_temp/data/ref.fa

bam_all=.DAJIN_temp/bam/
mkdir -p "${bam_all}"

#===============================================================================
#? Temporal
#===============================================================================

tmp_genome_location=.DAJIN_temp/data/tmp_genome_location

################################################################################
#! Server response
################################################################################

#===============================================================================
#? DAS server of UCSC Genome Browser
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

cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    awk -v genome=${genome} '{
        seq=substr($0, 1, 100)
    print "wget -q -O - \"https://genome.ucsc.edu/cgi-bin/hgBlat?db="genome"&type=BLAT%27s+guess&userSeq="seq"\""}' |
    sh - |
    grep "100.0%" | grep chr | awk '$6==100' |
    sed "s/.*chr/chr/g" |
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' |
cat > "${tmp_genome_location}"

#===============================================================================
#? Right flank
#===============================================================================

cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    awk -v genome=${genome} \
        '{seq=substr($0, length($0)-99, length($0));
    print "wget -q -O - \"https://genome.ucsc.edu/cgi-bin/hgBlat?db="genome"&type=BLAT%27s+guess&userSeq="seq"\""}' |
    sh - |
    grep "100.0%" | grep chr | awk '$6==100' |
    sed "s/.*chr/chr/g" |
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' |
cat >> "${tmp_genome_location}"

#===============================================================================
#? Assign genome location to variables
#===============================================================================

[ $(cat "${tmp_genome_location}" | wc -l) -ne 2 ] &&
error_exit 1 '
No matched sequence found in reference genome:
Check FASTA sequence and reference genome.'

chromosome=$(cat "${tmp_genome_location}" | head -n 1 | cut -f 1)
start=$(cat "${tmp_genome_location}" | sort -k 3,3n | head -n 1 | cut -f 3)
end=$(cat "${tmp_genome_location}" | sort -k 3,3nr | head -n 1 | cut -f 4)

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

rm ${tmp_genome_location}

exit 0
