#!/bin/sh

# ============================================================================
# Initialize shell environment
# ============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# For Windows Subsystem for Linux
uname -a | grep Microsoft 1>/dev/null 2>/dev/null &&
alias python="python.exe"


# ============================================================================
# Define the functions for printing usage and error message
# ============================================================================

print_usage_and_exit () {
cat <<-USAGE 1>&2

Usage     : ${0##*/} [options] \\
            -i [FASTA] \\
            -ont_dir [FASTA|FASTQ]  \\
            -ont_ref [FASTA|FASTQ] \\
            -genome [reference genome]
            -name [string] (optional) \\
            -seq [string] (optional) \\
            -t [integer] (optional)

Example   : ./DAJIN/${0##*/} \\
            -i DAJIN/example/cables2_flox.fa \\
            -ont DAJIN/example/demultiplex \\
            -ont_ref DAJIN/example/demultiplex/barcode21.fastq.gz \\
            -genome mm10 \\
            -seq ATAACTTCGTATAATGTATGCTATACGAAGTTAT \\
            -t 8

Options :   -i         : Multi-FASTA file. It must includes ">target" and ">wt"
            -ont_dir   : Directory of ONT demultiplexed reads
            -ont_ref   : Reference (or wild-type) reads from ONT MinION
            -genome    : Reference genome (e.g. hg19, mm10):
                         See https://gggenome.dbcls.jp/en/help.html#db_list
            -name         : Output file name :default=sequence
            -seq       : Target knock-in sequence
	            -t         : Number of threads: default=4

USAGE
    exit 1
}

error_exit() {
    ${2+:} false && echo "${0##*/}: $2" 1>&2
    exit $1
}

# ============================================================================
# Prerequisit
# ============================================================================
set +e
type python 1>/dev/null 2>/dev/null
[ "$?" -eq 1 ] && error_exit 1 'Please install Python'
type git 1>/dev/null 2>/dev/null
[ "$?" -eq 1 ] && error_exit 1 'Please install git'
type samtools 1>/dev/null 2>/dev/null
[ "$?" -eq 1 ] &&  error_exit 1 'Please install samtools'
type minimap2 1>/dev/null 2>/dev/null
[ "$?" -eq 1 ] &&  error_exit 1 'Please install minimap2'
type gzip 1>/dev/null 2>/dev/null
[ "$?" -eq 1 ] && error_exit 1 'Please install gzip'
#
python -c \
"from tensorflow.python.client import device_lib;
print(device_lib.list_local_devices())" \
1>/dev/null 2>/dev/null
[ "$?" -eq 1 ] &&  error_exit 1 'Please install Tensorflow'
set -e

# ============================================================================
# Parse arguments
# ============================================================================

[ $# -eq 0 ] && print_usage_and_exit

option_cnt=""
while [ $# -gt 0 ]; do
    option_cnt=${option_cnt}${1}
    case "$1" in
        --) shift
            break
        ;;
        -[hv]|--help|--version)
            print_usage_and_exit
        ;;
        -i)  [ ! -s "${2}" ] && error_exit 1 'Invalid -i option: Empty file'
            fasta=${2}
        ;;
        -ont_dir) [ ! -d "${2}" ] && error_exit 1 'Invalid -ont_dir option: Directory not found'
            ont_dir=${2}
        ;;
        -ont_ref) [ ! -s "${2}" ] && error_exit 1 'Invalid -ont_ref option: Empty file'
            ont_ref=${2}
        ;;
        -genome) [ ! -n "${2}" ] && error_exit 1 'Invalid -genome option'
            genome=${2}
        ;;
        -name) [ $(echo "${2}" | awk '{print NF}') -ne 1 ] && error_exit 1 'Invalid -name option'
            output_file=${2}
        ;;
        -t) [ ! "${2}" -eq "${2}" ] && error_exit 1 'Invalid -t option'
            threads=${2}
        ;;
    esac
    shift
done

set +e
echo $option_cnt | grep "\-i" > /dev/null ||
error_exit 1 "FASTA file is required: See ${0##*/} -h"

echo $option_cnt | grep "\-ont" > /dev/null ||
error_exit 1 "ONT read file (FASTA or FASTQ) is required: See ${0##*/} -h"

echo $option_cnt | grep "\-ont_ref" > /dev/null ||
error_exit 1 "ONT reads from wild-type is required: See ${0##*/} -h"

echo $option_cnt | grep "\-genome" > /dev/null ||
error_exit 1 "Reference genome should be specified: See ${0##*/} -h"

if [ "$(grep -c '>target' ${fasta})" -eq 0 ] || [ "$(grep -c '>wt' ${fasta})" -eq 0 ]; then
    error_exit 1 "FASTA requires including \">target\" and \">wt\" headers. See ./${0##*/} -h"
fi
set -e

# ##################################################
# Allocate threads
# ##################################################
set +u
# Linux and similar...
[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/1.5+0.5)}')
# FreeBSD and similar...
[ -z "$threads" ] && threads=$(getconf NPROCESSORS_ONLN | awk '{print int($0/1.5)+0.5}')
# Solaris and similar...
[ -z "$threads" ] && threads=$(ksh93 -c 'getconf NPROCESSORS_ONLN' | awk '{print int($0/1.5+0.5)}')
# Give up...
[ -z "$threads" ] && threads=1
set -u

# ============================================================================
# Setting Directory
# ============================================================================
rm -rf ".DAJIN_temp" 2>/dev/null
dirs="fasta fasta_conv fasta_ont NanoSim data \
    results/svg results/png results/igvjs"

echo "$dirs" | sed "s/ /\n/g" |
while read -r dir; do
    mkdir -p ".DAJIN_temp/$dir"
done

# mkdir -p DAJIN_results/bam 
# ============================================================================
# Format FASTA file
# ============================================================================

# CRLF to LF
cat "${fasta}" |
tr -d "\r" |
grep -v "^$" \
> .DAJIN_temp/fasta/fasta.fa
fasta_LF=".DAJIN_temp/fasta/fasta.fa"

# Separate multiple-FASTA into FASTA files
cat ${fasta_LF} |
sed "s/^/@/g" |
tr -d "\n" |
sed -e "s/@>/\n>/g" -e "s/$/\n/g" |
grep -v "^$" |
while read -r input; do
    output=$(echo "${input}" | sed -e "s/@.*//g" -e "s#>#.DAJIN_temp/fasta/#g" -e "s/$/.fa/g")
    echo "${input}" | sed "s/@/\n/g" > "${output}"
done

# When mutation point(s) are closer to もし変異部がFASTAファイルの5'側より3'側に近い場合、
# right flanking than left flanking,   reverse complementにする。
# convert a sequence into its reverse-complement

wt_seqlen=$(cat .DAJIN_temp/fasta/wt.fa | awk '!/[>|@]/ {print length($0)}')

if [ -e .DAJIN_temp/fasta/target.fa ]; then
    convert_revcomp=$(minimap2 -ax splice .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/target.fa --cs 2>/dev/null |
        awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
        sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
        tr -d "\~\*\-\+atgc" |
        awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
        awk -v wt_seqlen="$wt_seqlen" \
        '{if(wt_seqlen-$2>$1) print 0; else print 1}')

    if [ "$convert_revcomp" -eq 1 ] ; then
        ./DAJIN/src/revcomp.sh "${fasta_LF}" \
        > .DAJIN_temp/fasta/fasta_revcomp.fa &&
        fasta_LF=".DAJIN_temp/fasta/fasta_revcomp.fa"
    fi
fi

cat ${fasta_LF} | sed "s/^/@/g" | tr -d "\n" | sed -e "s/@>/\n>/g" -e "s/$/\n/g" | grep -v "^$" |
while read -r input; do
    output=$(echo "${input}" | sed -e "s/@.*//g" -e "s#>#.DAJIN_temp/fasta_conv/#g" -e "s/$/.fa/g")
    echo "${input}" | sed "s/@/\n/g" > "$output"
done

# ============================================================================
# Format ONT reads into FASTA file
# ============================================================================
# Check wheather the files are binary:
set +e
for input in ${ont_dir}/* ; do
    output=$(echo "${input}" | sed -e "s#.*/#.DAJIN_temp/fasta_ont/#g" -e "s#\.f.*#.fa#g")
    printf "${output} is now generating...\n"
    #
    if [ "$(file ${input} | grep -c compressed)" -eq 1 ]; then
        gzip -dc "${input}" |
        awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' \
        > ".DAJIN_temp"/tmp_$$
    else
        cat "${input}" |
        awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' \
        > ".DAJIN_temp"/tmp_$$
    fi
    mv ".DAJIN_temp"/tmp_$$ "${output}"
done
set -e
#
ont_ref_nanosim=$(echo "${ont_ref}" |
    sed -e "s#.*/#.DAJIN_temp/fasta_ont/#g" -e "s#\.f.*#.fa#g")
ont_ref_barcodeID=$(echo "${ont_ref}" |
    sed -e "s#.*/##g" -e "s#\.f.*##g")

# ============================================================================
# Setting NanoSim (v2.5.0)
# ============================================================================

printf "
+++++++++++++++++++++
NanoSim simulation starts
+++++++++++++++++++++\n"

printf "Read analysis...\n"
./DAJIN/utils/NanoSim/src/read_analysis.py genome \
    -i "$ont_ref_nanosim" \
    -rg .DAJIN_temp/fasta_conv/wt.fa \
    -t ${threads:-1} \
    -o .DAJIN_temp/NanoSim/training

wt_seqlen=$(cat .DAJIN_temp/fasta/wt.fa | awk '!/[>|@]/ {print length($0)}')
for input in .DAJIN_temp/fasta_conv/*; do
    printf "${input} is now simulating...\n"
    output=$(echo "$input" | sed -e "s#fasta_conv/#fasta_ont/#g" -e "s/.fasta$//g" -e "s/.fa$//g")
    #
    ## For deletion allele
    input_seqlength=$(cat "${input}" | awk '!/[>|@]/ {print length($0)-100}')
    if [ "$input_seqlength" -lt "$wt_seqlen" ]; then
        len=${input_seqlength}
    else
        len=${wt_seqlen}
    fi
    ##
    ./DAJIN/utils/NanoSim/src/simulator.py genome \
        -dna_type linear \
        -c .DAJIN_temp/NanoSim/training \
        -rg "${input}" \
        -n 10000 \
        -t "${threads:-1}" \
        -min "${len}" \
        -o "${output}_simulated"
    ##
    rm .DAJIN_temp/fasta_ont/*_error_* .DAJIN_temp/fasta_ont/*_unaligned_* 2>/dev/null
done

rm -rf DAJIN/utils/NanoSim/src/__pycache__

printf 'Success!!\nSimulation is finished\n'


# ============================================================================
# Mapping by minimap2 for IGV visualization
# ============================================================================
printf \
"+++++++++++++++++++++
Generate BAM files
+++++++++++++++++++++\n"

./DAJIN/src/igvjs.sh "${genome:-mm10}" "${threads:-1}"
rm .DAJIN_temp/tmp_* 2>/dev/null
printf "BAM files are saved at bam\n"
printf "Next converting BAM to MIDS format...\n"

# ============================================================================
# MIDS conversion
# ============================================================================
printf \
"++++++++++++
Converting ACGT into MIDS format
++++++++++++\n"

reference=".DAJIN_temp/fasta_conv/wt.fa"
query=".DAJIN_temp/fasta_conv/target.fa"

# Get mutation loci...
true > .DAJIN_temp/data/mutation_points
for query in .DAJIN_temp/fasta_conv/target*.fa; do
    minimap2 -ax splice "${reference}" "${query}" --cs 2>/dev/null |
    awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
    sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
    tr -d "\~\*\-\+atgc" |
    awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' \
    >> .DAJIN_temp/data/mutation_points
done

if [ "$(wc -l .DAJIN_temp/data/mutation_points | cut -d ' ' -f 1)" -gt 1 ]; then
    cat .DAJIN_temp/data/mutation_points |
    awk 'NR==1 {print $1}
    END{print $NF}' |
    tr "\n" " " |
    sed "s/ $/\n/g" \
    > .DAJIN_temp/data/mutation_points_$$
    mv .DAJIN_temp/data/mutation_points_$$ .DAJIN_temp/data/mutation_points
fi

# MIDS conversion...
find .DAJIN_temp/fasta_ont -type f | sort |
awk '{print "./DAJIN/src/mids_convertion.sh",$0, "wt", NR, "&"}' |
awk -v th=${threads:-1} '{
    if (NR%th==0) gsub("&","&\nwait",$0)
    print}
    END{print "wait"}' |
sh -
#
cat .DAJIN_temp/data/MIDS_* |
sed -e "s/_aligned_reads//g" |
sort -k 1,1 \
> ".DAJIN_temp/data/DAJIN_MIDS.txt"

rm .DAJIN_temp/data/MIDS_*
# ============================================================================
# Training
# ============================================================================

# Split DAJIN data
rm -rf .DAJIN_temp/data/split 2>/dev/null
mkdir -p .DAJIN_temp/data/split

# cat ".DAJIN_temp/data/DAJIN_MIDS.txt" |
# grep simulate \
# > ".DAJIN_temp/data/DAJIN_MIDS_sim.txt"
# |
# awk 'BEGIN{OFS="\t"} 
#     {print $0, ++names[$3]}' |
# sort -k 4,4n |
# cut -f 1,2,3 |
# split -l 30000 - \
# .DAJIN_temp/data/split/DAJIN_split_

# python DAJIN/src/ml_simulate_reads_kotaro.py ".DAJIN_temp/data/DAJIN_sim.txt"
# split_num=$(find .DAJIN_temp/data/split -name "DAJIN_split_*" | wc -l)
# num=1
# for file in .DAJIN_temp/data/split/DAJIN_split_*; do
#     echo "training model: ${num}/${split_num}"
#     python DAJIN/src/ml_simulate_reads_kotaro.py "${file}"
#     num=$((num+1))
# done

# ============================================================================
# Abnormal detection and prediction
# ============================================================================

# cat ".DAJIN_temp/data/DAJIN_MIDS.txt" |
# grep -v simulate \
# > ".DAJIN_temp/data/DAJIN_real.txt"

# ============================================================================
# Prediction
# ============================================================================
printf "Start allele prediction...\n"
#
python DAJIN/src/ml_allele_profile.py ".DAJIN_temp/data/DAJIN_MIDS.txt"
# Rscript DAJIN/src/ml_abnormal_detection.R ".DAJIN_temp"/data/DAJIN_MIDS_MIDS.txt
# python DAJIN/src/ml_abnormal_detection.py ".DAJIN_temp"/data/DAJIN_MIDS_MIDS.txt
# Rscript DAJIN/src/ml_prediction.R ".DAJIN_temp"/data/DAJIN_MIDS_MIDS.txt


# python DAJIN/src/anomaly_detection.py ".DAJIN_temp"/data/DAJIN_MIDS_trimmed.txt
#
# rm .DAJIN_temp/onehot_*
# mutation_type=$(
#     minimap2 -ax splice ${reference} ${query} --cs 2>/dev/null |
#     awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
#     grep "~" | wc -l)
#
# if [ $(echo ${mutation_type}) -eq 1 ]; then
#     ./DAJIN/src/anomaly_exondeletion.sh ${genome} ${threads}
# else
#     cp .DAJIN_temp/anomaly_classification.txt .DAJIN_temp/anomaly_classification_revised.txt
# fi
#
# cp .DAJIN_temp/anomaly_classification.txt .DAJIN_temp/anomaly_classification_revised.txt
#
# python DAJIN/src/prediction.py ".DAJIN_temp"/data/DAJIN_MIDS_trimmed.txt
#
printf "Prediction was finished...\n"
#
# ============================================================================
# Report allele percentage
# ============================================================================
input=".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt"
output=".DAJIN_temp/data/DAJIN_MIDS_prediction_allele_percentage.txt"
# ----------------------------------------------------------------------------
cat "${input}"  |
cut -f 1,3 |
sort |
uniq -c |
tee -a ".DAJIN_temp/tmp_prediction" |
awk '{barcode[$2]+=$1} END{for(key in barcode) print key,barcode[key]}' |
sort |
join -1 1 -2 2 - ".DAJIN_temp/tmp_prediction" |
awk '{print $1, int($3/$2*100+0.5), $4}' \
> ".DAJIN_temp/tmp_prediction_proportion"

# Filter low-percent alleles ---------
per_refab=$(cat ".DAJIN_temp"/tmp_prediction_proportion | 
    grep "${ont_ref_barcodeID:=barcode21}" | #! define "barcode30" by automate manner
    grep abnormal |
    cut -d " " -f 2)

cat .DAJIN_temp/tmp_prediction_proportion |
awk -v refab="${per_refab}" \
    '!($2<refab+5 && $3 == "abnormal")' |
# Less than 5% allele type is removed except for target allele
awk '($2 > 5 && $3 != "target") || ($2 > 0 && $3 == "target")' |
tee -a ".DAJIN_temp/tmp_prediction_filtered" |
# Report allele percentage -------------
awk '{array[$1]+=$2}
    END{for(key in array) print key, array[key]}' |
sort |
join - ".DAJIN_temp/tmp_prediction_filtered" |
awk '{print $1, int($3*100/$2+0.5),$4}' \
> "${output}"

rm .DAJIN_temp/tmp_*

# ============================================================================
# Clustering within each allele type
# ============================================================================

input=".DAJIN_temp/data/DAJIN_MIDS_prediction_allele_percentage.txt"
control="${ont_ref_barcodeID:=barcode26}" #! define "barcode30" by automate manner
# ./DAJIN/src/test_clustering.sh ${barcode} ${control} ${allele}
# cat .tmp_/clustering_results_*

#
# rm -rf .DAJIN_temp/clustering
# rm -rf DAJIN_Report/bam_clustering
# rm -rf DAJIN_Report/allele_type
# ============================================================================
# Temporal directory
# ============================================================================
temp_dir=".DAJIN_temp/clustering/"
mkdir -p "${temp_dir}"

cat "${input}" |
cut -d " " -f 1,3 |
awk -v cont=${control} \
    '{print "./DAJIN/src/clustering.sh",$1, cont, $2, "&"}' |
 #! ---------------------------------
# grep -e barcode18 -e barcode23 -e barcode26 |
grep -e barcode18 |
 #! ---------------------------------
 awk -v th=${threads:-1} '{
    if (NR%th==0) gsub("&","&\nwait",$0)
    print}
    END{print "wait"}' |
sh -
#
# rm -rf .DAJIN_temp/clustering

# ============================================================================
# Joint sequence logo in 2-cut Exon deletion
# ============================================================================

if [ $(echo ${mutation_type}) -eq 1 ]; then
    printf "Check the intactness of a joint sequence of deletion...\n"
    ./DAJIN/src/intact_2cutdeletion.sh ${threads:-1}
fi

# ============================================================================
# KI sequence intactness
# ============================================================================
reference=fasta/wt.fa
query=fasta/target.fa

set +e
minimap2 -a ${reference} ${query} --cs 2>/dev/null |
awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
sed -e "s/cs:Z:://g" -e "s/:/\t/g" |
grep -q -i -e "ATAACTTCGTATAATGTATGCTATACGAAGTTAT" \
    -e "ATAACTTCGTATAGCATACATTATACGAAGTTAT"
if [ $? -eq 0 ]; then
    printf "Check the intactness of loxP sequence loci...\n"
    ./DAJIN/src/intact_preparation.sh
    printf "Generate sequence logo at loxP sites...\n"
    ./DAJIN/src/intact_seqlogo.sh
    printf "Search loxP exactly matched reads...\n"
    ./DAJIN/src/intact_fullmatch.sh
    python ./DAJIN/src/intact_fullmatch.py
fi
set -e

# ============================================================================
# SVG allele types
# ============================================================================
# wt_seqlen=$(cat .tmp_/wt.fa | awk '!/[>|@]/ {print length($0)}')
# cat .tmp_/mutation_points |
# awk -v reflen=${wt_seqlen} \
#     '{print int($1/reflen*100) + 10 ,int($2/reflen*100) + 10 }'

# ============================================================================
# Alignment viewing
# ============================================================================

printf "Visualizing alignment reads...\n"
printf "Browser will be launched. Click 'igvjs.html'.\n"
{ npx live-server DAJIN_Report/igvjs/ & } 1>/dev/null 2>/dev/null

# rm -rf .tmp_
# rm .DAJIN_temp/tmp_* .DAJIN_temp/clustering/tmp_* 2>/dev/null

printf "Completed! \nCheck 'results/figures/' directory.\n"

exit 0

