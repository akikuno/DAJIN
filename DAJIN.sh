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
            -o [string] (optional) \\
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
            -o         : Output file name :default=sequence
            -seq       : Target knock-in sequence
	            -t         : Number of threads: default=4

USAGE
    exit 1
}

error_exit() {
    ${2+:} false && echo "${0##*/}: $2" 1>&2
    exit $1
}
warning() {
    ${1+:} false && echo "${0##*/}: $1" 1>&2
}

# ============================================================================
# Parse arguments
# ============================================================================

[ $# -eq 0 ] && print_usage_and_exit; fi

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
        -ont) [ ! -d "${2}" ] && error_exit 1 'Invalid -ont option: Directory not found'
            ont=${2}
        ;;
        -ont_ref) [ ! -s "${2}" ] && error_exit 1 'Invalid -ont_ref option: Empty file'
            ont_ref=${2}
        ;;
        -genome) [ ! -n "${2}" ] && error_exit 1 'Invalid -genome option'
            genome=${2}
        ;;
        -o) [ $(echo "${2}" | awk '{print NF}') -ne 1 ] && error_exit 1 'Invalid -o option'
            output_file=${2}
        ;;
        -seq) [ $(echo "${2}" | awk '{print NF}') -ne 1 ] && error_exit 1 'Invalid -seq option'
            mut_seq=${2}
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

if [ $(grep -i '>target' $fasta | wc -l) -ne 1 -a $(grep -i '>wt' $fasta | wc -l) -ne 1 ]; then
    error_exit 1 "FASTA requires including \">target\" and \">wt\" headers. See ./${0##*/} -h"
fi
set -e

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
[ "$?" -eq 1 ] &&  error_exit 1 'GPU is not recognized'
set -e

# ============================================================================
# Setting Directory
# ============================================================================

dirs="fasta fasta_ont bam data_for_ml results .tmp_ results/figures/svg results/figures/png results/igvjs"
rm -rf ${dirs} .tmp_/NanoSim
mkdir -p ${dirs}

# ============================================================================
# Format FASTA file
# ============================================================================

# CRLF to LF
cat ${fasta} |
tr -d "\r" |
grep -v "^$" \
> .tmp_/fasta.fa
fasta=".tmp_/fasta.fa"

# Separate multiple-FASTA into FASTA files
cat ${fasta} |
sed "s/^/@/g" |
tr -d "\n" |
sed -e "s/@>/\n>/g" -e "s/$/\n/g" |
grep -v "^$" |
while read input; do
    output=$(echo ${input} | sed -e "s/@.*//g" -e "s#>#fasta/#g" -e "s/$/.fa/g")
    echo ${input} | sed "s/@/\n/g" > ${output}
done
cp fasta/wt.fa .tmp_/
cp fasta/target.fa .tmp_/

# When mutation point(s) are closer to もし変異部がFASTAファイルの5'側より3'側に近い場合、
# right flanking than left flanking,   reverse complementにする。
# convert a sequence into its reverse-complement

ref_seqlength=$(cat .tmp_/wt.fa | awk '!/[>|@]/ {print length($0)}')

conv_revcomp=$(minimap2 -ax splice .tmp_/wt.fa .tmp_/target.fa --cs 2>/dev/null |
awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
tr -d "\~\*\-\+atgc" |
awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
awk -v ref_seqlength="$ref_seqlength" \
'{if(ref_seqlength-$2>$1) print 0; else print 1}')

if [ $(echo "$conv_revcomp") -eq 1 ] ; then
    ./DAJIN/src/revcomp.sh "${fasta}" \
    > .tmp_/fasta_revcomp.fa &&
    fasta=".tmp_/fasta_revcomp.fa"
    #
    cat ${fasta} | sed "s/^/@/g" | tr -d "\n" | sed -e "s/@>/\n>/g" -e "s/$/\n/g" | grep -v "^$" |
    while read input; do
        output=$(echo "${input}" | sed -e "s/@.*//g" -e "s#>#fasta/#g" -e "s/$/.fa/g")
        echo "${input}" | sed "s/@/\n/g" > $output
    done
fi

# ============================================================================
# Format ONT reads into FASTA file
# ============================================================================
# Check wheather the files are binary:
set +e
for input in ${ont}/* ; do
    output=$(echo ${input} | sed -e "s#.*/#fasta_ont/#g" -e "s#\..*#.fa#g")
    printf "${output} is now generating...\n"
    #
    if [ $(file ${input} | grep compressed | wc -l) -eq 1 ]; then
        cat ${input} | gzip -dc |
        awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' \
        > .tmp_/tmp_$$
    else
        cat ${input} |
        awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' \
        > .tmp_/tmp_$$
    fi
    mv .tmp_/tmp_$$ ${output}
done
set -e
#
ont_ref=$(echo ${ont_ref} |
    sed -e "s#.*/#fasta_ont/#g" -e "s#\..*#.fa#g")

# ============================================================================
# Setting NanoSim (v2.5.0)
# ============================================================================

printf "
+++++++++++++++++++++
NanoSim simulation starts
+++++++++++++++++++++\n"

printf "Read analysis...\n"
./DAJIN/utils/NanoSim/src/read_analysis.py genome \
    -i "$ont_ref" \
    -rg fasta/wt.fa \
    -t ${threads:-1} \
    -o .tmp_/NanoSim/training

ref_seqlength=$(cat fasta/wt.fa | awk '!/[>|@]/ {print length($0)}')
for input in fasta/*; do
    printf "${input} is now simulating...\n"
    output=$(echo $input | sed -e "s#fasta/#fasta_ont/#g" -e "s/.fasta$//g" -e "s/.fa$//g")
    input_seqlength=$(cat ${input} | sed 1d | awk '{print length($0)-100}')
    ## For deletion allele
    if [ "$input_seqlength" -lt "$ref_seqlength" ]; then
        len=${input_seqlength}
    else
        len=${ref_seqlength}
    fi
    ##
    ./DAJIN/utils/NanoSim/src/simulator.py genome \
        -dna_type linear -c .tmp_/NanoSim/training \
        -rg $input -n 3000 -t ${threads:-1} \
        -min ${len} \
        -o ${output}_simulated
    ##
    rm fasta_ont/*_error_* fasta_ont/*_unaligned_* 2>/dev/null
done

rm -rf .tmp_/NanoSim \
    DAJIN/utils/NanoSim/src/__pycache__

printf 'Success!!\nSimulation is finished\n'


# ============================================================================
# Mapping by minimap2 for IGV visualization
# ============================================================================
printf \
"+++++++++++++++++++++
Generate BAM files
+++++++++++++++++++++\n"

./DAJIN/src/igvjs.sh ${genome:-mm10} ${threads:-1}

printf "BAM files are saved at bam\n"
printf "Next converting BAM to MIDS format...\n"

# ============================================================================
# MIDS conversion
# ============================================================================
printf \
"++++++++++++
Converting ACGT into MIDS format
++++++++++++\n"

reference=fasta/wt.fa
query=fasta/target.fa

# Get mutation loci...
minimap2 -ax splice ${reference} ${query} --cs 2>/dev/null |
awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
tr -d "\~\*\-\+atgc" |
awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' \
> .tmp_/mutation_points

# MIDS conversion...
find fasta_ont -type f | sort |
awk '{print "./DAJIN/src/mids_convertion.sh",$0, "&"}' |
awk -v th=${threads:-1} '{
    if (NR%th==0) gsub("&","&\nwait",$0)
    print}
    END{print "wait"}' |
sh -
cat .tmp_/MIDS_* |
sort -k 1,1 \
> data_for_ml/${output_file:=DAJIN}.txt

# One-hot encording...
for i in M I D S; do
    { cat data_for_ml/${output_file}.txt |
    cut -f 2 |
    sed -e "s/^/MIDS=/" |
    sed -e "s/[^${i}]/0 /g" |
    sed -e "s/${i}/1 /g" |
    sed -e "s/ $//" \
    > .tmp_/onehot_${i}.txt & } 1>/dev/null 2>/dev/null
done
wait 2>/dev/null

cat data_for_ml/${output_file}.txt |
cut -f 1,3 \
> data_for_ml/${output_file}_trimmed.txt
#
# bgzip -f data_for_ml/${output_file}.txt & 1>/dev/null 2>/dev/null
# wait 2>/dev/null

printf "Finished.\n${output_file}.txt is generated.\n"

# ============================================================================
# Prediction
# ============================================================================
printf "Start allele prediction...\n"
#
python DAJIN/src/anomaly_detection.py data_for_ml/${output_file:-DAJIN}_trimmed.txt
#
rm .tmp_/onehot_*
# mutation_type=$(
#     minimap2 -ax splice ${reference} ${query} --cs 2>/dev/null |
#     awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
#     grep "~" | wc -l)
#
# if [ $(echo ${mutation_type}) -eq 1 ]; then
#     ./DAJIN/src/anomaly_exondeletion.sh ${genome} ${threads}
# else
#     cp .tmp_/anomaly_classification.txt .tmp_/anomaly_classification_revised.txt
# fi
#
# cp .tmp_/anomaly_classification.txt .tmp_/anomaly_classification_revised.txt
#
python DAJIN/src/prediction.py data_for_ml/${output_file:-DAJIN}.txt.gz
#
printf "Prediction was finished...\n"
#
# ============================================================================
# Report allele percentage
# ============================================================================
#
cat .tmp_/DAJIN_prediction_result.txt  | cut -f 1,3 | sort | uniq -c > .tmp_/tmp_prediction
#
cat .tmp_/tmp_prediction |
awk '{barcode[$2]+=$1} END{for(key in barcode) print key,barcode[key]}' |
sort |
join -1 1 -2 2 - .tmp_/tmp_prediction |
awk '{print $1, int($3/$2*100+0.5), $4}' \
> .tmp_/tmp1_prediction

per_refab=$(cat .tmp_/tmp1_prediction | 
    grep barcode30 | #! define "barcode30" by automate manner
    grep abnormal |
    cut -d " " -f 2)

# Filter low-percent alleles ------------------------------------------------
cat .tmp_/tmp1_prediction |
awk -v refab=${per_refab} \
    '{if( !($2<refab+5 && $3 == "abnormal") && ($2>5) ) print $0}' \
> .tmp_/tmp_prediction_filtered

# Report allele percentage ------------------------------------------------
cat .tmp_/tmp_prediction_filtered |
awk '{array[$1]+=$2}
    END{for(key in array) print key, array[key]}' |
sort |
join - .tmp_/tmp_prediction_filtered |
awk '{print $1, int($3*100/$2+0.5),$4}' \
> .tmp_/DAJIN_prediction_allele_percentage

# ============================================================================
# Clustering within each allele type
# ============================================================================

barcode=barcode30
allele=wt
cat .tmp_/DAJIN_prediction_allele_percentage |
cut -d " " -f 1,3 |
while read input; do
    barcode=$(echo ${input} | cut -d " " -f 1)
    allele=$(echo ${input} | cut -d " " -f 2)
    #

done
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
ref_seqlength=$(cat .tmp_/wt.fa | awk '!/[>|@]/ {print length($0)}')
cat .tmp_/mutation_points |
awk -v reflen=${ref_seqlength} \
    '{print int($1/reflen*100) + 10 ,int($2/reflen*100) + 10 }'

# ============================================================================
# Alignment viewing
# ============================================================================

printf "Visualizing alignment reads...\n"
printf "Browser will be launched. Click 'igvjs.html'.\n"
{ npx live-server results/igvjs/ & } 1>/dev/null 2>/dev/null

# rm -rf .tmp_
printf "Completed! \nCheck 'results/figures/' directory.\n"

exit 0

