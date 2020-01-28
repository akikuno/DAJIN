#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ======================================
# Define the functions for printing usage and error message
# ======================================

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

Example   : ./miniogenotype/allele_profiler.sh \\
            -i miniogenotype/example/cables2_flox.fa \\
            -ont miniogenotype/example/demultiplex \\
            -ont_ref miniogenotype/example/demultiplex/barcode21.fastq.gz \\
            -genome mm10 \\
            -seq ATAACTTCGTATAATGTATGCTATACGAAGTTAT
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

# ======================================
# Parse arguments
# ======================================

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

if [ $(grep '>target' $fasta | wc -l) -ne 1 -a $(grep '>wt' $fasta | wc -l) -ne 1 ]; then
    error_exit 1 "FASTA requires including \">target\" and \">wt\" headers. See ./${0##*/} -h"
fi
set -e

# ======================================
# Prerequisit
# ======================================
set +e
type python 1>/dev/null 2>/dev/null
if test "$?" -eq 1; then error_exit 1 'Please install Python'; fi
type git 1>/dev/null 2>/dev/null
if test "$?" -eq 1; then error_exit 1 'Please install git'; fi
type samtools 1>/dev/null 2>/dev/null
if test "$?" -eq 1; then error_exit 1 'Please install samtools'; fi
type minimap2 1>/dev/null 2>/dev/null
if test "$?" -eq 1; then error_exit 1 'Please install minimap2'; fi
type bgzip 1>/dev/null 2>/dev/null
if test "$?" -eq 1; then error_exit 1 'Please install bgzip'; fi
#
uname -a |
grep Microsoft 1>/dev/null 2>/dev/null &&
alias python="python.exe"
python -c \
"from tensorflow.python.client import device_lib;
print(device_lib.list_local_devices())" \
1>/dev/null 2>/dev/null
if test "$?" -eq 1; then error_exit 1 'GPU is not recognized'; fi
set -e

# Define threads
set +u
# Linux and similar...
[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/2)}')
# FreeBSD and similar...
[ -z "$threads" ] && threads=$(getconf NPROCESSORS_ONLN | awk '{print int($0/2)}')
# Solaris and similar...
[ -z "$threads" ] && threads=$(ksh93 -c 'getconf NPROCESSORS_ONLN' | awk '{print int($0/2)}')
# Give up...
[ -z "$threads" ] && threads=1
set -u

# ======================================
# Setting Directory
# ======================================

dirs="fasta fasta_ont bam data_for_ml results .tmp_ results/figures/svg results/figures/png results/igvjs"
rm -rf $dirs NanoSim
mkdir -p $dirs

# ======================================
# Format FASTA file
# ======================================

# CRLF to LF
cat ${fasta} |
sed -e "s/\r$//" |
grep -v "^$" \
> .tmp_/fasta
fasta=.tmp_/fasta

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

ref_seqlength=$(cat fasta/wt.fa | sed 1d | awk '{print length($0)}')
minimap2 -ax splice fasta/wt.fa fasta/target.fa --cs 2>/dev/null |
awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
sed -e "s/cs:Z:://g" -e "s/:/\t/g" |
tr -d "\*\-\+atgc" |
awk '{print $1,$1+$2}' |
awk -v ref_seqlength="$ref_seqlength" \
'{if(ref_seqlength-$2>$1) print 0; else print 1}' \
> .tmp_/revcomp

if [ $(cat .tmp_/revcomp) -eq 1 ] ; then
    ./miniogenotype/src/revcomp.sh "${fasta}" \
    > .tmp_/fasta_revcomp && fasta=".tmp_/fasta_revcomp"
    cat $fasta | sed "s/^/@/g" | tr -d "\n" | sed -e "s/@>/\n>/g" -e "s/$/\n/g" | grep -v "^$" |
    while read input; do
        output=$(echo $input | sed -e "s/@.*//g" -e "s#>#fasta/#g" -e "s/$/.fa/g")
        echo $input | sed "s/@/\n/g" > $output
    done
fi

# ======================================
# Format ONT reads into FASTA file
# ======================================
# Check wheather the files are binary:
set +e
for input in ${ont}/* ; do
    output=$(echo $input | sed -e "s#.*/#fasta_ont/#g" -e "s#\..*#.fa#g")
    printf "${output} is now generating...\n"
    #
    if [ $(file ${input} | grep compressed | wc -l) -eq 1 ]; then
        cat ${input} | bgzip -dc |
        awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' > .tmp_/tmp_$$
    else
        cat ${input} | awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' > .tmp_/tmp_$$
    fi
        cat .tmp_/tmp_$$ > ${output}
done
set -e
ont_ref=$(echo $ont_ref | sed -e "s#.*/#fasta_ont/#g" -e "s#\..*#.fa#g")

# ======================================
# Setting NanoSim (v2.5.0)
# ======================================

# git clone https://github.com/bcgsc/NanoSim # version 2.5.0
# cat NanoSim/src/simulator.py  |
# sed "/Simulate unaligned reads/,/Merging error subfiles/d" \
# > tmp
# chmod 755 tmp && mv tmp NanoSim/src/simulator.py

printf "
+++++++++++++++++++++
NanoSim simulation starts
+++++++++++++++++++++\n"

printf "Read analysis...\n"
./miniogenotype/NanoSim/src/read_analysis.py genome \
    -i "$ont_ref" \
    -rg fasta/wt.fa \
    -t ${threads:-4} \
    -o NanoSim/training

ref_seqlength=$(cat fasta/wt.fa | sed 1d | awk '{print length($0)}')
for input in $(ls fasta/* | grep -v igv); do
    printf "${input} is now simulating...\n"
    output=$(echo $input | sed -e "s#fasta/#fasta_ont/#g" -e "s/.fasta$//g" -e "s/.fa$//g")
    input_seqlength=$(cat $input | sed 1d | awk '{print length($0)-100}')
    ## For deletion allele
    if [ $input_seqlength -lt $ref_seqlength ]; then
        len=$input_seqlength
    else
        len=$ref_seqlength
    fi
    ##
    ./miniogenotype/NanoSim/src/simulator.py genome \
        -dna_type linear -c NanoSim/training \
        -rg $input -n 3000 -t ${threads:-4} \
        -min ${len} \
        -o ${output}_simulated
    rm fasta_ont/*_error_* # fasta_ont/*_unaligned_*
done

rm NanoSim/training*

printf 'Success!!\nSimulation is finished\n'


# ======================================
# Mapping by minimap2 for IGV visualization
# ======================================
printf \
"+++++++++++++++++++++
Generate BAM files
+++++++++++++++++++++\n"

./miniogenotype/src/igvjs.sh ${genome}

printf "BAM files are saved at bam\n"
printf "Next converting BAM to MIDS format...\n"

# ======================================
# MIDS conversion
# ======================================
printf \
"++++++++++++
Converting ACGT into MIDS format
++++++++++++\n"

reference=fasta/wt.fa
query=fasta/target.fa

minimap2 -a ${reference} ${query} --cs 2>/dev/null |
awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
sed -e "s/cs:Z:://g" -e "s/:/\t/g" |
tr -d "\*\-\+atgc" |
awk '{print $1,$1+$2}' > .tmp_/mutation_points

reflength=$(cat ${reference} | grep -v "^>" | awk '{print length($0)}')
ext=${ext:=100}
first_flank=$(cat .tmp_/mutation_points | awk -v ext=${ext} '{print $1-ext}')
second_flank=$(cat .tmp_/mutation_points | awk -v ext=${ext} '{if(NF==2) print $2+ext; else print $1+ext}')
if [ "$first_flank" -lt 1 ]; then first_flank=1; fi
if [ "$second_flank" -gt "$reflength" ]; then second_flank=$(($reflength)); fi
# echo $first_flank $second_flank

true > data_for_ml/${output_file:=sequence_MIDS}.txt

for input in ./fasta_ont/*; do
    output=$(echo ${input} | sed -e "s#.*/##g" -e "s#\..*##g" -e "s/_aligned_reads//g")
    #
    minimap2 -t ${threads:-4} --cs=long -ax splice ${reference} ${input} 2>/dev/null |
    awk '$3 == "wt"' > .tmp_/${output}
    #
    ./miniogenotype/src/mids_convertion.sh .tmp_/${output} ${first_flank} ${second_flank} \
    >> data_for_ml/${output_file}.txt
done

gzip -f data_for_ml/${output_file}.txt

printf "Finished.\n${output_file}.txt.gz is generated.\n"

printf "Start allele prediction...\n"
python miniogenotype/src/allele_prediction.py data_for_ml/${output_file}.txt.gz
printf "Prediction was finished...\n"

# ======================================
# KI sequence intactness
# ======================================
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
    ./miniogenotype/src/ki_seq_search.sh
    ./miniogenotype/src/fullmatch_search.sh
    printf "Finished!...\n"
fi
set -e

printf "Completed\n"

exit 0

