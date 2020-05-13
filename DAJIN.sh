#!/bin/sh

# ============================================================================
# Initialize shell environment
# ============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ============================================================================
# Define the functions for printing usage and error message
# ============================================================================
VERSION=1.0

usage(){
cat <<- USAGE 1>&2
Usage     : ./DAJIN.sh -f [text file](described at "Input")

Example   : ./DAJIN.sh -f DAJIN/example/example.txt

Input     : Input file should be formatted as below:
            # Example
            ------
            design=DAJIN/example/input.txt
            sequence=DAJIN/example/demultiplex
            control=barcode21
            genome=mm10
            grna=CCCTGCGGCCAGCTTTCAGGCAG
            threads=10
            ------
            - desing: a multi-FASTA file contains sequences of each genotype. ">wt" and ">target" must be included. 
            - sequence: a directory contains FASTA or FASTQ files of long-read sequencing
            - control: control barcode ID
            - genome: reference genome. e.g. mm10, hg38
            - grna: gRNA sequence(s). multiple gRNA sequences must be deliminated by comma. e.g. CCCTGCGGCCAGCTTTCAGGCAG,CCCTGCGGCCAGCTTTCAGGCAG
            - threads: optional. default is two-thirds of available CPU threads.
USAGE
}

usage_and_exit(){
    usage
    exit "$1"
}

error(){
    echo "$@" 1>&2
    usage_and_exit 1
}

error_exit() {
  ${2+:} false && echo "${0##*/}: $2" 1>&2
  exit $1
}

# ============================================================================
# Parse arguments
# ============================================================================
[ $# -eq 0 ] && usage_and_exit 1

while [ $# -gt 0 ]
do
    case "$1" in
        --help | --hel | --he | --h | '--?' | -help | -hel | -he | -h | '-?')
            usage_and_exit 0
            ;;
        --version | --versio | --versi | --vers | --ver | --ve | --v | \
        -version | -versio | -versi | -vers | -ver | -ve | -v )
            echo "DAJIN version: $VERSION"
            exit 0
            ;;
        --file | -f )
            fasta=$(cat "$2" | grep "design" | sed -e "s/ //g" -e "s/.*=//g")
            ont_dir=$(cat "$2" | grep "sequence" | sed -e "s/ //g" -e "s/.*=//g")
            ont_cont=$(cat "$2" | grep "control" | sed -e "s/ //g" -e "s/.*=//g")
            genome=$(cat "$2" | grep "genome" | sed -e "s/ //g" -e "s/.*=//g")
            grna=$(cat "$2" | grep "grna" | sed -e "s/ //g" -e "s/.*=//g")
            threads=$(cat "$2" | grep "threads" | sed -e "s/ //g" -e "s/.*=//g")
            ;;
        -* )
        error "Unrecognized option : $1"
            ;;
        *)
            break
            ;;
    esac
    shift
done

set +e

if [ -z "$fasta" ] || [ -z "$ont_dir" ] || [ -z "$ont_cont" ] || [ -z "$genome" ] || [ -z "$grna" ]
then
    error_exit 1 "Required arguments are not specified"
fi

if [ "$(grep -c '>target' ${fasta})" -eq 0 ] || [ "$(grep -c '>wt' ${fasta})" -eq 0 ]
then
    error_exit 1 "FASTA requires including \">target\" and \">wt\". "
fi

[ -f "$fasta" ] || error_exit 1 "No such file"
[ -d "$ont_dir" ] || error_exit 1 "No such directory"

set -e

# ============================================================================
# Define threads
# ============================================================================

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
# Required software
# ============================================================================

set +e

type python 1>/dev/null 2>/dev/null || error_exit 1 'Command "python" not found'
type samtools 1>/dev/null 2>/dev/null || error_exit 1 'Command "samtools" not found'
type minimap2 1>/dev/null 2>/dev/null || error_exit 1 'Command "minimap2" not found'
type gzip 1>/dev/null 2>/dev/null || error_exit 1 'Command "gzip" not found'

python -c "import tensorflow as tf" \
1>/dev/null 2>/dev/null ||  error_exit 1 '"Tensorflow" not found'
set -e

# ============================================================================
# For WSL (Windows Subsystem for Linux)
# ============================================================================

uname -a | 
grep Microsoft 1>/dev/null 2>/dev/null &&
alias python="python.exe"

# ============================================================================
# Setting Directory
# ============================================================================
rm -rf ".DAJIN_temp" 2>/dev/null
dirs="fasta fasta_conv fasta_ont NanoSim data \
    results/svg results/png results/igvjs \
    clustering/temp seqlogo/temp"

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
    grep -v "^$" |
cat - > .DAJIN_temp/fasta/fasta.fa
fasta_LF=".DAJIN_temp/fasta/fasta.fa"

# Separate multiple-FASTA into FASTA files
cat ${fasta_LF} |
sed "s/^/@/g" |
tr -d "\n" |
sed -e "s/@>/\n>/g" -e "s/$/\n/g" |
grep -v "^$" |
while read -r input; do
    output=$(echo "${input}" |
    sed -e "s/@.*//g" -e "s#>#.DAJIN_temp/fasta/#g" -e "s/$/.fa/g")
    #
    echo "${input}" |
        sed "s/@/\n/g" |
        awk '{if($1 !~ "^>") $0=toupper($0)
            print}' |
    cat - > "${output}"
done

# When mutation point(s) are closer to もし変異部がFASTAファイルの5'側より3'側に近い場合、
# right flanking than left flanking,   reverse complementにする。
# convert a sequence into its reverse-complement

wt_seqlen=$(cat .DAJIN_temp/fasta/wt.fa | awk '!/[>|@]/ {print length($0)}')

convert_revcomp=$(
    minimap2 -ax splice \
    .DAJIN_temp/fasta/wt.fa \
    .DAJIN_temp/fasta/target.fa --cs 2>/dev/null |
    awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
    sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
    tr -d "\~\*\-\+atgc" |
    awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
    awk -v wt_seqlen="$wt_seqlen" \
    '{if(wt_seqlen-$2>$1) print 0; else print 1}'
    )

if [ "$convert_revcomp" -eq 1 ] ; then
    ./DAJIN/src/revcomp.sh "${fasta_LF}" \
    > .DAJIN_temp/fasta/fasta_revcomp.fa &&
    fasta_LF=".DAJIN_temp/fasta/fasta_revcomp.fa"
fi

# 別々のFASTAファイルとして保存する
cat "${fasta_LF}" |
sed "s/^/@/g" |
tr -d "\n" |
sed -e "s/@>/\n>/g" -e "s/$/\n/g" |
grep -v "^$" |
while read -r input; do
    output=$(
        echo "${input}" |
        sed -e "s/@.*//g" \
            -e "s#>#.DAJIN_temp/fasta_conv/#g" \
            -e "s/$/.fa/g"
        )
    echo "${input}" |
        sed "s/@/\n/g" |
        awk '{if($1 !~ "^>") $0=toupper($0)
            print}' |
    cat - > "${output}"
done

# ----------------------------------------------------------------------------
# 変異のタイプ（deletion, knock-in, or point-mutation）を判定する
# ----------------------------------------------------------------------------
mutation_type=$(
    minimap2 -ax map-ont \
    .DAJIN_temp/fasta/wt.fa \
    .DAJIN_temp/fasta/target.fa \
    --cs 2>/dev/null |
    grep -v "^@" |
    awk '{
        cstag=$(NF-1)
        if(cstag ~ "-") print "D"
        else if(cstag ~ "+") print "I"
        else if(cstag ~ "*") print "P"
        }'
)

# ----------------------------------------------------------------------------
# Targetが一塩基変異の場合: 
# Cas9の切断部に対してgRNA部自体の欠損およびgRNA長分の塩基挿入したものを異常アレルとして作成する
# ----------------------------------------------------------------------------

if [ "$mutation_type" = "P" ]; then
    grna=CCCTGCGGCCAGCTTTCAGGCAG #!----------------------------------------------------------------------------
    grna_len=$(awk -v grna="$grna" 'BEGIN{print length(grna)}')
    grna_firsthalf=$(awk -v grna="$grna" 'BEGIN{print substr(grna, 1, int(length(grna)/2))}')
    grna_secondhalf=$(awk -v grna="$grna" 'BEGIN{print substr(grna, int(length(grna)/2)+1, length(grna))}')
    #
    ins_seq=$(
        seq_length="$grna_len" &&
        od -A n  -t u4 -N $(($seq_length*100)) /dev/urandom |
        tr -d "\n" |
        sed 's/[^0-9]//g' |
        sed "s/[4-9]//g" |
        sed -e "s/0/A/g" -e "s/1/G/g" -e "s/2/C/g" -e "s/3/T/g" |
        awk -v seq_length=$seq_length '{print substr($0, 1, seq_length)}'
    )
    # insertion
    cat .DAJIN_temp/fasta_conv/wt.fa |
        sed "s/$grna/$grna_firsthalf,$grna_secondhalf/g" |
        sed "s/,/$ins_seq/g" |
    cat - > .DAJIN_temp/fasta_conv/wt_ins.fa
    # deletion
    cat .DAJIN_temp/fasta_conv/wt.fa |
        sed "s/$grna//g" |
    cat - > .DAJIN_temp/fasta_conv/wt_del.fa
fi

# ============================================================================
# Format ONT reads into FASTA file
# ============================================================================
set +e
for input in ${ont_dir}/* ; do
    output=$(echo "${input}" | sed -e "s#.*/#.DAJIN_temp/fasta_ont/#g" -e "s#\.f.*#.fa#g")
    # printf "${output} is now generating...\n"
    #
    # Check wheather the files are binary:
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
#! --------------------------------------------------------
ont_cont=barcode32

# ont_cont_nanosim=$(echo "${ont_cont}" |
#     sed -e "s#.*/#.DAJIN_temp/fasta_ont/#g" -e "s#\.f.*#.fa#g")
# ont_cont_barcodeID=$(echo "${ont_cont}" |
#     sed -e "s#.*/##g" -e "s#\.f.*##g")

# ============================================================================
# Setting NanoSim (v2.5.0)
# ============================================================================

printf "
+++++++++++++++++++++
NanoSim simulation starts
+++++++++++++++++++++\n"

printf "Read analysis...\n"
./DAJIN/utils/NanoSim/src/read_analysis.py genome \
    -i ".DAJIN_temp/fasta_ont/${ont_cont}.fa" \
    -rg .DAJIN_temp/fasta_conv/wt.fa \
    -t ${threads:-1} \
    -o .DAJIN_temp/NanoSim/training

wt_seqlen=$(cat .DAJIN_temp/fasta/wt.fa | awk '!/[>|@]/ {print length($0)}')
for input in .DAJIN_temp/fasta_conv/*; do
    printf "${input} is now simulating...\n"
    output=$(
        echo "$input" |
        sed -e "s#fasta_conv/#fasta_ont/#g" \
            -e "s/.fasta$//g" -e "s/.fa$//g"
        )
    #
    ## For deletion allele
    input_seqlength=$(
        cat "${input}" |
        awk '!/[>|@]/ {print length($0)-100}'
        )
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
# MIDS conversion
# ============================================================================
printf \
"++++++++++++
Converting ACGT into MIDS format
++++++++++++\n"

reference=".DAJIN_temp/fasta_conv/wt.fa"
query=".DAJIN_temp/fasta_conv/target.fa"

# Get mutation loci...
# true > .DAJIN_temp/data/mutation_points
# for query in .DAJIN_temp/fasta_conv/target*.fa; do
cat "$reference" |
    minimap2 -ax splice - "${query}" --cs 2>/dev/null |
    awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
    sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
    tr -d "\~\*\-\+atgc" |
    awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
cat - > .DAJIN_temp/data/mutation_points
# done

# if [ "$(wc -l .DAJIN_temp/data/mutation_points | cut -d ' ' -f 1)" -gt 1 ]; then
#     cat .DAJIN_temp/data/mutation_points |
#         awk 'NR==1 {print $1}
#         END{print $NF}' |
#         tr "\n" " " |
#         sed "s/ $/\n/g" |
#     cat - > .DAJIN_temp/data/mutation_points_$$
#     mv .DAJIN_temp/data/mutation_points_$$ .DAJIN_temp/data/mutation_points
# fi

# MIDS conversion...
find .DAJIN_temp/fasta_ont -type f | sort |
    awk '{print "./DAJIN/src/mids_convertion.sh",$0, "wt", "&"}' |
    awk -v th=${threads:-1} '{
        if (NR%th==0) gsub("&","&\nwait",$0)
        print}
        END{print "wait"}' |
sh -
rm .DAJIN_temp/tmp_*
#
if [ "$mutation_type" = "P" ]; then
    rm .DAJIN_temp/data/MIDS_target*
fi

cat .DAJIN_temp/data/MIDS_* |
    sed -e "s/_aligned_reads//g" |
    sort -k 1,1 |
cat - > ".DAJIN_temp/data/DAJIN_MIDS.txt"

rm .DAJIN_temp/data/MIDS_*

# ============================================================================
# Prediction
# ============================================================================
printf "Start allele prediction...\n"
#
python DAJIN/src/ml_abnormal_detection.py ".DAJIN_temp"/data/DAJIN_MIDS.txt

# ----------------------------------------------------------------
# 2-cut deletionの場合は、大丈夫そうなabnormalを検出する
# ----------------------------------------------------------------

# if [ "$mutation_type" = "D" ]; then
#     ./DAJIN/src/anomaly_exondeletion.sh ${genome} ${threads}
# else
#     cp .DAJIN_temp/anomaly_classification.txt .DAJIN_temp/anomaly_classification_revised.txt
# fi
#
# cp .DAJIN_temp/anomaly_classification.txt .DAJIN_temp/anomaly_classification_revised.txt
#
if [ "$mutation_type" = "P" ]; then
    mv ".DAJIN_temp/data/DAJIN_anomaly_classification.txt" ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt"
else
    python DAJIN/src/ml_prediction.py ".DAJIN_temp/data/DAJIN_MIDS.txt"
fi
#
printf "Prediction was finished...\n"
#
# ============================================================================
# Report allele percentage
# ============================================================================
input=".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt"
output=".DAJIN_temp/data/DAJIN_MIDS_prediction_allele_percentage.txt"
# ----------------------------------------------------------------------------

# --------------------------------
# 各サンプルに含まれるアレルの割合を出す
# --------------------------------
cat "${input}"  |
    cut -f 1,3 |
    sort |
    uniq -c |
    tee ".DAJIN_temp/tmp_prediction" |
    awk '{barcode[$2]+=$1} END{for(key in barcode) print key,barcode[key]}' |
    sort |
    join -1 1 -2 2 - ".DAJIN_temp/tmp_prediction" |
    awk '{print $1, int($3/$2*100+0.5), $4}' |
cat - > ".DAJIN_temp/tmp_prediction_proportion"

# --------------------------------
# コントロールの異常アレルの割合を出す
# --------------------------------
per_refab=$(
    cat ".DAJIN_temp"/tmp_prediction_proportion | 
    grep "${ont_cont:=barcode32}" | #! define "control" by automate manner
    grep abnormal |
    cut -d " " -f 2)

# --------------------------------
# Filter low-percent alleles
# --------------------------------

cat .DAJIN_temp/tmp_prediction_proportion |
    # --------------------------------
    # 各サンプルの異常アレルの割合が
    # 「コントロールの異常アレルの割合＋5%以内」の場合、
    # その判定は偽陽性と判断し、取り除く
    # --------------------------------
    awk -v refab="${per_refab}" \
        '!($2<refab+5 && $3 == "abnormal")' |
    # --------------------------------
    # Less than 5% allele type is removed except for target allele
    # --------------------------------
    awk '($2 > 5 && $3 != "target") || ($2 > 0 && $3 == "target")' |
    tee -a ".DAJIN_temp/tmp_prediction_filtered" |
    # --------------------------------
    # Report allele percentage
    # --------------------------------
    awk '{array[$1]+=$2}
        END{for(key in array) print key, array[key]}' |
    sort |
    join - ".DAJIN_temp/tmp_prediction_filtered" |
    awk '{print $1, int($3*100/$2+0.5),$4}' |
cat - > "${output}"

rm .DAJIN_temp/tmp_*

# ============================================================================
# Clustering within each allele type
# ============================================================================
input=".DAJIN_temp/data/DAJIN_MIDS_prediction_allele_percentage.txt"
# ----------------------------------------------------------------------------

temp_dir=".DAJIN_temp/clustering/"
mkdir -p "${temp_dir}"

cat "${input}" |
    cut -d " " -f 1,3 |
    awk -v cont="${ont_cont}" \
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

rm .DAJIN_temp/tmp_*

# ============================================================================
# Sequence logo
# ============================================================================



# # ============================================================================
# # Joint sequence logo in 2-cut Exon deletion
# # ============================================================================

# if [ $(echo ${mutation_type}) -eq 1 ]; then
#     printf "Check the intactness of a joint sequence of deletion...\n"
#     ./DAJIN/src/intact_2cutdeletion.sh ${threads:-1}
# fi

# # ============================================================================
# # KI sequence intactness
# # ============================================================================
# reference=fasta/wt.fa
# query=fasta/target.fa

# set +e
# minimap2 -a ${reference} ${query} --cs 2>/dev/null |
# awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
# sed -e "s/cs:Z:://g" -e "s/:/\t/g" |
# grep -q -i -e "ATAACTTCGTATAATGTATGCTATACGAAGTTAT" \
#     -e "ATAACTTCGTATAGCATACATTATACGAAGTTAT"
# if [ $? -eq 0 ]; then
#     printf "Check the intactness of loxP sequence loci...\n"
#     ./DAJIN/src/intact_preparation.sh
#     printf "Generate sequence logo at loxP sites...\n"
#     ./DAJIN/src/intact_seqlogo.sh
#     printf "Search loxP exactly matched reads...\n"
#     ./DAJIN/src/intact_fullmatch.sh
#     python ./DAJIN/src/intact_fullmatch.py
# fi
# set -e

# # ============================================================================
# # SVG allele types
# # ============================================================================
# # wt_seqlen=$(cat .tmp_/wt.fa | awk '!/[>|@]/ {print length($0)}')
# # cat .tmp_/mutation_points |
# # awk -v reflen=${wt_seqlen} \
# #     '{print int($1/reflen*100) + 10 ,int($2/reflen*100) + 10 }'

# ============================================================================
# Mapping by minimap2 for IGV visualization
# ============================================================================
printf \
"+++++++++++++++++++++
Generate BAM files
+++++++++++++++++++++\n"
if [ "$mutation_type" = "P" ]; then
    rm .DAJIN_temp/fasta_ont/wt_ins* .DAJIN_temp/fasta_ont/wt_del*
fi

./DAJIN/src/igvjs.sh "${genome:-mm10}" "${threads:-1}"
rm .DAJIN_temp/tmp_* 2>/dev/null
printf "BAM files are saved at bam\n"
printf "Next converting BAM to MIDS format...\n"

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

