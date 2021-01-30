#!/bin/bash

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

VERSION=0.4

usage(){
cat <<- USAGE
Usage     : ./DAJIN/DAJIN.sh -i [text file] (described at "Input")

Example   : ./DAJIN/DAJIN.sh -i DAJIN/example/design.txt

Input     : Input file should be formatted as below:
            # Example
            ------
            design=DAJIN/example/example.fa
            input_dir=DAJIN/example/fastq
            control=barcode01
            genome=mm10
            grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
            output_dir=DAJIN_example
            threads=10
            filter=on
            ------
            - design: a multi-FASTA file contains sequences of each genotype. ">wt" and ">target" must be included.
            - input_dir: a directory contains FASTA or FASTQ files of long-read sequencing
            - control: control barcode ID
            - genome: reference genome. e.g. mm10, hg38
            - grna: gRNA sequence(s) including PAM. multiple gRNA sequences must be deliminated by comma.
            - output_dir (optional): output directory name. optional. Default is "DAJIN_results"
            - threads (optional; integer): Default is two-thirds of available CPU threads.
            - filter (optional; "on" or "off"): set filter to remove very minor allele (less than 3%). Default is "on"
USAGE
}

usage_and_exit(){
    usage
    exit 1
}

error_exit() {
    echo "$@" 1>&2
    exit 1
}

################################################################################
#! Parse arguments
################################################################################
[ $# -eq 0 ] && usage_and_exit

while [ $# -gt 0 ]
do
    case "$1" in
        --help | --hel | --he | --h | '--?' | -help | -hel | -he | -h | '-?')
            usage_and_exit
            ;;
        --version | --versio | --versi | --vers | --ver | --ve | --v | \
        -version | -versio | -versi | -vers | -ver | -ve | -v )
            echo "DAJIN version: $VERSION" && exit 0
            ;;
        --input | --in | --i | -i )
            if ! [ -r "$2" ]; then
                error_exit "$2: No such file"
            fi
            design=$(cat "$2" | grep "^design" | sed -e "s/ //g" -e "s/.*=//g")
            input_dir=$(cat "$2" | grep "^input_dir" | sed -e "s/ //g" -e "s/.*=//g")
            control=$(cat "$2" | grep "^control" | sed -e "s/ //g" -e "s/.*=//g")
            genome=$(cat "$2" | grep "^genome" | sed -e "s/ //g" -e "s/.*=//g")
            grna=$(cat "$2" | grep "^grna" | sed -e "s/ //g" -e "s/.*=//g")
            output_dir=$(cat "$2" | grep "^output_dir" | sed -e "s/ //g" -e "s/.*=//g")
            threads=$(cat "$2" | grep "^threads" | sed -e "s/ //g" -e "s/.*=//g")
            filter=$(cat "$2" | grep "^filter" | sed -e "s/ //g" -e "s/.*=//g")
            TEST=$(cat "$2" | grep "^TEST" | sed -e "s/ //g" -e "s/.*=//g")
            ;;
        -* )
        error_exit "Unrecognized option : $1"
            ;;
        *)
            break
            ;;
    esac
    shift
done

#===========================================================
#? Check required arguments
#===========================================================

[ -z "$design" ] && error_exit "design argument is not specified"
[ -z "$input_dir" ] && error_exit "input_dir argument is not specified"
[ -z "$control" ] && error_exit "control argument is not specified"
[ -z "$genome" ] && error_exit "genome argument is not specified"
[ -z "$grna" ] && error_exit "grna argument is not specified"

#===========================================================
#? Check fasta file
#===========================================================

[ -e "$design" ] || error_exit "$design: No such file"

[ "$(grep -c -e '>wt' -e '>target' ${design})" -ne 2 ] &&
    error_exit "$design: design must include '>target' and '>wt'. "

#===========================================================
#? Check directory
#===========================================================

[ -d "${input_dir}" ] || error_exit "$input_dir: No such directory"

fastq_num=$(find ${input_dir}/* -type f | grep -c -e ".fq" -e ".fastq")

[ "$fastq_num" -eq 0 ] && error_exit "$input_dir: No FASTQ file in directory"

#===========================================================
#? Check control
#===========================================================

if find ${input_dir}/ -type f | grep -q "${control}"; then
    :
else
    error_exit "$control: No control file in ${input_dir}"
fi

#===========================================================
#? Check genome
#===========================================================

genome_check=$(
    wget -q -O - "http://hgdownload.soe.ucsc.edu/downloads.html" |
    grep hgTracks |
    grep -c "${genome:-XXX}"
)

[ "$genome_check" -eq 0 ] &&
    error_exit "$genome: No such reference genome"

#===========================================================
#? Check grna
#===========================================================

set $(echo "${grna}" | sed "s/,/ /g")
x=1
while [ "${x}" -le $# ]
do
    [ "$(grep -c ${1} ${design})" -eq 0 ] && error_exit "No gRNA sites"
    x=$(( "${x}" + 1 ))
done

#===========================================================
#? Check output directory name
#===========================================================

[ $(echo "$output_dir" | sed "s/[_a-zA-Z0-9]*//g" | wc | awk '{print $2}') -ne 0 ] &&
    error_exit "$output_dir: invalid directory name"

#===========================================================
#? Check "filter"
#===========================================================

if [ -z "${filter}" ]; then
    filter=on
elif [ _"${filter}" = _"on" ] || [ _"${filter}" = _"off" ]; then
    :
else
    error_exit "${filter}: invalid filter name (on/off)"
fi

#===========================================================
#? Define threads
#===========================================================

{
unset max_threads tmp_threads
max_threads=$(getconf _NPROCESSORS_ONLN)
[ -z "$max_threads" ] && max_threads=$(getconf NPROCESSORS_ONLN)
[ -z "$max_threads" ] && max_threads=$(ksh -c 'getconf NPROCESSORS_ONLN')
[ -z "$max_threads" ] && max_threads=1
tmp_threads=$(("${threads}" + 0))
}  2>/dev/null || true

if [ "${tmp_threads:-0}" -gt 1 ] && [ "${tmp_threads}" -lt "${max_threads}" ]
then
    :
else
    threads=$(echo "${max_threads}" | awk '{print int($0*2/3+0.5)}')
fi

################################################################################
#! Setting Conda environment
################################################################################

. DAJIN/src/conda_setting.sh

################################################################################
#! Formatting environments
################################################################################

#===========================================================
#? Make temporal directory
#===========================================================

dirs="fasta fasta_conv fasta_ont NanoSim data"
echo "${dirs}" |
    sed "s:^:.DAJIN_temp/:g" |
    sed "s: : .DAJIN_temp/:g" |
xargs mkdir -p

./DAJIN/src/format_fasta.sh "$design" "$input_dir" "$grna"

################################################################################
#! NanoSim (v2.5.0)
################################################################################

cat << EOF >&2
--------------------------------------------------------------------------------
NanoSim read simulation
--------------------------------------------------------------------------------
EOF

set +u
conda activate DAJIN_nanosim
set -u

if [ "$(find .DAJIN_temp/fasta_ont | grep -c simulated)" -eq 0 ]; then
    ./DAJIN/src/nanosim.sh "${control}" "${threads}"
fi

################################################################################
#! MIDS conversion
################################################################################

cat << EOF >&2
--------------------------------------------------------------------------------
Preprocessing
--------------------------------------------------------------------------------
EOF

set +u
conda activate DAJIN
set -u

# #===========================================================
# #? Get mutation loci
# #===========================================================
ref=".DAJIN_temp/fasta_conv/wt.fa"
que=".DAJIN_temp/fasta_conv/target.fa"

minimap2 -ax splice "$ref" "$que"  --cs 2>/dev/null |
  awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
  sed -e "s/cs:Z:://g" -e "s/:/ /g" -e "s/~/ /g" |
  tr -d "\~\*\-\+atgc" |
  awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
cat > .DAJIN_temp/data/mutation_points

unset ref que

#===========================================================
#? MIDS conversion
#===========================================================

find .DAJIN_temp/fasta_ont -type f |
  sort |
  awk '{print "./DAJIN/src/mids_classification.sh", $0, "wt", "&"}' |
  awk -v th=${threads:-1} 'NR%th==0 {sub("&", "&\nwait")}1;END{print "wait"}' |
sh - 2>/dev/null

################################################################################
#! Prediction
################################################################################

cat << EOF >&2
--------------------------------------------------------------------------------
Predict allele types
--------------------------------------------------------------------------------
EOF

./DAJIN/src/ml_prediction.sh "${control}" "${threads}" > .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt ||
  exit 1


################################################################################
#! Clustering
################################################################################

cat << EOF >&2
--------------------------------------------------------------------------------
Clustering alleles
--------------------------------------------------------------------------------
EOF

rm -rf .DAJIN_temp/clustering 2>/dev/null || true
mkdir -p .DAJIN_temp/clustering/temp

#===========================================================
#? Prepare control's score to define sequencing error
#===========================================================

./DAJIN/src/clustering_control_scoring.sh "${control}" "${threads}"

#===========================================================
#? Clustering
#===========================================================

cat .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
  cut -f 2,3 |
  sort -u |
  awk -v ctrl="$control" '$1 $2 != ctrl "wt"' |
  awk -v th="${threads:-1}" '{print "./DAJIN/src/clustering.sh", $1, $2, th}' |
sh -

cat .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v ctrl="$control" '$2 $3 == ctrl "wt" {print $1"\t"1}' |
cat > ".DAJIN_temp/clustering/temp/hdbscan_${control}_wt"
cat ".DAJIN_temp/clustering/temp/hdbscan_${control}_wt" > .DAJIN_temp/clustering/temp/query_seq_${control}_wt
true > ".DAJIN_temp/clustering/temp/possible_true_mut_${control}_wt"

#===========================================================
#? Allele percentage
#===========================================================

rm -rf ".DAJIN_temp/clustering/allele_per/" 2>/dev/null
mkdir -p ".DAJIN_temp/clustering/allele_per/"

cat .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
  cut -f 2 |
  sort -u |
  awk -v filter="${filter:-on}" '{print "./DAJIN/src/clustering_allele_percentage.sh", $1, filter, "&"}' |
  awk -v th=${threads:-1} 'NR%th==0 {sub("&", "&\nwait")}1;END{print "wait"}' |
sh -

################################################################################
#! Get consensus sequence in each cluster
################################################################################

cat << EOF >&2
--------------------------------------------------------------------------------
Report consensus sequence
--------------------------------------------------------------------------------
EOF

#===========================================================
#? Setting directory
#===========================================================

rm -rf .DAJIN_temp/consensus 2>/dev/null || true
mkdir -p .DAJIN_temp/consensus/temp .DAJIN_temp/consensus/sam

#===========================================================
#? Generate temporal SAM files
#===========================================================

cat .DAJIN_temp/clustering/allele_per/label* |
  cut -d " " -f 1,2 |
  sort -u |
  grep -v abnormal |  #TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
while read -r input; do
  barcode="${input%% *}"
  mapping_alleletype="$(echo "${input##* }" | sed "s/abnormal/wt/g" | sed "s/normal/wt/g")"

  cat .DAJIN_temp/clustering/allele_per/readid_cl_mids_"${barcode}"_"${mapping_alleletype}" |
    awk '{print ">"$1}' |
    sort |
  cat > .DAJIN_temp/consensus/tmp_id

  cat .DAJIN_temp/fasta_ont/"${barcode}".fa |
    awk 'BEGIN{RS=">"}{print ">"$1" "$NF}' |
    sed 1d |
    sort |
    join - .DAJIN_temp/consensus/tmp_id |
    awk '{gsub(" ", "\n")}1' |
    minimap2 -ax map-ont -t "${threads}" ".DAJIN_temp/fasta/${mapping_alleletype}.fa" - --cs=long 2>/dev/null |
    sort |
  cat > .DAJIN_temp/consensus/sam/"${barcode}"_"${mapping_alleletype}".sam
rm .DAJIN_temp/consensus/tmp_id
done

#===========================================================
#? Execute consensus.sh
#===========================================================

cat .DAJIN_temp/clustering/allele_per/label* |
    awk '{nr[$1]++; print $0, nr[$1]}' |
    grep -v abnormal |  #TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    awk '{print "./DAJIN/src/consensus.sh", $0, "&"}' |
  awk -v th=${threads:-1} 'NR%th==0 {sub("&", "&\nwait")}1;END{print "wait"}' |
sh -

################################################################################
#! Summarize to Details.csv and Details.pdf
################################################################################

./DAJIN/src/details.sh

################################################################################
#! Mapping by minimap2 for IGV visualization
################################################################################

cat << EOF >&2
--------------------------------------------------------------------------------
Generate BAM files
--------------------------------------------------------------------------------
EOF

./DAJIN/src/generate_bam.sh "${genome}" "${threads}"

################################################################################
#! Move output files
################################################################################

rm -rf "${output_dir:=DAJIN_results}" 2>/dev/null || true
mkdir -p "${output_dir:-DAJIN_results}"/BAM
mkdir -p "${output_dir:-DAJIN_results}"/Consensus

#===========================================================
#? BAM
#===========================================================

rm -rf .DAJIN_temp/bam/temp 2>/dev/null || true
cp -r .DAJIN_temp/bam/* "${output_dir:-DAJIN_results}"/BAM/ 2>/dev/null

#===========================================================
#? Consensus
#===========================================================

(   find .DAJIN_temp/consensus/* -type d |
    grep -v -e "consensus/temp" -e "sam" |
    xargs -I @ cp -f -r @ "${output_dir:-DAJIN_results}"/Consensus/
) 2>/dev/null || true

#===========================================================
#? Details
#===========================================================

cp .DAJIN_temp/details/* "${output_dir:-DAJIN_results}"/

################################################################################
#! Finish call
################################################################################

[ -z "${TEST}" ] && rm -rf .DAJIN_temp/

set +u
conda deactivate

cat << EOF >&2
--------------------------------------------------------------------------------
Completed!
Check ${output_dir:-DAJIN_results} directory
--------------------------------------------------------------------------------
EOF

exit 0