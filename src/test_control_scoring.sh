#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -eu
umask 0022
export LC_ALL=C
export UNIX_STD=2003  # to make HP-UX conform to POSIX

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# control=barcode30
# threads=12

#===========================================================
#? Auguments
#===========================================================

control="${1}"
threads="${2}"

#===========================================================
#? Input
#===========================================================

alleletype="wt"

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/"

#===========================================================
#? Temporal
#===========================================================

MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${control}_${alleletype}"
tmp_MIDS=".DAJIN_temp/clustering/temp/tmp_MIDS_${control}_${alleletype}"
tmp_control=".DAJIN_temp/clustering/temp/tmp_control_${control}_${alleletype}"
tmp_strecher=".DAJIN_temp/clustering/temp/tmp_strecher_${control}_${alleletype}"
tmp_strecher_wt=".DAJIN_temp/clustering/temp/tmp_strecher_wt_${control}_${alleletype}"
tmp_strecher_label=".DAJIN_temp/clustering/temp/tmp_strecher_label${control}_${alleletype}"

################################################################################
#! Generate mutation scores of WT alleles
################################################################################

#===========================================================
#? MIDS conversion for clustering
#===========================================================

./DAJIN/src/mids_clustering.sh "${control}" "${alleletype}" > "${MIDS_ref}"

#===========================================================
#? Select WT allele
#===========================================================

cat "${MIDS_ref}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v alelle="$alleletype" '$NF==alelle' |
    cut -d " " -f 2 |
    awk -F "" 'BEGIN{OFS=","}{$1=$1}1' |
cat > "${tmp_MIDS}"

#===========================================================
#? Generate mutation scoring of wt alleles
#===========================================================

Rscript DAJIN/src/test_control_scoring_wt.R "${tmp_MIDS}" "${threads}" #! RENAME <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
#! Generate mutation scores of other alleles
################################################################################

#===========================================================
#? Annotate mutation loci of other possible alleles
#? (mutation loci = 1)
#===========================================================

cat ".DAJIN_temp/fasta/${alleletype}.fa" |
    sed 1d |
    awk -F "" '{for(i=1;i<=NF;i++) print $i, 0}' |
    awk '{print NR"_"$1, $2}' |
    sort -t " " -k 1,1 |
cat > "${tmp_control}"

find .DAJIN_temp/fasta/ -type f |
    grep -v wt.fa |
    grep -v fasta.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
while read -r label; do

    stretcher \
        -asequence .DAJIN_temp/fasta/wt.fa \
        -bsequence .DAJIN_temp/fasta/"${label}".fa \
        -aformat markx2 \
        -outfile "${tmp_strecher}" 2>/dev/null

    cat "${tmp_strecher}" |
        awk 'NF==2' |
        awk '$0 ~ /[ACGT.-]/' |
        awk 'NR%2==1 {printf $2}' |
        awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    cat > "${tmp_strecher_wt}"

    cat "${tmp_strecher}" |
        awk 'NF==2' |
        awk '$0 ~ /[ACGT.-]/' |
        awk 'NR%2==0 {printf $2}' |
        awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    cat > "${tmp_strecher_label}"

    splice_num="$(
        minimap2 -ax map-ont .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/${label}.fa 2>/dev/null |
        awk '$2==0 || $2==16 || $2==2048 || $2==2064' |
        grep -c -v "^@"
        )"

    if [ "${splice_num}" -eq 3 ]; then
    # inversion
    seq_len="$(cat .DAJIN_temp/fasta/"$label".fa | sed 1d | awk '{print length}')"
    minimap2 -ax map-ont .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/${label}.fa 2>/dev/null |
        awk '$2==0 || $2==16 || $2==2048 || $2==2064 {print $4}' |
        sort -n |
        awk -v seq_len="$seq_len" '{print} END {print seq_len}' |
        tr "\n" " " |
        awk '{for(i=$1; i<=$2; i++) print 0
            for(i=$3; i>$2; i--) print 2
            for(i=$3+1; i <=$4; i++) print 0}' |
        awk '{print NR","$1}'
    else
    # deletion/knockin/point mutation
    paste "${tmp_strecher_wt}" "${tmp_strecher_label}" |
        awk '$1=="-" {print $0, NR; next}
            {refrow++; querow=NR; print refrow"_"$0, querow}' |
        sort |
        join -a 1 - "${tmp_control}" |
        awk '$2=="-" {$4=1}1' | # deletion
        awk '$2~/[ACGT]/ {$4=0}1' | # point mutation
        awk '$1=="-" {$4=1}1' | # knockin
        sort -t " " -k 3,3n |
        awk '{print $3","$4}'
    fi |
    cat > ".DAJIN_temp/clustering/temp/control_score_${label}"
done

# #!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# # このcontrol_scoreはいらなくなるかも. rsdで全部代用できる可能性がある. 
# [ _"${mutation_type}" = "_S" ] &&
# mv ".DAJIN_temp/clustering/temp/control_score_target" "${control_score}"
# #!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#===========================================================
#? Generate mutation scoring of other alleles
#===========================================================

find .DAJIN_temp/fasta/ -type f |
    grep -v wt.fa |
    grep -v fasta.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
while read -r label; do
    Rscript DAJIN/src/test_control_scoring_others.R \
    ".DAJIN_temp/clustering/temp/control_score_${label}" "${threads}" #! RENAME <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
done

exit 0