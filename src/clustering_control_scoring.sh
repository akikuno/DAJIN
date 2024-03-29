#!/bin/sh

################################################################################
# Initialize shell environment
################################################################################

set -eu
umask 0022
export LC_ALL=C

################################################################################
# I/O naming
################################################################################

#===========================================================
# Auguments
#===========================================================

control="${1}"
threads="${2}"

#===========================================================
# Input
#===========================================================

alleletype="wt"

#===========================================================
# Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/"

#===========================================================
# Temporal
#===========================================================

MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${control}_${alleletype}".csv
tmp_MIDS=".DAJIN_temp/clustering/temp/tmp_MIDS_${control}_${alleletype}".csv
tmp_control=".DAJIN_temp/clustering/temp/tmp_control_${control}_${alleletype}".csv
tmp_strecher=".DAJIN_temp/clustering/temp/tmp_strecher_${control}_${alleletype}".csv
tmp_strecher_wt=".DAJIN_temp/clustering/temp/tmp_strecher_wt_${control}_${alleletype}".csv
tmp_strecher_label=".DAJIN_temp/clustering/temp/tmp_strecher_label${control}_${alleletype}".csv

################################################################################
# Generate mutation scores of WT alleles
################################################################################

#===========================================================
# MIDS conversion for clustering
#===========================================================

./DAJIN/src/mids_clustering.sh "${control}" "${alleletype}" >"${MIDS_ref}"

#===========================================================
# Select WT allele
#===========================================================

cat "${MIDS_ref}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v alelle="$alleletype" '$NF==alelle' |
    cut -d " " -f 2 |
    sed "s/[1-9a-z]/I/g" |
    awk '{n=split($0,array,""); for(i=1;i<=n;i++) printf array[i]","; print ""}' |
    sed "s/,$//" |
    cat >"${tmp_MIDS}"

#===========================================================
# Generate mutation scoring of wt alleles
#===========================================================

Rscript DAJIN/src/clustering_control_scoring_wt.R "${tmp_MIDS}" "${threads}" 2>/dev/null

################################################################################
# Generate mutation scores of other alleles
################################################################################

#===========================================================
# Annotate mutation loci of other possible alleles
# (mutation loci = 1)
#===========================================================

cat ".DAJIN_temp/fasta/${alleletype}.fa" |
    sed 1d |
    awk '{n=split($0,array,""); for(i=1;i<=n;i++) print array[i], 0}' |
    awk '{print NR"_"$1, $2}' |
    sort -t " " -k 1,1 |
    cat >"${tmp_control}"

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
            awk '{n=split($0,array,""); for(i=1;i<=n;i++) print array[i]}' |
            cat >"${tmp_strecher_wt}"

        cat "${tmp_strecher}" |
            awk 'NF==2' |
            awk '$0 ~ /[ACGT.-]/' |
            awk 'NR%2==0 {printf $2}' |
            awk '{n=split($0,array,""); for(i=1;i<=n;i++) print array[i]}' |
            cat >"${tmp_strecher_label}"

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
                awk '$2=="-" {$4=1}1' |     # deletion
                awk '$2~/[ACGT]/ {$4=1}1' | # point mutation
                awk '$1=="-" {$4=1}1' |     # knockin
                sort -t " " -k 3,3n |
                awk '{print $3","$4}'
        fi |
            cat >".DAJIN_temp/clustering/temp/control_score_${label}".csv
    done

#===========================================================
# Generate mutation scoring of other alleles
#===========================================================

find .DAJIN_temp/fasta/ -type f |
    grep -v wt.fa |
    grep -v fasta.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
    while read -r label; do
        Rscript DAJIN/src/clustering_control_scoring_others.R \
            ".DAJIN_temp/clustering/temp/control_score_${label}".csv "${threads}" 2>/dev/null
    done

rm .DAJIN_temp/clustering/temp/tmp_*

exit 0
