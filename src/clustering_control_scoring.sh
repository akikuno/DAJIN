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

# control=barcode41
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

mutation_type=$(
    minimap2 -ax splice \
        .DAJIN_temp/fasta/wt.fa \
        .DAJIN_temp/fasta/target.fa \
        --cs 2>/dev/null |
    grep -v "^@" |
    awk '{
        cstag=$(NF-1)
        if(cstag ~ "~") print "D"
        else if(cstag ~ "\+") print "I"
        else if(cstag ~ "\*") print "S"
        }' 2>/dev/null
)

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/"
control_score=".DAJIN_temp/clustering/temp/control_score_${alleletype}"

#===========================================================
#? Temporal
#===========================================================

MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${control}_${alleletype}"
tmp_MIDS=".DAJIN_temp/clustering/temp/tmp_MIDS_${control}_${alleletype}"
tmp_mask=".DAJIN_temp/clustering/temp/tmp_mask_${control}_${alleletype}"
tmp_control=".DAJIN_temp/clustering/temp/tmp_control_${control}_${alleletype}"
tmp_strecher=".DAJIN_temp/clustering/temp/tmp_strecher_${control}_${alleletype}"
tmp_strecher_wt=".DAJIN_temp/clustering/temp/tmp_strecher_wt_${control}_${alleletype}"
tmp_strecher_label=".DAJIN_temp/clustering/temp/tmp_strecher_label${control}_${alleletype}"

################################################################################
#! MIDS conversion for clustering
################################################################################

#===========================================================
#? MIDS conversion
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

################################################################################
#! Define sequence error
################################################################################

#===========================================================
#? Mask repeat sequences
#===========================================================

cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    sed "s/A[CGT]A/AAA/g" |
    sed "s/C[AGT]C/CCC/g" |
    sed "s/G[ACT]G/GGG/g" |
    sed "s/T[ACG]T/TTT/g" |
    #
    sed "s/\(AA*\)/\1 /g" |
    sed "s/\(CC*\)/\1 /g" |
    sed "s/\(TT*\)/\1 /g" |
    sed "s/\(GG*\)/\1 /g" |
    #
    awk '{for(i=1;i<=NF;i++) if(length($i) > 5) $i=tolower($i)}1' |
    sed "s/ //g" |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
cat > "${tmp_mask}"_

cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    paste - "${tmp_mask}"_ |
    awk '$2~/[acgt]/ {$1=tolower($1)} {print $1}' |
cat > "${tmp_mask}"

#----------------------------------------
#? Sequence error detection
# Sequence error is defined by 90% or less for M or 5% or more for each IDS.
#----------------------------------------

Rscript DAJIN/src/clustering_control_scoring.R "${tmp_MIDS}" "${tmp_mask}" "${control_score}" "${threads}"

################################################################################
#! Generate scores of other possible alleles
################################################################################

cat ".DAJIN_temp/fasta/${alleletype}.fa" |
    sed 1d |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    paste - "${control_score}" |
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
        minimap2 -ax splice .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/${label}.fa 2>/dev/null |
        grep -c -v "^@"
        )"

    if [ "${splice_num}" -eq 3 ]; then
    # inversion
    set $(
        paste "${tmp_strecher_wt}" "${tmp_strecher_label}" |
        awk '$1=="-" {print "-"; next}
            {refrow++; print $1, $NF}' |
        awk '$2=="-" {$1="-"}1' |
        awk '{printf $1}' |
        sed "s/-\([ACTG].*\)-/-\n\1\n/g" | # flox-ki
        sed "s/--*$//g" | # flox-ki
        awk 'NR==2{gsub(".", "-", $0)} {printf $0}' | # flox-ki
        awk '{mutlen=gsub("-", " ", $0)
            print length($1)+1,length($1)+mutlen+2}'
        )

    cat "${tmp_control}" |
        sort -k 1,1n |
        awk '{print $NF}' |
        sed -e "$1i @" |
        sed -e "$2i @" |
        tr -d "\n" |
        sed "s/@/\n/g" |
        awk -F "" 'NR==2{
            seq=""
            for(i=NF; i>0; i--) seq=seq""$i
            print seq; next}1' |
        tr -d "\n" |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}'
    else
    # deletion/knockin/point mutation
    paste "${tmp_strecher_wt}" "${tmp_strecher_label}" |
        awk '$1=="-" {print $0, NR; next}
            {refrow++; querow=NR; print refrow"_"$0, querow}' |
        sort |
        join -a 1 - "${tmp_control}" |
        awk '$2!="-"' | # deletion
        awk '$2~/[ACGT]/ {$4=1}1' | # point mutation
        awk '$1=="-" {$4=100}1' | # knockin
        sort -t " " -k 3,3n |
        awk '{print $NF}'
    fi |
    cat > ".DAJIN_temp/clustering/temp/control_score_${label}"

done

[ _"${mutation_type}" = "_S" ] &&
mv ".DAJIN_temp/clustering/temp/control_score_target" "${control_score}"


exit 0