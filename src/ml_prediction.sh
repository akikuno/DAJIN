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

#===============================================================================
#? TEST Aurguments
#===============================================================================
# barcode="barcode32"
# alleletype="abnormal"

#===========================================================
#? Auguments
#===========================================================

control="${1}"
threads="${2}"

#===========================================================
#? Input
#===========================================================

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

#===========================================================
#? Temporal
#===========================================================
tmp_prediction=".DAJIN_temp/data/tmp_DAJIN_MIDS_prediction_result.txt"

################################################################################
#! Train models
################################################################################

cat .DAJIN_temp/data/MIDS_* |
    grep "_sim" |
    sed -e "s/_aligned_reads//g" |
cat > ".DAJIN_temp/data/DAJIN_MIDS_sim.txt"

cat .DAJIN_temp/data/MIDS_"${control}"_wt |
    grep -v "IIIIIIIIII" |
    grep -v "DDDDDDDDDD" |
    grep -v "SSSSSSSSSS" |
    head -n 10000 |
    sed "s/${control}$/wt_simulated/g" |
cat >> ".DAJIN_temp/data/DAJIN_MIDS_sim.txt"

python ./DAJIN/src/ml_simulated.py \
    ".DAJIN_temp/data/DAJIN_MIDS_sim.txt" "${threads}" >&2

################################################################################
#! Predict allele type
################################################################################
true > "${tmp_prediction}"

find .DAJIN_temp/data/MIDS* |
    grep -v sim |
    sort |
while read -r input; do
    barcode=$(echo "${input}" | cut -d "_" -f 3)
    echo "Prediction of ${barcode} is now processing..." >&2

    python ./DAJIN/src/ml_real.py \
        "${input}" "${mutation_type}" "${threads}" ||
    exit 1
done

cat "${tmp_prediction}" |
    sort |
cat

rm "${tmp_prediction}"