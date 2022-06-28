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

target_mutation_type=$(cat .DAJIN_temp/target_mutation_type)

################################################################################
# Train models
################################################################################

[ "_${target_mutation_type}" = "_S" ] && rm .DAJIN_temp/data/MIDS_target* 2>/dev/null || true

cat .DAJIN_temp/data/MIDS_* |
  grep "_sim" |
  sed -e "s/_aligned_reads//g" |
  cat >".DAJIN_temp/data/DAJIN_MIDS_sim.txt"

cat .DAJIN_temp/data/MIDS_"${control}"_wt |
  grep -v "IIIIIIIIII" |
  grep -v "DDDDDDDDDD" |
  grep -v "SSSSSSSSSS" |
  head -n 10000 |
  sed "s/${control}$/wt_simulated/g" |
  cat >>".DAJIN_temp/data/DAJIN_MIDS_sim.txt"

echo "Model training..." >&2
python ./DAJIN/src/ml_simulated.py \
  ".DAJIN_temp/data/DAJIN_MIDS_sim.txt" "${threads}" >&2

################################################################################
# Predict allele type
################################################################################
tmp_prediction=".DAJIN_temp/data/tmp_DAJIN_MIDS_prediction_result.txt"

true >"${tmp_prediction}"

find .DAJIN_temp/data/MIDS* |
  grep -v sim |
  sort |
  while read -r input; do
    barcode=$(echo "${input%_*}" | sed "s/.*MIDS_//")
    echo "Prediction of ${barcode} is now processing..." >&2
    python ./DAJIN/src/ml_real.py "${input}" "${target_mutation_type}" "${threads}" ||
      exit 1
  done

cat "${tmp_prediction}" |
  sort |
  cat

rm "${tmp_prediction}"
