#!/bin/sh

set -eu

control=$1
threads=$2

error_exit() {
  echo "$@" 1>&2
  exit 1
}

#===========================================================
# NanoSim
#===========================================================

if ! read_analysis.py -v 2>/dev/null | grep -q "NanoSim 3.0.0"; then
  mamba install -y nanosim=3.0.0 >/dev/null 2>&1
fi

read_analysis.py genome \
  -i ".DAJIN_temp/fasta_ont/${control}.fa" \
  -rg .DAJIN_temp/fasta_conv/wt.fa \
  -t ${threads:-1} \
  -o .DAJIN_temp/NanoSim/training 1>&2 ||
  error_exit 'nanosim error...'

wt_seqlen=$(awk '!/^>/ {print length}' .DAJIN_temp/fasta_conv/wt.fa)

for input in .DAJIN_temp/fasta_conv/*; do
  echo "${input} is now simulating..." 1>&2
  output=$(echo "${input%.*}" | sed "s;fasta_conv;fasta_ont;g")
  ## For deletion allele
  input_seqlength=$(awk '!/[>|@]/ {print length-100}' "${input}")
  if [ "$input_seqlength" -lt "$wt_seqlen" ]; then
    len="${input_seqlength}"
  else
    len="${wt_seqlen}"
  fi
  ##
  simulator.py genome \
    -dna_type linear \
    -c .DAJIN_temp/NanoSim/training \
    -r "${input}" \
    -n 10000 \
    -t "${threads:-1}" \
    -min "${len}" \
    -o "${output}_simulated" 1>&2
  ##
  rm .DAJIN_temp/fasta_ont/*_error_* .DAJIN_temp/fasta_ont/*_unaligned_* 2>/dev/null || true
done

rm -rf DAJIN/utils/NanoSim/src/__pycache__ || true
