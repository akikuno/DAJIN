#!/bin/bash

set -eu

mkdir -p .DAJIN_dev
threads=28

#===========================================================
# Split fasta
#===========================================================

fasta="DAJIN/tests/sv-alleles/cables2-wt.fa"
grep ^ "$fasta" |
    paste - - |
    while read -r line; do
        output="$(echo ${line#>} | cut -d " " -f 1)".fa
        echo "$line" | tr "\t" "\n" >DAJIN/tests/sv-alleles/fasta/"$output"
    done

#===========================================================
# NanoSim
#===========================================================

gzip -dc DAJIN/example/fastq/barcode01.fq.gz |
    cut -d " " -f 1 |
    paste - - - - |
    cut -f 1,2 |
    tr "\t" "\n" >DAJIN/tests/sv-alleles/fasta/barcode01.fa

read_analysis.py genome \
    -i DAJIN/tests/sv-alleles/fasta/barcode01.fa \
    -rg DAJIN/tests/sv-alleles/fasta/wt.fa \
    -t ${threads:-1} \
    -o DAJIN/tests/sv-alleles/nanosim/ 1>&2

wt_seqlen=$(grep -v ">" DAJIN/tests/sv-alleles/fasta/wt.fa | wc -c)

for input in DAJIN/tests/sv-alleles/fasta/wt*.fa*; do
    echo "${input} is now simulating..." 1>&2
    output=$(echo "${input%.*}" | sed "s|/fasta/|/fasta-simulated/|")
    ## For deletion allele
    input_seqlength=$(awk '!/[>|@]/ {print length-100}' "${input}")
    if [ "$input_seqlength" -lt "$wt_seqlen" ]; then
        len="${input_seqlength}"
    else
        len="${wt_seqlen}"
    fi
    ## Simulate!!
    ./DAJIN/utils/NanoSim/src/simulator.py genome \
        -dna_type linear \
        -c DAJIN/tests/sv-alleles/nanosim/ \
        -rg "${input}" \
        -n 1000 \
        -t "${threads:-1}" \
        -min "${len}" \
        -o "${output}" 1>&2
done
rm DAJIN/tests/sv-alleles/fasta-simulated/*_error_* 2>/dev/null || :
rm DAJIN/tests/sv-alleles/fasta-simulated/*_unaligned_* 2>/dev/null || :
rm -rf DAJIN/tests/sv-alleles/nanosim || :
rm -rf DAJIN/utils/NanoSim/src/__pycache__ || :
