#!/bin/sh

################################################################################
# Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C

################################################################################
# I/O naming
################################################################################

#===========================================================
# Arguments
#===========================================================

genome="${1}"
threads="${2}"

#===========================================================
# Variable
#===========================================================

target_mutation_type=$(cat .DAJIN_temp/target_mutation_type)

#===========================================================
# Setting directory
#===========================================================

rm -rf .DAJIN_temp/bam/ 2>/dev/null
mkdir -p .DAJIN_temp/bam/temp .DAJIN_temp/bam/reads100

################################################################################
# Generate BAM files
################################################################################

if [ "_$target_mutation_type" = "_S" ]; then
    mv .DAJIN_temp/fasta_ont/wt_ins* .DAJIN_temp/
    mv .DAJIN_temp/fasta_ont/wt_del* .DAJIN_temp/
fi

if [ "${genome}" = "no_ref" ]; then
    for input in .DAJIN_temp/fasta_ont/*; do
        output=$(
            basename "${input}" |
                sed -e "s|^|.DAJIN_temp/bam/|" -e "s|\.f.*$|.bam|g"
        )
        minimap2 -t "${threads:-1}" -ax map-ont --cs=long .DAJIN_temp/fasta/wt.fa "${input}" 2>/dev/null |
            samtools sort -@ "${threads:-1}" - 2>/dev/null |
            cat >"${output}"
        samtools index -@ "${threads:-1}" "${output}"
    done
    cp .DAJIN_temp/fasta/wt.fa .DAJIN_temp/bam/reference.fa
else
    ./DAJIN/src/mapping.sh "${genome:-mm10}" "${threads:-1}" || exit 1
fi

if [ "_$target_mutation_type" = "_S" ]; then
    mv .DAJIN_temp/wt_ins* .DAJIN_temp/fasta_ont/
    mv .DAJIN_temp/wt_del* .DAJIN_temp/fasta_ont/
fi

#===========================================================
# Generate BAM files on each cluster
#===========================================================

cat .DAJIN_temp/clustering/allele_per/label* |
    awk '{nr[$1]++; print $0, nr[$1]}' |
    while read -r allele; do
        barcode=$(echo ${allele} | cut -d " " -f 1)
        alleletype=$(echo ${allele} | cut -d " " -f 2)
        cluster=$(echo ${allele} | cut -d " " -f 3)
        alleleid=$(echo ${allele} | cut -d " " -f 5)
        #
        input_bam="${barcode}_${alleletype}"
        output_bam="${barcode}_allele${alleleid}"
        #
        find .DAJIN_temp/clustering/allele_per/readid_cl_mids* |
            grep "${input_bam}" |
            xargs cat |
            awk -v cl="${cluster}" '$2==cl' |
            cut -f 1 |
            sort |
            cat >".DAJIN_temp/bam/temp/tmp_id_$$"

        samtools view -h ".DAJIN_temp/bam/${barcode}".bam |
            awk '/^@/{print}
            NR==FNR{a[$1];next}
            $1 in a' \
                .DAJIN_temp/bam/temp/tmp_id_$$ - |
            samtools sort -@ "${threads:-1}" 2>/dev/null |
            cat >.DAJIN_temp/bam/"${output_bam}".bam
        samtools index .DAJIN_temp/bam/"${output_bam}".bam
    done

#===========================================================
# Generate BAM files with 100 reads
#===========================================================

find .DAJIN_temp/bam/ -name "*bam" -type f |
    grep -v "reads100" |
    while read -r input_bam; do

        output_bam=$(
            echo "${input_bam}" |
                sed "s%.DAJIN_temp/bam/%.DAJIN_temp/bam/reads100/%g"
        )

        header_num=$(samtools view -H "${input_bam}" | wc -l)
        bam_num=$((100 + header_num))

        samtools view -h "${input_bam}" |
            awk '$1 ~ /^@/ || $6 != "*"' |
            head -n "${bam_num}" |
            samtools sort -@ "${threads:-1}" 2>/dev/null |
            cat >"${output_bam}"
        samtools index "${output_bam}"
    done
