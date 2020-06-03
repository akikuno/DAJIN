#!/bin/sh
# Original:"test_seqlogo.sh"

barcode=${1}
[ -z "${barcode}" ] && exit 1
# printf "${barcode} is now processing...\n"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract 200 bp where may include cutting sites
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input="bam/${barcode}.bam"
output_left=".tmp_/cutting_${barcode}_left"
output_right=".tmp_/cutting_${barcode}_right"
# ----------------------------------------
start=$(cat .tmp_/gggenome_location | cut -f 2)
cut_left=$(cat .tmp_/cutting_sites | cut -d " " -f 1 | head -n 1)
cut_right=$(cat .tmp_/cutting_sites | cut -d " " -f 1 | tail -n 1)
#echo $cut_left $cut_right
#
samtools view ${input} |
sort |
join -1 1 - -2 2 .tmp_/sorted_prediction_result |
# #!--------------------
# grep -e "ACATAAGGGCAGATGTGATGCTCGGCTGTGGCTTGGTGACAGAGCCCTTG" \
#     -e "CAAGGGCTCTGTCACCAAGCCACAGCCGAGCATCACATCTGCCCTTATGT" \
#     -e "GACTCAGACATAAAGTGGTTGCGCTCTTGGGAACAACGCTGTTTTTAAGG" \
#     -e "CCTTAAAAACAGCGTTGTTCCCAAGAGCGCAACCACTTTATGTCTGAGTC" |
# #!--------------------
awk '$2==0 || $2==16 {print $4,$6,$10}' |
# Remove softclipping  ------------------------------------
awk '{sub("S.*","",$2)
    print $1,substr($3,$2+1)}' |
# add "N" to compensate starting position ------------------------------------
awk -v start=${start} \
    '{seq=""
    st=$1-start+1
    for(i=1;i<st;i++) seq=seq"N"
    print seq$2
    }' |
# Extract joint sequence ------------------------------------
awk -v left=${cut_left} -v right=${cut_right} \
    -v out_l=${output_left} -v out_r=${output_right} \
    '{left_seq=substr($1,left-100,200)
    right_seq=substr($1,right-100,200)
    print ">left"NR"@" "\n" left_seq > out_l
    print ">right"NR"@" "\n" right_seq > out_r}'

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split fasta files for the following alignment
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
for site in left right; do
    # ----------------------------------------
    input=".tmp_/cutting_${barcode}_${site}"
    output=".tmp_/split_${barcode}_${site}"
    # ----------------------------------------
    rm -rf ${output} 2>/dev/null
    mkdir -p ${output}
    split -l 2 ${input} ${output}"/split_"
done
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Align reads to joint sequence
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for site in left right; do
    # ----------------------------------------
    input=".tmp_/split_${barcode}_${site}"
    output=".tmp_/lalign_${barcode}_${site}"
    # ----------------------------------------
    find ${input}/ -name split_* -type f |
    xargs -I @ ./DAJIN/src/test_intact_lalign.sh \
        .tmp_/cutting_sites_${site}.fa @ \
        > ${output}
done
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract reads with more than 50% of alignment score
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for site in left right; do
    # ----------------------------------------
    input=".tmp_/lalign_${barcode}_${site}"
    output=".tmp_/lalign_${barcode}_${site}.fa"
    percentile=0.5
    # ----------------------------------------
    #
    per=$(cat ${input} |
        awk -v per=${percentile} '{
            percentile[NR]=$1}
            END{asort(percentile)
            print percentile[int(NR*per)]
        }')
    #
    mut_length=$(cat .tmp_/cutting_sites_${site}.fa |
    sed 1d |
    awk '{print length($0)}')
    #
    cat ${input} |
    awk -v per=${per} '$1>=per' |
    cut -d " " -f 2- |
    sed "s/ /\n/g" |
    awk -v mut_len=${mut_length} \
    '{print substr($0,0,mut_len)}' \
    > ${output} # 1>/dev/null
done
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract probable intact or non-intact reads
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for site in left right; do
    # ----------------------------------------
    input=".tmp_/lalign_${barcode}_${site}.fa"
    output_intact=".tmp_/lalign_intact_${barcode}_${site}.fa"
    output_nonintact=".tmp_/lalign_nonintact_${barcode}_${site}.fa"
    # ----------------------------------------
    output_alignment=".tmp_/numseq_alignment_${barcode}_${site}"
    output_intact_ratio=".tmp_/numseq_intact_${barcode}_${site}"
    # ----------------------------------------
    #
    joint_seq=$(cat .tmp_/cutting_sites_${site}.fa | 
    sed 1d |
    awk -F "" '{print substr($0,int(NF)/2-2,6)}')
    #
    cat ${input} |
    sed -e "s/>/HOGE>/g" -e "s/$/FUGA/g" |
    tr -d "\n" |
    sed -e "s/HOGE/\n/g" -e "s/FUGA/\t/g" | 
    grep ${joint_seq} |
    sed "s/\t/\n/g" |
    grep -v "^$" \
    > ${output_intact}
    #
    cat ${input} |
    sed -e "s/>/HOGE>/g" -e "s/$/FUGA/g" |
    tr -d "\n" |
    sed -e "s/HOGE/\n/g" -e "s/FUGA/\t/g" | 
    grep -v ${joint_seq} |
    sed "s/\t/\n/g" |
    grep -v "^$" \
    > ${output_nonintact}
    #
done
#
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output read alignment and ratio of intact/nonintact
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for site in left right; do
    # ----------------------------------------
    output_intact=".tmp_/lalign_intact_${barcode}_${site}.fa"
    output_nonintact=".tmp_/lalign_nonintact_${barcode}_${site}.fa"
    # ----------------------------------------
    output_alignment=".tmp_/numseq_alignment_${barcode}_${site}"
    output_intact_ratio=".tmp_/numseq_intact_${barcode}_${site}"
    # ----------------------------------------
    #
    true > ${output_alignment}
    true > ${output_intact_ratio}
    #
    cat ".tmp_/cutting_${barcode}_${site}" | grep -v "^>" | wc -l \
    >> ${output_alignment}
    #
    cat ".tmp_/lalign_${barcode}_${site}" | grep -v "^>" | wc -l \
    >> ${output_alignment}
    #
    cat ${output_intact} | grep -v "^>" | wc -l \
    >> ${output_intact_ratio}
    #
    cat ${output_nonintact} | grep -v "^>" | wc -l \
    >> ${output_intact_ratio}
done
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Report the results of Sequence logo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
arg1=".tmp_/cutting_sites_${site}.fa"
arg2=${output_intact}
arg3=${output_nonintact}
arg4=${output_alignment}
arg5=${output_intact_ratio}
# ----------------------------------------
python DAJIN/src/test_logomaker.py ${arg1} ${arg2} ${arg3} ${arg4} ${arg5}

#exit 0
