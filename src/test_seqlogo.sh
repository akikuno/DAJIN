#!/bin/sh

barcode=${1}
[ -z "${barcode}" ] && exit 1
# printf "${barcode} is now processing...\n"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract 200 bp where includes joint sequence
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input="bam/${barcode}.bam"
# output=".tmp_/mutation_locus_${barcode}"
output=".tmp_/mutsite_split_${barcode}"
# ----------------------------------------
samtools view ${input} |
sort |
join -1 1 - -2 2 .tmp_/sorted_prediction_result |
awk '$2==0 || $2==16' |
awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i,$10}' |
# Remove deletion
sed "s/\([-|+|=|*]\)/ \1/g" |
awk '{seq="";
    for(i=1; i<=NF-1;i++) if($i !~ /-/) seq=seq$i
    print seq, $NF}' |
# Extract joint sequence
awk '{match($1, "~"); seq=substr($1,RSTART-50,50)
    gsub("[-|+|=]", "", seq)
    gsub("*[a-z]", "", seq)
    seq=toupper(seq)
    match($2,seq); seq=substr($2,RSTART+length(seq)-100,200)
    print ">"NR"\n"seq}' \
> ${output}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split fasta files for the following alignment
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input=".tmp_/mutsite_split_${barcode}"
output=".tmp_/split_${barcode}"
# ----------------------------------------
rm -rf ${output} 2>/dev/null
mkdir -p ${output}
split -l 2 ${input} ${output}"/split_"
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Align reads to joint sequence
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input=".tmp_/split_${barcode}"
output=".tmp_/tmp_lalign_${barcode}"
# ----------------------------------------
#
find ${input}/ -name split_* -type f |
xargs -I @ ./DAJIN/src/test_intact_lalign.sh \
    .tmp_/mutation.fa @ > ${output}
#
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract reads with more than 50% of alignment score
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input=".tmp_/tmp_lalign_${barcode}"
output=".tmp_/tmp_lalign_${barcode}.fa"
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
mut_length=$(cat .tmp_/mutation.fa |
sed 1d |
awk '{print length($0)}')
#
cat ${input} |
awk -v per=${per} '$1>per' |
cut -d " " -f 2- |
sed "s/ /\n/g" |
awk -v mut_len=${mut_length} '{print substr($0,0,mut_len)}' \
> ${output} # 1>/dev/null
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract probable intact or non-intact reads
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input=".tmp_/tmp_lalign_${barcode}.fa"
output_intact=".tmp_/lalign_intact_${barcode}.fa"
output_nonintact=".tmp_/lalign_nonintact_${barcode}.fa"
# ----------------------------------------
#
joint_seq=$(cat .tmp_/mutation.fa | 
sed 1d |
awk -F "" '{print substr($0,int(NF)/2-4,8)}')
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
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output read alignment and ratio of intact/nonintact
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
output_alignment=".tmp_/numseq_alignment_${barcode}"
output_intact_ratio=".tmp_/numseq_intact_${barcode}"
# ----------------------------------------
true > ${output_alignment}
true > ${output_intact_ratio}
#
cat ".tmp_/mutsite_split_${barcode}" | grep -v "^>" | wc -l >> ${output_alignment}
cat ".tmp_/tmp_lalign_${barcode}" | grep -v "^>" | wc -l >> ${output_alignment}
cat ${output_intact} | grep -v "^>" | wc -l >> ${output_intact_ratio}
cat ${output_nonintact} | grep -v "^>" | wc -l >> ${output_intact_ratio}
#
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Report the results of Sequence logo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
arg1=".tmp_/mutation.fa"
arg2=${output_intact}
arg3=${output_nonintact}
arg4=${output_alignment}
arg5=${output_intact_ratio}
# ----------------------------------------
python DAJIN/src/test_logomaker.py ${arg1} ${arg2} ${arg3} ${arg4} ${arg5}

exit 0
