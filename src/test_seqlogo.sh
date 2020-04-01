#!/bin/sh

barcode=barcode14
alleletype=target
suffix="${barcode}_${alleletype}"

mkdir -p .DAJIN_temp/seqlogo/
fasta=.DAJIN_temp/fasta_ont/barcode14.fa
clustering_id=.DAJIN_temp/clustering/allele_id_barcode14_target_HOGE
output_fa=.DAJIN_temp/seqlogo/"${suffix}".fa

cat "${fasta}" |
awk '{if(NR % 2 != 0) $1="HOGE"$1
    else $1="FUGA"$1
    print}' |
tr -d "\n" |
sed -e "s/HOGE/\n/g" -e "s/FUGA/ /g" |
awk '{print $1,$NF}' |
sed "s/^[@|>]//g" |
grep -v "^$" |
sort -t " " |
join - "${clustering_id}" |
awk '{print ">"$NF"@@@"$1"\n"$2}' \
> ${output_fa}


# ============================================================================
# Split fasta files for the following alignment
# ============================================================================
# ----------------------------------------
# input=".tmp_/mutsite_split_${barcode}"
output_split_dir=".DAJIN_temp/seqlogo/split_${suffix}"
# ----------------------------------------
rm -rf ${output_split_dir} 2>/dev/null
mkdir -p ${output_split_dir}
split -l 2 ${output_fa} ${output_split_dir}"/split_"
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# ============================================================================
# Align reads to target sequence (gRNA and Target mutation)
# ============================================================================
# input=".tmp_/split_${barcode}"
output_targetseq=".DAJIN_temp/seqlogo/tmp_targetseq_${suffix}"
output_lalign=".DAJIN_temp/seqlogo/tmp_lalign_${suffix}"
# ----------------------------------------

#! gRNAのデータを保存する
# grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
cat << EOF > tmp_grna.fa
CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
EOF
cat tmp_grna.fa | DAJIN/src/revcomp.sh - > tmp_grna_conv.fa

grna=$(cat tmp_grna_conv.fa)

true > "${output_targetseq}"
cat .DAJIN_temp/fasta_conv/${alleletype}.fa |
# cat .DAJIN_temp/fasta_conv/wt.fa |
sed 1d |
awk -v grna="${grna}" \
'{  seq=""
    split(grna,array,",")
    for(key in array) {
        match($0,array[key])
        if(RLENGTH == length(array[key]))
            seq=seq" "substr($0, RSTART-10, RLENGTH+20)
    }
    print seq
}' |
sed -e "s/ /\n/g" |
grep -v "^$" |
awk '{print ">grna_"NR,$0}' \
>> "${output_targetseq}"


if [ "${alleletype}" = "target" ]; then
    minimap2 -ax map-ont .DAJIN_temp/fasta_conv/wt.fa .DAJIN_temp/fasta_conv/target.fa --cs 2>/dev/null |
    awk '{print toupper($(NF-1)),$10}' |
    grep "CS:Z::" |
    sed -e "s/CS:Z:://g" \
        -e "s/[-|+|:|*]/ /g" \
        -e "s/[0-9]*//g" |
    awk '{seq=""
        for(i=1; i<NF; i++){
            match($NF,$i)
            seq=seq" "substr($NF, RSTART-10, RLENGTH+20)
        }
        print seq}' |
    sed -e "s/ /\n/g" |
    grep -v "^$" |
    awk '{print ">target_"NR,$0}' \
    >> "${output_targetseq}"
fi


target_seqlen=$(awk '!/[>|@]/ {print length($0)}' .DAJIN_temp/fasta_conv/target.fa)

true > test.sam
cat "${output_targetseq}" |
while read -r input; do
    label=$(echo "${input}" | cut -d " " -f 1)
    echo "$input" |
    awk -v tar="${target_seqlen}" \
    '{  s=sprintf("%."tar"d","0")
        gsub("0","N",s)
        print $1,$2s}' |
    sed "s/ /\n/g" \
    > .DAJIN_temp/seqlogo/tmp_"${suffix}"
    #
    minimap2 -ax map-ont .DAJIN_temp/seqlogo/tmp_"${suffix}" \
    ${output_fa} --cs=long 2>/dev/null \
    >> test.sam
done

#cat test.sam | grep "cs:" | cut -f 2 | sort | uniq -c
cat test.sam |
awk '$2==0 || $2==16 {print $1, $2, $3, $(NF-1)}' |
awk '{cstag=$NF
    sub("cs:Z:", "", cstag)
    gsub(/\*[a-z]/, "", cstag)
    gsub(/\-[a-z]*[+|*|=]/, "", cstag)
    gsub("[=|+]", "", cstag)
    #
    gsub(/@@@.*/, "", $1)
    print $1, $2, $3, toupper(cstag)
    }' \
> test2.sam

cat test2.sam |
cut -d " " -f 1,3 | sort -u |
while read -r input; do
    cl=$(echo "$input" | cut -d " " -f 1)
    target=$(echo "$input" | cut -d " " -f 2)
    awk -v cl="${cl}" -v target="${target}"\
    '$1==cl && $3==target' test2.sam |
    head -n 5
done

#* ========================================


true > "${output_lalign}"
cat "${output_targetseq}" |
while read -r input; do
    label=$(echo "${input}" | cut -d " " -f 1)
    echo "${input}" | sed "s/ /\n/g" \
    > .DAJIN_temp/seqlogo/tmp_"${suffix}"
    #
    find ${output_split_dir}/ -name split_* -type f |
    xargs -I @ ./DAJIN/src/test_intact_lalign.sh \
        .DAJIN_temp/seqlogo/tmp_"${suffix}" @ |
    sed "s/^/${label} /g" \
    >> "${output_lalign}"
done
#



#? =======================================================

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
