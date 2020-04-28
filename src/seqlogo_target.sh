#!/bin/sh

# ============================================================================
# I/O and Arguments
# ============================================================================
# mkdir -p .DAJIN_temp/seqlogo/
# gRNAのデータを保存する
# grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
# grna=$(echo "${grna}" |
#     DAJIN/src/revcomp.sh - |
#     sed "s/^/${grna},/g")

barcode=barcode14
alleletype=target
suffix="${barcode}_${alleletype}"

fasta=.DAJIN_temp/fasta_ont/"${barcode}".fa
clustering_id=$(find .DAJIN_temp/clustering/result_allele_id* | grep "$suffix")

tmp_fa=".DAJIN_temp/seqlogo/temp/${suffix}.fa"
tmp_targetseq=".DAJIN_temp/seqlogo/temp/targetseq_${suffix}"
tmp_lalign=".DAJIN_temp/seqlogo/temp/lalign_${suffix}"

# ============================================================================
# FASTAリードにクラスタ番号をラベルづけする
# ============================================================================

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
> ${tmp_fa}

# ============================================================================
# Align reads to target sequence
# ============================================================================

mutation_type=$(
    minimap2 -ax map-ont \
    .DAJIN_temp/fasta/target.fa \
    .DAJIN_temp/fasta/wt.fa 2>/dev/null |
    grep -v "^@" |
    cut -f 6 |
    awk '{if($0~"I") print "D"
        else if($0~"D") print "I"
        else if($0~"S") print "P"
        }'
)
# ------------------------------------------------------------
# 2cut deletionの場合
# ------------------------------------------------------------

# 結合部から±25塩基を抽出する
# if [ "${alleletype}" = "target" ]; then
minimap2 -ax map-ont \
    .DAJIN_temp/fasta_conv/wt.fa \
    .DAJIN_temp/fasta_conv/target.fa --cs=long 2>/dev/null |
grep -v "^@" |
awk '{print $(NF-1)}' |
sed "s/[a-z]*=//g" |
awk '{match($0, "-")
print substr($0, RSTART-25, 51)
}' |
sed "s/-//g" \
> "$tmp_targetseq"

seqlen=$(awk '!/[>|@]/ {print length($0)}' .DAJIN_temp/fasta_conv/target.fa)

true > test_cstag
cat "${tmp_targetseq}" |
while read -r input; do
    label=$(echo "${input}" | cut -d " " -f 1)
    echo "$input" |
    awk -v tar="${seqlen}" \
    '{  s=sprintf("%."tar"d","0")
        gsub("0","N",s)
        print $1,$2s}' |
    sed "s/ /\n/g" \
    > .DAJIN_temp/seqlogo/tmp_"${suffix}"
    #
    minimap2 -ax map-ont .DAJIN_temp/seqlogo/tmp_"${suffix}" \
    ${tmp_fa} --cs=long 2>/dev/null |
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
    >> test_cstag
done

true > "${tmp_lalign}"
cat test_cstag |
cut -d " " -f 1,3 | sort -u |
while read -r input; do
    cl=$(echo "$input" | cut -d " " -f 1)
    target=$(echo "$input" | cut -d " " -f 2)
    label="${cl}"@"${target}"
    cat "${tmp_targetseq}" |
    grep "${target}" - |
    sed "s/ /\n/g" \
    > test_targetseq.fa
    #
    cat test_targetseq.fa | ./DAJIN/src/revcomp.sh - \
    > test_revcomp_targetseq.fa
    #
    awk -v cl="${cl}" -v target="${target}"\
    '$1==cl && $3==target' test_cstag |
    # head -n 3 |
    awk '{print ">"$1"@"$3"\n"$NF}' \
    > test.fa
    #
    cat test.fa | ./DAJIN/src/revcomp.sh - \
    > test_revcomp.fa
    # ============================================================================
    # Split fasta files for the following alignment
    # ============================================================================
    tmp_split_dir=".DAJIN_temp/seqlogo/split_${suffix}"
    # ----------------------------------------
    rm -rf ${tmp_split_dir} 2>/dev/null
    mkdir -p ${tmp_split_dir}
    split -l 2 test.fa ${tmp_split_dir}"/split_"
    split -l 2 test_revcomp.fa ${tmp_split_dir}"/revcomp_split_"
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #
    time find ${tmp_split_dir}/ -name revcomp_split_* -type f | head -n 100 |
    xargs -I @ ./DAJIN/src/test_intact_lalign.sh \
        test_targetseq.fa @ |
    awk -v label="${label}" '{$2=label; print}' \
    > test_lalign
    #
    time find ${tmp_split_dir}/ -name split_* -type f | head -n 100 |
    xargs -I @ ./DAJIN/src/test_intact_lalign.sh \
        test_revcomp_targetseq.fa @ |
    awk -v label="${label}" '{$2=label; print ">"$1"#"$2"\n"$3}' |
    ./DAJIN/src/revcomp.sh - |
    sed -e "s/^/#/g" |
    tr -d "\n" |
    sed "s/#>/\n>/g" | 
    sed "s/#/ /g" \
    > test_revcomp_lalign

    #
    cat test*_lalign > test_lalign_merge
    cat test_lalign > test_lalign_merge
    cat test_revcomp_lalign > test_lalign_merge
    percentile=0.5
    per=$(cat test_lalign_merge |
    awk -v per=${percentile} '{
        percentile[NR]=$1}
        END{asort(percentile)
        print percentile[int(NR*per)]
    }')
    #
    mut_length=$(awk '!/[>|@]/ {print length($0)}' test_targetseq.fa)
    #
    cat test_lalign_merge |
    awk -v per=${per} '$1>per' |
    cut -d " " -f 3 |
    awk -v mut_len=${mut_length} '{print substr($0,0,mut_len)}' \
    > ${tmp_lalign} # 1>/dev/null
done

#* ========================================


true > "${tmp_lalign}"
cat "${tmp_targetseq}" |
while read -r input; do
    label=$(echo "${input}" | cut -d " " -f 1)
    echo "${input}" | sed "s/ /\n/g" \
    > .DAJIN_temp/seqlogo/tmp_"${suffix}"
    #
    find ${tmp_split_dir}/ -name split_* -type f |
    xargs -I @ ./DAJIN/src/test_intact_lalign.sh \
        .DAJIN_temp/seqlogo/tmp_"${suffix}" @ |
    sed "s/^/${label} /g" \
    >> "${tmp_lalign}"
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
tmp_intact=".tmp_/lalign_intact_${barcode}.fa"
tmp_nonintact=".tmp_/lalign_nonintact_${barcode}.fa"
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
> ${tmp_intact}
#
cat ${input} |
sed -e "s/>/HOGE>/g" -e "s/$/FUGA/g" |
tr -d "\n" |
sed -e "s/HOGE/\n/g" -e "s/FUGA/\t/g" | 
grep -v ${joint_seq} |
sed "s/\t/\n/g" |
grep -v "^$" \
> ${tmp_nonintact}
#
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output read alignment and ratio of intact/nonintact
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
tmp_alignment=".tmp_/numseq_alignment_${barcode}"
tmp_intact_ratio=".tmp_/numseq_intact_${barcode}"
# ----------------------------------------
true > ${tmp_alignment}
true > ${tmp_intact_ratio}
#
cat ".tmp_/mutsite_split_${barcode}" | grep -v "^>" | wc -l >> ${tmp_alignment}
cat ".tmp_/tmp_lalign_${barcode}" | grep -v "^>" | wc -l >> ${tmp_alignment}
cat ${tmp_intact} | grep -v "^>" | wc -l >> ${tmp_intact_ratio}
cat ${tmp_nonintact} | grep -v "^>" | wc -l >> ${tmp_intact_ratio}
#
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Report the results of Sequence logo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
arg1=".tmp_/mutation.fa"
arg2=${tmp_intact}
arg3=${tmp_nonintact}
arg4=${tmp_alignment}
arg5=${tmp_intact_ratio}
# ----------------------------------------
python DAJIN/src/test_logomaker.py ${arg1} ${arg2} ${arg3} ${arg4} ${arg5}

exit 0
