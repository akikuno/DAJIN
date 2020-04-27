#!/bin/sh

# ============================================================================
# Initialize shell environment
# ============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX


error_exit() {
    ${2+:} false && echo "${0##*/}: $2" 1>&2
    exit $1
}

# ============================================================================
# I/O naming
# ============================================================================
# ----------------------------------------
# Input
# ----------------------------------------
# barcode="barcode18"
# alleletype="target"

# control="barcode26"
# alleletype_original=${alleletype}
# suffix="${barcode}"_"${alleletype}"
# echo $suffix
# [ "$alleletype" = "abnormal" ] && alleletype="wt"

barcode="${1}"
control="${2}"
alleletype="${3}"
alleletype_original=${3}
suffix="${barcode}"_"${alleletype}"
[ "$alleletype" = "abnormal" ] && alleletype="wt"

# ----------------------------------------
# Output
# ----------------------------------------
# MIDS conversion
MIDS_que=".DAJIN_temp/clustering/temp/tmp_MIDS_${suffix}"
MIDS_ref=".DAJIN_temp/clustering/temp/tmp_MIDS_${control}_${alleletype}"

# Mutation scoring of samples
output_label=".DAJIN_temp/clustering/temp/query_labels_${suffix}"
output_query_seq=".DAJIN_temp/clustering/temp/query_seq_${suffix}"

# Output Genomic coodinates (Control)
output_ref_score=".DAJIN_temp/clustering/temp/control_score_${suffix}"

# Output Genomic coodinates (Query)
output_query_score=".DAJIN_temp/clustering/temp/query_score_${suffix}"

# Output Plot
hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"
output_plot=".DAJIN_temp/clustering/temp/plot_${suffix}"
plot_mutsites=.DAJIN_temp/clustering/temp/tmp_mutation_"${suffix}"

# -------------------------------
# Report allele mutation info
# -------------------------------
output_alleleper=".DAJIN_temp/clustering/allele_percentage_${suffix}"
output_result=".DAJIN_temp/clustering/result_alleleinfo_${suffix}"

# Get max sequence length
seq_maxnum=$(
    cat .DAJIN_temp/fasta/fasta.fa |
    grep -v "^>" |
    awk '{if(max<length($0)) max=length($0)}
    END{print max}'
    )

# ============================================================================
# MIDS conversion
# ============================================================================

find .DAJIN_temp/fasta_ont/ -type f | grep "${barcode}" |
xargs -I @ ./DAJIN/src/mids_convertion.sh @ "${alleletype}" &&
mv ".DAJIN_temp/data/MIDS_${barcode}_${alleletype}" "${MIDS_que}"

# If no control MIDS files, output... 
if [ ! -s "${MIDS_ref}" ]; then
    find .DAJIN_temp/fasta_ont/ -type f | grep "${control}" |
    xargs -I @ ./DAJIN/src/mids_convertion.sh @ "${alleletype}" "control" &&
    mv ".DAJIN_temp/data/MIDS_${control}_${alleletype}" "${MIDS_ref}"
fi

rm .DAJIN_temp/tmp_${barcode}*${alleletype}
rm .DAJIN_temp/tmp_${control}*${alleletype}

# ============================================================================
# Mutation scoring of samples
# ============================================================================
# ----------------------------------------
# 配列IDとラベルをつくる
# ----------------------------------------

cat "${MIDS_que}" |
grep "${barcode}" |
sort -k 1,1 |
join -1 1 -2 2 - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
awk -v atype="${alleletype_original}" \
'$NF==atype' |
cut -d " " -f 1,3 |
sed "s/ /\t/g" \
> "${output_label}"


# ----------------------------------------
# 挿入塩基を1つの挿入塩基数にまとめて配列のズレを無くす
# ----------------------------------------
cat "${MIDS_que}" |
grep "${barcode}" |
sort -k 1,1 |
join -1 1 -2 2 - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
awk -v atype="${alleletype_original}" \
'$NF==atype' |
cut -d " " -f 2 |
awk -F "" '{
    for(i=1; i<=NF; i++){
        if($i=="I") num=num+1
        if($i=="I" && $(i+1)!="I") {
            ### e.g) if num=10, num becomes "a"
            if(num>=10 && num<=35) num=sprintf("%c", num+87)
            else if(num>=36) num="z"
            ###
            $(i+1)=num; num=0}
        }
    print $0}' |
sed -e "s/I//g" -e "s/ //g" \
> "${output_query_seq}"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output Genomic coodinates (Control)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ----------------------------------------
# 挿入塩基を1つの挿入塩基数にまとめて配列のズレを無くす
# ----------------------------------------
cat "${MIDS_ref}" |
grep "${control}" |
sort -k 1,1 |
join -1 1 -2 2 - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
awk '$NF=="wt"' |
cut -d " " -f 2 |
# Insertion annotation
awk -F "" '{
    for(i=1; i<=NF; i++){
        if($i=="I") num=num+1
        if($i=="I" && $(i+1)!="I") {
            ### e.g) if num=10, num becomes "a"
            if(num>=10 && num<=35) {num=sprintf("%c", num+87)}
            else if(num>=36) num="z"
            ###
            $(i+1)=num; num=0}
        }
    print $0
    }' |
sed -e "s/I//g" -e "s/ //g" |
# ----------------------------------------
# 短い配列をPaddingする
# ----------------------------------------
awk -F "" -v seqnum="${seq_maxnum}" \
    '{for(i=1;i<=seqnum;i++) {
    if($i=="") $i="="
    seq[i]=seq[i]$i
}}
END{for(key in seq) print seq[key]}' |
# ----------------------------------------
# シークエンスエラーを描出する
# ----------------------------------------
awk -F "" '{
    sum[1]=gsub("=","=",$0)
    sum[2]=gsub("M","M",$0)
    sum[3]=gsub(/[1-9]|[a-z]/,"@",$0)
    sum[4]=gsub("D","D",$0)
    sum[5]=gsub("S","S",$0)
    # ----------------------------------------
    ### Controlにおいて変異塩基が10%を超える塩基部位をシークエンスエラーとする
    # ----------------------------------------
    per=10
    if(sum[3]+sum[4]+sum[5] > NF*per/100) num = 2
    else num=1
    #
    print NR, "@", sum[1], sum[2], sum[3], sum[4], sum[5], "@", \
        (sum[1]+sum[2])/NF, sum[3]/NF,sum[4]/NF,sum[5]/NF, num
}' \
> "${output_ref_score}"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output Genomic coodinates (Query)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cat "${output_query_seq}" |
awk -F "" -v seqnum="${seq_maxnum}" \
    '{for(i=1;i<=seqnum;i++) {
        if($i=="") $i="="
        seq[i]=seq[i]$i
    }}
END{for(key in seq) print seq[key]}' |
awk -F "" 'BEGIN{OFS=","}{
    totalI=gsub(/[1-9]|[a-z]/,"@",$0)
    totalD=gsub("D","D",$0)
    totalS=gsub("S","S",$0)
    for(i=1; i<=NF; i++){
        if($i=="=") $i=0
        else if($i=="M") $i=0
        else if($i=="@") $i=totalI
        else if($i=="D") $i=totalD*(-1)
        else if($i=="S") $i=totalS
    }
    print $0
    }' |
# ----------------------------------------
# シークエンスエラーはMatchとしてあつかつ
# ----------------------------------------
paste - "${output_ref_score}" |
awk '{if($NF==2) $1=0
    print $1}' \
> "${output_query_score}"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#* Clustering by HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Rscript DAJIN/src/clustering.R \
"${output_query_score}" "${output_label}" 2>/dev/null
[ "$?" -eq 1 ] && error_exit 1 'Clustering error...'

# ============================================================================
# Remove minor allele (< 10%)
# ============================================================================

hdbscan_id_NR=$(cat "${hdbscan_id}" | wc -l)

cat "${hdbscan_id}" | awk '{print $NF}' | sort | uniq -c |
awk -v nr="${hdbscan_id_NR}" \
'{if($1/nr>0.1) print $2,int($1/nr*100+0.5)}' \
> .DAJIN_temp/clustering/temp/tmp_"${suffix}"

per=$(awk '{sum+=$2} END{print sum}' .DAJIN_temp/clustering/temp/tmp_"${suffix}" )

cat .DAJIN_temp/clustering/temp/tmp_"${suffix}" |
awk -v per="${per}" '{print $1, NR, int($2*100/per+0.5)}' \
> "${output_alleleper}"


# true > .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
# cat "${output_alleleper}" |
# while read -r input; do
#     before=$(echo $input | cut -d " " -f 1)
#     after=$(echo $input | cut -d " " -f 2)
#     #
#     cat "${hdbscan_id}" |
#     awk -v bf="${before}" -v af="${after}" \
#     '$2==bf {$2=af; print}' \
#     >> .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
# done
# cat .DAJIN_temp/clustering/temp/tmp_id_"${suffix}" |
# sed "s/ /\t/g" |
# sort \
# > ${output_id}

# ============================================================================
# Generate BAM files on each cluster
# ============================================================================

output_bamdir="DAJIN_Report/bam_clustering"
mkdir -p "${output_bamdir}"
set +e
rm ${output_bamdir}/${barcode}_${alleletype_original}* 2>/dev/null
set -e

for i in $(cat "${output_alleleper}" | cut -d " " -f 2 | sort -u); do
    index=$(cat "${output_alleleper}" | sed -n "${i}"p | cut -d " " -f 1)
    #
    cat "${hdbscan_id}" | grep "${index}$" | cut -f 1 | sort \
    > ".DAJIN_temp/clustering/temp/tmp_id_${suffix}"
    #
    samtools view -h DAJIN_Report/bam/"${barcode}".bam |
    grep "^@" > ".DAJIN_temp/clustering/temp/tmp_header_${suffix}" 
    #
    samtools view DAJIN_Report/bam/"${barcode}".bam |
    sort |
    join - ".DAJIN_temp/clustering/temp/tmp_id_${suffix}" 2>/dev/null |
    sed "s/ /\t/g" 2>/dev/null |
    head -n 100 \
    >> ".DAJIN_temp/clustering/temp/tmp_header_${suffix}"
    #
    samtools sort ".DAJIN_temp/clustering/temp/tmp_header_${suffix}" \
    > "${output_bamdir}/${barcode}_${alleletype_original}_${i}.bam"
    samtools index "${output_bamdir}/${barcode}_${alleletype_original}_${i}.bam"
    #
done

# ============================================================================
# Plot mutation loci
# ============================================================================

minimap2 -ax map-ont \
    .DAJIN_temp/fasta_conv/target.fa \
    .DAJIN_temp/fasta_conv/wt.fa --cs 2>/dev/null |
grep -v "^@" |
awk '{print $(NF-1)}' |
sed -e "s/cs:Z:://g" | 
sed -e "s/:/ /g" |
sed -e "s/\([-|+|*]\)/ \1 /g" |
awk '{$NF=""
    for(i=1; i<NF; i++){if($i~/[a|t|g|c]/) $i=num+length($i)}
    print $0}' |
awk '{num=0
    for(i=1; i<=NF; i++){ if($i!~/[-|+|*]/) {num=num+$i; $i=num} }
    print $0}' |
sed -e "s/[-|+|*|=]/,/g" |
sed -e "s/ , /,/g" -e "s/ /,/g" \
> "${plot_mutsites}"

# =================================================
#* Subtract Control from Query
# =================================================

# ------------------------------------------
# annotate Deletion(D), Knock-in(I), or Point mutation(P)
# ------------------------------------------
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

cut_start=$(cut -d "," -f 1 "${plot_mutsites}")
del_size=$(awk -F "," '{print $2-$1}' "${plot_mutsites}")

true > "${output_plot}"
for cluster in $(cat "${output_alleleper}" | cut -d " " -f 2 | sort -u); do
    index=$(cat "${output_alleleper}" | sed -n "${cluster}"p | cut -d " " -f 1)
    #
    paste "${output_query_seq}" "${hdbscan_id}" |
    awk -v cl="${index}" '$NF==cl' |
    cut -f 1 |
    awk -F "" -v seqnum="${seq_maxnum}" \
        '{for(i=1;i<=seqnum;i++) {
        if($i=="") $i="="
        seq[i]=seq[i]$i
        }}
        END{for(key in seq) print seq[key]}' |
    #
    # head -n 1979 tmp_test | #!------------------------------------ 
    awk -F "" '{sequence=$0
        sum[1]=gsub("=","=",sequence)
        sum[2]=gsub("M","M",sequence)
        sum[3]=gsub(/[1-9]|[a-z]/,"@", sequence)
        sum[4]=gsub("D","D",sequence)
        sum[5]=gsub("S","S",sequence)
        # ------------------------------------------
        # 各塩基部位において最多の変異(MIDS=)をレポートする
        # ------------------------------------------
        max=sum[1]; num=1
        for(i=2; i<=5;i++){if(max<sum[i]){max=sum[i]; num=i}}
        # print max, num
        # ------------------------------------------
        # Insertion数をレポートする
        # ------------------------------------------
        max=0; ins_num=0
        if(num==3) {
            for(i=1; i<=NF; i++){array[$i]++}
            for(key in array){if(max<array[key]) {max=array[key]; ins_num=key}} 
        }
        
        print num, NR, ins_num, "@", (sum[1]+sum[2])/NF,sum[3]/NF,sum[4]/NF,sum[5]/NF
        }' |
    #
    paste - "${output_ref_score}" |
    # head -n 1979 | #! -------------------------------------------
    # ------------------------------------------
    # 各塩基部位にたいして「Mの頻度、Iの頻度、Dの頻度、Sの頻度、Iの個数」を表示する
    # ------------------------------------------
    awk 'function abs(v) {return v < 0 ? -v : v}
        $NF==1 {
            I=abs($6-$(NF-3))
            D=abs($7-$(NF-2))
            S=abs($8-$(NF-1))
            M=abs(1-I-D-S)
            print NR, M, I, D, S, $3
        }
        # Sequence error annontated as Match
        $NF==2 {
            print NR, 1, 0, 0, 0, 0
        }' |
    # ------------------------------------------
    # 各塩基部位にたいして「最大頻度の変異と挿入塩基数」を表示する
    # ------------------------------------------
    awk '{max=0; num=0
        for(i=2; i<=5;i++){if(max<$i){max=$i; num=i}}
        print num, $1, $NF}' |
    awk -v cl="${cluster}" \
    '{if($1==1) print $2, "M", cl, $NF
    else if($1==2) print $2, "M", cl, $NF
    else if($1==3) print $2, "I", cl, $NF
    else if($1==4) print $2, "D", cl, $NF
    else if($1==5) print $2, "S", cl, $NF}' |
    # ------------------------------------------
    # 「2cut-deletionかつアレルタイプがTarget」のとき、
    # 変異箇所の行番号に変異サイズを追加して、
    # seq_maxより長い配列をトリミングします。
    # ------------------------------------------
    if [ "${mutation_type}" = "D" ] && [ "${alleletype_original}" = "target" ] ; then    
        cat - |
        awk -v cut="${cut_start}" -v del="${del_size}" \
        '{if($1>cut) $1=$1+del
        print}' |
        awk -v seqnum="${seq_maxnum}" '$1 <= seqnum'
    else
        cat -
    fi \
    >> "${output_plot}"
done


# =================================================
#* Plot mutation loci
# =================================================

# printf "Plot mutation loci... \n"
mkdir -p DAJIN_Report/alleletypes
for cluster in $(cat "${output_alleleper}" | cut -d " " -f 2 | sort -u)
do
    cat "${output_plot}" |
    awk -v cl="${cluster}" '$3==cl' \
    > .DAJIN_temp/clustering/temp/tmp_"${suffix}"_"${cluster}"
    #
    Rscript DAJIN/src/clustering_alleleplot.R \
    .DAJIN_temp/clustering/temp/tmp_"${suffix}"_"${cluster}" "${plot_mutsites}" 2>/dev/null
done

# ============================================================================
# Report allele mutation info
# ============================================================================

cat "${output_plot}" |
awk '{
    num=1
    cl_mut[$3]=cl_mut[$3]$2
    if($2!="M"){
        loc[NR] = $1
        mut[NR] = $2
        cl[NR] = $3
        ins[NR] = $4
    }}
END{
    # -------------------------------
    # もし変異がなければintactと表示する
    # -------------------------------
    for(j in cl_mut) {
        if (cl_mut[j] !~/[I|D|S]/) {print j, 0,"intact"}
    }
    # -------------------------------
    # 同じ変異が5つ飛ばし以内で続いている場合は連続した変異とみなす
    # また、挿入塩基の場合は挿入塩基数を直接表記する
    # -------------------------------
    for(i in loc){
        if(loc[i+1] - loc[i] == 1 || \
            loc[i+2] - loc[i] == 2 || \
            loc[i+3] - loc[i] == 3 || \
            loc[i+4] - loc[i] == 4 || \
            loc[i+5] - loc[i] == 5) {num++}
        # if( for(j=1; j<=5; j++){loc[j+1]-log[j] == 1}) num++
        else if(ins[i]>0) {print cl[i], i, ins[i]""mut[i], loc[i]}
        else {print cl[i], i, num""mut[i], loc[i]-num; num=1}
}}' |
# -------------------------------
# 挿入塩基数が10以上の場合に数値情報に逆変換する
# -------------------------------
awk '{
    if($3 ~ /[a-z]I/) {
        for (i=10; i<=36; i++) {
            num=i+87
            ins=sprintf("%c", num)
            if($3==ins"I") $3=i"I"
        }
    }
    print $0}' |
sed "s/35I/>35I/g" |
sort -t " " -k 1,1 -k 2,2n \
> ".DAJIN_temp/clustering/temp/tmp_${suffix}"

# -------------------------------
# barcode, alleletype, クラスター番号, アレル頻度、変異、変異部位、ソート番号を出力する
# -------------------------------
cat "$output_alleleper" | cut -d " " -f 2- |
sort |
join - ".DAJIN_temp/clustering/temp/tmp_${suffix}" |
sort -t " " -k 1,1n  -k 3,3n |
awk '{print $1,$2,$4,$5, NR}' |
sed "s/^/${barcode} ${alleletype_original} /g" \
> "${output_result}"
cat "${output_result}"
cat $output_plot | grep -v M | awk '$1 < 2800'

echo "${suffix}" |
sed "s/_/ /g" |
cut -d " " -f 1,2 |
sed "s/$/ is finished.../g"

set +eu
# rm .DAJIN_temp/clustering/temp/tmp_*
# exit 0
