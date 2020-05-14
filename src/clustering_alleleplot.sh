#!/bin/sh

# ==============================================================================
# Initialize shell environment
# ==============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX


error_exit() {
    ${2+:} false && echo "${0##*/}: $2" 1>&2
    exit "$1"
}

# ==============================================================================
# I/O naming
# ==============================================================================
# ----------------------------------------
# Input
# ----------------------------------------
# barcode="barcode12"
# alleletype="normal"

# control="barcode21" # cables2
# control="barcode26" # prdm14
# control="barcode32" # tyr point mutation
# alleletype_original=${alleletype}
# suffix="${barcode}"_"${alleletype}"
# echo $suffix
# [ "$alleletype" = "normal" ] && alleletype="wt"
# [ "$alleletype" = "abnormal" ] && alleletype="wt"

barcode="${1}"
control="${2}"
alleletype="${3}"
alleletype_original=${3}
suffix="${barcode}"_"${alleletype}"
[ "$alleletype" = "normal" ] && alleletype="wt"
[ "$alleletype" = "abnormal" ] && alleletype="wt"

mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため

# ----------------------------------------
# Temporal Output
# ----------------------------------------

# Mutation scoring of samples
output_query_seq=".DAJIN_temp/clustering/temp/query_seq_${suffix}"

# Output Genomic coodinates (Control)
output_ref_score=".DAJIN_temp/clustering/temp/control_score_${suffix}"

# Output Plot
hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"
output_plot=".DAJIN_temp/clustering/temp/plot_${suffix}"
mutation_sites=.DAJIN_temp/clustering/temp/tmp_mutation_"${suffix}"

# allele percentage on each cluster
output_alleleper=".DAJIN_temp/clustering/temp/allele_percentage_${suffix}".txt

# ----------------------------------------------------------
# Output results
# ----------------------------------------------------------
output_id=".DAJIN_temp/clustering/result_allele_id_${suffix}".txt
output_result=".DAJIN_temp/clustering/result_allele_mutinfo_${suffix}".txt

# ----------------------------------------------------------
# Get max sequence length
# ----------------------------------------------------------
seq_maxnum=$(
    cat .DAJIN_temp/fasta/fasta.fa |
    grep -v "^>" |
    awk '{if(max<length($0)) max=length($0)}
    END{print max}'
)

# ==============================================================================
# Summarize and plot mutation loci
# ==============================================================================
# ----------------------------------------------------------
# Remove minor allele (< 10%) 
# ----------------------------------------------------------

cat "${hdbscan_id}" |
    awk '{print $NF}' |
    sort |
    uniq -c |
    awk -v nr="$(cat "${hdbscan_id}" | wc -l))" \
    '{if($1/nr>0.1) print $2,int($1/nr*100+0.5)}' |
cat - > .DAJIN_temp/clustering/temp/tmp_"${suffix}"

# ----------------------------------------------------------------
# 取り除かれたぶんの割合を調整して、合計の割合を100％とする
# ----------------------------------------------------------------
cat .DAJIN_temp/clustering/temp/tmp_"${suffix}" |
    awk -v per="$(awk '{sum+=$2} END{print sum}' .DAJIN_temp/clustering/temp/tmp_"${suffix}")" \
    '{print $1, NR, int($2*100/per+0.5)}' |
cat - > "${output_alleleper}"

# ----------------------------------------------------------
# Extract mutation sites
# ----------------------------------------------------------

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
    sed -e "s/ , /,/g" -e "s/ /,/g" |
cat - > "${mutation_sites}"

# ------------------------------------------
# annotate Deletion(D), Knock-in(I), or Point mutation(P)
# ------------------------------------------
mutation_type=$(
    minimap2 -ax map-ont \
    .DAJIN_temp/fasta/wt.fa \
    .DAJIN_temp/fasta/target.fa \
    --cs 2>/dev/null |
    grep -v "^@" |
    awk '{
        cstag=$(NF-1)
        if(cstag ~ "-") print "D"
        else if(cstag ~ "+") print "I"
        else if(cstag ~ "*") print "P"
        }'
)

cut_start=$(cut -d "," -f 1 "${mutation_sites}")
del_size=$(awk -F "," '{print $2-$1}' "${mutation_sites}")

true > "${output_plot}"
cluster=3

cat "${output_alleleper}" |
cut -d " " -f 2 |
sort -u |
while read -r cluster
do
    index=$(cat "${output_alleleper}" |
        sed -n "${cluster}"p |
        cut -d " " -f 1)
    #
    paste "${output_query_seq}" "${hdbscan_id}" |
    awk -v cl="${index}" '$NF==cl' |
    cut -f 1 |
    # ----------------------------------------
    # 行を「リード指向」から「塩基部位指向」に変換する
    # ----------------------------------------
    awk -F "" \
    '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } }
    END {    
        for(j=1; j<=NF; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){ str=str""a[i,j] }
            print str }
    }' |
    # head -n 740 | #! -----------------------------
    # ------------------------------------------
    # 各塩基部位において最多の変異をレポートする
    # ------------------------------------------
    awk -F "" '{sequence=$0
        sum[1]=gsub("=","=",sequence)
        sum[2]=gsub("M","M",sequence)
        sum[3]=gsub(/[1-9]|[a-z]/,"@", sequence)
        sum[4]=gsub("D","D",sequence)
        sum[5]=gsub("S","S",sequence)
        max=sum[1]; num=1
        for(i=2; i<=5;i++){if(max<sum[i]){max=sum[i]; num=i}}

        # print max, num

        # ------------------------------------------
        # Insertion数をレポートする
        # ------------------------------------------
        max=0; ins_num=0
        if(num==3) {
            for(i=1; i<=NF; i++) { if($i ~ /[0-9]|[a-z]/) array[$i]++ }
            for(key in array){if(max<array[key]) {max=array[key]; ins_num=key}} 
        }
        
        print num, NR, ins_num, "@", (sum[1]+sum[2])/NF,sum[3]/NF,sum[4]/NF,sum[5]/NF
        }' |
    #
    paste - "${output_ref_score}" |
    #head -n 740 | tail -n 5 | #! -----------------------------
    # ------------------------------------------
    # 各塩基部位にたいして「Mの頻度、Iの頻度、Dの頻度、Sの頻度、Iの個数」を表示する
    # ------------------------------------------
    # head test |
    awk 'function abs(v) {return v < 0 ? -v : v}
        $NF==1 {
            I=abs($6-$(NF-3))
            D=abs($7-$(NF-2))
            S=abs($8-$(NF-1))
            M=abs(1-I-D-S)
            print $2, M, I, D, S, $3
        }
        # Sequence error annontated as Match
        $NF==2 {
            print $2, 1, 0, 0, 0, 0
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
    fi |
    # ------------------------------------------
    # 「knock-inかつアレルタイプがTarget」かつ「KI箇所がM」のとき、
    # KI箇所の配列情報を”T”に置換します
    # ------------------------------------------
    if [ "${mutation_type}" = "I" ] && [ "${alleletype_original}" = "target" ] ; then    
        cat - |
        awk -v mut=$(cat "$mutation_sites") \
        '{split(mut, array, ",")
        for(i=1; i<=length(array); i=i+2){
            if($1 >= array[i] && $1 <= array[i+1] && $2 == "M") $2="T"
            }
        print}' |
        awk -v seqnum="${seq_maxnum}" '$1 <= seqnum'
    else
        cat -
    fi |
    cat - >> "${output_plot}"
done


# ----------------------------------------------------------
# Plot mutation loci
# ----------------------------------------------------------

# printf "Plot mutation loci... \n"
mkdir -p DAJIN_Report/alleletypes

cat "${output_alleleper}" |
cut -d " " -f 2 |
sort -u |
while read -r cluster 
do
    cat "${output_plot}" |
        awk -v cl="${cluster}" '$3==cl' |
    cat - > .DAJIN_temp/clustering/temp/tmp_"${suffix}"_"${cluster}"
    #
    Rscript DAJIN/src/clustering_alleleplot.R \
        .DAJIN_temp/clustering/temp/tmp_"${suffix}"_"${cluster}" \
        "${mutation_sites}" 2>/dev/null
done
