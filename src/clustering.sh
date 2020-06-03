#!/bin/sh

# ==============================================================================
# Initialize shell environment
# ==============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX


# ==============================================================================
# I/O naming
# ==============================================================================
# ----------------------------------------
# Input
# ----------------------------------------
# barcode="barcode03"
# alleletype="flox_deletion"
# original_percentage=56
# suffix="${barcode}"_"${alleletype}"

# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

barcode="${1}"
alleletype="${2}"
# original_percentage="${3}"
suffix="${barcode}"_"${alleletype}"

mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"


# ----------------------------------------
# Input files
# ----------------------------------------
# control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"

# ----------------------------------------------------------
# Output files
# ----------------------------------------------------------
mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため
# temporal -----------
MIDS_que=".DAJIN_temp/clustering/temp/MIDS_${suffix}"
query_seq=".DAJIN_temp/clustering/temp/query_seq_${suffix}"
# hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"
# tmp_allele_percentage=".DAJIN_temp/clustering/temp/allele_percentage_${suffix}".txt

# resuts -----------
query_score=".DAJIN_temp/clustering/temp/query_score_${suffix}"
query_label=".DAJIN_temp/clustering/temp/query_labels_${suffix}"
# allele_id=".DAJIN_temp/clustering/result_allele_id_${suffix}".txt
# allele_percentage=".DAJIN_temp/clustering/result_allele_percentage_${suffix}".txt

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
# Clustering
# ==============================================================================

# # ----------------------------------------------------------
# # MIDS conversion
# # ----------------------------------------------------------

find .DAJIN_temp/fasta_ont/ -type f |
    grep "${barcode}" |
    xargs -I @ ./DAJIN/src/mids_convertion.sh @ "${mapping_alleletype}"
cp ".DAJIN_temp/data/MIDS_${barcode}_${mapping_alleletype}" "${MIDS_que}"

# ----------------------------------------------------------
# Mutation scoring of samples
# ----------------------------------------------------------
# ----------------------------------------
# 配列IDとラベルをつくる
# ----------------------------------------

cat "${MIDS_que}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v atype="${alleletype}" '$NF==atype' |
    cut -d " " -f 1,3 |
    sed "s/ /,/g" |
cat - > "${query_label}"

# ----------------------------------------
# 挿入塩基を1つの挿入塩基数にまとめて配列のズレを無くす
# ----------------------------------------
cat "${MIDS_que}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v atype="${alleletype}" '$NF==atype' |
    cut -d " " -f 2 |
    awk -F "" '{
        for(i=1; i<=NF; i++){
            if($i=="I") num=num+1
            if($i=="I" && $(i+1)!="I") {
                # -----------------------------------
                # e.g) if num=10, num becomes "a"
                # -----------------------------------
                if(num>=10 && num<=35) num=sprintf("%c", num+87)
                else if(num>=36) num="z"
                #
                $(i+1)=num; num=0}
            }
        print $0}' |
    # ----------------------------------------
    # MIDS変換で末尾がDになった配列を=に変換する
    # ----------------------------------------
    sed -e "s/I//g" -e "s/ //g" |
    sed "s/\(D*$\)/ \1/g" |
    awk '{
        for(i=1; i<=NF; i++) if($i~/^D*$/) gsub(/./, "=", $i)
    }1' |
    sed "s/ //g" |
    # ----------------------------------------
    # 短い配列を"="でPaddingする
    # ----------------------------------------
    awk -v seqnum="${seq_maxnum}" \
        'BEGIN{OFS=""}
        { if(length($0) < seqnum){
            seq="="
            for(i=length($0)+1; i<=seqnum; i++) $i=seq
            print $0}
        }' |
cat - > "${query_seq}"

# ----------------------------------------------------------
# Output Genomic coodinates (Query)
# ----------------------------------------------------------
cat "${query_seq}" |
    awk -F '' 'BEGIN{OFS=","} {$1=$1;print $0}' |
    sed "s/=/M/g" |
    sed "s/[0-9]/I/g" |
    sed "s/[a-z]/I/g" |
cat > "${query_score}"

# cp "${control_score}" test_control
# echo "${query_score}" "${query_label}" "${control_score}"
# time Rscript DAJIN/src/test_clustering.R "${query_score}" "${query_label}" "${control_score}"
# ls -lh "${hdbscan_id}" #! =============================


#? --------------------------------------------------------


# cat "${query_seq}" | head -n 5 > test_seq

# cat test_seq |
#     sed "s/[M=]/0 /g" |
#     sed "s/[1-9]/1 /g" |
#     sed "s/[a-z]/1 /g" |
#     sed "s/D/-1 /g" |
#     sed "s/S/1 /g" |
#     sed "s/ $//g" |
# cat > tmp_test


# cat "${control_score}" |
#     awk -F "" \
#     '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } }
#     END {    
#         for(j=1; j<=NF; j++) {
#             str=a[1,j]
#             for(i=2; i<=NR; i++){ str=str" "a[i,j] }
#             print str }
#     }' |
# cat - tmp_test |
# awk '{for(i=1; i<=NF; i++){if(NR==1 && $i==2) array[i] = i}
#         {for(key in array) $array[key] = 0; print $0}
#     }' |
# sed 1d |
# sed "s/ /,/g" |
# cat > tmp_score

# cat tmp_score |
#     awk -F "," '{for(i=1; i<=NF; i++) sum[i]+=$i}
#         END{for(key in sum) printf sum[key]","}' |
#     sed "s/,$/\n/g" |
# cat > tmp_weight

# cp "${query_label}" tmp_label
# #? --------------------------------------------------------

# cat "${query_seq}" |
# # cat test_seq |
# awk -F '' 'BEGIN{OFS=","} {$1=$1;print $0}' |
# sed "s/=/M/g" |
# sed "s/[0-9]/I/g" |
# sed "s/[a-z]/I/g" |
# cat > test_MIDS.csv
# cp "${control_score}" test_control

#     # ----------------------------------------
#     # 行を「リード指向」から「塩基部位指向」に変換する
#     # ----------------------------------------
#     awk -F "" \
#     '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } }
#     END {    
#         for(j=1; j<=NF; j++) {
#             str=a[1,j]
#             for(i=2; i<=NR; i++){ str=str""a[i,j] }
#             print str }
#     }' |
#     # ----------------------------------------
#     # 変異部の数を数える
#     # ----------------------------------------
#     awk -F "" 'BEGIN{OFS=","}{
#         totalI=gsub(/[1-9]|[a-z]/,"@",$0)
#         totalD=gsub("D","D",$0)
#         totalS=gsub("S","S",$0)
#         for(i=1; i<=NF; i++){
#             if($i=="=") $i=0
#             else if($i=="M") $i=0
#             else if($i=="@") $i=totalI
#             else if($i=="D") $i=totalD*(-1)
#             else if($i=="S") $i=totalS
#             }
#         }1' |
#     # ----------------------------------------
#     # シークエンスエラーはMatchとしてあつかつ
#     # ----------------------------------------
#     paste - "${control_score}" |
#     awk '{if($NF==2) $1=0
#         print $1}' |
# cat - > "${query_score}"

# ls -lh "${query_score}" #!-----------------------------------
# # ----------------------------------------------------------
# # Clustering by HDBSCAN
# # ----------------------------------------------------------

# Rscript DAJIN/src/clustering.R \
#     "${query_score}" "${query_label}" \
#     2>/dev/null
# if [ "$?" -eq 1 ]; then
#     echo "Clustering error..." 1>&2
# 	exit 1
# fi
#
# ls -lh ".DAJIN_temp/clustering/temp/hdbscan_${suffix}" #!-----------------------------------

# # ==============================================================================
# # Summarize and plot mutation loci
# # ==============================================================================
# # ----------------------------------------------------------
# # Remove minor allele (< 10%) 
# # 全体の10%以下のアレルは削除する
# # ----------------------------------------------------------

# cat "${hdbscan_id}" |
#     awk '{print $NF}' |
#     sort |
#     uniq -c |
#     awk -v per="${original_percentage}" -v nr="$(cat "${hdbscan_id}" | wc -l))" \
#     '{allele_per=$1/nr*per
#     if(allele_per>10) {
#         total+=allele_per
#         allele[NR]=$2" "allele_per}}
#     END{for(key in allele) print allele[key],total, per}' |
#     awk '{print $1, NR, int($2/$3*$4+0.5)}' |
# cat - > "${tmp_allele_percentage}"
# # cat - > .DAJIN_temp/clustering/temp/tmp_"${suffix}"

# # # ----------------------------------------------------------------
# # # 取り除かれたぶんの割合を調整して、合計の割合を100％とする
# # # ----------------------------------------------------------------

# # cat .DAJIN_temp/clustering/temp/tmp_"${suffix}" |
# #     awk -v per="$(awk '{sum+=$2} END{print sum}' .DAJIN_temp/clustering/temp/tmp_"${suffix}")" \
# #     '{print $1, NR, int($2*100/per+0.5)}' |
# # cat - > "${allele_percentage}"

# # ============================================================================
# # Report allele mutation info
# # 各リードとクラスターの対応付を行う
# #（次のVCF作製とSequence logo描出のために必要）
# # ============================================================================

# before=$(cat "${tmp_allele_percentage}" | cut -d " " -f 1 | xargs echo)
# after=$(cat "${tmp_allele_percentage}" | cut -d " " -f 2 | xargs echo)

# paste "${hdbscan_id}" "${query_seq}" |
#     awk -v bf="${before}" -v af="${after}" \
#     'BEGIN{OFS="\t"
#         split(bf,bf_," ")
#         split(af,af_," ")}
#     {for(i in bf_){if($2==bf_[i]){$2=af_[i]; print}}
#     }' |
#     sed "s/ /\t/g" |
#     sort |
# cat - > "${allele_id}"

# cat "${tmp_allele_percentage}" |
#     cut -d " " -f 2- |
#     sed "s/^/${suffix} /g" |
# cat - > "${allele_percentage}"

#! >>>>>>>>>>>>>>>>>>>>>>>> clustering_variantcall.sh


# # ============================================================================
# # Generate BAM files on each cluster
# # ============================================================================
# output_bam="" #!----------------------------------------------------------
# output_bamdir=".DAJIN_temp/clustering/bam_clustering"
# mkdir -p "${output_bamdir}"

# cat "${allele_percentage}" |
# cut -d " " -f 3 |
# sort -u |
# while read -r cluster
# do
#     cat "${allele_id}" |
#         awk -v cl="${cluster}" '$2==cl' |
#         cut -f 1 |
#         sort |
#     cat - > ".DAJIN_temp/clustering/temp/tmp_id_${suffix}"
#     #
#     samtools view -h DAJIN_results/BAM/"${barcode}".bam |
#         awk '/^@/{print}
#             NR==FNR{a[$1];next}
#             $1 in a' \
#             ".DAJIN_temp/clustering/temp/tmp_id_${suffix}" - |
#     samtools sort > DAJIN_results/BAM/"${output_bam}".bam
#     samtools index "${output_bam}".bam
# done


# # ----------------------------------------------------------
# # Extract mutation sites
# # 各クラスタに含まれる変異の種類と位置
# # ----------------------------------------------------------

# minimap2 -ax map-ont \
#     .DAJIN_temp/fasta_conv/target.fa \
#     .DAJIN_temp/fasta_conv/wt.fa --cs 2>/dev/null |
#     grep -v "^@" |
#     awk '{print $(NF-1)}' |
#     sed -e "s/cs:Z:://g" | 
#     sed -e "s/:/ /g" |
#     sed -e "s/\([-|+|*]\)/ \1 /g" |
#     awk '{$NF=""
#         for(i=1; i<NF; i++){if($i~/[a|t|g|c]/) $i=num+length($i)}
#         print $0}' |
#     awk '{num=0
#         for(i=1; i<=NF; i++){ if($i!~/[-|+|*]/) {num=num+$i; $i=num} }
#         print $0}' |
#     sed -e "s/[-|+|*|=]/,/g" |
#     sed -e "s/ , /,/g" -e "s/ /,/g" |
# cat - > "${plot_mutsites}"

# # ------------------------------------------
# # annotate Deletion(D), Knock-in(I), or Point mutation(P)
# # ------------------------------------------
# mutation_type=$(
#     minimap2 -ax map-ont \
#     .DAJIN_temp/fasta/wt.fa \
#     .DAJIN_temp/fasta/target.fa \
#     --cs 2>/dev/null |
#     grep -v "^@" |
#     awk '{
#         cstag=$(NF-1)
#         if(cstag ~ "-") print "D"
#         else if(cstag ~ "+") print "I"
#         else if(cstag ~ "*") print "P"
#         }'
# )

# cut_start=$(cut -d "," -f 1 "${plot_mutsites}")
# del_size=$(awk -F "," '{print $2-$1}' "${plot_mutsites}")

# true > "${output_plot}"
# cluster=3

# cat "${allele_percentage}" |
# cut -d " " -f 2 |
# sort -u |
# while read -r cluster
# do
#     index=$(cat "${allele_percentage}" |
#         sed -n "${cluster}"p |
#         cut -d " " -f 1)
#     #
#     paste "${query_seq}" "${hdbscan_id}" |
#     awk -v cl="${index}" '$NF==cl' |
#     cut -f 1 |
#     # ----------------------------------------
#     # 行を「リード指向」から「塩基部位指向」に変換する
#     # ----------------------------------------
#     awk -F "" \
#     '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } }
#     END {    
#         for(j=1; j<=NF; j++) {
#             str=a[1,j]
#             for(i=2; i<=NR; i++){ str=str""a[i,j] }
#             print str }
#     }' |
#     # head -n 740 | #! -----------------------------
#     # ------------------------------------------
#     # 各塩基部位において最多の変異をレポートする
#     # ------------------------------------------
#     awk -F "" '{sequence=$0
#         sum[1]=gsub("=","=",sequence)
#         sum[2]=gsub("M","M",sequence)
#         sum[3]=gsub(/[1-9]|[a-z]/,"@", sequence)
#         sum[4]=gsub("D","D",sequence)
#         sum[5]=gsub("S","S",sequence)
#         max=sum[1]; num=1
#         for(i=2; i<=5;i++){if(max<sum[i]){max=sum[i]; num=i}}

#         # print max, num

#         # ------------------------------------------
#         # Insertion数をレポートする
#         # ------------------------------------------
#         max=0; ins_num=0
#         if(num==3) {
#             for(i=1; i<=NF; i++) { if($i ~ /[0-9]|[a-z]/) array[$i]++ }
#             for(key in array){if(max<array[key]) {max=array[key]; ins_num=key}} 
#         }
        
#         print num, NR, ins_num, "@", (sum[1]+sum[2])/NF,sum[3]/NF,sum[4]/NF,sum[5]/NF
#         }' |
#     #
#     paste - "${control_score}" |
#     #head -n 740 | tail -n 5 | #! -----------------------------
#     # ------------------------------------------
#     # 各塩基部位にたいして「Mの頻度、Iの頻度、Dの頻度、Sの頻度、Iの個数」を表示する
#     # ------------------------------------------
#     # head test |
#     awk 'function abs(v) {return v < 0 ? -v : v}
#         $NF==1 {
#             I=abs($6-$(NF-3))
#             D=abs($7-$(NF-2))
#             S=abs($8-$(NF-1))
#             M=abs(1-I-D-S)
#             print $2, M, I, D, S, $3
#         }
#         # Sequence error annontated as Match
#         $NF==2 {
#             print $2, 1, 0, 0, 0, 0
#         }' |
#     # ------------------------------------------
#     # 各塩基部位にたいして「最大頻度の変異と挿入塩基数」を表示する
#     # ------------------------------------------
#     awk '{max=0; num=0
#         for(i=2; i<=5;i++){if(max<$i){max=$i; num=i}}
#         print num, $1, $NF}' |
#     awk -v cl="${cluster}" \
#     '{if($1==1) print $2, "M", cl, $NF
#     else if($1==2) print $2, "M", cl, $NF
#     else if($1==3) print $2, "I", cl, $NF
#     else if($1==4) print $2, "D", cl, $NF
#     else if($1==5) print $2, "S", cl, $NF}' |
#     # ------------------------------------------
#     # 「2cut-deletionかつアレルタイプがTarget」のとき、
#     # 変異箇所の行番号に変異サイズを追加して、
#     # seq_maxより長い配列をトリミングします。
#     # ------------------------------------------
#     if [ "${mutation_type}" = "D" ] && [ "${alleletype}" = "target" ] ; then    
#         cat - |
#         awk -v cut="${cut_start}" -v del="${del_size}" \
#         '{if($1>cut) $1=$1+del
#         print}' |
#         awk -v seqnum="${seq_maxnum}" '$1 <= seqnum'
#     else
#         cat -
#     fi |
#     # ------------------------------------------
#     # 「knock-inかつアレルタイプがTarget」かつ「KI箇所がM」のとき、
#     # KI箇所の配列情報を”T”に置換します
#     # ------------------------------------------
#     if [ "${mutation_type}" = "I" ] && [ "${alleletype}" = "target" ] ; then    
#         cat - |
#         awk -v mut=$(cat "$plot_mutsites") \
#         '{split(mut, array, ",")
#         for(i=1; i<=length(array); i=i+2){
#             if($1 >= array[i] && $1 <= array[i+1] && $2 == "M") $2="T"
#             }
#         print}' |
#         awk -v seqnum="${seq_maxnum}" '$1 <= seqnum'
#     else
#         cat -
#     fi |
#     cat - >> "${output_plot}"
# done


# # ----------------------------------------------------------
# # Plot mutation loci
# # ----------------------------------------------------------

# # printf "Plot mutation loci... \n"
# mkdir -p DAJIN_Report/alleletypes

# cat "${allele_percentage}" |
# cut -d " " -f 2 |
# sort -u |
# while read -r cluster 
# do
#     cat "${output_plot}" |
#         awk -v cl="${cluster}" '$3==cl' |
#     cat - > .DAJIN_temp/clustering/temp/tmp_"${suffix}"_"${cluster}"
#     #
#     Rscript DAJIN/src/clustering_alleleplot.R \
#         .DAJIN_temp/clustering/temp/tmp_"${suffix}"_"${cluster}" \
#         "${plot_mutsites}" 2>/dev/null
# done

# # ============================================================================
# # Report allele mutation info
# # 各リードとクラスターの対応付を行う（次のSequence logo描出のために必要）
# # ============================================================================

# true > .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
# cat "${allele_percentage}" |
# while read -r input; do
#     before=$(echo "$input" | cut -d " " -f 1)
#     after=$(echo "$input" | cut -d " " -f 2)
#     #
#     cat "${hdbscan_id}" |
#     awk -v bf="${before}" -v af="${after}" \
#     '$2==bf {$2=af; print}' \
#     >> .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
# done

# cat .DAJIN_temp/clustering/temp/tmp_id_"${suffix}" |
#     sed "s/ /\t/g" |
#     sort |
# cat - > "${allele_id}"

# # ============================================================================
# # Report allele mutation info
# # ============================================================================

# cat "${output_plot}" |
#     awk '{
#         num=1
#         cl_mut[$3]=cl_mut[$3]$2
#         if($2!="M"){
#             loc[NR] = $1
#             mut[NR] = $2
#             cl[NR] = $3
#             ins[NR] = $4
#         }}
#     END{
#         # ----------------------------------------------------------
#         # もし変異がなければintactと表示する
#         # ----------------------------------------------------------
#         for(j in cl_mut) {
#             if (cl_mut[j] !~/[I|D|S]/) {print j, 0,"intact"}
#         }
#         # ----------------------------------------------------------
#         # 同じ変異が5つ飛ばし以内で続いている場合は連続した変異とみなす
#         # また、挿入塩基の場合は挿入塩基数を直接表記する
#         # ----------------------------------------------------------
#         for(i in loc){
#             if(loc[i+1] - loc[i] == 1 || \
#                 loc[i+2] - loc[i] == 2 || \
#                 loc[i+3] - loc[i] == 3 || \
#                 loc[i+4] - loc[i] == 4 || \
#                 loc[i+5] - loc[i] == 5) {num++}
#             # if( for(j=1; j<=5; j++){loc[j+1]-log[j] == 1}) num++
#             else if(ins[i]>0) {print cl[i], i, ins[i]""mut[i], loc[i]}
#             else {print cl[i], i, num""mut[i], loc[i]-num; num=1}
#     }}' |
#     # ----------------------------------------------------------
#     # 挿入塩基数が10以上の場合に数値情報に逆変換する
#     # ----------------------------------------------------------
#     awk '{
#         if($3 ~ /[a-z]I/) {
#             for (i=10; i<=36; i++) {
#                 num=i+87
#                 ins=sprintf("%c", num)
#                 if($3==ins"I") $3=i"I"
#             }
#         }
#         print $0}' |
#     sed "s/35I/>35I/g" |
#     sort -t " " -k 1,1 -k 2,2n |
# cat - > ".DAJIN_temp/clustering/temp/tmp_${suffix}"

# # ----------------------------------------------------------
# # barcode, alleletype, クラスター番号, アレル頻度、変異、変異部位、ソート番号を出力する
# # ----------------------------------------------------------
# cat "$allele_percentage" |
#     cut -d " " -f 2- |
#     sort |
#     join - ".DAJIN_temp/clustering/temp/tmp_${suffix}" |
#     sort -t " " -k 1,1n  -k 3,3n |
#     awk '{print $1,$2,$4,$5, NR}' |
#     sed "s/^/${barcode} ${alleletype} /g" |
# cat - > "${output_result}"


# # ============================================================================
# # Finish call
# # ============================================================================

# echo "${suffix}" |
# sed "s/_/ /g" |
# cut -d " " -f 1,2 |
# sed "s/$/ is finished.../g"

# rm .DAJIN_temp/clustering/temp/*${suffix}*

# set +eu
# exit 0
