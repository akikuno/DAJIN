#!/bin/sh

# ============================================================================
# Initialize shell environment
# ============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ============================================================================
# Parse auguments
# ============================================================================

# barcode="barcode02"
# control="barcode30"
# alleletype="wt"
# alleletype_original=${alleletype}
# pid=$$
# suffix="${barcode}"_"${alleletype}"_"${pid}"
# [ "$alleletype" = "abnormal" ] && alleletype="wt"

barcode="${1}"
control="${2}"
alleletype="${3}"
alleletype_original=${alleletype}
pid=$$
suffix="${barcode}"_"${alleletype}"_"${pid}"
[ "$alleletype" = "abnormal" ] && alleletype="wt"

temp_dir=".DAJIN_temp/clustering/"
mkdir -p "${temp_dir}"

# ============================================================================
# MIDS conversion
# ============================================================================
MIDS_que=".DAJIN_temp/clustering/1_MIDS_${suffix}"
MIDS_ref=".DAJIN_temp/clustering/1_MIDS_${control}_${alleletype}"
# ----------------------------------------

find .DAJIN_temp/fasta_ont/ -type f | grep "${barcode}" |
xargs -I @ ./DAJIN/src/mids_convertion.sh @ "${alleletype}" "${pid}" &&
mv .DAJIN_temp/data/MIDS_"${barcode}"_"${pid}" "${MIDS_que}"

# If no control MIDS files, output... 
if [ ! -s "${MIDS_ref}" ]; then
    find .DAJIN_temp/fasta_ont/ -type f | grep "${control}" |
    xargs -I @ ./DAJIN/src/mids_convertion.sh @ "${alleletype}" "${pid}" "control" &&
    mv ".DAJIN_temp/data/MIDS_${control}_${pid}" "${MIDS_ref}"
fi

# ============================================================================
# Mutation scoring of samples
# ============================================================================
output_label="${temp_dir}/2_labels_${suffix}"
output_seq="${temp_dir}/2_seq_${suffix}"
# ----------------------------------------
cat "${MIDS_que}" |
grep "${barcode}" |
sort -k 1,1 |
join -1 1 -2 2 - .DAJIN_temp/data/DAJIN_prediction_result.txt |
awk -v atype=${alleletype_original} \
'$NF==atype' |
cut -d " " -f 1,3 |
sed "s/ /\t/g" \
> "${output_label}"
#
cat "${MIDS_que}" |
grep "${barcode}" |
sort -k 1,1 |
join -1 1 -2 2 - .DAJIN_temp/data/DAJIN_prediction_result.txt |
awk -v atype=${alleletype_original} \
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
> "${output_seq}"

# Get max sequence length
seq_maxnum=$(cat "${output_seq}" |
    grep -v "^$" |
    awk -F "" 'BEGIN{max=0}
    {if(max<length($0)) max=length($0)}
    END{print max}')

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output Genomic coodinates (Control)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# コントロールとサンプルで**同じ箇所に同じ変異がある場合**、
# その塩基はMatchとして扱うようにする
output_ref="${temp_dir}/3_score_control_${suffix}"
# ----------------------------------------
cat "${MIDS_ref}" |
grep "${control}" |
sort -k 1,1 |
join -1 1 -2 2 - .DAJIN_temp/data/DAJIN_prediction_result.txt |
# join -1 1 -2 2 - ".DAJIN_temp/prediction_result.txt" |
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
#
sed -e "s/I//g" -e "s/ //g" |
awk -F "" -v seqnum="${seq_maxnum}" \
    '{for(i=1;i<=seqnum;i++) {
    if($i=="") $i="="
    seq[i]=seq[i]$i
}}
END{for(key in seq) print seq[key]}' |

awk -F "" '{
    sum[1]=gsub("=","=",$0)
    sum[2]=gsub("M","M",$0)
    sum[3]=gsub(/[1-9]|[a-z]/,"@",$0)
    sum[4]=gsub("D","D",$0)
    sum[5]=gsub("S","S",$0)
    #
    # max=sum[1]; num=1
    # for(i=2; i<5;i++){if(max<sum[i]){max=sum[i]; num=i}}
    # ##
    ### ControlにおいてIDSがtotalの10%を超える場合をシークエンスエラーありとする。
    per=10
    if(sum[3]+sum[4]+sum[5] > NF*per/100) \
        # print NR, int(sum[3]/NF*100+0.5),int(sum[4]/NF*100+0.5),int(sum[5]/NF*100+0.5)
        num = 2
    else num=1
    #
    #print NR, num, "@", sum[1], sum[2], sum[3], sum[4], sum[5], "@", int(sum[3]/NF*100+0.5),int(sum[4]/NF*100+0.5),int(sum[5]/NF*100+0.5)
    # print NR, num, "@", sum[1], sum[2], sum[3], sum[4], sum[5], "@", sum[1]/NF*100, sum[2]/NF*100, sum[3]/NF*100,sum[4]/NF*100,sum[5]/NF*100
    print num
}' \
> "${output_ref}"
#? ===========================================================
# per=0.95
# thresI=$(cat $output_ref |
# # awk '{print $(NF-4)+$(NF-3)}' |
# awk '{print $6}' | # Insertion
# # awk '{print $7}' | # Deletion
# # awk '{print $8}' | # Substitution
# # awk '{print $(NF-2)+$(NF-1)+$(NF)}' |
# awk -v per="${per}" \
#     '{percentile[NR]=$1}
#     END{asort(percentile)
#     print percentile[int(NR*per)]
# }')

# thresD=$(cat $output_ref |
# # awk '{print $(NF-4)+$(NF-3)}' |
# # awk '{print $6}' | # Insertion
# awk '{print $7}' | # Deletion
# # awk '{print $8}' | # Substitution
# # awk '{print $(NF-2)+$(NF-1)+$(NF)}' |
# awk -v per="${per}" \
#     '{percentile[NR]=$1}
#     END{asort(percentile)
#     print percentile[int(NR*per)]
# }')

# thresS=$(cat $output_ref |
# # awk '{print $(NF-4)+$(NF-3)}' |
# # awk '{print $6}' | # Insertion
# # awk '{print $7}' | # Deletion
# awk '{print $8}' | # Substitution
# # awk '{print $(NF-2)+$(NF-1)+$(NF)}' |
# awk -v per="${per}" \
#     '{percentile[NR]=$1}
#     END{asort(percentile)
#     print percentile[int(NR*per)]
# }')
# echo $thresI $thresD $thresS
# cat  $output_ref |
# awk -v I="${thresI}" -v D="${thresD}" -v S="${thresS}" \
# '{if($6>I || $7>D || $8>S) $2=2
#     print}' \
# > tmp_ref

# cat  $output_ref | head -n 410 | tail -n20 # Deletion
# cat  $output_ref | head -n 2370 | tail -n10 # Substitution=2365
# cat  $output_ref | head -n 2745 | tail -n10 # 2740 insertion #02
# cat  $output_ref | head -n 2930 | tail -n20 # Deletion 2919  (6080-3152

# cat $tmp_plot | grep -v Match #Stx2 #02 cluster2
# # 2012 genome Deletion 2
# # 2013 genome Deletion 2
# # 2740 genome Insertion 2
# # 2919 genome Deletion 2

# per=0.50
# cat $output_ref |
# awk '{print $(NF-4)+$(NF-3)}' |
# # awk '{print $6}' | # Insertion
# # awk '{print $7}' | # Deletion
# # awk '{print $8}' | # Substitution
# # awk '{print $(NF-2)+$(NF-1)+$(NF)}' |
# awk -v per="${per}" \
#     '{percentile[NR]=$1}
#     END{asort(percentile)
#     print percentile[int(NR*per)]
# }'
#? ===========================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output Genomic coodinates (Query)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
# マッチが8割以上の塩基部位は解析から除く
# Controlと同じ塩基部位に同じ変異がある場合は除く
input="${output_seq}"
output_que="${temp_dir}/4_score_${suffix}"
# ----------------------------------------

cat "${input}" |
# grep a | head -n 1 | grep a |
awk -F "" -v seqnum="${seq_maxnum}" \
    '{for(i=1;i<=seqnum;i++) {
        if($i=="") $i="="
        seq[i]=seq[i]$i
    }}
END{for(key in seq) print seq[key]}' |
awk -F "" 'BEGIN{OFS=","}{
    len=length($0)
    totalGap=gsub("=","=",$0)
    totalM=gsub("M","M",$0)
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
    #
    # perG=int(totalGap/len*100+0.5)
    # perM=int(totalM/len*100+0.5)
    # perI=int(totalI/len*100+0.5)
    # perD=int(totalD/len*100+0.5)
    # perS=int(totalS/len*100+0.5)
    # print NR" "perI" "perD" "perS
    }' |
paste - ${output_ref} |
awk '{if($NF==2) $1=0
    print $1}' \
> "${output_que}"
# join -a 2 "${output_ref}" - | 
# # awk '$1>0 && $1<500' |
# awk '{
#     if(NF==5) {print $NF}
#     else if(NF>5){
#         num_ref=0; num_que=0
#         for(i=2; i<=4;i++) {if(max_ref<$i) {max_ref=$i; num_ref++}}
#         for(i=5; i<=7;i++) {if(max_que<$i) {max_que=$i; num_que++}}
#         if(num_ref!=num_que) print $NF}
#     }' \
# > "${output_que}"
# #     wc -l



# echo "${output_ref} and ${output_que} are finished!"
# exit 0

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Clustering by HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Clustering... \n"
Rscript DAJIN/src/test_clustering.R \
"${output_que}" "${output_label}" 2>/dev/null

# echo $?
printf "Finish clustering... \n"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract high-impact nucreotides
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Output figures... \n"
input_id="${temp_dir}/hdbscan_${suffix}"
output_plot="${temp_dir}/5_plot_${suffix}"

plot_mutsites=.DAJIN_temp/clustering/tmp_mutation_"${suffix}"
tmp_plot=.DAJIN_temp/clustering/tmp_plot_"${suffix}"
# ----------------------------------------

# true > "${output_plot}"

wt_seqlen=$(cat .DAJIN_temp/fasta/wt.fa | awk '!/[>|@]/ {print length($0)}')
minimap2 -ax map-ont .DAJIN_temp/fasta/target.fa .DAJIN_temp/fasta/wt.fa --cs 2>/dev/null |
grep -v "^@" |
awk '{print $(NF-1)}' |
sed -e "s/cs:Z:://g" | 
sed -e "s/:/ /g" |
sed -e "s/\([-|+|*]\)/ \1 /g" |
awk '{for(i=1; i<NF; i++){if($i~/[a|t|g|c]/) $i=length($i)}
    $NF=""
    print $0}' |
awk -v wt="${wt_seqlen}" \
    '{for(i=1; i<NF; i++){ if($i~/[-|+|*]/) $(i+1)=$(i+1)+$(i-1) }
    print $0}' |
sed -e "s/[-|+|*|=]/,/g" \
> "${plot_mutsites}"
cat "${plot_mutsites}"
# for cluster in $(cat "${input_id}" | cut -f 2 | sort -u); do
#     cat "${tmp_mutation}" |
#     awk -v cl="${cluster}" 'BEGIN{OFS="\t"}
#         {for(i=1;i<=NF; i++) {
#         if($i~/[-|+|*]/) {
#             for(j=$(i-2);j<=$(i-1);j++) print j, "genome", "Match", cl
#             for(j=$(i-1)+1;j<$(i+1);j++) print j, "genome", "Target", cl
#             for(j=$(i+1);j<=$(i+2);j++) print j, "genome", "Match", cl }
#         }}' |
#     sort -un \
#     >> "${output_plot}"
# done
# rm "${tmp_mutation}"

cut_start=$(cut -d " " -f 1 .DAJIN_temp/data/mutation_points)
del_size=$(cat .DAJIN_temp/data/mutation_points | awk '{print $2-$1}')

true > "${output_plot}"
for cluster in $(cat "${input_id}" | cut -f 2 | sort -u); do
    paste "${output_seq}" "${input_id}" |
    awk -v cl="${cluster}" '$NF==cl' |
    cut -f 1 |
    awk -F "" -v seqnum=${seq_maxnum} \
        '{for(i=1;i<=seqnum;i++) {
        if($i=="") $i="="
        seq[i]=seq[i]$i
        }}
        END{for(key in seq) print seq[key]}' |
    awk -F "" '{
        sum[1]=gsub("=","=",$0)
        sum[2]=gsub("M","M",$0)
        sum[3]=gsub(/[1-9]|[a-z]/,"@", $0)
        sum[4]=gsub("D","D",$0)
        sum[5]=gsub("S","S",$0)
        ###
        max=sum[1]; num=1
        for(i=2; i<5;i++){if(max<sum[i]){max=sum[i]; num=i}}
        # print NR, num, sum[1],sum[2],sum[3],sum[4],sum[5]
        print num
        }' |
    paste - ${output_ref} |
    awk '{if($NF==2) $1=1
        print $1}' |
    #
    awk -v cl="${cluster}" \
    '{if($1==1) print NR,"genome","Match", cl
    else if($1==2) print NR,"genome","Match", cl
    else if($1==3) print NR,"genome","Insertion", cl
    else if($1==4) print NR,"genome","Deletion", cl
    else if($1==5) print NR,"genome","Substitusion", cl}' \
    >> "${output_plot}"
    # # もしアレルタイプが2cut-deletionならば、変異箇所の行番号に変異サイズを追加する。
    # if [ "${alleletype_original}" = "target" ]; then    
    #     cat "${tmp_plot}" |
    #     awk -v cut="${cut_start}" -v del="${del_size}" \
    #     '{if($1>cut) $1=$1+del
    #     print }' |
    #     sed "s/ /\t/g" \
    #     > .DAJIN_temp/tmp_2_$$
    # else
    #     cp "${tmp_plot}" .DAJIN_temp/tmp_2_$$
    # fi
    # #
    # # もし最終行番号がWTの配列長に満たなかった場合、WTの配列長に合わせるように行を追加する。
    # tmp_NR=$(cat .DAJIN_temp/tmp_2_$$ | awk 'END{print $1}')
    # if [ $(echo "${tmp_NR}") -lt $(echo "${wt_seqlen}") ]; then
    #     addrow=$(awk -v nr="${tmp_NR}" -v wt="${wt_seqlen}" BEGIN'{print wt-nr}')
    #     cat $output_plot |
    #     awk -v add="${addrow}" \
    #     '{array[NR]=$0}
    #     END{for(i=NR-add; i<=NR;i++) print array[i]}' |
    #     # tail -n "${tmp_addrow}" |
    #     awk -v cl="${cluster}" 'BEGIN{OFS="\t"} {$NF=cl; print}' \
    #     >> "${output_plot}"
    # else
    #     cat .DAJIN_temp/tmp_2_$$ >> "${output_plot}"
    # fi
done
# rm "${tmp_plot}" .DAJIN_temp/clustering/tmp_2_$$

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot mutation loci
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Plot mutation loci... \n"
Rscript DAJIN/src/test_2ndclustering.R \
"${output_plot}" "${plot_mutsites}" 2>/dev/null

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Separate BAM files
# input_bamdir="DAJIN_Report/bam/"
output_bamdir="DAJIN_Report/bam_clustering/"
mkdir -p "${output_bamdir}"
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Output BAM files... \n"

# rm ${barcode}*.bam* 2>/dev/null
for i in $(cat "${input_id}" | cut -f 2 | sort -u);do
    cat "${input_id}" | grep "${i}$" | cut -f 1 | sort \
    > .DAJIN_temp/tmp_id_${suffix}
    #
    samtools view -h DAJIN_Report/bam/"${barcode}".bam |
    grep "^@" > .DAJIN_temp/tmp_header_${suffix} 
    #
    samtools view DAJIN_Report/bam/"${barcode}".bam |
    sort |
    join - .DAJIN_temp/tmp_id_${suffix} |
    sed "s/ /\t/g" |
    head -n 100 \
    >> .DAJIN_temp/tmp_header_${suffix}
    #
    samtools sort .DAJIN_temp/tmp_header_${suffix} \
    > "${output_bamdir}/${barcode}_${alleletype_original}_${i}.bam"
    samtools index "${output_bamdir}/${barcode}_${alleletype_original}_${i}.bam"
    #
done

mkdir -p DAJIN_Report/allele_type
mv .DAJIN_temp/clustering/*png DAJIN_Report/allele_type/
echo "${suffix} is successfully finished!"
# rm .DAJIN_temp/*${suffix}
rm .DAJIN_temp/tmp_*${suffix}

exit 0
