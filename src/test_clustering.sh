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

barcode="barcode02"
control="barcode30"
alleletype="wt"
pid=1
suffix="${barcode}"_"${alleletype}"_"${pid}"

barcode="${1}"
control="${2}"
alleletype="${3}"
pid="${4}"
suffix="${barcode}"_"${alleletype}"_"${pid}"
echo "$suffix"
# ============================================================================
# Sort prediction results
# ============================================================================

cat .tmp_/DAJIN_prediction_result.txt |
grep -e "${barcode}" -e "${control}" |
awk -v barcode="${barcode}" -v control="${control}" \
    -v que="${alleletype}" -v ref="wt" \
    '($1==barcode && $3==que) || ($1==control && $3==ref)' |
sort -k 2,2 \
> .tmp_/clustering_prediction_result_"${suffix}"

# ============================================================================
# MIDS conversion
# ============================================================================
[ "$alleletype" = "abnormal" ] && alleletype="wt"

find fasta_ont/ -type f | grep "${barcode}" |
xargs -I @ ./DAJIN/src/mids_convertion.sh @ "$alleletype" "${pid}" &&
mv .tmp_/MIDS_"${barcode}"_"${pid}" .tmp_/MIDS_"${suffix}"
#
MIDS_que=.tmp_/MIDS_"${suffix}"
#
# If no control MIDS files, output... 
if [ ! -s .tmp_/MIDS_"${control}"_"${alleletype}" ]; then
    find fasta_ont/ -type f | grep "${control}" |
    xargs -I @ ./DAJIN/src/mids_convertion.sh @ "$alleletype" "${pid}" "control" &&
    mv .tmp_/MIDS_"${control}" .tmp_/MIDS_"${control}"_"${alleletype}"
fi
# ./DAJIN/src/mids_convertion.sh fasta/wt.fa "$alleletype"

MIDS_ref=.tmp_/MIDS_"${control}"_"${alleletype}"
#
# echo "MIDS conversion of ${suffix} is finished!"
# exit 0
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Mutation scoring of samples
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
output_label=".tmp_/clustering_labels_${suffix}"
output_seq=".tmp_/clustering_seq_${suffix}"
# ----------------------------------------
cat "${MIDS_que}" |
grep "${barcode}" |
sort -k 1,1 |
join -1 1 -2 2 - ".tmp_/clustering_prediction_result_${suffix}" |
cut -d " " -f 1,3 |
sed "s/ /\t/g" \
> "${output_label}"
#
cat "${MIDS_que}" |
grep "${barcode}" |
sort -k 1,1 |
join -1 1 -2 2 - ".tmp_/clustering_prediction_result_${suffix}" |
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
output_ref=".tmp_/clustering_score_control_${suffix}"
# ----------------------------------------
cat "${MIDS_ref}" |
grep ${control} |
sort -k 1,1 |
join -1 1 -2 2 - ".tmp_/clustering_prediction_result_${suffix}" |
# join -1 1 -2 2 - ".tmp_/prediction_result.txt" |
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
awk -F "" -v seqnum=${seq_maxnum} \
    '{for(i=1;i<=seqnum;i++) {
    if($i=="") $i="="
    seq[i]=seq[i]$i
}}
END{for(key in seq) print seq[key]}' |
awk -F "" '{
    # sum[1]=gsub("=","=",$0)
    # sum[2]=gsub("M","M",$0)
    sum[3]=gsub(/[1-9]|[a-z]/,"@",$0)
    sum[4]=gsub("D","D",$0)
    sum[5]=gsub("S","S",$0)
    # max=sum[1]; num=1
    # for(i=2; i<5;i++){if(max<sum[i]){max=sum[i]; num=i}}
    # print num, sum[1],sum[2],sum[3],sum[4],sum[5]
    # print num
    # ##
    # ControlにおいてIDSがtotalの5%を超える場合をシークエンスエラーありとする。
    per=5
    # ##
    if(sum[3] > NF*per/100) num=2
    else if(sum[4] > NF*per/100) num=2
    else if(sum[5] > NF*per/100) num=2
    else num=1
    print num
}' \
> ${output_ref}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output Genomic coodinates (Query)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input=".tmp_/clustering_seq_${suffix}"
output_que=".tmp_/clustering_score_${suffix}"
# ----------------------------------------

cat ${input} |
# grep a | head -n 1 | grep a |
awk -F "" -v seqnum=${seq_maxnum} \
    '{for(i=1;i<=seqnum;i++) {
        if($i=="") $i="="
        seq[i]=seq[i]$i
    }}
END{for(key in seq) print seq[key]}' |
awk -F "" 'BEGIN{OFS=","}{
    totalGap=gsub("=","=",$0)
    totalM=gsub("M","M",$0)
    totalI=gsub(/[1-9]|[a-z]/,"@",$0)
    totalD=gsub("D","D",$0)
    totalS=gsub("S","S",$0)
    for(i=1; i<=NF; i++){
        if($i=="=") $i=0
        else if($i=="M") $i=0
        else if($i=="@") $i=totalI
        else if($i=="D") $i=totalD
        else if($i=="S") $i=totalS
    }
    # print $0,totalGap,totalM,totalI,totalD,totalS
    print $0}' |
paste - ${output_ref} |
awk '{if($NF==2) $1=0
    print $1}' \
> ${output_que}

# echo "${output_ref} and ${output_que} are finished!"
# exit 0

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Clustering by HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Clustering... \n"
Rscript DAJIN/src/test_clustering.R \
${output_que} ${output_label} 2>/dev/null

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract high-impact nucreotides
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Output figures... \n"
input_id=".tmp_/clustering_id_${suffix}"
output_plot=".tmp_/clustering_plot_${suffix}"
# ----------------------------------------
sum_input_id=$(cat "$input_id" | wc -l)
cat "$input_id" | cut -f 2 | sort | uniq -c | awk -v sum=${sum_input_id} \
'{print int($1/sum*100+0.5)}'


if [ "$alleletype" = "target" ]; then
target_loci=$(
    minimap2 -ax map-ont fasta/target.fa fasta/wt.fa --cs 2>/dev/null |
    grep -v "^@" |
    awk '{print $(NF-1)}' |
    sed -e "s/cs:Z:://g" | 
    sed -e "s/:/ /g" |
    sed -e "s/\([-|+|*]\)/ \1 /g" |
    awk '{for(i=1; i<NF; i++){if($i~/[a|t|g|c]/) $i=length($i)}
        $NF=""
        print}'
)
fi
# awk -v target="${target_loci}" 'BEGIN{cnt=split(target, path, " ")
#   for(i=1; i<=cnt; i++) {sum+=path[i]; array[i]=sum}
#   for(i=1; i<=cnt; i=i+3) {print array[i], array[i+3]}
#   }'

true > "${output_plot}"
for cluster in $(cat "${input_id}" | cut -f 2 | sort -u); do
    paste ".tmp_/clustering_seq_${suffix}" "${input_id}" |
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
        # print num, sum[1],sum[2],sum[3],sum[4],sum[5]
        print num
        }' |
    # Sequence errorのポジションはMatchとする。
    paste - ${output_ref} |
    awk '{if($2==2) $1=1
        print $1}' |
    #
    awk -v cl="${cluster}" \
    '{if($1==1) print NR,"genome","Match", cl
    else if($1==2) print NR,"genome","Match", cl
    else if($1==3) print NR,"genome","Insertion", cl
    else if($1==4) print NR,"genome","Deletion", cl
    else if($1==5) print NR,"genome","Substitusion", cl
    else if($1==6) print NR,"genome","Target", cl
    }' \
    >> "${output_plot}"
done

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot mutation loci
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript DAJIN/src/test_2ndclustering.R "${output_plot}"  2>/dev/null

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Separate BAM files
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Output BAM files... \n"

rm ${barcode}*.bam* 2>/dev/null
for i in $(cat ${input_id} | cut -f 2 | sort -u);do
    cat ${input_id} | grep ${i}$ | cut -f 1 | sort > tmp_id
    samtools view -h bam/${barcode}.bam | grep "^@" > tmp_header 
    #
    samtools view bam/${barcode}.bam |
    sort |
    join - tmp_id |
    sed "s/ /\t/g" |
    head -n 100 \
    >> tmp_header
    #
    samtools sort tmp_header > ${barcode}_${i}.bam
    samtools index ${barcode}_${i}.bam
done

# echo "${suffix} is successfully finished!"
# rm .tmp_/*${suffix}

exit 0
