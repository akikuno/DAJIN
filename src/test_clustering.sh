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

barcode="${1}"
control="${2}"
alleletype="${3}"
# MIDS_file="${4}"
suffix="${barcode}"_"${alleletype}"

# barcode="target_merged"
# control="barcode30"
# alleletype="target"
# MIDS_file="data_for_ml/DAJIN.txt"
# suffix="${barcode}"_"${alleletype}"

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

# cat .tmp_/DAJIN_prediction_result.txt |
# grep -e "${barcode}" -e "${control}" |
# awk -v atype="${alleletype}" '$3==atype' |
# sort -k 2,2 \
# > .tmp_/clustering_prediction_result_"${suffix}"


# ============================================================================
# MIDS conversion
# ============================================================================
find fasta_ont/ -type f | grep "${barcode}" |
# find fasta/ -type f | grep target_1 |
xargs -I @ ./DAJIN/src/mids_convertion.sh @ "$alleletype" &&
mv .tmp_/MIDS_"${barcode}" .tmp_/MIDS_"${suffix}"
#
MIDS_que=.tmp_/MIDS_"${suffix}"
#
# If no control MIDS files, output... 
if [ ! -s .tmp_/MIDS_"${control}"_"${alleletype}" ]; then
    find fasta_ont/ -type f | grep "${control}" |
    xargs -I @ ./DAJIN/src/mids_convertion.sh @ "$alleletype" "control" &&
    mv .tmp_/MIDS_"${control}" .tmp_/MIDS_"${control}"_"${alleletype}"
fi
# ./DAJIN/src/mids_convertion.sh fasta/wt.fa "$alleletype"

MIDS_ref=.tmp_/MIDS_"${control}"_"${alleletype}"
#
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

# Output \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ----------------------------------------
input=".tmp_/clustering_seq_${suffix}"
output_que=".tmp_/clustering_score_${suffix}"
# ----------------------------------------
#
seqnum=$(cat ${input} |
    grep -v "^$" |
    awk -F "" 'BEGIN{min="inf"} \
    {if(min>length($0)) min=length($0)} \
    END{print min}')
#
cat ${input} |
# grep a | head -n 1 | grep a |
awk -F "" -v seqnum=${seqnum} \
    '{for(i=1;i<=seqnum;i++) {
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
        else if($i=="S") $i=totalD
    }
    # print $0,totalGap,totalM,totalI,totalD,totalS
    print $0}' \
> "${output_que}"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Mutation scoring of Control
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
output_ref=".tmp_/clustering_score_control_${suffix}"
# ----------------------------------------
cat "${MIDS_ref}" |
grep ${control} |
sort -k 1,1 |
join -1 1 -2 2 - ".tmp_/clustering_prediction_result_${suffix}" |
# #? ===================================================
# head |
# #? ===================================================
cut -d " " -f 2 |
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
awk -F "" -v seqnum=${seqnum} \
    '{for(i=1;i<=seqnum;i++) {
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
        else if($i=="S") $i=totalD
    }
    # print $0,totalGap,totalM,totalI,totalD,totalS
    print $0}' \
> "${output_ref}"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Clustering... \n"
Rscript DAJIN/src/test_clustering.R \
${output_ref} ${output_que} ${output_label} 2>/dev/null

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Extract nucreotide
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_feat=".tmp_/clustering_features_${suffix}"
input_id=".tmp_/clustering_id_${suffix}"
#
output=".tmp_/clustering_results_${suffix}"
output_tmp=".tmp_/clustering_tmpresult_${suffix}"
# -----------------------------------------
left=$(cat .tmp_/mutation_points | cut -d " " -f 1)
right=$(cat .tmp_/mutation_points | cut -d " " -f 2)

genome=mm10
chr=$(cat .tmp_/gggenome_location | cut -f 1)
start=$(cat .tmp_/gggenome_location | cut -f 2)
#
# cluster=allele2
for cluster in $(cat ${input_id} | cut -f 2 | sort -u); do
    percent=$(paste ".tmp_/clustering_seq_${suffix}" "${input_id}" |
    grep "${cluster}" |
    wc -l ".tmp_/clustering_seq_${suffix}"  - |
    grep -v total |
    tr -d "\n" |
    awk '{printf "%.1f\n", ($3/$1*100)}')
    #
    paste ".tmp_/clustering_seq_${suffix}" "${input_id}" |
    grep "${cluster}" |
    sort -k 3,3 |
    join -1 3 -2 2 - "${input_feat}" |
    awk '{ID[$NF]=ID[$NF]","substr($2,$NF,1)}
    END{for(key in ID) print key, ID[key]}' |
    #? ---------------------------
    # head -n 1 |
    #? ---------------------------
    awk '{ID=$1; seq=$2
        totalGap=gsub("=","",$2)
        totalM=gsub("M","",$2)
        totalI=gsub(/[1-9]|[a-z]/,"@",seq)
        totalD=gsub("D","",$2)
        totalS=gsub("S","",$2)
        if ( (totalM < totalI || totalM < totalD || totalM < totalS) && (totalD < totalI && totalS < totalI) ){
            print ID, totalGap, totalM, totalI, totalD, totalS, "@", $2}
        else if(totalM < totalI || totalM < totalD || totalM < totalS){
            print ID, totalGap, totalM, totalI, totalD, totalS}
        }' |
    # Extract insertion numbers
    awk -F "@" '{max=0; maxI=0
        for(i=1; i<=10; i++){
            numI=gsub(i,i,$2)
            if(max<numI) {max=numI; maxI=i}
            }
        for (i=10; i<36; i++){
            var=sprintf("%c", i+87)
            numI=gsub(var,var,$2)
            if(max<numI) {max=numI; maxI=i}
            }
        array[$1]=maxI}
        END{for(key in array) print key,array[key]}' |
    #
    awk '{max=0; col=0
        for(i=2;i<=NF-1;i++) if(max<$i) {max=$i; col=i}
        if(col == 2) col="gap"
        if(col == 3) col="M"
        if(col == 4) col="I"
        if(col == 5) col="D"
        if(col == 6) col="S"
        print $1, col, $NF}' |
    sort -k 1,1n |
    awk -v left=${left} -v right=${right} 'BEGIN{seq=1}
        {seqnum[$1] = $1
        seqmut[$1] = $2
        if(seqnum[$1] == seqnum[$1-1] + 1 && seqmut[$1] == seqmut[$1-1]) {seq++; loc=$1; mut=$2}
        else {loc=$1; mut=$2; seq=1}
        #
        if(loc > left-25 && loc < left + 25) {position=loc; loc="Left"}
        else if(loc > right-25 && loc < right + 25) {position=loc; loc="Right"}
        else {position=loc; loc="other"}
        #
        print $NF, loc, mut, seq, position
        }' |
    #
    awk '{if(max[$1" "$2" "$3]<$4){max[$1" "$2" "$3]=$4; position[$1" "$2" "$3]=$NF}}
        END{for(key in max) print key, max[key], position[key]}' |
    # OUTPUT
    awk -v genome=${genome:-mm10} -v chr=${chr:=chr5} -v start=${start} \
    '{if($3=="I") {s=start+$NF-$1; e=start+$NF
        print $2"-"$1""$3, genome, chr":"e"-"s}
    else {s=start+$NF-$4; e=start+$NF
        print $2"-"$4""$3, genome, chr":"s"-"e}}' |
    sed "s/^/${barcode} ${cluster} ${percent} /g" |
    sort -k 1,1 -k 6,6 |
    sed "s/ /,/g" > ${output_tmp}_${cluster}
    #
    if [ -s ${output_tmp}_${cluster} ]; then
      cat ${output_tmp}_${cluster} > ${output}_${cluster}
    else
      printf "${barcode},${cluster},${percent},intact ${alleletype}\n" > ${output}_${cluster}
    fi
done

# cat ${output}*
printf "${barcode} is finished...\n"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Separate BAM files
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# cat fasta_ont/* > fasta_ont/all_merged.fa
# cat fasta_ont/wt* > fasta_ont/wt_merged.fa
# ./DAJIN/src/igvjs.sh ${genome:-mm10} ${threads:-1}

#barcode=barcode02
#cat ${input_id} | cut -f 2 | sort | uniq -c
#
# rm ${barcode}*.bam* 2>/dev/null
# for i in $(cat ${input_id} | cut -f 2 | sort -u);do
#     cat ${input_id} | grep ${i}$ | cut -f 1 | sort > tmp_id
#     samtools view -h bam/${barcode}.bam | grep "^@" > tmp_header 
#     #
#     samtools view bam/${barcode}.bam |
#     sort |
#     join - tmp_id |
#     sed "s/ /\t/g" |
#     head -n 10 \
#     >> tmp_header
#     #
#     samtools sort tmp_header > ${barcode}_${i}.bam
#     samtools index ${barcode}_${i}.bam
# done

# samtools view -h bam/${barcode}.bam | grep "^@" > tmp_header 
# #
# samtools view bam/${barcode}.bam |
# sort |
# join - tmp_id |
# sed "s/ /\t/g" |
# head -n 100 \
# >> tmp_header
# #
# samtools sort tmp_header > ${barcode}_all.bam
# samtools index ${barcode}_all.bam

# rm tmp_*
