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
# alleletype="wt"
# barcode="barcode04"
# alleletype="target"
# barcode="barcode26"
# alleletype="abnormal"
#
# control="barcode30"
# alleletype_original=${alleletype}
# pid=$$
# suffix="${barcode}"_"${alleletype}"_"${pid}"
# [ "$alleletype" = "abnormal" ] && alleletype="wt"
# echo $suffix

barcode="${1}"
control="${2}"
alleletype="${3}"; alleletype_original=${3}
[ "$alleletype_original" = "target" ] && pid="@"
[ "$alleletype_original" = "wt" ] && pid="@@"
[ "$alleletype_original" = "abnormal" ] && pid="@@@"
suffix="${barcode}"_"${alleletype}"_"${pid}"
[ "$alleletype" = "abnormal" ] && alleletype="wt"


temp_dir=".DAJIN_temp/clustering/"
mkdir -p "${temp_dir}"
# ============================================================================
# MIDS conversion
# ============================================================================
MIDS_que=".DAJIN_temp/clustering/tmp_MIDS_${suffix}"
MIDS_ref=".DAJIN_temp/clustering/tmp_MIDS_${control}_${alleletype}"
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
output_label="${temp_dir}/query_labels_${suffix}"
output_query_seq="${temp_dir}/query_seq_${suffix}"
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
# cat $MIDS_que | grep -e 7e0b197a -e a9fb4ea1 > tmp 
# cat $MIDS_que | grep -e 2c6bb00c -e 96c238ac > tmp 

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
> "${output_query_seq}"

# Get max sequence length
seq_maxnum=$(cat "${output_query_seq}" |
    grep -v "^$" |
    awk -F "" 'BEGIN{max=0}
    {if(max<length($0)) max=length($0)}
    END{print max}')

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output Genomic coodinates (Control)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# コントロールとサンプルで**同じ箇所に同じ変異がある場合**、
# その塩基はMatchとして扱うようにする
output_ref_seq="${temp_dir}/control_seq_${suffix}"
output_ref_score="${temp_dir}/control_score_${suffix}"
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
sed -e "s/I//g" -e "s/ //g" \
> "${output_ref_seq}"

cat "${output_ref_seq}" |
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
    # print NR, "@", sum[1], sum[2], sum[3], sum[4], sum[5], "@", int(sum[3]/NF*100+0.5),int(sum[4]/NF*100+0.5),int(sum[5]/NF*100+0.5), num
    print NR, "@", sum[1], sum[2], sum[3], sum[4], sum[5], "@", sum[1]/NF*100, sum[2]/NF*100, sum[3]/NF*100,sum[4]/NF*100,sum[5]/NF*100, num
    # print num
}' \
> "${output_ref_score}"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output Genomic coodinates (Query)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ----------------------------------------
input="${output_query_seq}"
output_query_score="${temp_dir}/query_score_${suffix}"
# ----------------------------------------

cat "${input}" |
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
    }' |
paste - ${output_ref_score} |
awk '{if($NF==2) $1=0
    print $1}' \
> "${output_query_score}"

# echo "${output_ref_score} and ${output_query_score} are finished!"
# exit 0

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#* Clustering by HDBSCAN
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Clustering... \n"
Rscript DAJIN/src/test_clustering.R \
"${output_query_score}" "${output_label}" 2>/dev/null

# echo $?
printf "Finish clustering... \n"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#* Remove minor allele (< 10%)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
input_id="${temp_dir}/hdbscan_${suffix}"
output_alleleper="${temp_dir}/allele_percentage_${suffix}"
# ----------------------------------------

input_id_NR=$(cat "${input_id}" | wc -l)

cat "${input_id}" | awk '{print $NF}' | sort | uniq -c |
awk -v nr="${input_id_NR}" \
'{if($1/nr>0.1) print $2,int($1/nr*100+0.5)}' \
> .DAJIN_temp/tmp_"${suffix}"

per=$(cat .DAJIN_temp/tmp_"${suffix}" | awk '{sum+=$2} END{print sum}')

cat .DAJIN_temp/tmp_"${suffix}" |
awk -v per="${per}" '{print $1, int($2*100/per+0.5)}' \
> "${output_alleleper}"


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#* # Separate BAM files
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# input_bamdir="DAJIN_Report/bam/"
output_bamdir="DAJIN_Report/bam_clustering/"
mkdir -p "${output_bamdir}"
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Output BAM files... \n"

rm "${output_bamdir}"/${barcode}_${alleletype_original}*.bam* 2>/dev/null
for i in $(cat "${output_alleleper}" | cut -d " " -f 1 | sort -u);do
    cat "${input_id}" | grep "${i}$" | cut -f 1 | sort \
    > .DAJIN_temp/tmp_id_${suffix}
    #
    samtools view -h DAJIN_Report/bam/"${barcode}".bam |
    grep "^@" > .DAJIN_temp/tmp_header_${suffix} 
    #
    samtools view DAJIN_Report/bam/"${barcode}".bam |
    sort |
    join - .DAJIN_temp/tmp_id_${suffix} 2>/dev/null |
    sed "s/ /\t/g" 2>/dev/null |
    head -n 100 \
    >> .DAJIN_temp/tmp_header_${suffix}
    #
    samtools sort .DAJIN_temp/tmp_header_${suffix} \
    > "${output_bamdir}/${barcode}_${alleletype_original}_${i}.bam"
    samtools index "${output_bamdir}/${barcode}_${alleletype_original}_${i}.bam"
    #
done

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#* Plot mutation loci
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
printf "Output figures... \n"
input_id="${temp_dir}/hdbscan_${suffix}"
output_plot="${temp_dir}/plot_${suffix}"

plot_mutsites=.DAJIN_temp/clustering/tmp_mutation_"${suffix}"
# ----------------------------------------

# --------------------------------------------------------------------------------
#* Subtract Control from Query
# --------------------------------------------------------------------------------
# Control
cat "${output_ref_seq}" |
cut -f 1 |
awk -F "" -v seqnum=${seq_maxnum} \
    '{for(i=1;i<=seqnum;i++) {
    if($i=="") $i="="
    seq[i]=seq[i]$i
    }}
    END{for(key in seq) print seq[key]}' |
awk -F "" '{sequence=$0
    sum[1]=gsub("=","=",sequence)
    sum[2]=gsub("M","M",sequence)
    sum[3]=gsub(/[1-9]|[a-z]/,"@", sequence)
    sum[4]=gsub("D","D",sequence)
    sum[5]=gsub("S","S",sequence)
    print (sum[1]+sum[2])/NF, sum[3]/NF, sum[4]/NF, sum[5]/NF
    }' \
> .DAJIN_temp/clustering/tmp_control_"${suffix}"


# annotate Deletion(D), Knock-in(I), or Point mutation(P)
mutation_type=$(
    minimap2 -ax map-ont .DAJIN_temp/fasta/target.fa .DAJIN_temp/fasta/wt.fa 2>/dev/null |
    grep -v "^@" |
    cut -f 6 |
    awk '{if($0~"I") print "D"
        else if($0~"D") print "I"
        else if($0~"S") print "P"
        }'
)

cut_start=$(cut -d " " -f 1 .DAJIN_temp/data/mutation_points)
del_size=$(cat .DAJIN_temp/data/mutation_points | awk '{print $2-$1}')

true > "${output_plot}"
for cluster in $(cat "${output_alleleper}" | cut -d " " -f 1 | sort -u); do
    paste "${output_query_seq}" "${input_id}" |
    awk -v cl="${cluster}" '$NF==cl' |
    cut -f 1 |
    awk -F "" -v seqnum=${seq_maxnum} \
        '{for(i=1;i<=seqnum;i++) {
        if($i=="") $i="="
        seq[i]=seq[i]$i
        }}
        END{for(key in seq) print seq[key]}' |
    awk -F "" '{sequence=$0
        sum[1]=gsub("=","=",sequence)
        sum[2]=gsub("M","M",sequence)
        sum[3]=gsub(/[1-9]|[a-z]/,"@", sequence)
        sum[4]=gsub("D","D",sequence)
        sum[5]=gsub("S","S",sequence)
        ### 各塩基部位において最多の変異(MIDS=)をレポートする
        max=sum[1]; num=1
        for(i=2; i<5;i++){if(max<sum[i]){max=sum[i]; num=i}}
        ### もしInsertionだった場合にはInsertion数をレポートする
        max=0; ins_num=0
        if(num==3) {
            for(i=1; i<=NF; i++){array[$i]++}
            for(key in array){if(max<array[key]) {max=array[key]; ins_num=key}} 
        }
        
        print num, NR, ins_num, "@", (sum[1]+sum[2])/NF,sum[3]/NF,sum[4]/NF,sum[5]/NF
        # print num, ins_num
        }' |
    #
    paste - .DAJIN_temp/clustering/tmp_control_"${suffix}" |
    awk 'function abs(v) {return v < 0 ? -v : v}
        {# num=""; for(i=4;i<8;i++) num=num" "abs($i-$(i+4))
        I=abs($6-$10)
        D=abs($7-$11)
        S=abs($8-$12)
        M=abs(1-I-D-S)
        print NR, M, I, D, S}' |
    # もしアレルタイプが2cut-deletionならば、変異箇所の行番号に変異サイズを追加する。
    if [ "${mutation_type}" = "D" ] && [ "${alleletype_original}" = "target" ] ; then    
        cat - |
        awk -v cut="${cut_start}" -v del="${del_size}" \
        '{if($1>cut) $1=$1+del
        print }' |
        sed "s/ /\t/g"
    else cat -
    fi |
    sed "s/$/\t"${cluster}"/g" \
    >> "${output_plot}"
done

minimap2 -ax map-ont .DAJIN_temp/fasta/target.fa .DAJIN_temp/fasta/wt.fa --cs 2>/dev/null |
grep -v "^@" |
awk '{print $(NF-1)}' |
sed -e "s/cs:Z:://g" | 
sed -e "s/:/ /g" |
sed -e "s/\([-|+|*]\)/ \1 /g" |
awk '{for(i=1; i<NF; i++){if($i~/[a|t|g|c]/) $i=length($i)}
    $NF=""
    print $0}' |
awk '{for(i=1; i<NF; i++){ if($i~/[-|+|*]/) $(i+1)=$(i+1)+$(i-1) }
    print $0}' |
sed -e "s/[-|+|*|=]/,/g" \
> "${plot_mutsites}"
cat "${plot_mutsites}"

# --------------------------------------------------------------------------------
#* Plot
# --------------------------------------------------------------------------------

printf "Plot mutation loci... \n"

for cluster in $(cat "${output_alleleper}" | cut -d " " -f 1 | sort -u); do
    cat "${output_plot}" |
    grep "${cluster}$" > tmp_"${suffix}"_"${cluster}"
    Rscript DAJIN/src/clustering_alleleplot.R \
    tmp_"${suffix}"_"${cluster}" "${plot_mutsites}" 2>/dev/null
done

# mkdir -p DAJIN_Report/allele_type
# mv .DAJIN_temp/clustering/*png DAJIN_Report/allele_type/
# echo "${suffix} is successfully finished!"
# # rm .DAJIN_temp/*${suffix}
# rm .DAJIN_temp/tmp_*${suffix}


# exit 0

cat $output_alleleper
cat $output_plot | grep -v Match

# rm .DAJIN_temp/tmp_"${suffix}"
# rm .DAJIN_temp/tmp_* .DAJIN_temp/clustering/tmp_*


    # head -n 1575 | tail # ins=1570
    # head -n 1382 | tail # del
    # cat - > tmp_hoge
    # paste - ${output_ref_score} |
    # awk '{if($NF==2) $1=1
    #     print $1,$2,$3,$4,$5,$6}' > tmp2 #! ===================================================
    # cat tmp2 > tmp_${cluster}
    # awk -v cl="${cluster}" \
    # '{if($1==1) print NR,"Match", cl,$2
    # else if($1==2) print NR,"Match", cl,$2
    # else if($1==3) print NR,"Insertion", cl,$2
    # else if($1==4) print NR,"Deletion", cl,$2
    # else if($1==5) print NR,"Substitusion", cl,$2}' |
    ###

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