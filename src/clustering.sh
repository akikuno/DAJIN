#!/bin/sh

################################################################################
# Initialize shell environment
################################################################################

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX


################################################################################
#! I/O naming
################################################################################

#===========================================================
#? Auguments
#===========================================================

# barcode="barcode02"
# alleletype="target"
# suffix="${barcode}"_"${alleletype}"

# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

barcode="${1}"
alleletype="${2}"
suffix="${barcode}"_"${alleletype}"

mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"


#===========================================================
#? Input
#===========================================================

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため
query_score=".DAJIN_temp/clustering/temp/query_score_${suffix}"
query_label=".DAJIN_temp/clustering/temp/query_labels_${suffix}"

#===========================================================
#? Temporal
#===========================================================
MIDS_que=".DAJIN_temp/clustering/temp/MIDS_${suffix}"
query_seq=".DAJIN_temp/clustering/temp/query_seq_${suffix}"


################################################################################
#! Clustering
################################################################################

#===========================================================
#? MIDS conversion
#===========================================================

find .DAJIN_temp/fasta/ -type f |
    grep "${barcode}" |
    xargs -I @ ./DAJIN/src/mids_convertion.sh @ "${mapping_alleletype}"
cp ".DAJIN_temp/data/MIDS_${barcode}_${mapping_alleletype}" "${MIDS_que}"

#===========================================================
#? Output Sequence ID and Lable
#===========================================================

cat "${MIDS_que}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v atype="${alleletype}" '$NF==atype' |
    cut -d " " -f 1,3 |
    sed "s/ /,/g" |
cat - > "${query_label}"

#===========================================================
#? Mutation scoring
#===========================================================

#---------------------------------------
#* Get max sequence length
#---------------------------------------
seq_length=$(
    cat .DAJIN_temp/fasta/"${mapping_alleletype}".fa |
    grep -v "^>" |
    awk '{print length($0)}'
)

#---------------------------------------
#* 挿入塩基を1つの挿入塩基数にまとめて配列のズレを無くす
#---------------------------------------
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
    awk -v seqnum="${seq_length}" \
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

exit 0