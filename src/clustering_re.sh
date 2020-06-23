#!/bin/sh

################################################################################
#! Initialize shell environment
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
#? TEST Aurguments
#===========================================================
# barcode="barcode12"
# alleletype="wt"
# suffix="${barcode}"_"${alleletype}"

# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

#===========================================================
#? Auguments
#===========================================================

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
query_seq=".DAJIN_temp/clustering/temp/query_seq_${suffix}"
query_label=".DAJIN_temp/clustering/temp/query_labels_${suffix}"

#===========================================================
#? Temporal
#===========================================================
MIDS_que=".DAJIN_temp/clustering/temp/MIDS_${suffix}"

################################################################################
#! Clustering
################################################################################

#===========================================================
#? MIDS conversion
#===========================================================

./DAJIN/src/mids_clustering.sh "${barcode}" "${mapping_alleletype}" > "${MIDS_que}"

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
cat > "${query_label}"

#===========================================================
#? Mutation scoring
#===========================================================

#---------------------------------------
#* 挿入塩基を1つの挿入塩基数にまとめて配列のズレを無くす
#---------------------------------------
cat "${MIDS_que}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v atype="${alleletype}" '$NF==atype' |
    cut -d " " -f 2 |
cat > "${query_seq}"

# ----------------------------------------------------------
# Output Genomic coodinates (Query)
# ----------------------------------------------------------
cat "${query_seq}" |
    awk -F '' 'BEGIN{OFS=","} {$1=$1;print $0}' |
    sed "s/=/M/g" |
    sed "s/[0-9]/I/g" |
    sed "s/[a-z]/I/g" |
cat > "${query_score}"

rm "${MIDS_que}"

exit 0