#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -eu
umask 0022
export LC_ALL=C
export UNIX_STD=2003  # to make HP-UX conform to POSIX


################################################################################
#! I/O naming
################################################################################

#===============================================================================
#? TEST Aurguments
#===============================================================================

# barcode="barcode43"
# alleletype="wt"
# threads=14

#===========================================================
#? Auguments
#===========================================================

barcode="${1}"
alleletype="${2}"
threads="${3}"

#===========================================================
#? Input
#===========================================================

suffix="${barcode}"_"${alleletype}"
mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"

#===========================================================
#? Output
#===========================================================

mkdir -p ".DAJIN_temp/clustering/temp/"
# hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"

#===========================================================
#? Temporal
#===========================================================

MIDS_que=".DAJIN_temp/clustering/temp/MIDS_${suffix}"
query_score=".DAJIN_temp/clustering/temp/query_score_${suffix}"
query_seq=".DAJIN_temp/clustering/temp/query_seq_${suffix}"
query_label=".DAJIN_temp/clustering/temp/query_labels_${suffix}"

################################################################################
#! MIDS conversion
################################################################################

./DAJIN/src/mids_clustering.sh "${barcode}" "${alleletype}" > "${MIDS_que}"

################################################################################
#! Query seq (compressed MIDS) and Query score (comma-sep MIDS)
################################################################################

#===========================================================
#? Output query seq for `clustering_allele_percentage.sh`
#===========================================================

cat "${MIDS_que}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v atype="${alleletype}" '$NF==atype' |
    cut -d " " -f 1,2 |
cat > "${query_seq}"

#===========================================================
#? Output query score
#===========================================================

cat "${query_seq}" |
    cut -d " " -f 2 |
    awk -F '' 'BEGIN{OFS=","}{$1=$1}1' |
    sed "s/[0-9]/I/g" |
    sed "s/[a-z]/I/g" |
cat > "${query_score}"

################################################################################
#! Query label (seqID,barcodeID)
################################################################################

cat "${MIDS_que}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v atype="${alleletype}" '$NF==atype' |
    cut -d " " -f 1,3,4 |
    sed "s/ /,/g" |
cat > "${query_label}"

################################################################################
#! Clustering
################################################################################

echo "Clustering ${barcode} ${alleletype}..." >&2
Rscript DAJIN/src/clustering.R "${query_score}" "${query_label}" "${control_score}" "${threads}" 2>/dev/null || true
ps -au | grep -e "clustering.R" -e "joblib" | awk '{print $2}'| xargs kill 2>/dev/null || true