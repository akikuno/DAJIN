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
# I/O naming
################################################################################
#===========================================================
#? TEST Aurguments
#===========================================================
# barcode="barcode02"
# alleletype="target"
# suffix="${barcode}"_"${alleletype}"

# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

#===========================================================
#? Aurguments
#===========================================================

barcode="${1}"
alleletype="${2}"
suffix="${barcode}"_"${alleletype}"

# mapping_alleletype="wt"
mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

#===========================================================
#? Input
#===========================================================
control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"
query_score=".DAJIN_temp/clustering/temp/query_score_${suffix}"
query_label=".DAJIN_temp/clustering/temp/query_labels_${suffix}"

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため
# hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"

################################################################################
# HDBSCAN
################################################################################

Rscript DAJIN/src/test_clustering.R "${query_score}" "${query_label}" "${control_score}"

echo "Clustering ${barcode} ${alleletype} was finished..."

exit 0