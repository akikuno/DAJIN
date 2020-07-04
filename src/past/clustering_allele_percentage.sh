#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
export PATH="$(command -p getconf PATH 2>/dev/null)${PATH+:}${PATH-}"
case $PATH in :*) PATH=${PATH#?};; esac
export UNIX_STD=2003  # to make HP-UX conform to POSIX


################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================
# barcode="barcode21"
# alleletype="inversion"
# original_percentage=0.00465658
# suffix="${barcode}"_"${alleletype}"

#===========================================================
#? Auguments
#===========================================================

barcode="${1}"
alleletype="${2}"
original_percentage="${3}"
suffix="${barcode}"_"${alleletype}"

#===========================================================
#? Input
#===========================================================
query_seq=".DAJIN_temp/clustering/temp/query_seq_${suffix}"
hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"

#===========================================================
#? Output
#===========================================================
allele_id=".DAJIN_temp/clustering/result_allele_id_${suffix}".txt
allele_percentage=".DAJIN_temp/clustering/result_allele_percentage_${suffix}".txt

#===========================================================
#? Temporal
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため
tmp_allele_percentage=".DAJIN_temp/clustering/temp/allele_percentage_${suffix}".txt

################################################################################
#! Remove minor allele (< 5%) 
################################################################################

cat "${hdbscan_id}" |
    awk '{print $NF}' |
    sort |
    uniq -c |
    awk -v per="${original_percentage}" -v nr="$(cat "${hdbscan_id}" | wc -l))" \
    '{allele_per=$1/nr*per
    if(allele_per>5) {
        total+=allele_per
        allele[NR]=$2" "allele_per}}
    END{for(key in allele) print allele[key],total, per}' |
    awk '{print $1, NR, int($2/$3*$4+0.5)}' |
cat > "${tmp_allele_percentage}"

################################################################################
#! Report allele mutation info
################################################################################

before=$(cat "${tmp_allele_percentage}" | cut -d " " -f 1 | xargs echo)
after=$(cat "${tmp_allele_percentage}" | cut -d " " -f 2 | xargs echo)

paste "${hdbscan_id}" "${query_seq}" |
    awk -v bf="${before}" -v af="${after}" \
    'BEGIN{OFS="\t"
        split(bf,bf_," ")
        split(af,af_," ")}
    {for(i in bf_){if($2==bf_[i]){$2=af_[i]; print}}
    }' |
    sed "s/ /\t/g" |
    sort |
cat > "${allele_id}"

cat "${tmp_allele_percentage}" |
    cut -d " " -f 2- |
    sed "s/^/${suffix} /g" |
cat > "${allele_percentage}"

exit 0