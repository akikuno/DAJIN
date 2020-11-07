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

#===========================================================
#? TEST Auguments
#===========================================================

# barcode=barcode14
# filter=on

#===========================================================
#? Auguments
#===========================================================

barcode="${1}"
filter="${2}"

#===========================================================
#? Input
#===========================================================
# query_seq
# hdbscan_id

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/allele_per/"

# allele_id=.DAJIN_temp/clustering/allele_per/readid_cl_mids_
allele_percentage=.DAJIN_temp/clustering/allele_per/label_cl_percentage_"${barcode}"

#===========================================================
#? Temporal
#===========================================================

tmp_clusterid=.DAJIN_temp/clustering/temp/clusterid_"${barcode}"
tmp_alleleper_before=.DAJIN_temp/clustering/temp/alleleper_before_"${barcode}"
tmp_alleleper_after=.DAJIN_temp/clustering/temp/alleleper_after_"${barcode}"

################################################################################
#! Retain Target (if filter=on)
################################################################################

#===========================================================
#? Count total read numbers
#===========================================================

true > "${tmp_clusterid}"

find .DAJIN_temp/clustering/temp/hdbscan_* |
    grep "${barcode}" |
while read -r input; do
    label=$(echo $input | sed "s/.*hdbscan_//g")
    #
    cat "${input}" |
    sed "s/^/${label}\t/g" |
    cat >> "${tmp_clusterid}"
done

total_reads=$(
    cat "${tmp_clusterid}" |
    cut -f 1,3 |
    sort |
    uniq -c |
    awk '{sum+=$1} END{print sum}'
)

#===========================================================
#? Retain Target > 3%; others > 3% (if filter=on)
#===========================================================

cat "${tmp_clusterid}" |
    cut -f 1,3 |
    sort |
    uniq -c |
    sed "s/$/\t${total_reads}/g" |
    awk '{$NF=$1/$NF*100}1' |
    if [ "_${filter}" = "_on" ]; then
        awk '($2!~"target" && $NF > 3) ||
        ($2~"target" && $NF > 3)'
    else
        cat -
    fi |
    cut -d " " -f 2- |
    awk '{nr[$1]++; print $1, $2, nr[$1], $3}' |
cat > "${tmp_alleleper_before}"


#===========================================================
#? Adjust to total 100%
#===========================================================

total_percentage=$(
    cat "${tmp_alleleper_before}" |
    awk '{sum+=$4} END{print sum}'
)

cat "${tmp_alleleper_before}" |
    sed "s/$/ ${total_percentage}/g" |
    if [ "_${filter}" = "_on" ]; then
        awk '{$4=sprintf("%.1f", $4*100/$5)
            print $1,$2,$3,$4}'
    else
        awk '{printf $1" "$2" "$3" "; printf "%.3f\n", $4}'
    fi |
cat > "${tmp_alleleper_after}"

################################################################################
#! Extract reads
################################################################################

rm .DAJIN_temp/clustering/readid_cl_mids_"${barcode}"* 2>/dev/null || true

cat "${tmp_alleleper_after}" |
while read -r input; do

    id=$(echo "${input}" | cut -d " " -f 1 | xargs echo)
    before=$(echo "${input}" | cut -d " " -f 2 | xargs echo)
    after=$(echo "${input}" | cut -d " " -f 3 | xargs echo)

    hdbscan_id=.DAJIN_temp/clustering/temp/hdbscan_"${id}"
    query_seq=.DAJIN_temp/clustering/temp/query_seq_"${id}"
    allele_id=.DAJIN_temp/clustering/allele_per/readid_cl_mids_"${id}"

    join "${hdbscan_id}" "${query_seq}" |
        awk -v bf="${before}" -v af="${after}" \
        'BEGIN{OFS="\t"
            split(bf,bf_," ")
            split(af,af_," ")}
        {for(i in bf_){if($2==bf_[i]){$2=af_[i]; print}}
        }' |
        sed "s/ /\t/g" |
        sort |
    cat >> "${allele_id}"
done

################################################################################
#! Report cluster number and percentage after filtration
################################################################################

cat "${tmp_alleleper_after}" |
    awk '{print $1,$3,$4}' |
    sed "s/_/ /" |
cat > "${allele_percentage}"


################################################################################
#! remove temporal files
################################################################################

rm "${tmp_clusterid}" "${tmp_alleleper_before}" "${tmp_alleleper_after}"

exit 0