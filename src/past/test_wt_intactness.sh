#!/bin/sh

# Original:"intactness_2cutdeletion.sh"
# ======================================
# Initialize shell environment
# ======================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ======================================
# Extract barcodes which have more than 20% wt reads.
# ======================================

cat .tmp_/prediction_result.txt |
sed 1d |
awk '{barcode[$1]++
    if($3 == "wt") barcode_target[$1]++}
    END{for(key in barcode)
    print key, barcode_target[key]/barcode[key]*100}' |
awk '{if($2 > 20) print $1"@@@"}' |
sort -t " " \
> .tmp_/prediction_barcodelist

# ======================================
# Extract wt reads
# ======================================

cat .tmp_/prediction_result.txt |
awk '$3 == "wt"' |
awk '{$1=$1"@@@"; print}' |
grep -f .tmp_/prediction_barcodelist |
sed -e "s/@@@//g" -e "s/ /\t/g" |
sort -k 2,2 \
> .tmp_/sorted_prediction_result

# ======================================
# Positive controls
# ======================================

minimap2 -ax map-ont .tmp_/ref.fa .tmp_/target.fa --cs=long 2>/dev/null |
awk '$1 !~ "@" {sub("cs:Z:=","",$(NF-1)); print $(NF-1)}' |
awk -F "" '{
    for(i=1;i<=NF;i++) if($i ~ /[-|+|*|=]/) cutsites[$i] = i}
    END { for(key in cutsites) print cutsites[key], substr($0,cutsites[key]-25,51)}' |
awk '{gsub(/[-|+|*|=]/,"",$2); print toupper($0)}' \
> .tmp_/cutting_sites

cat .tmp_/cutting_sites | head -n 1 | awk '{print ">left" "\n" $2}' > .tmp_/cutting_sites_left.fa
cat .tmp_/cutting_sites | tail -n 1 | awk '{print ">right" "\n" $2}' > .tmp_/cutting_sites_right.fa

# ============================================================================
# Extract Joint sequence 
# ============================================================================

cat .tmp_/prediction_barcodelist | sed "s/@@@//g" |
awk '{print "./DAJIN/src/test_wt_seqlogo.sh",$0, "&"}' |
awk -v th=${threads:-1} '{
    if (NR%th==0) gsub("&","&\nwait",$0)
    print}' |
sed -e "$ a wait" |
sh -

exit 0
