#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ======================================
# Extract barcodes which have more than 20% target reads.
# ======================================

cat .tmp_/prediction_result.txt |
sed 1d |
awk '{barcode[$1]++
    if($3 == "target") barcode_target[$1]++}
    END{for(key in barcode)
    print key, barcode_target[key]/barcode[key]*100}' |
awk '{if($2 > 20) print $1"@@@"}' |
sort -t " " \
> .tmp_/prediction_barcodelist

# ======================================
# Extract target reads
# ======================================

cat .tmp_/prediction_result.txt |
awk '$3 == "target"' |
awk '{$1=$1"@@@"; print}' |
grep -f .tmp_/prediction_barcodelist |
sed -e "s/@@@//g" -e "s/ /\t/g" |
sort -k 2,2 \
> .tmp_/sorted_prediction_result

# ======================================
# Positive controls
# ======================================

minimap2 -ax map-ont fasta/wt.fa fasta/target.fa 2>/dev/null |
awk '$1 !~ "@"' |
awk '{sub("M.*","",$6)
    print substr($10,$6-25,50)}' |
sed "s/^/>mut\n/g" \
> .tmp_/mutation.fa

# ============================================================================
# Extract Joint sequence 
# ============================================================================

cat .tmp_/prediction_barcodelist | sed "s/@@@//g" |
awk '{print "./DAJIN/src/test_seqlogo.sh",$0, "&"}' |
awk -v th=${threads:-1} '{
    if (NR%th==0) gsub("&","&\nwait",$0)
    print}' |
sed -e "$ a wait" |
sh -E -

exit 0
