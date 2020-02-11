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
# Extract barcodes which have more than 1% target reads.
# ======================================

cat .tmp_/prediction_result.txt |
sed 1d |
awk '{barcode[$1]++
    if($3 == "target") barcode_target[$1]++}
    END{for(key in barcode)
    print key, barcode_target[key]/barcode[key]*100}' |
awk '{if($2 > 1) print $1"@@@"}' |
sort -t " " \
> .tmp_/prediction_barcodelist

# ======================================
# Extract target reads
# ======================================

cat .tmp_/prediction_result.txt |
grep target |
awk '{$1=$1"@@@"; print}' |
grep -f .tmp_/prediction_barcodelist |
sed -e "s/@@@//g" -e "s/ /\t/g" |
sort -k 2,2 > .tmp_/sorted_prediction_result

# ======================================
# Detect Mutation location
# ======================================
mutation_profile="ATAACTTCGTATAATGTATGCTATACGAAGTTAT"

printf ">mut\n${mutation_profile}\n" \
> .tmp_/mutation.fa

reference=.tmp_/mutation.fa
query=.tmp_/target.fa
if grep -q '-' .tmp_/gggenome_location; then # 2>/dev/null 1>/dev/null
    ./DAJIN/src/revcomp.sh .tmp_/target.fa > .tmp_/target_rev.fa
    query=.tmp_/target_rev.fa
fi
lalign36 -m 3 ${reference} ${query} |
grep "100.0%" |
cut -d ":" -f 2 |
sed -e "s/^/@ /g" -e "s/-/ /g" -e "s/)//g" |
sort -t " " -k2,2n |
awk '{print int(($3+$2)/2)}' \
> .tmp_/lalign_mut_center
