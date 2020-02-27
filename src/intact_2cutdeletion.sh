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

#seqlogo barcode04
#barcode=barcode04

cat .tmp_/prediction_barcodelist | sed "s/@@@//g" |
awk '{print "./DAJIN/src/test_seqlogo.sh",$0, "&"}' |
awk -v th=${threads:-1} '{
    if (NR%th==0) gsub("&","&\nwait",$0)
    print}' |
sed -e "$ a wait" |
sh -E -

# mkdir -p results/figures/png/seqlogo/ results/figures/svg/seqlogo/
# for barcode in $(cat .tmp_/prediction_barcodelist | sed "s/@@@//g"); do
# echo $barcode
# done



#
# for direction in Fw Rv; do
#     for fasta in tmp1 tmp2; do
#         #input=".tmp_/tmp_lalign.fa" &&
#         input=.tmp_/${fasta}.fa
#         output="${barcode}.png"
#         #
#         weblogo --title "${barcode}: ${direction} joint sequence" \
#         --scale-width yes -n 100 --errorbars no -c classic --format pdf \
#         < ${input} > hoge_${fasta}.pdf #results/figures/png/seqlogo/${output} &
#     done
# done #} 1>/dev/null 2>/dev/null
#wait 1>/dev/null 2>/dev/null
#

# ======================================
## Positive control
# ======================================
#
# for direction in Fw Rv; do
#     i=0
#     true > .tmp_/seqlogo_postion.fa
#     while [ $i -lt 100 ] ;do
#         cat ".tmp_/mutation.fa" \
#         >> .tmp_/seqlogo_postion.fa
#         i=$((i+1))
#     done
# done
# #
# for direction in Fw Rv; do
#     ## PNG
#     { weblogo \
#         --title "${direction} Expected Joint sequence" \
#         -n 50 \
#         --errorbars no -c classic \
#         --format png_print \
#         < .tmp_/seqlogo_postion.fa \
#     > results/figures/png/seqlogo/Expected.png & } \
#     1>/dev/null 2>/dev/null
#     ## SVG
#     { weblogo \
#         --title "${direction} Expected Joint sequence" \
#         -n 50 \
#         --errorbars no \
#         -c classic \
#         --format svg \
#         < .tmp_/seqlogo_postion.fa \
#         > results/figures/svg/seqlogo/Expected.svg & } \
#     1>/dev/null 2>/dev/null
#     wait 1>/dev/null 2>/dev/null
# done

exit 0
