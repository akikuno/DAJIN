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
joint_site=$(cat .tmp_/mutation_points | cut -d " " -f 1)
echo ${joint_site}

cat fasta/target.fa |
sed 1d |
awk -v site=${joint_site} '{print ">joint\n"substr($0, site+1-15,30)}' \
> .tmp_/mutation.fa

# ======================================
# 
# ======================================
joint_site=$(cat .tmp_/mutation_points | cut -d " " -f 1)
echo ${joint_site}

barcode=barcode24
mkdir -p results/figures/png/seqlogo/ results/figures/svg/seqlogo/
for barcode in $(cat .tmp_/prediction_barcodelist | sed "s/@@@//g"); do
    bam=bam/${barcode}.bam
    printf " ----------------------- \n ${bam} is processing... \n ----------------------- \n"
    # -----------------------------------------------------------------
    samtools view ${bam} |
    sort |
    join -1 1 - -2 2 .tmp_/sorted_prediction_result |
    cut -d " " -f 1-5,10 \
    > .tmp_/sorted_bam
    #
    cat .tmp_/sorted_bam |
    awk -v joint=${joint_site} '{
        print ">"$1"\n"substr($6,joint-100, 200)}' \
    > .tmp_/mutsite_split
    #
    rm -rf .tmp_/split 2>/dev/null 1>/dev/null
    mkdir -p .tmp_/split
    split -l 2 .tmp_/mutsite_split .tmp_/split/split_
    #
    printf "Align reads to joint sequence...\n"
    find .tmp_/split/ -name split_* -type f |
        xargs -I {} ./DAJIN/src/intact_lalign.sh .tmp_/mutation.fa {} \
    > .tmp_/lalign.fa # 1>/dev/null 2>/dev/null
    #
    # -----------------------------------------------------------------
    # Multiple alignment by clustal omega
    # -----------------------------------------------------------------
    printf "Output sequence logo at loxP loci...\n"
    clustalo --threads=${threads:-1} -t DNA --auto -i .tmp_/lalign.fa \
    > .tmp_/clustalo.fa 2>/dev/null
    # -----------------------------------------------------------------
    # REMOVE GAP
    # -----------------------------------------------------------------
    output_rmgap=$(echo .tmp_/clustalo.fa | sed -e "s/.fa/_rmgap.fa/g")
    # Extract gap-enriched nucreotide location
    true > .tmp_/remove_gaprow
    seqnum=$(cat .tmp_/clustalo.fa | awk -F "" '{if(NR==2) print length($0)}')
    for i in $(awk -v num=${seqnum} 'BEGIN{for(i=1;i<=num;i++) print i}'); do
        # echo "$i ==============="
        cat .tmp_/clustalo.fa |
        awk -F "" -v i=${i} '{if(NR%2==0) print $i}' |
        sort |
        uniq -c |
        awk -v i=${i} '{sum+=$1; if(max<$1) {max=$1; nuc=$2}}
        END{print i,nuc,max/sum*100}' |
        #Extract nucleotide position with gap "-" > 20%
        awk '$2 == "-" && $3>20' |
        cut -d " " -f 1 >> .tmp_/remove_gaprow
    done
    # Remove gap-enriched nucreotide location
    cat .tmp_/remove_gaprow |
    sed -e "s/^/\$/g" \
    -e 's/$/="";@/g' |
    tr -d "\n" |
    sed -e "s/@/ /g" \
    -e "s/^/{if(NR%2==0){/g" \
    -e "s/$/ print} else print}/g" \
    > .tmp_/remove_gap.awk
    #
    cat .tmp_/clustalo.fa |
    awk -F "" -f .tmp_/remove_gap.awk |
    sed "s/ //g" \
    > ${output_rmgap}
    #
    # Output sequence logo
    ## PNG
    { weblogo --title "${barcode}: Joint sequence" --scale-width no -n 50 --errorbars no -c classic --format png_print \
    < ${output_rmgap} > results/figures/png/seqlogo/${barcode}.png & } 1>/dev/null 2>/dev/null
    ## SVG
    { weblogo --title "${barcode}: Joint sequence" --scale-width no -n 50 --errorbars no -c classic --format svg \
    < ${output_rmgap} > results/figures/svg/seqlogo/${barcode}.svg & } 1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null
done
## Positive control PNG
{ weblogo --title "Expected Joint sequence" --scale-width no -n 50 --errorbars no -c classic --format png_print \
< .tmp_/mutation.fa > results/figures/png/seqlogo/expected.png & } 1>/dev/null 2>/dev/null
## Positive control SVG
{ weblogo --title "Expected Joint sequence" --scale-width no -n 50 --errorbars no -c classic --format svg \
< .tmp_/mutation.fa > results/figures/svg/seqlogo/expected.svg & } 1>/dev/null 2>/dev/null
wait 1>/dev/null 2>/dev/null

exit 0
