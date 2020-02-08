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
# Pairwise alignment between KI sequence and reads
# ======================================

flank1=$(cat .tmp_/lalign_mut_center | head -n 1)
flank2=$(cat .tmp_/lalign_mut_center | tail -n 1)

# barcode=barcode20
mkdir -p results/figures/png/seqlogo/ results/figures/svg/seqlogo/
for barcode in $(cat .tmp_/prediction_barcodelist | sed "s/@@@//g"); do
    bam=bam/${barcode}.bam
    printf " ----------------------- \n ${bam} is processing... \n ----------------------- \n"
    #
    samtools view ${bam} |
    sort |
    join -1 1 - -2 2 .tmp_/sorted_prediction_result |
    cut -d " " -f 1-5,10 \
    > .tmp_/sorted_bam
    #
    cat .tmp_/sorted_bam |
    awk -v f1=${flank1} -v f2=${flank2} '{
        print ">"$1"\n"substr($6,f1-100, 200) > ".tmp_/mutsite_split1"
    print ">"$1"\n"substr($6,f2-100, 200) > ".tmp_/mutsite_split2"}'
    #
    rm -rf .tmp_/split 2>/dev/null 1>/dev/null
    mkdir -p .tmp_/split
    #
    split -l 2 .tmp_/mutsite_split1 .tmp_/split/split1_
    split -l 2 .tmp_/mutsite_split2 .tmp_/split/split2_
    #
    printf "Align reads to loxP sequence...\n"
    { find .tmp_/split |
        grep split1_ |
        xargs -I {} ./DAJIN/src/intact_lalign.sh .tmp_/mutation.fa {} \
    > .tmp_/lalign1.fa & } 1>/dev/null 2>/dev/null
    #
    { find .tmp_/split |
        grep split2_ |
        xargs -I {} ./DAJIN/src/intact_lalign.sh .tmp_/mutation.fa {} \
    > .tmp_/lalign2.fa & } 1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null # 1 min...
    #
    printf "Output sequence logo at loxP loci...\n"
    # Multiple alignment by clustal omega
    { clustalo -t DNA --auto -i .tmp_/lalign1.fa \
    > .tmp_/clustalo_left.fa & } 1>/dev/null 2>/dev/null
    #
    { clustalo -t DNA --auto -i .tmp_/lalign2.fa \
    > .tmp_/clustalo_right.fa & } 1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null
    # -----------------------------------------------------------------
    # REMOVE GAP
    # -----------------------------------------------------------------
    for input in .tmp_/clustalo_left.fa .tmp_/clustalo_right.fa; do
        output_rmgap=$(echo ${input} | sed -e "s/.fa/_rmgap.fa/g")
        # Extract gap-enriched nucreotide location
        true > .tmp_/remove_gaprow
        seqnum=$(cat ${input} | awk -F "" '{if(NR==2) print length($0)}')
        for i in $(awk -v num=${seqnum} 'BEGIN{for(i=1;i<=num;i++) print i}'); do
            # echo "$i ==============="
            cat ${input} |
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
        cat ${input} |
        awk -F "" -f .tmp_/remove_gap.awk |
        sed "s/ //g" \
        > ${output_rmgap}
        #
        # Output sequence logo
        ## PNG
        output_logo=$(echo ${output_rmgap} | sed -e "s/.*clustalo_//g" -e "s/_rmgap.fa//g")
        { weblogo --title "${barcode} ${output_logo}" --scale-width no -n 50 --errorbars no -c classic --format png_print \
        < ${output_rmgap} > results/figures/png/seqlogo/${barcode}_${output_logo}.png & } 1>/dev/null 2>/dev/null
        ## SVG
        { weblogo --title "${barcode} ${output_logo}" --scale-width no -n 50 --errorbars no -c classic --format svg \
        < ${output_rmgap} > results/figures/svg/seqlogo/${barcode}_${output_logo}.svg & } 1>/dev/null 2>/dev/null
        wait 1>/dev/null 2>/dev/null
    done
done