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
#./DAJIN/src/revcomp.sh .tmp_/target.fa > .tmp_/target_rev.fa
# STRECHER
## compare original and revcomp
# score1=$(stretcher -asequence .tmp_/ref.fa -bsequence .tmp_/target.fa \
# -outfile stdout 2>/dev/null | grep Score | awk '{print $NF}')
# score2=$(stretcher -asequence .tmp_/ref.fa -bsequence .tmp_/target_rev.fa \
# -outfile stdout 2>/dev/null | grep Score | awk '{print $NF}')
# query=.tmp_/target.fa
# [ "$score2" -gt "$score1" ] && query=.tmp_/target_rev.fa

# stretcher -asequence .tmp_/ref.fa -bsequence ${query} \
# -outfile stdout -aformat sam 2>/dev/null|
# grep -v "^@" |
# cut -f 6,10 |
# awk '{gsub(/[A-Z].*/,"",$1)
#     print $1 > ".tmp_/joint_site"
#     print ">mut\n"substr($2,$1-24,50) > ".tmp_/mutation.fa"}'

minimap2 -ax map-ont fasta/wt.fa fasta/target.fa 2>/dev/null |
awk '$1 !~ "@"' |
awk '{sub("M.*","",$6)
    print substr($10,$6-24,50)}' |
sed "s/^/>mut\n/g" \
> .tmp_/mutation_Fw.fa

./DAJIN/src/revcomp.sh .tmp_/mutation_Fw.fa \
> .tmp_/mutation_Rv.fa

# ======================================
# Extract Joint sequence 
# ======================================

barcode=barcode03
mkdir -p results/figures/png/seqlogo/ results/figures/svg/seqlogo/
for barcode in $(cat .tmp_/prediction_barcodelist | sed "s/@@@//g"); do
    bam=bam/${barcode}.bam
    printf " ----------------------- \n ${bam} is processing... \n ----------------------- \n"
    # -----------------------------------------------------------------
    samtools view ${bam} |
    sort |
    join -1 1 - -2 2 .tmp_/sorted_prediction_result |
    awk '$2==0 || $2==16' |
    awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i,$10}' |
    # Remove deletion
    sed "s/\([-|+|=|*]\)/ \1/g" |
    awk '{seq="";
        for(i=1; i<=NF-1;i++) if($i !~ /-/) seq=seq$i
        print seq, $NF}' |
    #
    awk '{
        match($1, "~"); seq=substr($1,RSTART-50,50)
        gsub("[-|+|=]", "", seq)
        gsub("*[a-z]", "", seq)
        seq=toupper(seq)
        match($2,seq); seq=substr($2,RSTART+length(seq)-100,200)
        print ">"NR"\n"seq}' \
    > .tmp_/mutsite_split
    #
    rm -rf .tmp_/split 2>/dev/null 1>/dev/null
    mkdir -p .tmp_/split
    split -l 2 .tmp_/mutsite_split .tmp_/split/split_
    #
    printf "Align reads to joint sequence...\n"
    #
    { for direction in Fw Rv; do
    true > .tmp_/tmp_lalign_${direction} &&
        for file in $(find .tmp_/split/ -name split_* -type f); do
            ./DAJIN/src/test_intact_lalign.sh \
            .tmp_/mutation_${direction}.fa ${file} \
            >> .tmp_/tmp_lalign_${direction}
        done &&
    echo $direction &
    done } 1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null
    #
    { for direction in Fw Rv; do
        output=".tmp_/lalign_test_${direction}.fa" &&
        #
        per=$(cat .tmp_/tmp_lalign_${direction} |
            awk '{percentile[NR]=$1}
            END{asort(percentile)
            print percentile[int(NR*0.5)]}'
        ) &&
        #
        cat .tmp_/tmp_lalign_${direction} |
        awk -v per=${per} '$1>per' |
        cut -d " " -f 2- |
        sed "s/ /\n/g" \
        #awk '{print substr($0,0,50)}' \ #?--------------
        > ${output} & # 1>/dev/null
    done } 2>/dev/null
    wait 1>/dev/null 2>/dev/null
    #
    { for direction in Fw Rv; do
        input=".tmp_/lalign_test_${direction}.fa" &&
        output="${barcode}_test_${direction}.png" &&
        #
        weblogo --title "${barcode}: ${direction} joint sequence" \
        --scale-width yes -n 100 --errorbars no -c classic --format png_print \
        < ${input} > results/figures/png/seqlogo/${output} &
    done } 1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null
done

# ======================================
## Positive control
# ======================================

for direction in Fw Rv; do
    i=0
    true > .tmp_/seqlogo_postion_${direction}.fa
    while [ $i -lt 100 ] ;do
        cat ".tmp_/mutation_${direction}.fa" \
        >> .tmp_/seqlogo_postion_${direction}.fa
        i=$((i+1))
    done
done

for direction in Fw Rv; do
    ## PNG
    { weblogo \
        --title "${barcode}: ${direction} Expected Joint sequence" \
        -n 50 \
        --errorbars no -c classic \
        --format png_print \
        < .tmp_/seqlogo_postion_${direction}.fa \
    > results/figures/png/seqlogo/${barcode}_${direction}_expected.png & } \
    1>/dev/null 2>/dev/null
    ## SVG
    { weblogo \
        --title "${barcode}: ${direction} Expected Joint sequence" \
        -n 50 \
        --errorbars no \
        -c classic \
        --format svg \
        < .tmp_/seqlogo_postion_${direction}.fa \
        > results/figures/svg/seqlogo/${barcode}_${direction}_expected.svg & } \
    1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null
done

exit 0
