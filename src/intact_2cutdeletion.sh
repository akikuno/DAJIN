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
./DAJIN/src/revcomp.sh .tmp_/target.fa > .tmp_/target_rev.fa
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
# 
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
    # grep ATGTGATGCTCGGCTTGGGAACAACGCTGT |
    # grep 9a53fc27 |
    #head -n5 |
    awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i,$10}' |
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
    find .tmp_/split/ -name split_* -type f |
        xargs -I {} ./DAJIN/src/intact_lalign.sh .tmp_/mutation.fa {} \
    > .tmp_/lalign.fa # 1>/dev/null 2>/dev/null
    # cp .tmp_/lalign_test.fa .tmp_/clustalo.fa 
    # -----------------------------------------------------------------
    # Multiple alignment by clustal omega
    # -----------------------------------------------------------------
    printf "Output sequence logo at loxP loci...\n"
    clustalo --threads=${threads:-1} -t DNA --outfmt=vie --auto -i .tmp_/lalign.fa \
    > .tmp_/clustalo.fa 2>/dev/null
    # -----------------------------------------------------------------
    # REMOVE GAP
    # -----------------------------------------------------------------
    output_rmgap=$(echo .tmp_/clustalo.fa | sed -e "s/.fa/_rmgap.fa/g")
    # Extract gap-enriched nucreotide location
    true > .tmp_/remove_gaprow
    seqnum=$(cat .tmp_/clustalo.fa | awk -F "" '{if(NR==2) print length($0)}')
    for i in $(awk -v num=${seqnum} 'BEGIN{for(i=1;i<=num;i++) print i}'); do
        #echo "$i ==============="  #! ===================================
        cat .tmp_/clustalo.fa | #! ===================================
        #cat .tmp_/lalign_test.fa | #! ===================================
        awk -F "" -v i=${i} '{if($1!~/^>/) print $i}' |
        sort |
        uniq -c |
        awk '$2!=""' |
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
    -e 's/^/{if($1 !~ "^>"){/g' \
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
    { weblogo --title "${barcode}: Joint sequence" \
    --revcomp \
    --scale-width no -n 100 --errorbars no -c classic --format png_print \
    < ${output_rmgap} > results/figures/png/seqlogo/${barcode}.png & } 1>/dev/null 2>/dev/null
    ## SVG
    { weblogo --title "${barcode}: Joint sequence" \
    --revcomp \
    --scale-width no -n 100 --errorbars no -c classic --format svg \
    < ${output_rmgap} > results/figures/svg/seqlogo/${barcode}.svg & } 1>/dev/null 2>/dev/null
    wait 1>/dev/null 2>/dev/null
done

i=0
true > .tmp_/seqlogo_postion.fa
while [ $i -lt 1000 ] ;do
    cat ".tmp_/mutation.fa" >> .tmp_/seqlogo_postion.fa
    i=$((i+1))
done

## Positive control PNG
{ weblogo --title "Expected Joint sequence" \
--revcomp \
--scale-width no -n 50 --errorbars no -c classic --format png_print \
< .tmp_/seqlogo_postion.fa > results/figures/png/seqlogo/expected.png & } 1>/dev/null 2>/dev/null
## Positive control SVG
{ weblogo --title "Expected Joint sequence" \
--revcomp \ 
--scale-width no -n 50 --errorbars no -c classic --format svg \
< .tmp_/seqlogo_postion.fa > results/figures/svg/seqlogo/expected.svg & } 1>/dev/null 2>/dev/null
wait 1>/dev/null 2>/dev/null

exit 0
