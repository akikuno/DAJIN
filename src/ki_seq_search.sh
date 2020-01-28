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
# Extract barcodes which have more than 10% target reads.
# ======================================

cat .tmp_/prediction_result.txt |
sed 1d |
awk '{barcode[$1]++
    if($3 ~ "target") barcode_target[$1]++}
    END{for(key in barcode)
print key, barcode_target[key]/barcode[key]*100}' |
awk '{if($2 > 10) print $1"@@@"}' |
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
# for: ATAACTTCGTATAATGTATGCTATACGAAGTTAT
# rev: ATAACTTCGTATAGCATACATTATACGAAGTTAT
# grep -e ATAACTTCGTATAATGTATGCTATACGAAGTTAT -e ATAACTTCGTATAGCATACATTATACGAAGTTAT | sort | uniq -c
#  mutation_profile_for=$1
mutation_profile="ATAACTTCGTATAATGTATGCTATACGAAGTTAT"

printf ">mut
${mutation_profile}
" > .tmp_/mutation.fa

reference=.tmp_/mutation.fa
query=.tmp_/target.fa
if grep -q '-' .tmp_/gggenome_location; then # 2>/dev/null 1>/dev/null
    ./miniogenotype/src/revcomp.sh .tmp_/target.fa > .tmp_/target_rev.fa
    query=.tmp_/target_rev.fa
fi
lalign36 -m 3 ${reference} ${query} |
grep "100.0%" |
cut -d ":" -f 2 |
sed -e "s/^/@ /g" -e "s/-/ /g" -e "s/)//g" \
> .tmp_/lalign_mut_location

cat .tmp_/gggenome_location | sed "s/^/@ /g" |
join - .tmp_/lalign_mut_location |
awk 'BEGIN{OFS="\t"}{print $2,$3+$6-20,$3+$7+20}' |
sort -k1,1 -k2,2n |
awk 'BEGIN{OFS="\t"}{print $0, NR}' \
> .tmp_/mut_location.bed


# ======================================
# Pairwise alignment between KI sequence and reads
# ======================================

# barcode=barcode14
for barcode in $(cat .tmp_/prediction_barcodelist | sed "s/@@@//g"); do
    bam=bam/${barcode}.bam
    printf "${bam} is processing..."
    #
    samtools view ${bam} |
    sort |
    join -1 1 - -2 2 .tmp_/sorted_prediction_result |
    cut -d " " -f 1-5,10 \
    > .tmp_/sorted_bam
    #
    rm -rf .tmp_/tmp_bam 2>/dev/null 1>/dev/null
    mkdir -p .tmp_/tmp_bam
    #
    cat .tmp_/sorted_bam |
    awk '{$6="HOGE"$6; print}' |
    sed -e "s/^/>/g" -e "s/HOGE/\n/g" |
    split -l 2 - .tmp_/tmp_bam/split_
    #
    time find .tmp_/tmp_bam | grep split_ |
    xargs -I {} \
    ./miniogenotype/src/format_laligh.sh .tmp_/mutation.fa {} > .tmp_/lalign.txt # 15 min...
    #
    progress_total=$(find .tmp_/tmp_bam | grep split_ | wc -l)
    true > .tmp_/progress
    for ten_break in $(awk BEGIN'{for (i=10; i<=90; i=i+10) print i}'); do
        awk -v prog=${progress_total} \
        -v ten=${ten_break} \
        'BEGIN{print "@"int(prog*ten*0.01)"@", ten}' >> .tmp_/progress
    done
    #
    printf "Detect loxP intactness...\n"
    true > .tmp_/lalign.txt
    progress_i=1
    time for fa in $(find .tmp_/tmp_bam | grep split_); do
        ./miniogenotype/src/format_laligh.sh .tmp_/mutation.fa ${fa} >> .tmp_/lalign.txt
        #
        grep @${progress_i}@ .tmp_/progress |
        awk '{print "Approx "$2"% complete"}'
        progress_i=$((progress_i + 1))
    done
    #
    cat .tmp_/lalign.txt |
    sort -t " " |
    join - .tmp_/sorted_bam |
    awk '{$NF=""; print $0}' |
    awk 'BEGIN{OFS="\t"} {print $7,$8+$2,$8+$3,$5,$1"@@@"$4"___"$6}' |
    bedtools intersect -a - -b .tmp_/mut_location.bed -wb |
    awk 'BEGIN{OFS="\t"} {print $NF"HOGE"$5,$4}' |
    sed "s/@@@/ /g" |
    awk '{if(max[$1]<$3) {max[$1]=$3; seq[$1]=$0}}
    END{for(key in max) print seq[key]}' |
    sed "s/___/ /g" |
    cut -d " " -f 1-2 |
    sed -e "s/^/>/g" -e "s/ /\n/g" \
    > .tmp_/mut_bestscore.fa
    
    # separate left and right insertion site
    cat .tmp_/mut_bestscore.fa |
    awk '{if($1 !~ /^>/) print "@@@"$0; else print}' |
    tr -d "\n" | sed -e "s/>/\n>/g" -e "s/@@@/ @@@/g" |
    grep "1HOGE" | sed -e "s/1HOGE//g" -e "s/@@@/\n/g" |
    sed "s/-//g" \
    > .tmp_/mut_bestscore_1.fa
    
    cat .tmp_/mut_bestscore.fa |
    awk '{if($1 !~ /^>/) print "@@@"$0; else print}' |
    tr -d "\n" | sed -e "s/>/\n>/g" -e "s/@@@/ @@@/g" |
    grep "2HOGE" | sed -e "s/2HOGE//g" -e "s/@@@/\n/g" \
    > .tmp_/mut_bestscore_2.fa
    
    printf "Output sequence logo at loxP loci..."
    # Multiple alignment by clustal omega
    clustalo -i .tmp_/mut_bestscore_1.fa \
    --dealign --iter 5 --threads=${threads:-4} \
    > .tmp_/clustalo_1.fa 2>/dev/null
    #
    clustalo -i .tmp_/mut_bestscore_2.fa \
    --dealign --iter 5 --threads=${threads:-4} \
    > .tmp_/clustalo_2.fa 2>/dev/null
    # Sequence logo
    ## PNG
    weblogo -n 50 --errorbars no -c classic --format png_print \
    < .tmp_/clustalo_1.fa > results/figures/png/seqlogo_left.png
    weblogo -n 50 --errorbars no -c classic --format png_print \
    < .tmp_/clustalo_2.fa > results/figures/png/seqlogo_right.png
    ## SVG
    weblogo -n 50 --errorbars no -c classic --format svg \
    < .tmp_/clustalo_1.fa > results/figures/svg/seqlogo_left.svg
    weblogo -n 50 --errorbars no -c classic --format svg \
    < .tmp_/clustalo_2.fa > results/figures/svg/seqlogo_right.svg
done