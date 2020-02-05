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
    if($3 ~ "target") barcode_target[$1]++}
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
# for: ATAACTTCGTATAATGTATGCTATACGAAGTTAT
# rev: ATAACTTCGTATAGCATACATTATACGAAGTTAT
# grep -e ATAACTTCGTATAATGTATGCTATACGAAGTTAT -e ATAACTTCGTATAGCATACATTATACGAAGTTAT | sort | uniq -c
#  mutation_profile_for=$1
mutation_profile="ATAACTTCGTATAATGTATGCTATACGAAGTTAT"

printf ">mut\n${mutation_profile}\n" \
> .tmp_/mutation.fa

reference=.tmp_/mutation.fa
query=.tmp_/target.fa
if grep -q '-' .tmp_/gggenome_location; then # 2>/dev/null 1>/dev/null
    ./miniogenotype/src/revcomp.sh .tmp_/target.fa > .tmp_/target_rev.fa
    query=.tmp_/target_rev.fa
fi
lalign36 -m 3 ${reference} ${query} |
grep "100.0%" |
cut -d ":" -f 2 |
sed -e "s/^/@ /g" -e "s/-/ /g" -e "s/)//g" |
sort -t " " -k2,2n
> .tmp_/lalign_mut_location

cat .tmp_/lalign_mut_location |
awk '{print int(($3+$2)/2)}' \
> .tmp_/lalign_mut_center

flank1=$(cat .tmp_/lalign_mut_center | head -n 1)
flank2=$(cat .tmp_/lalign_mut_center | tail -n 1)

# ======================================
# Pairwise alignment between KI sequence and reads
# ======================================

# barcode=barcode14
mkdir -p results/figures/png/seqlogo/ results/figures/svg/seqlogo/
for barcode in $(cat .tmp_/prediction_barcodelist | sed "s/@@@//g"); do
    bam=bam/${barcode}.bam
    printf "##########\n${bam} is processing...\n##########\n"
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
    xargs -I {} ./miniogenotype/src/format_laligh.sh .tmp_/mutation.fa {} \
    > .tmp_/lalign1.fa & } 1>/dev/null 2>/dev/null
    #
    { find .tmp_/split |
    grep split2_ |
    xargs -I {} ./miniogenotype/src/format_laligh.sh .tmp_/mutation.fa {} \
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
            # echo "$i --------" #! ==============================
            cat ${input} |
            awk -F "" -v i=${i} '{if(NR%2==0) print $i}' |
            sort |
            uniq -c |
            # awk -v i=${i} '{sum+=$1; if($2=="-") gapnum=$1}
            # END{print i,$1/sum*100}' |
            # awk '$2>50' |
            awk -v i=${i} '{sum+=$1; if(max<$1) {max=$1; nuc=$2}}
            END{print i,nuc,max/sum*100}' |
            #Extract nucleotide position with gap "-" > 20%
            awk '$2 == "-" && $3>50' |
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
        { weblogo -n 50 --errorbars no -c classic --format png_print \
            < ${output_rmgap} > results/figures/png/seqlogo/${barcode}_${output_logo}.png & } 1>/dev/null 2>/dev/null
        ## SVG
        { weblogo -n 50 --errorbars no -c classic --format svg \
            < ${output_rmgap} > results/figures/svg/seqlogo/${barcode}_${output_logo}.svg & } 1>/dev/null 2>/dev/null
        wait 1>/dev/null 2>/dev/null
    done
done

output_logo=$(echo ${output_rmgap} | sed -e "s/.*clustalo_//g" -e "s/_rmgap.fa//g")
{ weblogo -n 50 --errorbars no -c classic --format png_print \
    < .tmp_/clustalo_right.fa > results/figures/png/seqlogo/negacon.png & } 1>/dev/null 2>/dev/null

# Sequence logo
## PNG
# { weblogo -n 50 --errorbars no -c classic --format png_print \
#     < .tmp_/clustalo1.fa > results/figures/png/seqlogo/${barcode}_left.png & } 1>/dev/null 2>/dev/null
# { weblogo -n 50 --errorbars no -c classic --format png_print \
#     < .tmp_/clustalo2.fa > results/figures/png/seqlogo/${barcode}_right.png & } 1>/dev/null 2>/dev/null
# ## SVG
# { weblogo -n 50 --errorbars no -c classic --format svg \
#     < .tmp_/clustalo1.fa > results/figures/svg/seqlogo/${barcode}_left.svg & } 1>/dev/null 2>/dev/null
# { weblogo -n 50 --errorbars no -c classic --format svg \
#     < .tmp_/clustalo2.fa > results/figures/svg/seqlogo/${barcode}_right.svg & } 1>/dev/null 2>/dev/null
# wait 1>/dev/null 2>/dev/null
