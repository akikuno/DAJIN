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
# Extract mutation matched sequence
# ======================================

mutation_profile_for="ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
printf "${mutation_profile_for}\n" > .tmp_/tmp_mutation_profile_for
mutation_profile_rev=$(./miniogenotype/src/revcomp.sh .tmp_/tmp_mutation_profile_for)

# output_file=.tmp_/target_perfectmatch.csv
# true > ${output_file}
# for fasta in $(ls fasta_ont/*fa | grep -v -e simulated); do
#     label=$(echo ${fasta} | sed -e "s#.*/##g" -e "s/\..*//g")
#     cat ${fasta} |
#     awk '{if(NR%2 != 0) print "1STFOO"$0; else print "2NDBAR"$0}' |
#     tr -d "\n" |
#     sed -e "s/1STFOO@/\n/g" -e "s/2NDBAR/\t/g" |
#     awk '{print $1,$NF}' |
#     sort |
#     join -1 1 - -2 2 .tmp_/sorted_prediction_result |
#     cut -d " " -f 2 |
#     sed -e "s/${mutation_profile_for}/ /g" \
#     -e "s/${mutation_profile_rev}/ /g" > .tmp_/tmp_mut
#     #
#     if [ $(cat .tmp_/tmp_mut | wc -l) -ne 0 ]; then
#     cat .tmp_/tmp_mut |
#     awk '{print NF}' |
#     awk -v label=${label} '{
#         if($1==1 || length($1)==0) {print label",no-match"}
#         else {print label",match x"($1-1)}}' \
#     >> ${output_file}
#     else
#     printf "${label},no-match\n" >> ${output_file}
#     fi
# done


# ======================================
# Search reads with loxP exact matching
# ======================================

flank1=$(cat .tmp_/lalign_mut_center | head -n 1)
flank2=$(cat .tmp_/lalign_mut_center | tail -n 1)

true > .tmp_/target_perfectmatch.csv

# barcode=barcode14
mkdir -p results/figures/png/seqlogo/ results/figures/svg/seqlogo/
for barcode in $(cat .tmp_/prediction_barcodelist | sed "s/@@@//g"); do
    bam=bam/${barcode}.bam
    # printf "##########\n${bam} is processing...\n##########\n"
    #
    samtools view ${bam} |
    sort |
    join -1 1 - -2 2 .tmp_/sorted_prediction_result |
    cut -d " " -f 1-5,10 \
    > .tmp_/sorted_bam
    #
    cat .tmp_/sorted_bam |
    awk -v f1=${flank1} -v f2=${flank2} '{
        print $1,substr($6,f1-100, 200) > ".tmp_/mutsite_split1"
        print $1,substr($6,f2-100, 200) > ".tmp_/mutsite_split2"}'
    #
    cat .tmp_/mutsite_split1 |
    grep -e ${mutation_profile_for} -e ${mutation_profile_rev} |
    cut -d " " -f 1 |
    sort \
    > .tmp_/mutsite_targetpositive1
    cat .tmp_/mutsite_split2 |
    grep -e ${mutation_profile_for} -e ${mutation_profile_rev} |
    cut -d " " -f 1 |
    sort \
    > .tmp_/mutsite_targetpositive2
    # Output results...
    join .tmp_/mutsite_targetpositive1 .tmp_/mutsite_targetpositive2 |
    sed "s/^/${barcode},flox\t/g" |
    cut -f 1 \
    >> .tmp_/target_perfectmatch.csv
    join -v 1 .tmp_/mutsite_targetpositive1 .tmp_/mutsite_targetpositive2 |
    sed "s/^/${barcode},left loxP\t/g" |
    cut -f 1 \
    >> .tmp_/target_perfectmatch.csv
    join -v 2 .tmp_/mutsite_targetpositive1 .tmp_/mutsite_targetpositive2 |
    sed "s/^/${barcode},right loxP\t/g" |
    cut -f 1 \
    >> .tmp_/target_perfectmatch.csv
    cat .tmp_/mutsite_targetpositive1 .tmp_/mutsite_targetpositive2 |
    sort -u |
    join -v 1 .tmp_/sorted_bam - |
    sed "s/^/${barcode},No exact match\t/g" |
    cut -f 1 >> .tmp_/target_perfectmatch.csv
done
