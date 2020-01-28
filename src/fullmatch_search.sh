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
# Prepare sorted file
# ======================================

cat .tmp_/prediction_result.txt |
grep target |
sort -k 2,2 \
> .tmp_/sorted_prediction_result

# ======================================
# Detect Mutation type
# ======================================

reference=fasta/wt.fa
query=fasta/target.fa

minimap2 -a ${reference} ${query} --cs 2>/dev/null |
awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
sed -e "s/cs:Z:://g" -e "s/:/\t/g" > .tmp_/mutation_profile

grep "+" .tmp_/mutation_profile 1>/dev/null 2>/dev/null
if test "$?" -eq 1; then
    printf "The target mutation is not knock-in.\nfinished.\n"
    exit 0
fi

# ======================================
# Extract mutation matched sequence
# ======================================

mutation_profile_for="ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
printf "${mutation_profile_for}\n" > .tmp_/tmp_mutation_profile_for
mutation_profile_rev=$(./miniogenotype/src/revcomp.sh .tmp_/tmp_mutation_profile_for)

output_file=.tmp_/target_perfectmatch.csv
true > ${output_file}
for fasta in $(ls fasta_ont/*fa | grep -v -e simulated); do
    label=$(echo ${fasta} | sed -e "s#.*/##g" -e "s/\..*//g")
    cat ${fasta} |
    awk '{if(NR%2 != 0) print "1STFOO"$0; else print "2NDBAR"$0}' |
    tr -d "\n" |
    sed -e "s/1STFOO@/\n/g" -e "s/2NDBAR/\t/g" |
    awk '{print $1,$NF}' |
    sort |
    join -1 1 - -2 2 .tmp_/sorted_prediction_result |
    cut -d " " -f 2 |
    sed -e "s/${mutation_profile_for}/ /g" \
    -e "s/${mutation_profile_rev}/ /g" > .tmp_/tmp_mut
    #
    if [ $(cat .tmp_/tmp_mut | wc -l) -ne 0 ]; then
    cat .tmp_/tmp_mut |
    awk '{print NF}' |
    awk -v label=${label} '{
        if($1==1 || length($1)==0) {print label",no-match"}
        else {print label",match x"($1-1)}}' \
    >> ${output_file}
    else
    printf "${label},no-match\n" >> ${output_file}
    fi
done
