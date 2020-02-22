#!/bin/sh
# TGATGCTCGGC ref >2655 zczc >4703 zfzw
# TGATGCTCGGT mut >249 jo >5147 zgqy
# grep ">8546" .tmp_/split/*
# grep ">7737" .tmp_/split/*

mut_length=$(cat "${1}" | tail -n 1 | awk '{print length($0)}')
file_name=$(echo "${2}" | sed "s#.*/##g")
# lalign36 -m 3 .tmp_/mutation.fa .tmp_/split/split_zevq |
lalign36 -m 3 "${1}" "${2}" |
sed -n "/identity/,/^$/p" |
grep -v "^$" |
sed "s/:.*//g" |
sed "s/ ..$/HOGE/g" |
sed "s/.*identity (/FUGA/g" |
tr -d "\n" |
sed -e "s/HOGE/ /g" -e "s/FUGA/\n/g" -e "s/>/\n>/g" |
#grep -v "mut" |
sed "s/\(^[0-9].*%\)/FUGA\1/g" |
tr -d "\n" |
sed -e "s/FUGA/\n/g" -e "s/>/ >/g" |
sed -e "s/% similar).*overlap (/ /g" -e "s/-/ /" |
grep -v "^$" |
sort -n |
tail -n 1 |
awk -v mut_len=${mut_length} 'length($NF)>mut_len-10' \
> .tmp_/split/tmp_${file_name}
#
# Reverse complement
#
 [ ! -s .tmp_/split/tmp_${file_name} ] && exit 1

start=$(cat .tmp_/split/tmp_${file_name} | cut -d " " -f 2)
end=$(cat .tmp_/split/tmp_${file_name} | cut -d " " -f 3)
seq=$(cat .tmp_/split/tmp_${file_name} | awk '{print $NF}')
query=".tmp_/split/tmp_${file_name}"

if [ "$start" -gt "$end" ]; then
    cat .tmp_/split/tmp_${file_name} | awk '{print $NF}' \
    > .tmp_/split/tmp_${file_name}_seq
    #
    revseq=$(./DAJIN/src/revcomp.sh .tmp_/split/tmp_${file_name}_seq)
    #
    cat .tmp_/split/tmp_${file_name} | sed "s/${seq}/${revseq}/g" |
    awk -v s=${start} -v e=${end} '{$2=e;$3=s; print}' \
    > .tmp_/split/tmp_${file_name}_rev
    #
    query=".tmp_/split/tmp_${file_name}_rev"
fi
#
cat ${query} |
awk -v mut_len=${mut_length} '{
    score=$1
    start=$2-1
    ## design end with gap consideration
    match($5,"-+")
    if(RLENGTH != -1) end=mut_len-$3-RLENGTH
    else end=mut_len-$3
    # print start,end
    s_seq=""; e_seq=""
    for(i=1;i<=int(start);i++) s_seq=s_seq"-"
    for(i=1;i<=int(end);i++) e_seq=e_seq"-"
    #print s_seq, e_seq
    print $1,$6,s_seq""$NF""e_seq
}'