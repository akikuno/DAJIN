mut_length=$(cat "${1}" | tail -n 1 | awk '{print length($0)}')
label=$(cat "${2}" | head -n 1 | cut -d " " -f 1 | sed "s/^>//g")

#lalign36 -m 3 mutation.fa split_aj |
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
awk -v mut_len=${mut_length} '{
    score=$1
    start=$2-1
    ## design end with gap consideration
    match($5,"-+")
    if(RLENGTH != -1) end=mut_len-$3-RLENGTH
    else end=0
    # print start,end, mut_len, $3, RLENGTH
    s_seq=""; e_seq=""
    for(i=1;i<=int(start);i++) s_seq=s_seq"-"
    for(i=1;i<=int(end);i++) e_seq=e_seq"-"
    #print s_seq, e_seq
    print $6"_"$1"\n"s_seq,$NF,e_seq
}' |
sed "s/ //g"
# Extract 5(50-6+1) and 12(50-38)
awk -v mut_len=${mut_length} '{
    if($0 ~ "similar") {
        # Score store
        score=$0
        sub("% similar.*", "", score)
        # Sequence imputation
        s_seq=""; e_seq=""
        start=$NF
        sub("\\(", "", start)
        sub("-.*", "", start)
        start=start-1
        for(i=1;i<=int(start);i++) s_seq=s_seq"-"
        #
        end=$NF
        sub("^.*-", "", end)
        for(i=1;i<=int(mut_len-end);i++) e_seq=e_seq"-"
        # print start, int(mut_len-end), s_seq, e_seq
    }
    else if($0 ~ "^>mut") {match($2, /\-+/); for(i=1;i<=RLENGTH;i++) e_seq=e_seq"-"; print $1, s_seq$2e_seq, score }
    else {print $1, s_seq$2e_seq, score}
}'






# # Extract 5(50-6+1) and 12(50-38)
# awk -v mut_len=${mut_length} '{
#     if($0 ~ "identity") {
#         # Score store
#         score=$0
#         sub("% identity.*", "", score)
#         # Sequence imputation
#         s_seq=""; e_seq=""
#         start=$NF
#         sub("\\(", "", start)
#         sub("-.*", "", start)
#         start=start-1
#         for(i=1;i<=int(start);i++) s_seq=s_seq"-"
#         #
#         end=$NF
#         sub("^.*-", "", end)
#         for(i=1;i<=int(mut_len-end);i++) e_seq=e_seq"-"
#         print start, int(mut_len-end), s_seq, e_seq
#     }
#     else if($0 !~ "^>") {allseq=s_seq$0e_seq; print allseq}
#     else {gsub(" ..", "", $0); print $0"_"score}
#     }'


