#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================

set -u
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ======================================

label=$(cat "${2}" | head -n 1 | cut -d " " -f 1 | sed "s/^>//g")
mut_length=$(cat .tmp_/mutation.fa | tail -n 1 | awk '{print length($0)}')

#lalign36 -m 3 .tmp_/mutation.fa .tmp_/split/split2_zbed |
lalign36 -m 3 "${1}" "${2}" |
grep -v -e ">--" |
sed "s/ \.\./ /g" |
sed -e "s/Scan time.*$/SCANTIME/g" -e "s/^$/EMPTY_ROW/g" |
tr -d "\n" |
sed -e "s/SCANTIME/\n/g" -e "s/EMPTY_ROW/\n/g" |
grep -v -e "^$" -e "Scomplib" -e "LALIGN" |
grep "identity" |
sed "s/>/\n/g" |
grep -v "mut" |
awk '{if($0 ~ ")$") {gsub("-"," ",$0); print $0"@@@"}
    else {print $0"XXX"}}' |
tr -d "\n" |
sed -e "s/@@@/ /g" -e "s/XXX/\n/g" |
sed -e "s/% identity.*in / /" \
    -e "s/(//g" \
    -e "s/)//g" \
    -e "s/:/ /g" |
# Extract reads with the highest score
sort -k 1,1nr |
head -n 1 |
# Extract reads with sufficient seq length
awk -v mutlen=${mut_length} 'length($NF) > mutlen-3 && length($NF) < mutlen+3' |
# awk '{if(array[$3]<$1) {array[$3]=$1; seq[$3]=$0}}
# END{for(key in array) print seq[key]}' |
awk -v label=${label} '{print ">"label"_"$1"\n"$NF}'
