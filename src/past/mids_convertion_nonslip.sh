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
# Parse arguments
# ======================================

que="$1"
ref=.DAJIN_temp/fasta/${2}.fa

# ============================================================
# label
# ============================================================
label=$(echo "$que" |
    sed -e "s|.*/||g" \
        -e "s|_aligned.*||g" \
        -e "s|.fa*||g")
echo "$label"
wt_seqlen=$(awk '!/[>|@]/ {print length($0)}' "$ref")

# ============================================================
# mappinng
# ============================================================

minimap2 -ax map-ont "$ref" "$que" --cs=long 2>/dev/null |
grep -v "^@" |
awk -v ref="${2}" '$3 == ref' |
awk '$2==0 || $2==16' |
# cat - > test.sam
# ============================================================
# MIDS conv
# ============================================================
# cat test.sam | grep target#c140G-C_1408_aligned_35356_F_1118_178_1500
awk '{print $1, $2, $4, $(NF-1)}' |
awk -v seqlen="$wt_seqlen" -v label="$label" \
'{
    # annotate
    id=$1; flag=$2; start=$3; cs=$4
    # -------------------------------
    # CS tag MIDS formatting...
    # -------------------------------
    sub("cs:Z:","",cs)
    gsub("=", "", cs)
    # Deletion
    cs_del=""
    split(cs,array,"-")
    for(key in array){
        del="D"
        del_num=match(array[key], "^[a|c|g|t]+")
        for(i=1; i<RLENGTH; i++) del=del"D"
        sub("^[a|c|g|t]+",del,array[key])
        cs_del=cs_del""array[key]
    }
    cs=cs_del
    # Substitution
    gsub("*[a|c|g|t][a|c|g|t]", "S", cs)
    # Insertion
    gsub("+[a|c|g|t]+.", "I", cs)
    # Match
    gsub("[A|C|G|T]", "M", cs)
    # -------------------------------
    # Padding start sites...
    # -------------------------------
    seq=""
    for(i=1; i<start; i++) seq=seq"="
    cs=seq""cs

    # print cs

    # -------------------------------
    # Padding or trimming end sites...
    # -------------------------------
    # もし短い場合は"="を足す。長い場合は切る：
    if(length(cs)<seqlen){
        seq=""
        for(i=1; i<=seqlen-length(cs); i++) seq=seq"="
        cs=cs""seq
    } else {
        cs=substr(cs,1,seqlen)
    }

    print id, cs, label
}' |
sed "s/ /\t/g"

exit 0