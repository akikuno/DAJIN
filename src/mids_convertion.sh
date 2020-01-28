#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================

set -u
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export PATH="$(command -p getconf PATH)${PATH+:}${PATH-}"
export UNIX_STD=2003  # to make HP-UX conform to POSIX

reflength=$(cat fasta/wt.fa | grep -v "^>" | awk '{print length($0)}')

label=$(echo $1 | sed -e "s#.*/##g" -e "s#\..*##g")
printf "$label is now processing...\n" 1>&2

#---------------------------------
cat $1 |
# fetch sequence start and end sites
awk 'BEGIN{OFS="\t"}{
    cigar=$6;
    gsub("[0-9]*S","",cigar);
    gsub("[0-9]*H","",cigar);
    gsub("M|D|I|N","\t",cigar);
    gsub("+$","",cigar);
    print $1, $4, cigar}' |
awk '{sum=0; for(i=3; i<=NF; i++){ sum+=$i }; print $1,$2,$2+sum}' |
sort -t " " -n |
awk '{if(length(min[$1])==0) min[$1]="inf";
    if(min[$1]>$2) min[$1]=$2;
    if(max[$1]<$3) max[$1]=$3}
    END{for(key in min) print key, min[key], max[key]}' |
sort -t " " -n |
awk -v first=${2} -v second=${3} '{
    if($2<=first && $3>=second) print $1}' |
sort > .tmp_/tmp_sequenceID

cat $1 |
sort |
join - .tmp_/tmp_sequenceID |
# append alignment info
awk '{
    if($2==0 || $2==16) {alignment="primary"} else {alignment="secondary"};
    for(i=1;i<=NF;i++) if($i ~ /cs:Z/) print $1,$4,alignment, $i
    }' |
# cut -d " " -f 1-3 |
sort -t " " -k 1,1 -k 2,2n |
# concatenate primary and secondary (secondary is converted to "I")
awk '{s=$2; alignment=$3; cstag=$4;
    gsub("cs:Z:","",cstag);
    if(alignment=="secondary"){
        ins="I"; for(i=1; i<length(cstag); i++){ins=ins"I"}
        cstag=ins
        }
    ID[$1]=ID[$1]""cstag;
    if(length(start[$1])==0) start[$1]="inf"
    if(start[$1]>s) start[$1]=s;
    }
    END {for (key in ID) print key, start[key],ID[key]}' |
# replace matched nuc to "M"
awk '{cstag=$3;
    gsub(/[A|C|G|T]/, "M", cstag);
    print $1,$2,cstag}' |
# replace insertion to "I"
awk '{seq=cstag=$3
    insertion_num=gsub("+","+",cstag);
    if(insertion_num > 0){
    gsub(/[=M]/,"",$3); cnt=split($3,ins,"+");
    ins_seq=""
    for(i=1;i<=cnt;i++){
        gsub("[-*][acgt]+", "", ins[i])
        ins_seq=ins_seq"+"ins[i]
    }
    sub("++","+",ins_seq)
    cnt=split(ins_seq,ins,"+");
    for(i=1;i<=cnt;i++){
        ins_seq="+"ins[i]
        sub(ins_seq, toupper(ins_seq), seq)
    }
    gsub(/[A|T|G|C]/,"I", seq)
    print $1,$2,seq
    }
    else {print $0}
    }' |
# replace single-nuc substitution to "S"
awk '{gsub(/\*[a|t|g|c]*/, "S", $3); print $0}' |
# replace Deletion to "D"
awk '{gsub(/[a|t|g|c]/, "D", $3); print $0}' |
# erase remained symbols
awk '{gsub(/[-|+|=]/, "", $3); print $0}' |
# replace splice sites to "D"
awk '{seq=cstag=$3;
    gsub("~","",cstag);
    gsub("[^0-9~]","",seq);
    sub("^~","",seq);
    count=split(seq, splice_num, "~");
    for (i=1; i<=count; i++){
        SPLICE="D"
        for (j=1; j<splice_num[i]; j++){
            SPLICE=SPLICE"D"
        }
        sub(splice_num[i],SPLICE,cstag)
    }
    print $1,$2,toupper(cstag)
    }' |
# complement seqences to match sequence length (insert "=")
## start
awk '{start=""; for(i=1; i < $2; i++) start=start"="; print $1,$2,start""$3}' |
## end
awk -v reflen=${reflength} '{
    seqlen=length($3);
    end="";
    if(seqlen>reflen) {print substr($3,1,reflen),$1}
    else {for(i=seqlen; i < reflen; i++) end=end"="; print $3""end,$1}
}' |
sed -e "s/^/${label}\t/g" -e "s/ /\t/g"