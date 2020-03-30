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

# ======================================
# Parse arguments
# ======================================

# input=".DAJIN_temp/fasta_ont/barcode30.fa"
# input="barcode30.fa"
# genotype="target"
# alleletype_original=${alleletype}
# [ "$alleletype_original" = "target" ] && pid="HOGE"
# [ "$alleletype_original" = "wt" ] && pid="FUGA"
# [ "$alleletype_original" = "abnormal" ] && pid="FOO"
# suffix="${barcode}"_"${alleletype}"_"${pid}"
# [ "$alleletype" = "abnormal" ] && alleletype="wt"
# echo $suffix

# insertion_skip="control"

input=${1}
genotype=${2}
pid=${3}

set +u
if [ "${4}" = "" ]; then
    insertion_skip=""
else
    insertion_skip=${4}
fi
set -u

parent_dir=".DAJIN_temp"
output=$(echo "${input}" | sed -e "s#.*/##g" -e "s#\..*##g" -e "s/_aligned_reads//g")
output_MIDS="${parent_dir}/data/MIDS_${output}_${pid}"
tmp_mapping="${parent_dir}/tmp_${output}_mapping_${pid}"
tmp_seqID="${parent_dir}/tmp_${output}_seqID_${pid}"

# ======================================
# Mapping
# ======================================

reference=".DAJIN_temp/fasta_conv/${genotype}.fa"
ref=$(cat "${reference}" | grep "^>" | sed "s/>//g")
#
minimap2 -t 1 --cs=long -ax map-ont "${reference}" "${input}" 2>/dev/null |
awk -v ref="${ref}" '$3 == ref' > "${tmp_mapping}"

# ======================================
# Identify mutation sites
# ======================================
reflength=$(cat "${reference}" | grep -v "^>" | awk '{print length($0)}')
ext=${ext:=100}
#
first_flank=$(cat ${parent_dir}/data/mutation_points | awk -v ext=${ext} '{print $1-ext}')
second_flank=$(cat ${parent_dir}/data/mutation_points | awk -v ext=${ext} '{if(NF==2) print $2+ext; else print $1+ext}')
if [ "$first_flank" -lt 1 ]; then first_flank=1; fi
if [ "$second_flank" -gt "$reflength" ]; then second_flank=$(($reflength)); fi

# ======================================
# MIDS conversion
# ======================================

label=$(echo "${input}" | sed -e "s#.*/##g" -e "s#\..*##g")
# printf "${label} is now processing...\n" 1>&2

cat "${tmp_mapping}" |
grep -v "^@" |
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
awk -v first=${first_flank} -v second=${second_flank} '{
    if($2<=first && $3>=second) print $1}' |
sort \
> "${tmp_seqID}"
#
cat "${tmp_mapping}" |
sort |
join - "${tmp_seqID}" |
# join "${tmp_mapping}" - |
# # Remove unpredictable long reads あまりにも長いリードは除去する
# awk -v reflen=${reflength} \
#     'length($10) < reflen*1.1' |
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
# もしSecondaryしか存在せず、全てがIになったリードがあれば除去する。
awk '$3 !~ /^I+$/' |
# replace matched nuc to "M"
awk '{cstag=$3;
    gsub(/[A|C|G|T]/, "M", cstag);
print $1,$2,cstag}' |
# replace insertion to "I"
# WT配列をDeletion配列にマップした場合（Clusteringのときのみの限定的条件）、
# IIIIIの連続配列が邪魔をして後方の配列が取り除かれてしまう。
# そのためIIIIの連続配列をトリムして配列長を保つ。
if [ -z "${insertion_skip}" ]; then
awk '{seq=cstag=$3
    insertion_num=gsub("+","+",cstag);
    if(insertion_num > 0){
    gsub(/[=M]/,"",$3)
    cnt=split($3,ins,"+");
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
}' 
else
    # Delete insertion
    sed "s/+[a-z]*\([-|=|*]\)/\1/g"
fi |

# replace single-nuc substitution to "S"
awk '{gsub(/\*[a|t|g|c]*/, "S", $3); print $0}' |
# replace Deletion to "D"
awk '{gsub(/[a|t|g|c]/, "D", $3); print $0}' |
# erase remained symbols
awk '{gsub(/[-|+|=]/, "", $3); print $0}' |
# # replace splice sites to "D"
# awk '{seq=cstag=$3;
#     gsub("~","",cstag);
#     gsub("[^0-9~]","",seq);
#     sub("^~","",seq);
#     count=split(seq, splice_num, "~");
#     for (i=1; i<=count; i++){
#         SPLICE="D"
#         for (j=1; j<splice_num[i]; j++){
#             SPLICE=SPLICE"D"
#         }
#         sub(splice_num[i],SPLICE,cstag)
#     }
#     print $1,$2,toupper(cstag)
# }' |
# complement seqences to match sequence length (insert "=")
## start
awk '{start=""; for(i=1; i < $2; i++) start=start"="; print $1,$2,start""$3}' |
## end
awk -v reflen="${reflength}" '{
    seqlen=length($3);
    end="";
    if(seqlen>reflen) {print $1,substr($3,1,reflen)}
    else {for(i=seqlen; i < reflen; i++) end=end"="; print $1,$3""end}
}' |
sed -e "s/$/\t${label}/g" -e "s/ /\t/g" \
> "${output_MIDS}"
#

# printf "${output} is now converted...\n"

exit 0