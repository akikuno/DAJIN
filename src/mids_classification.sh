#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

################################################################################
#! I/O naming
################################################################################
#===========================================================
#? Auguments
#===========================================================
# input_fa=".DAJIN_temp/fasta_ont/wt_simulated_aligned_reads.fasta"
# genotype="wt"
# label=$(echo "${input_fa}" | sed -e "s#.*/##g" -e "s#\..*##g" -e "s/_aligned_reads//g")
# suffix="${label}_${genotype}"

input_fa=${1}
genotype=${2}
label=$(echo "${input_fa}" | sed -e "s#.*/##g" -e "s#\..*##g" -e "s/_aligned_reads//g")
suffix="${label}_${genotype}"

#===========================================================
#? Input
#===========================================================

#===========================================================
#? Output
#===========================================================
output_MIDS=".DAJIN_temp/data/MIDS_${suffix}"

#===========================================================
#? Temporal
#===========================================================
tmp_mapping=".DAJIN_temp/tmp_mapping_${suffix}"_$$
tmp_seqID=".DAJIN_temp/tmp_seqID_${suffix}"_$$

tmp_all=".DAJIN_temp/tmp_all_${suffix}"_$$
tmp_primary=".DAJIN_temp/tmp_primary_${suffix}"_$$
tmp_secondary=".DAJIN_temp/tmp_secondary_${suffix}"_$$


################################################################################
#! Function
################################################################################

mids_conv(){
    set /dev/stdin
    cat "${1}" |
        gawk '{id=$1; strand=$3; loc=$4; $0=$5
        sub("cs:Z:","",$0)
        gsub(/[ACGT]/, "M", $0)
        gsub(/\*[acgt][acgt]/, " S", $0)
        gsub("=", " ", $0)
        gsub("\+", " +", $0)
        gsub("-", " -", $0)
        for(i=1; i<=NF; i++){
            if($i ~ /^\+/){
                len=length($i)-1
                str=sprintf("%"len"s","")
                gsub(/ /,"I",str)
                $i=str }
            else if($i ~ "^-"){
                len=length($i)-1
                str=sprintf("%"len"s","")
                gsub(/ /,"D",str)
                $i=str }
            }
        gsub(" ", "", $0)
        print id, loc, $0}' |
    sort -t " " |
    cat
}

################################################################################
#! 変異部から±100塩基を含むリードのみを取り出す
################################################################################

#===========================================================
#? 変数の定義
# 切断面から含むべき塩基数
#===========================================================

reference=".DAJIN_temp/fasta_conv/wt.fa"
reflength=$(cat "${reference}" | grep -v "^>" | awk '{print length($0)}')
ref=$(cat "${reference}" | grep "^>" | sed "s/>//g")

ext=${ext:=100}
first_flank=$(
    cat .DAJIN_temp/data/mutation_points |
    awk -v ext=${ext} '{print $1-ext}'
    )
[ "${first_flank}" -lt 1 ] && first_flank=1

second_flank=$(
    cat .DAJIN_temp/data/mutation_points |
    awk -v ext=${ext} '{if(NF==2) print $2+ext; else print $1+ext}'
    )
[ "${second_flank}" -gt "${reflength}" ] && second_flank="${reflength}"

#===========================================================
#? 変異部から±100塩基を含むリードのみを取り出す
#===========================================================

minimap2 -ax map-ont "${reference}" "${input_fa}" --cs=long 2>/dev/null |
    awk -v ref="${ref}" -v reflen="${reflength}" '$3 == ref && length($10) < reflen * 1.1' |
    tee "${tmp_mapping}" |
    grep -v "^@" |
    # fetch sequence start and end sites
    awk 'BEGIN{OFS="\t"}{
        cigar=$6;
        gsub("[0-9]*S","",cigar);
        gsub("[0-9]*H","",cigar);
        gsub("M|D|I|N","\t",cigar);
        gsub(/\+$/,"",cigar);
        print $1, $4, cigar}' |
    awk '{sum=0; for(i=3; i<=NF; i++){ sum+=$i }
        print $1,$2,$2+sum}' |
    sort -t " " -n |
    awk '{if(length(min[$1])==0) min[$1]="inf";
        if(min[$1]>$2) min[$1]=$2;
        if(max[$1]<$3) max[$1]=$3}
        END{for(key in min) print key, min[key], max[key]}' |
    sort -t " " -n |
    awk -v first=${first_flank} -v second=${second_flank} '{
        if($2<=first && $3>=second) print $1}' |
    sort |
cat > "${tmp_seqID}"

################################################################################
#! MIDS conversion
################################################################################

#===========================================================
#? primaryとsecondaryを分けてMIDS変換を行う
#===========================================================

cat "${tmp_mapping}" |
    sort |
    join - "${tmp_seqID}" |
    # append alignment info
    awk '{
        if($2==0 || $2==16) {alignment="primary"} else {alignment="secondary"};
        if($2==0 || $2==2048) {strand="plus"} else {strand="minus"};
        for(i=1;i<=NF;i++) if($i ~ /cs:Z/) print $1,alignment,strand,$4,$i
    }' |
cat > "${tmp_all}"

cat "${tmp_all}" |
    grep "primary" |
    mids_conv |
cat > "${tmp_primary}"

cat "${tmp_all}" |
    grep "secondary" |
    mids_conv |
cat > "${tmp_secondary}"

#===========================================================
#? primaryとsecondaryを結合する
#===========================================================

cat "${tmp_primary}" "${tmp_secondary}" |
    sort -t " " -k 2,2n |
    awk '{seq_[$1]=seq_[$1]" "$2" "$3" "$4}
        END{for(key in seq_) print key,seq_[key]}' |
    # grep inversion_100_aligned_7213_F_54_2709_73 |
    awk 'NF==3
    #---------------------------------------
    #* Flox deletion
    #---------------------------------------
    NF==5{id=$1; minloc=$2; maxloc=$2 ; seq=""; len=0; str=""
        for(i=2; i<=NF; i=i+2){
            seq_[$i]=$(i+1)
            if(minloc>$i){ minloc=$i; len=length($(i+1))}
            if(maxloc<$i) maxloc=$i
        }
        len=maxloc-minloc-len
        str=sprintf("%"len"s","")
        gsub(/ /,"D",str) 
        for(key in seq_) seq=seq seq_[key] str
        print id, minloc, seq
        delete seq_
    }
    #---------------------------------------
    #* Inversion
    #---------------------------------------
    NF==7{id=$1; minloc=$2; maxloc=$2 ; seq=""; len=0; str=""
        for(i=2; i<=NF; i=i+2){
            if(minloc>$i){ minloc=$i; len=length($(i+1))}
            if(maxloc<$i) maxloc=$i
        }
        len=length($5)
        str=sprintf("%"len"s","")
        gsub(/ /,"=",str) 
        $5=str
        print id, minloc, $3 $5 $7
    }' |
    sed "s/D*$//g" |
    # complement seqences to match sequence length (insert "M")
    ## start
    awk '{start=""; for(i=1; i < $2; i++) start=start"M"; print $1,$2,start""$3}' |
    ## end
    awk -v reflen="${reflength}" '{
        seqlen=length($3);
        end="";
        if(seqlen>reflen) {print $1,substr($3,1,reflen)}
        else {for(i=seqlen; i < reflen; i++) end=end"M"; print $1,$3""end}
    }' |
    # 全てが変異になったリードがあれば除去する。
    awk '$2 !~ /^[I|D|S]+$/' |
    sed -e "s/$/\t${label}/g" -e "s/ /\t/g" | 
cat > "${output_MIDS}"

rm "${tmp_mapping}" "${tmp_seqID}" "${tmp_all}" "${tmp_primary}" "${tmp_secondary}"
exit 0