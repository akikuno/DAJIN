#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -eu
umask 0022
export LC_ALL=C
export UNIX_STD=2003 # to make HP-UX conform to POSIX

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? Auguments
#===========================================================

que_fa=${1}
genotype=${2}

label=${que_fa##*/}
label=${label%%_aligned*}
label=${label%.fa*}
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
tmp_mapping=".DAJIN_temp/data/tmp_mapping_${suffix}"
tmp_seqID=".DAJIN_temp/data/tmp_seqID_${suffix}"

tmp_all=".DAJIN_temp/data/tmp_all_${suffix}"
tmp_primary=".DAJIN_temp/data/tmp_primary_${suffix}"
tmp_secondary=".DAJIN_temp/data/tmp_secondary_${suffix}"

################################################################################
#! Function
################################################################################

mids_conv() {
    set /dev/stdin
    cat "${1}" |
        # long deletion
        sed "s/~[acgt][acgt]\([0-9][0-9]*\)[acgt][acgt]/ ~\1 /g" |
        awk '{for(i=5;i<=NF;i++){
            if($i ~ /\~/){
                sub("~","",$i)
                len=int($i)
                for(j=1; j<=len; j++) str = "D" str
                $i=str
                str=""
            }
        }}1' 2>/dev/null |
        awk '{printf $1" "$2" "$3" "$4" "
            for(i=5;i<=NF;i++) printf $i
            printf "\n"}' |
        # insertion/point mutation/inversion
        awk '{id=$1; strand=$3; loc=$4; $0=$5
            sub("cs:Z:","",$0)
            sub("D"," D",$0)
            gsub(/[ACGT]/, "M", $0)
            gsub(/\*[acgt][acgt]/, " S", $0)
            gsub("=", " ", $0)
            gsub("\\+", " +", $0)
            gsub("\\-", " -", $0)
            for(i=1; i<=NF; i++){
                if($i ~ /^\+/){
                    len=length($i)-1
                    for(len_=1; len_ <= len; len_++) str = "I" str
                    $i=str
                    str=""}
                else if($i ~ "^-"){
                    len=length($i)-1
                    for(len_=1; len_ <= len; len_++) str = "D" str
                    $i=str
                    str=""}
                }
            gsub(" ", "", $0)
        print id, loc, $0}' 2>/dev/null |
        sort -t " " |
        cat
}

################################################################################
#! Extract reads containing ±100 bases from the mutation
################################################################################

#===========================================================
#? Define variables
#===========================================================

ref_fa=".DAJIN_temp/fasta_conv/wt.fa"
ref_len=$(awk '$1!~/^>/ {print length}' "${ref_fa}")
ref_label="${ref_fa##*/}" && ref_label="${ref_label%.*}"

ext=${ext:=100}
first_flank=$(awk -v ext=${ext} '{print $1-ext}' .DAJIN_temp/data/mutation_points)
[ "${first_flank}" -lt 1 ] && first_flank=1

second_flank=$(awk -v ext=${ext} '{if(NF==2) $0=$2+ext; else $0=$1+ext}1' .DAJIN_temp/data/mutation_points)
[ "${second_flank}" -gt "${ref_len}" ] && second_flank="${ref_len}"

#===========================================================
#? Extract reads
#===========================================================

max_len=$(awk '$1!~/^>/ {if(max<length) max=length} END{print max}' .DAJIN_temp/fasta/fasta.fa)

minimap2 -ax splice "${ref_fa}" "${que_fa}" --cs=long 2>/dev/null |
    awk -v ref_label="${ref_label}" -v max_len="${max_len}" \
        '$3 == ref_label && length($10) < max_len * 1.1' |
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
    awk -v first="${first_flank}" -v second="${second_flank}" '{
        if($2<=first && $3>=second) print $1}' |
    sort |
    cat >"${tmp_seqID}"

################################################################################
#! MIDS conversion
################################################################################

#===========================================================
#? Separate primary and secondary reads
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
    cat >"${tmp_all}"

cat "${tmp_all}" |
    grep "primary" |
    mids_conv |
    cat >"${tmp_primary}"

cat "${tmp_all}" |
    grep "secondary" |
    mids_conv |
    cat >"${tmp_secondary}"

#===========================================================
#? Concat primary secondary
#===========================================================

cat "${tmp_primary}" "${tmp_secondary}" |
    sort -t " " -k 2,2n |
    awk '{seq_[$1]=seq_[$1]" "$2" "$3" "$4}
        END{for(key in seq_) print key,seq_[key]}' |
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
        for(len_=1; len_ <= len; len_++) str = "D" str
        for(key in seq_) seq=seq seq_[key] str
        print id, minloc, seq
        delete seq_
    }
    #---------------------------------------
    #* Inversion
    #---------------------------------------
    NF==7{id=$1; minloc=$2; maxloc=$2 ; seq=""; len=0; str=""
        for(i=2; i<=NF; i=i+2){
            if(minloc>$i){minloc=$i; len=length($(i+1))}
            if(maxloc<$i) maxloc=$i
        }
        len=length($5)
        for(len_=1; len_ <= len; len_++) str = "=" str
        $5=str
        print id, minloc, $3 $5 $7
    }' |
    sed "s/D*$//g" |
    # complement seqences to match sequence length (insert "M")
    ## start
    awk '{start=""; for(i=1; i < $2; i++) start=start"M"; print $1,$2,start""$3}' |
    ## end
    awk -v max_len="${max_len}" '{
        seqlen=length($3);
        end="";
        if(seqlen>max_len) {print $1,substr($3,1,max_len)}
        else {for(i=seqlen; i < max_len; i++) end=end"M"; print $1,$3""end}
    }' |
    # Remove all-mutated reads
    awk '$2 !~ /^[I|D|S]+$/' |
    sed "s/$/ ${label}/g" |
    awk '{gsub(" ", "\t")}1' |
    cat >"${output_MIDS}"

rm "${tmp_mapping}" "${tmp_seqID}" "${tmp_all}" "${tmp_primary}" "${tmp_secondary}"
exit 0
