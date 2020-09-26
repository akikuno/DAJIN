#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
export UNIX_STD=2003  # to make HP-UX conform to POSIX


################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# barcode=barcode02
# alleletype="abnormal"

#===========================================================
#? Auguments
#===========================================================

barcode=${1}
alleletype=${2}

#===========================================================
#? Input
#===========================================================

suffix="${barcode}_${alleletype}"
mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

ref_fa=".DAJIN_temp/fasta/${mapping_alleletype}.fa"
que_fa=".DAJIN_temp/fasta_ont/${barcode}.fa"

#===========================================================
#? Output
#===========================================================

#===========================================================
#? Temporal
#===========================================================

tmp_que_fa=".DAJIN_temp/clustering/tmp_que_fa_${suffix}"_$$
tmp_seqID=".DAJIN_temp/clustering/tmp_seqID_${suffix}"_$$

tmp_all=".DAJIN_temp/tmp_all_${suffix}"_$$
tmp_primary=".DAJIN_temp/tmp_primary_${suffix}"_$$
tmp_secondary=".DAJIN_temp/tmp_secondary_${suffix}"_$$

################################################################################
#! Function definitions
################################################################################

mids_compressed(){
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
        awk '{id=$1; strand=$3; loc=$4; $0=$5;
        sub("cs:Z:","",$0)
        sub("D"," D",$0)
        gsub(/[ACGT]/, "M", $0)
        gsub(/\*[acgt][acgt]/, " S", $0)
        gsub("=", " ", $0)
        gsub("\+", " +", $0)
        gsub("\-", " -", $0)
        for(i=1; i<=NF; i++){
            if($i ~ /^\+/){
                len=length($i)-1
                if(len>=10 && len<=35) len=sprintf("%c", len+87)
                else if(len>=36) len="z"
                $i=" "
                $(i+1)=len substr($(i+1),2) }
            else if($i ~ /^\-/){
                len=length($i)-1
                for(len_=1; len_ <= len; len_++){
                    str= "D" str
                }
                $i=str
                str=""}
            }
        gsub(" ", "", $0)
        print id, loc, $0, strand}' 2>/dev/null |
    sort -t " " |
    cat
}

################################################################################
#! Extract fasta file including "mapping_alleletype"
################################################################################

cat .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    grep "${barcode}" |
    grep "${alleletype}" |
    cut -f 1 |
    sort |
cat > "${tmp_seqID}"

cat "${que_fa}" |
    awk '{printf $1"\t"}' |
    sed "s/>/\n/g" |
    sort |
    join - "${tmp_seqID}" |
    sed "s/^/>/g" |
    sed "s/ /\n/g" |
    grep -v "^$" |
cat > "${tmp_que_fa}"

################################################################################
#! Mapping
################################################################################

ref_len=$(cat "${ref_fa}" | grep -v "^>" | awk '{print length}')

minimap2 -ax splice "${ref_fa}" "${tmp_que_fa}" --cs=long 2>/dev/null |
    awk -v allele="${mapping_alleletype}" -v ref_len="${ref_len}" \
        '$3 == allele && length($10) < ref_len * 1.1' |
    sort |
    # append flag info
    awk '{
        if($2==0 || $2==16) {flag="primary"} else {flag="secondary"};
        if($2==0 || $2==2048) {strand="plus"} else {strand="minus"};
        for(i=1;i<=NF;i++) if($i ~ /cs:Z/) print $1,flag,strand,$4,$i
    }' |
    # remove all secondary reads
    awk '{flag[$1]=flag[$1]" "$2
        num_pri=sub("primary","",flag[$1])
        if(num_pri>0) print $0
    }' |
cat > "${tmp_all}"

################################################################################
#! MIDS conversion
################################################################################

cat "${tmp_all}" |
    grep "primary" |
    mids_compressed |
cat > "${tmp_primary}"

cat "${tmp_all}" |
    grep "secondary" |
    mids_compressed |
cat > "${tmp_secondary}"

#===========================================================
#? Concatenate "primary" and "secondary"
#===========================================================

cat "${tmp_primary}" "${tmp_secondary}" |
    sort -t " " -k 2,2n |
    awk '{seq_[$1]=seq_[$1]" "$2" "$3" "$4}
        END{for(key in seq_) print key,seq_[key]}' |
    awk 'NF==4
    #---------------------------------------
    #* Flox deletion
    #---------------------------------------
    NF==7{id=$1; minloc=$2; maxloc=$2 ; seq=""; len=0; str=""
        for(i=2; i<=NF; i=i+2){
            seq_[$i]=$(i+1)
            if(minloc>$i){ minloc=$i; len=length($(i+1))}
            if(maxloc<$i) maxloc=$i
        }
        len=maxloc-minloc-len
        for(len_=1; len_<=len; len_++) str="D" str
        for(key in seq_) seq=seq seq_[key] str
        print id, minloc, seq
        delete seq_
    }
    #---------------------------------------
    #* Inversion
    #---------------------------------------
    NF==10{id=$1; minloc=$2; maxloc=$2 ; seq=""; len=0; str=""
        for(i=2; i<=NF; i=i+2){
            if(minloc>$i){ minloc=$i; len=length($(i+1))}
            if(maxloc<$i) maxloc=$i
        }
        len=length($5)
        for(len_=1; len_<=len; len_++) str="=" str
        $5=str
        print id, minloc, $3 $5 $7
    }' 2>/dev/null |
    sed "s/D*$//g" |
    # complement seqences to match sequence length (insert "M")
    ## start
    awk '{start=""
        for(i=1; i < $2; i++) start=start"M"
        $3=start""$3}1' |
    ## remove start site info
    awk '{$2=""}1' |
    ## end
    awk -v reflen="${ref_len}" '{
        seqlen=length($2);
        end="";
        if(seqlen>reflen) $2=substr($2,1,reflen)
        else {for(i=seqlen; i < reflen; i++) end=end"M"; $2=$2""end}
    }1' |
    # Remove all-mutated reads
    awk '$2 !~ /^[I|D|S]+$/' |
    sed -e "s/$/\t${barcode}/g" -e "s/ /\t/g" |
cat

rm "${tmp_que_fa}" "${tmp_seqID}" "${tmp_all}" "${tmp_primary}" "${tmp_secondary}"

exit 0