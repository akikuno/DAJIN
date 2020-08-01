#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -eu
umask 0022
export LC_ALL=C
export UNIX_STD=2003  # to make HP-UX conform to POSIX

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# control=barcode41
# threads=12

#===========================================================
#? Auguments
#===========================================================

control="${1}"
threads="${2}"

#===========================================================
#? Input
#===========================================================

alleletype="wt"

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/"
control_score=".DAJIN_temp/clustering/temp/control_score_${alleletype}"

#===========================================================
#? Temporal
#===========================================================
MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${control}_${alleletype}"
tmp_MIDS=".DAJIN_temp/clustering/temp/tmp_MIDS_${control}_${alleletype}"
tmp_control=".DAJIN_temp/clustering/temp/tmp_control_${control}_${alleletype}"
tmp_label_mapping=".DAJIN_temp/clustering/temp/tmp_label_mapping_${control}_${alleletype}"

################################################################################
#! Control scoring
################################################################################

#===========================================================
#? MIDS conversion
#===========================================================

./DAJIN/src/mids_clustering.sh "${control}" "${alleletype}" > "${MIDS_ref}"

#===========================================================
#? Mutation scoring of control samples
#===========================================================

cat "${MIDS_ref}" |
    grep "${control}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v alelle="$alleletype" '$NF==alelle' |
    cut -d " " -f 2 |
cat > "${tmp_MIDS}"

# ----------------------------------------
# Transpose matrix
# ----------------------------------------
nr=$(cat "${tmp_MIDS}" | wc -l)

cat "${tmp_MIDS}" |
    awk -v nr="${nr}" -F "" \
    '{for(i=1; i<=NF; i++){
            row[i]=row[i] $i
            if(NR==nr) print row[i]
        }
    }' |
    awk -F "" '{
        INS=gsub(/[1-9]|[a-z]/,"@",$0)
        DEL=gsub("D","D",$0)
        SUB=gsub("S","S",$0)
        # ----------------------------------------
        #* Define sequence error when control sample has more than 5% mutations
        # ----------------------------------------
        per=5
        if(INS > NF*per/100 || DEL > NF*per/100 || SUB > NF*per/100)
            num=2
        else
            num=1

        print num
    }' |
cat > "${control_score}"


################################################################################
#! Generate scores of other possible alleles
################################################################################

cat ".DAJIN_temp/fasta/${alleletype}.fa" |
    sed 1d |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    paste - "${control_score}" |
    awk '{print NR"_"$1, $2}' |
    sort -t " " -k 1,1 |
cat > "${tmp_control}"

find .DAJIN_temp/fasta/ -type f |
    grep -v wt.fa |
    grep -v fasta.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
while read -r label; do
    minimap2 -t "${threads}" -ax splice \
        .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/"${label}".fa \
        --cs=long 2>/dev/null |
        grep -v "^@" |
        awk '{print $4, $(NF-1)}' |
        sed "s/cs:Z://g" |
        sed "s/[=+]//g" |
        sed "s/\~[actg][actg]//g" |
        sed "s/\([0-9][0-9]*\)[actg][actg]/\n\1 /g" |
        sed "s/-[acgt][acgt]*//g" |
        sort -t " " -k 1,1n |
        awk '{
            ref_position=1
            start=$1
            split($2, seq_array, "")
            for(i=1; i<=length(seq_array); i++){
                que_position=i+start-1
                if(seq_array[i] ~ /[ACGT]/) {
                    loc=start+ref_position-1
                    print loc"_"seq_array[i], que_position
                    ref_position++}
                else{
                    loc=start+i-1
                    print "mut_"seq_array[i], que_position}
            }
        }' |
        sort |
        join -a 1 - "${tmp_control}" |
        awk 'NF==2 {print $0,1; next}1' |
        sort -t " " -k 2,2n |
        awk '{print $NF}' |
    cat > ".DAJIN_temp/clustering/temp/control_score_${label}"
done

rm "${MIDS_ref}" "${tmp_MIDS}" "${tmp_control}" "${tmp_label_mapping}"

exit 0