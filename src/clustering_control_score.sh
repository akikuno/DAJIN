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

cat "${control_score}" |
    awk '{print $1,NR}' |
    sort -t " " -k 2,2 |
cat > "${tmp_control}"

find .DAJIN_temp/fasta/ -type f |
    grep -v wt.fa |
    grep -v fasta.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
while read -r label; do
    minimap2 -t "${threads}" -ax map-ont \
        .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/"${label}".fa \
        --cs=long 2>/dev/null |
        grep -v "^@" |
        awk '{print $4, $(NF-1)}' |
        sed "s/cs:Z://g" |
        sed "s/=//g" |
        sort -t " " -k 1,1n |
        cut -d " " -f 2 |
    cat > "${tmp_label_mapping}"
    #
    nrow_tmp_label_mapping=$(cat "${tmp_label_mapping}" | wc -l)
    #
    if [ "${nrow_tmp_label_mapping}" -eq 3 ]; then # inversion
            awk 'NR==2{printf tolower($1); next} {printf $1}' "${tmp_label_mapping}"
    elif [ "${nrow_tmp_label_mapping}" -eq 2 ]; then # flox deletion
        wt_length=$(sed 1d .DAJIN_temp/fasta/wt.fa | awk '{print length}')
        que_length=$(sed 1d .DAJIN_temp/fasta/"${label}".fa | awk '{print length}')
        del_length=$((${wt_length}-${que_length}+1))

        cat "${tmp_label_mapping}" |
            awk -v del_len="${del_length}" 'NR==2{
                seq=""
                for(i=1;i<del_len;i++){seq=seq "a"}
                $0=substr($0,2)
                printf seq, $0}{printf $0}'
    else
        cat "${tmp_label_mapping}"
    fi |
    sed "s/*[acgt]//g" |
    sed "s/[=+]//g" |
    awk -F "" '{
        refnr=0
        if($0 !~ "-"){
            for(i=1; i<=NF; i++){
                if($i ~ /[A|C|G|T]/){
                    refnr++
                    $i="ref"
                    print $i, refnr, i
                }
                else {
                $i="mut"
                print $i, "NA", i
                }
            }
        }
        else{
            gsub("-","",$0)
            for(i=1; i<=NF; i++){
                if($i ~ /[A|C|G|T]/){refnr++; $i="ref"; print $i, refnr, i}
            }
        }}' |
    sort -t " " -k 2,2 |
    join -a 1 -1 2 -2 2 - "${tmp_control}" |
    awk 'NF==3{$4=1}{print $0}' |
    sort -k 3,3n |
    if [ "${nrow_tmp_label_mapping}" -eq 2 ]; then # flox deletion
        cat - | grep -v "^NA"
    else
        cat -
    fi |
    awk '{print $NF}' > ".DAJIN_temp/clustering/temp/control_score_${label}"
done

rm "${MIDS_ref}" "${tmp_MIDS}" "${tmp_control}" "${tmp_label_mapping}"

exit 0