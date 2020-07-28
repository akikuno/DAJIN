#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================

# barcode=barcode34
# alleletype=flox_deletion
# threads=12
# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

#===========================================================
#? Auguments
#===========================================================

barcode="${1}"
alleletype="${2}"
threads="${3}"
mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

#===========================================================
#? Input
#===========================================================

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/"
control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"

#===========================================================
#? Temporal
#===========================================================
MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${barcode}_${alleletype}"
MIDS_tmp=".DAJIN_temp/clustering/temp/MIDS_tmp_${barcode}_${alleletype}"
control_tmp=".DAJIN_temp/clustering/temp/control_tmp_${barcode}_${alleletype}"
label_nuc_tmp=".DAJIN_temp/clustering/temp/label_nuc_tmp_${barcode}_${alleletype}"

################################################################################
#! Control scoring
################################################################################

#===========================================================
#? MIDS conversion
#===========================================================

./DAJIN/src/mids_clustering.sh "${barcode}" "${mapping_alleletype}" > "${MIDS_ref}"

#===========================================================
#? Mutation scoring of control samples
#===========================================================

cat "${MIDS_ref}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v alelle="$alleletype" '$NF==alelle' |
    cut -d " " -f 2 |
cat > "${MIDS_tmp}"

# ----------------------------------------
# Transpose matrix
# ----------------------------------------
nr=$(cat "${MIDS_tmp}" | wc -l)

cat "${MIDS_tmp}" |
    awk -v nr="${nr}" -F "" \
    '{for(i=1; i<=NF; i++){
            row[i]=row[i] $i
            if(NR==nr) print row[i]
        }
    }' |
    # cat tmp | head -n 3113 | tail -n 1 |
    # head -n 3208 tmp | tail |
    awk -F "" '{
        INS=gsub(/[1-9]|[a-z]/,"@",$0)
        DEL=gsub("D","D",$0)
        SUB=gsub("S","S",$0)
        # ----------------------------------------
        #* Define sequence error when control sample has more than 5% mutations
        # ----------------------------------------
        per=5
        # print INS, DEL, SUB, NF*per/100 #!<<<<
        if(INS > NF*per/100 || DEL > NF*per/100 || SUB > NF*per/100) num = 2
        else num=1
        print num
    }' |
cat > "${control_score}"


################################################################################
#! Generate scores of other possible alleles
################################################################################

cat "${control_score}" |
    awk '{print $1,NR}' |
    sort -t " " -k 2,2 |
cat > "${control_tmp}"

find .DAJIN_temp/fasta/* |
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
    cat > "${label_nuc_tmp}"
    #
    if [ "$(cat "${label_nuc_tmp}" | wc -l)" -eq 3 ]; then # inversion
        cat "${label_nuc_tmp}" |
            awk 'NR==2{printf tolower($1); next} {printf $1}'
    elif [ "$(cat "${label_nuc_tmp}" | wc -l)" -eq 2 ]; then # flox deletion
        wt_length=$(sed 1d .DAJIN_temp/fasta/wt.fa | awk '{print length}')
        que_length=$(sed 1d .DAJIN_temp/fasta/"${label}".fa | awk '{print length}')
        del_length=$((${wt_length}-${que_length}+1))

        cat "${label_nuc_tmp}" |
            awk -v del_len="${del_length}" 'NR==2{
                seq=""
                for(i=1;i<del_len;i++){seq=seq "a"}
                $0=substr($0,2)
                printf seq, $0}{printf $0}'
    else
        cat "${label_nuc_tmp}"
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
    join -a 1 -1 2 -2 2 - "${control_tmp}" |
    awk 'NF==3{$4=1}{print $0}' |
    sort -k 3,3n |
    if [ "$(cat "${label_nuc_tmp}" | wc -l)" -eq 2 ]; then # flox deletion
        cat - | grep -v "^NA"
    else
        cat -
    fi |
    cat > ".DAJIN_temp/clustering/temp/control_score_${label}"
done

rm "${MIDS_ref}" "${MIDS_tmp}" "${control_tmp}" "${label_nuc_tmp}"

exit 0