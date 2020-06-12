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
#? Auguments
#===========================================================
# barcode=barcode26
# alleletype=wt

barcode="${1}"
alleletype="${2}"

mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

#===========================================================
#? Input
#===========================================================

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため
control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"

#===========================================================
#? Temporal
#===========================================================
MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${barcode}_${alleletype}"
MIDS_tmp=".DAJIN_temp/clustering/temp/MIDS_tmp_${barcode}_${alleletype}"
control_tmp=".DAJIN_temp/clustering/temp/control_tmp_${barcode}_${alleletype}"

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
# 行を「リード指向」から「塩基部位指向」に変換する
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
    awk -F "" '{
        INS=gsub(/[1-9]|[a-z]/,"@",$0)
        DEL=gsub("D","D",$0)
        SUB=gsub("S","S",$0)
        # ----------------------------------------
        ### Controlにおいて系統的な変異が10%を超える部位をシークエンスエラーとする
        # ----------------------------------------
        per=10
        # print INS, DEL, SUB, NF*per/100 
        if(INS > NF*per/100 || DEL > NF*per/100 || SUB > NF*per/100) num = 2
        else num=1
        print num
    }' |
cat > "${control_score}"


# ==============================================================================
# Generate  "${control_score}"
# ==============================================================================


cat "${control_score}" |
    awk '{print $1,NR}' |
    sort -t " " -k 2,2 |
cat > "${control_tmp}"

find .DAJIN_temp/fasta/* |
    grep -v wt.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
while read -r label; do
    minimap2 -ax map-ont .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/"${label}".fa --cs=long 2>/dev/null |
    awk '$1 !~ /^@/ {print $(NF-1)}' |
    sed "s/cs:Z://g" |
    sed "s/*[acgt]//g" |
    sed "s/[=+]//g" |
    awk -F "" '{
        refnr=0
        if($0 !~ "-"){
            for(i=1; i<=NF; i++){
                if($i ~ /[A|C|G|T]/){refnr++; $i="ref"; print $i, refnr, i}
                else {$i="mut"; print $i, "empty", i}
                }}
        else{
            gsub("-","",$0)
            for(i=1; i<=NF; i++){
                if($i ~ /[A|C|G|T]/){refnr++; $i="ref"; print $i, refnr, i}
            }
        }}' |
    sort -t " " -k 3,3 |
    join -a 1 -1 3 -2 2 -  "${control_tmp}" |
    sort -k 3,3n |
    awk 'NF==3{$4=1}{print $NF}' |
    cat > ".DAJIN_temp/clustering/temp/control_score_${label}"
done

rm "${MIDS_ref}" "${MIDS_tmp}" "${control_tmp}"

exit 0