#!/bin/sh
# ==============================================================================
# Because the controls are commonly used, they should be created in advance.
# コントロールは共通して用いるため、先立って作成しておく
# ==============================================================================

# ==============================================================================
# Initialize shell environment
# ==============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ==============================================================================
# I/O naming
# ==============================================================================
# ----------------------------------------
# Input
# ----------------------------------------

barcode="${1}"
alleletype="${2}"

mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

# ----------------------------------------
# Output
# ----------------------------------------
mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため

# temporal --------------------------------
MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${barcode}_${alleletype}"

# result --------------------------------
control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"

# ----------------------------------------------------------
# Get max sequence length
# ----------------------------------------------------------
seq_maxnum=$(
    cat .DAJIN_temp/fasta/fasta.fa |
    grep -v "^>" |
    awk '{if(max<length($0)) max=length($0)}
    END{print max}'
)

# ==============================================================================
# Generate  "${control_score}"
# ==============================================================================

# # ----------------------------------------------------------
# # MIDS conversion
# # ----------------------------------------------------------

find .DAJIN_temp/fasta_ont/ -type f |
    grep "${barcode}" |
    xargs -I @ ./DAJIN/src/mids_convertion.sh @ "${mapping_alleletype}" "control"
cp ".DAJIN_temp/data/MIDS_${barcode}_${mapping_alleletype}" "${MIDS_ref}"

# rm .DAJIN_temp/tmp_${barcode}*${alleletype}
# rm .DAJIN_temp/tmp_${barcode}*${alleletype}

# ----------------------------------------------------------
# Mutation scoring of control samples
# ----------------------------------------------------------
# ----------------------------------------
# 挿入塩基を1つの挿入塩基数にまとめて配列のズレを無くす
# ----------------------------------------
cat "${MIDS_ref}" |
    grep "${barcode}" |
    sort -k 1,1 |
    join -1 1 -2 2 - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v alelle="$alleletype" '$NF==alelle' |
    cut -d " " -f 2 |
    # Insertion annotation
    awk -F "" '{
        for(i=1; i<=NF; i++){
            if($i=="I") num=num+1
            if($i=="I" && $(i+1)!="I") {
                # -----------------------------------
                ### e.g) if num=10, num becomes "a"
                # -----------------------------------
                if(num>=10 && num<=35) {num=sprintf("%c", num+87)}
                else if(num>=36) num="z"
                ###
                $(i+1)=num; num=0}
            }
        print $0
        }' |
    sed -e "s/I//g" -e "s/ //g" |
    # ----------------------------------------
    # MIDS変換で末尾がDになった配列を=に変換する
    # ----------------------------------------
    sed -e "s/I//g" -e "s/ //g" |
    sed "s/\(D*$\)/ \1/g" |
    awk '{
        for(i=1; i<=NF; i++) if($i~/^D*$/) gsub(/./, "=", $i)
    }1' |
    sed "s/ //g" |
    # ----------------------------------------
    # 短い配列をPaddingする
    # ----------------------------------------
    awk -v seqnum="${seq_maxnum}" \
        'BEGIN{OFS=""}
        { if(length($0) < seqnum){
            seq="="
            for(i=length($0)+1; i<=seqnum; i++) $i=seq
            print $0}
        }' |
    # ----------------------------------------
    # 行を「リード指向」から「塩基部位指向」に変換する
    # 例：
    # MMM
    # MII
    # ↓
    # MM
    # MI
    # MI
    # ----------------------------------------
    awk -F "" \
        '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } }
        END {    
            for(j=1; j<=NF; j++) {
                str=a[1,j]
                for(i=2; i<=NR; i++){ str=str""a[i,j] }
                print str }
        }' |
    # ----------------------------------------
    # シークエンスエラーを描出する
    # ----------------------------------------
    awk -F "" '{
        sum[1]=gsub("=","=",$0)
        sum[2]=gsub("M","M",$0)
        sum[3]=gsub(/[1-9]|[a-z]/,"@",$0)
        sum[4]=gsub("D","D",$0)
        sum[5]=gsub("S","S",$0)
        # ----------------------------------------
        ### Controlにおいて系統的な変異が10%を超える部位をシークエンスエラーとする
        # ----------------------------------------
        per=10
        if(sum[3] > NF*per/100) num = 2
        else if(sum[4] > NF*per/100) num = 2
        else if(sum[5] > NF*per/100) num = 2
        else num=1
        #
        # print NR, "@", sum[1], sum[2], sum[3], sum[4], sum[5], "@", \
        #     (sum[1]+sum[2])/NF, sum[3]/NF,sum[4]/NF,sum[5]/NF, num
        print num
    }' |
cat - > "${control_score}"

