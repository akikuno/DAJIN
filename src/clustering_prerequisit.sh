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

seq_length=$(
    cat .DAJIN_temp/fasta/wt.fa |
    grep -v "^>" |
    awk '{print length($0)}'
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
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
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
    awk -v seqnum="${seq_length}" \
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
cat > "${control_score}"


# ==============================================================================
# Generate  "${control_score}"
# ==============================================================================


cat "${control_score}" |
    awk '{print $1,NR}' |
    sort -t " " -k 2,2 |
cat > .DAJIN_temp/clustering/temp/tmp_control_score

find .DAJIN_temp/fasta_conv/* |
    grep -v wt.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
    grep flox_deletion |
while read -r label; do
    minimap2 -ax map-ont .DAJIN_temp/fasta_conv/wt.fa .DAJIN_temp/fasta_conv/"${label}".fa --cs=long 2>/dev/null |
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
    join -a 1 -1 3 -2 2 - .DAJIN_temp/clustering/temp/tmp_control_score |
    sort -k 3,3n |
    awk 'NF==3{$4=1}{print $NF}' |
    cat > ".DAJIN_temp/clustering/temp/control_score_${label}"
done

rm .DAJIN_temp/clustering/temp/tmp_control_score
# find .DAJIN_temp/fasta_conv/* |
# grep -v wt.fa |
# sed "s:.*/::g" |
# sed "s/.fa.*$//g" |
# xargs -I @ \
#     minimap2 -ax map-ont .DAJIN_temp/fasta_conv/wt.fa .DAJIN_temp/fasta_conv/@.fa --cs=long 2>/dev/null |
#     awk '$1 !~ /^@/ {
#         id=$1
#         cstag=$(NF-1)
#         gsub("cs:Z:","",cstag)
#         gsub(/*[acgt]/,"",cstag)
#         gsub(/[+=-]/,"",cstag)
#         print id, cstag}' |
#     awk -F "" '{
#         refnr=0
#         for(i=1; i<=NF; i++){
#             if($i ~ /[A|C|G|T]/){$i="ref"; refnr++; print $1, $i, refnr, i}
#             else {$i="mut"; print $1, $i, "empty", i}
#     }}' |
#     sort -t " " -k 2,2 |
#     join -a 1 -1 2 -2 2 - test2 |
#     sort -k 3,3n |
#     awk 'NF==3{$4=1}{print $NF}' |
#     wc -l

# while read -r label; do
#     minimap2 -ax map-ont .DAJIN_temp/fasta_conv/wt.fa .DAJIN_temp/fasta_conv/"${label}".fa --cs=long 2>/dev/null |
#     awk '$1 !~ /^@/ {print $(NF-1)}' |
#     sed "s/cs:Z://g" |
#     sed "s/*[acgt]//g" |
#     sed "s/[=+-]//g" |
#     awk -F "" '{
#         refnr=0
#         for(i=1; i<=NF; i++){
#             if($i ~ /[A|C|G|T]/){$i="ref"; refnr++; print $i, refnr, i}
#             else {$i="mut"; print $i, "empty", i}
#     }}' |
#     sort -t " " -k 2,2 |
#     join -a 1 -1 2 -2 2 - test2 |
#     sort -k 3,3n |
#     awk 'NF==3{$4=1}{print $NF}' |
#     wc -l
# done
