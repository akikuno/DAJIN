#!/bin/sh

# ============================================================================
# I/O and Arguments
# ============================================================================
# mkdir -p .DAJIN_temp/seqlogo/
# gRNAのデータを保存する
# grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
# grna=$(echo "${grna}" |
#     DAJIN/src/revcomp.sh - |
#     sed "s/^/${grna},/g")

barcode=barcode14
alleletype=target
suffix="${barcode}_${alleletype}"

fasta=.DAJIN_temp/fasta_ont/"${barcode}".fa
clustering_id=$(find .DAJIN_temp/clustering/result_allele_id* | grep "$suffix")

mkdir -p .DAJIN_temp/seqlogo/temp/
tmp_fa=".DAJIN_temp/seqlogo/temp/${suffix}.fa"
tmp_targetseq=".DAJIN_temp/seqlogo/temp/targetseq_${suffix}"
tmp_targetseq_fa=".DAJIN_temp/seqlogo/temp/targetseq_${suffix}.fa"
tmp_targetseq_trimmed_fa=".DAJIN_temp/seqlogo/temp/targetseq_trimmed_${suffix}.fa"
tmp_mappedseq=".DAJIN_temp/seqlogo/temp/mapped_seq_${suffix}"
tmp_lalign=".DAJIN_temp/seqlogo/temp/lalign_${suffix}"
tmp_lalignseq=".DAJIN_temp/seqlogo/temp/lalignseq_${suffix}"

# ============================================================================
# FASTAリードにクラスタ番号をラベルづけする
# ============================================================================

cat "${fasta}" |
    awk '{if(NR % 2 != 0) $1="HOGE"$1
        else $1="FUGA"$1
        print}' |
    tr -d "\n" |
    sed -e "s/HOGE/\n/g" -e "s/FUGA/ /g" |
    awk '{print $1,$NF}' |
    sed "s/^[@|>]//g" |
    grep -v "^$" |
    sort -t " " |
    join - "${clustering_id}" |
    awk '{print ">"$NF"@@@"$1"\n"$2}' |
cat - > ${tmp_fa}

# ============================================================================
# Align reads to target sequence
# ============================================================================

mutation_type=$(
    minimap2 -ax map-ont \
    .DAJIN_temp/fasta/target.fa \
    .DAJIN_temp/fasta/wt.fa 2>/dev/null |
    grep -v "^@" |
    cut -f 6 |
    awk '{if($0~"I") print "D"
        else if($0~"D") print "I"
        else if($0~"S") print "P"
        }'
)
if [ "$mutation_type" = "D" ]; then
    # ------------------------------------------------------------
    # 2cut deletionの場合：
    # 結合部から±50塩基を抽出する
    # ------------------------------------------------------------
    minimap2 -ax map-ont \
            .DAJIN_temp/fasta_conv/wt.fa \
            .DAJIN_temp/fasta_conv/target.fa --cs=long 2>/dev/null |
        grep -v "^@" |
        awk '{print $(NF-1)}' |
        sed "s/[a-z]*=//g" |
        awk '{
            match($0, "-")
            print toupper(substr($0, RSTART-50, 101))
        }' |
        sed -e "s/^/>target /g" -e "s/-//g" |
    cat - > "$tmp_targetseq"
elif [ "$mutation_type" = "I" ]; then
    # ------------------------------------------------------------
    # Knock-inの場合：
    # Knock-in配列全長を抽出する
    # ------------------------------------------------------------
    minimap2 -ax map-ont \
            .DAJIN_temp/fasta_conv/wt.fa \
            .DAJIN_temp/fasta_conv/target.fa \
            --cs=long 2>/dev/null |
        grep -v "^@" |
        awk '{print $(NF-1), $10}' |
        awk '{seq=$2
            gsub("[A|C|G|T]", "",$1)
            gsub("cs:Z:=", "", $1)
            gsub("=", "", $1)
            split($1, array, "+")
            for(key in array) print toupper(array[key]), seq}' |
        sed 1d |
        awk '{
            seq_center=int(length($1)/2)
            match($2, $1)
        print toupper(substr($2, RSTART+seq_center-10, 20))
        }' |        
        # sed "s///g" |
        # awk '{; print}' |
        # sed "s/-/\n/g" |
        # sed "s/=//g" |
        # grep -v "^$" |
        awk '{print ">target"NR,$0}' |
    cat - > "$tmp_targetseq"
fi

# ------------------------------------------------------------
# クラスタごとのリードをターゲット配列にalignment
# ------------------------------------------------------------

seqlen=$(awk '!/[>|@]/ {print length($0)}' .DAJIN_temp/fasta_conv/target.fa)

true > "$tmp_mappedseq"
input=$(head -n1 "$tmp_targetseq")
cat "${tmp_targetseq}" |
while read -r input; do
    echo "$input" |
        # -----------------------------
        # 足りない配列をNで補う
        # -----------------------------
        awk -v seqlen="${seqlen}" \
        '{  s=sprintf("%."seqlen"d","0")
            gsub("0","N",s)
            s=substr(s,$2,seqlen)
            print $1,$2s}' |
        sed "s/ /\n/g" |
        # -----------------------------
        # Nで保管したTargetseqにmappingする
        # -----------------------------
        minimap2 -ax map-ont - \
            ${tmp_fa} --cs=long 2>/dev/null |
        awk '$2==0 || $2==16 {print $1, $2, $3, $(NF-1)}' |
        awk '{cstag=$NF
            sub("cs:Z:", "", cstag)
            gsub(/\*[a-z]/, "", cstag)
            gsub(/\-[a-z]*[+|*|=]/, "", cstag)
            gsub("[=|+]", "", cstag)
            #
            gsub(/@@@.*/, "", $1)
            print $1, $2, $3, toupper(cstag)
        }' |
    cat - >> "$tmp_mappedseq"
done

# ============================================================================
# Sequence logo
# ============================================================================

input=$(head -n1 "$tmp_mappedseq" |
cut -d " " -f 1,3 | sort -u)

cat "$tmp_mappedseq" |
cut -d " " -f 1,3 | sort -u |
while read -r input; do
    cl=$(echo "$input" | cut -d " " -f 1)
    target=$(echo "$input" | cut -d " " -f 2)
    label="${cl}"@"${target}"
    # ----------------------------------------
    # Split fasta files for the following alignment
    # ----------------------------------------
    tmp_split_dir=".DAJIN_temp/seqlogo/temp/split_${suffix}"
    # ----------------------------------------
    rm -rf "${tmp_split_dir}" 2>/dev/null
    mkdir -p "${tmp_split_dir}"
    # ----------------------------------------
    # マップしたリードを細分化します
    # ----------------------------------------
    cat "$tmp_mappedseq" |
        awk -v cl="${cl}" -v target="${target}"\
        '$1==cl && $3==target' |
        awk '{print ">"$1"@"$3"\n"$NF}' |
    split -l 2 - ${tmp_split_dir}"/split_"
    # ----------------------------------------
    # lalignでローカルアライメントします
    # ----------------------------------------
    cat "${tmp_targetseq}" |
        grep "${target}" - |
        sed "s/ /\n/g" |
    cat - > "$tmp_targetseq_fa"

    find "${tmp_split_dir}/" -name "split_*" -type f |
        xargs -I @ ./DAJIN/src/test_intact_lalign.sh \
            "$tmp_targetseq_fa" @ |
        awk -v label="${label}" '{$2=label; print}' |
    cat - > "$tmp_lalign"
    # ----------------------------------------
    # 低クオリティのリードを除きます
    # ----------------------------------------
    percentile=0.5
    per=$(
        cat "$tmp_lalign" |
        awk -v per=${percentile} '{
        percentile[NR]=$1}
        END{asort(percentile)
        print percentile[int(NR*per)]
    }')
    #
    # ----------------------------------------
    # 長さを同じにします
    # ----------------------------------------
    mut_length=$(awk '!/[>|@]/ {print length($0)}' "$tmp_targetseq_fa")
    #
    cat "$tmp_lalign" |
        awk -v per=${per} '$1>=per' |
        cut -d " " -f 3 |
        # 最初と最後の5塩基（合計10塩基）を取り除きます
        awk -v mut_len=${mut_length} '{print substr($0,5,mut_len-10)}' |
    cat - > ${tmp_lalignseq}
    #
    # 最初と最後の5塩基（合計10塩基）を取り除きます
    cat "$tmp_targetseq_fa" |
        awk '{if($1!~/^>/) $0=substr($0,5, length($0)-10)
            print}' |
    cat - > "$tmp_targetseq_trimmed_fa"
    # ----------------------------------------
    # 配列ロゴ
    # ----------------------------------------
    # echo $tmp_targetseq_fa $tmp_lalignseq
    python DAJIN/src/test_logomaker.py "$tmp_targetseq_trimmed_fa" "$tmp_lalignseq"
done

exit 0