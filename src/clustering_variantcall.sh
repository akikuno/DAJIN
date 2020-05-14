#!/bin/sh

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
# barcode="barcode08"
# alleletype="normal"
# alleletype_original=${alleletype}
# suffix="${barcode}"_"${alleletype}"
# echo $suffix
# [ "$alleletype" = "normal" ] && alleletype="wt"
# [ "$alleletype" = "abnormal" ] && alleletype="wt"

barcode="${1}"
alleletype="${3}"
alleletype_original=${3}
suffix="${barcode}"_"${alleletype}"
[ "$alleletype" = "normal" ] && alleletype="wt"
[ "$alleletype" = "abnormal" ] && alleletype="wt"

mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため

# ----------------------------------------
# Temporal Output
# ----------------------------------------

# Output Plot
hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"
output_plot=".DAJIN_temp/clustering/temp/plot_${suffix}"

# allele percentage on each cluster
output_alleleper=".DAJIN_temp/clustering/temp/allele_percentage_${suffix}".txt

# ----------------------------------------------------------
# Output results
# ----------------------------------------------------------
output_id=".DAJIN_temp/clustering/result_allele_id_${suffix}".txt
output_result=".DAJIN_temp/clustering/result_allele_mutinfo_${suffix}".txt

# ============================================================================
# Report allele mutation info
# ============================================================================

cat "${output_plot}" |
    awk '{
        num=1
        cl_mut[$3]=cl_mut[$3]$2
        if($2!="M"){
            loc[NR] = $1
            mut[NR] = $2
            cl[NR] = $3
            ins[NR] = $4
        }}
    END{
        # ----------------------------------------------------------
        # もし変異がなければintactと表示する
        # ----------------------------------------------------------
        for(j in cl_mut) {
            if (cl_mut[j] !~/[I|D|S]/) {print j, 0,"intact"}
        }
        # ----------------------------------------------------------
        # 同じ変異が5つ飛ばし以内で続いている場合は連続した変異とみなす
        # また、挿入塩基の場合は挿入塩基数を直接表記する
        # ----------------------------------------------------------
        for(i in loc){
            if(loc[i+1] - loc[i] == 1 || \
                loc[i+2] - loc[i] == 2 || \
                loc[i+3] - loc[i] == 3 || \
                loc[i+4] - loc[i] == 4 || \
                loc[i+5] - loc[i] == 5) {num++}
            # if( for(j=1; j<=5; j++){loc[j+1]-log[j] == 1}) num++
            else if(ins[i]>0) {print cl[i], i, ins[i]""mut[i], loc[i]}
            else {print cl[i], i, num""mut[i], loc[i]-num; num=1}
    }}' |
    # ----------------------------------------------------------
    # 挿入塩基数が10以上の場合に数値情報に逆変換する
    # ----------------------------------------------------------
    awk '{
        if($3 ~ /[a-z]I/) {
            for (i=10; i<=36; i++) {
                num=i+87
                ins=sprintf("%c", num)
                if($3==ins"I") $3=i"I"
            }
        }
        print $0}' |
    sed "s/35I/>35I/g" |
    sort -t " " -k 1,1 -k 2,2n |
cat - > ".DAJIN_temp/clustering/temp/tmp_${suffix}"

# ----------------------------------------------------------
# barcode, alleletype, クラスター番号, アレル頻度、変異、変異部位、ソート番号を出力する
# ----------------------------------------------------------
cat "$output_alleleper" |
    cut -d " " -f 2- |
    sort |
    join - ".DAJIN_temp/clustering/temp/tmp_${suffix}" |
    sort -t " " -k 1,1n  -k 3,3n |
    awk '{print $1,$2,$4,$5, NR}' |
    sed "s/^/${barcode} ${alleletype_original} /g" |
cat - > "${output_result}"


# ============================================================================
# Report allele mutation info
# 各リードとクラスターの対応付を行う
#（次のVCF作製とSequence logo描出のために必要）
# ============================================================================

true > .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
cat "${output_alleleper}" |
while read -r input; do
    before=$(echo "$input" | cut -d " " -f 1)
    after=$(echo "$input" | cut -d " " -f 2)
    #
    cat "${hdbscan_id}" |
    awk -v bf="${before}" -v af="${after}" \
    '$2==bf {$2=af; print}' \
    >> .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
done

cat .DAJIN_temp/clustering/temp/tmp_id_"${suffix}" |
    sed "s/ /\t/g" |
    sort |
cat - > "${output_id}"

# ============================================================================
# 変異情報の同定
# Variant call
# ============================================================================

# start=$(cat .DAJIN_temp/data/gggenome_location | cut -f 2)
cluster=2
cat $output_result |
awk -v cl="${cluster}" '$3==cl' |
while read -r input
do
mutation_type=$(echo $input | cut -d " " -f 5)
mutation_site=$(echo $input | cut -d " " -f 6)
echo $mutation_type $mutation_site "${cluster}"
done

cat "${output_id}" |
    grep "${cluster}$" |
    cut -f 1 |
    sort -u |
cat - > .DAJIN_temp/clustering/temp/tmp_sorted_id_"${suffix}"_"${cluster}"

ref=.DAJIN_temp/fasta/wt.fa
cat .DAJIN_temp/fasta_ont/"${barcode}".fa |
minimap2 -ax map-ont "${ref}" - --cs=long 2>/dev/null |
    sort |
    join - .DAJIN_temp/clustering/temp/tmp_sorted_id_"${suffix}"_"${cluster}" |
    awk '{print $4, $(NF-1)}' |
cat - > tmp_test

cat tmp_test |
sed "s/cs:Z://g" |
awk '{seq=""
    padding=$1-1
    for(i=1; i<=padding; i++) seq=seq"-"
    print seq""$2}' |
sed "s/*[a-z]/ /g" |
sed "s/[-|=]/ /g" |
sed "s/+/ +/g" |
awk -v mutsite="${mutation_site}" \
'{len=0
for(i=1; i<=NF;i++) {
    if(len >= mutsite-1) {print len, $i; break}
    else {
        if($i !~ /+/) len=len+length($i)}
    }
}' |
sort |
uniq -c |
awk '{if(max<$1){max=$1; loc=$2; max_mut=$3}}
    END{print loc, max_mut}'

# cat "$ref" |
# sed 1d |
# awk '{mut=substr($0, 739,1)
#     print mut}'



# ref:    TTGCCAGATA
# target: TTGCCAGATTA
# sample: TTAGCCGATTG
# 8塩基のAがATになる

mutation_sites=20
echo "=TT+attt=GC+acccc=C-a=GA+tt=T*ag" |
sed "s/*[a-z]/ /g" |
sed "s/[-|=]/ /g" |
sed "s/+/ +/g" |
awk -v mutsite="${mutation_sites}" \
'{len=0
for(i=1; i<=NF;i++) {
    if(len >= mutsite) {print len, $i; break}
    else {
        if($i !~ /+/) len=len+length($i)}
    }
}'



ref=.DAJIN_temp/fasta/wt.fa
que=test.fa
minimap2 -ax map-ont "${ref}" "${que}" --cs=long 2>/dev/null |
awk '{print $4, $(NF-1)}' |
cat - > tmp_test
cat tmp_test |
sed "s/cs:Z://g" |
awk '{seq=""
    padding=$1-1
    for(i=1; i<=padding; i++) seq=seq"-"
    print seq""$2}' |
sed "s/*[a-z]/ /g" |
sed "s/[-|=]/ /g" |
sed "s/+/ +/g" |
awk -v mutsite="${mutation_sites}" \
'{len=0
for(i=1; i<=NF;i++) {
    if(len >= mutsite-1) {print len, $i; break}
    else {
        if($i !~ /+/) len=len+length($i)}
    }
}'


head -n 1000 .DAJIN_temp/fasta_ont/target_simulated_aligned_reads.fasta |
minimap2 -ax map-ont "${ref}" - --cs=long 2>/dev/null |
awk '{print $4, $(NF-1)}' |
cat - > tmp_test



