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
control="${2}"
alleletype="${3}"
alleletype_original="${3}"
suffix="${barcode}"_"${alleletype}"
[ "$alleletype" = "normal" ] && alleletype="wt"
[ "$alleletype" = "abnormal" ] && alleletype="wt"

mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため

# ----------------------------------------
# Temporal Output
# ----------------------------------------

# Output Plot
hdbscan_id=".DAJIN_temp/clustering/temp/hdbscan_${suffix}"
plot_mutation=".DAJIN_temp/clustering/temp/plot_${suffix}"

# allele percentage on each cluster
allepe_percentage=".DAJIN_temp/clustering/temp/allele_percentage_${suffix}".txt

# ----------------------------------------------------------
# Output results
# ----------------------------------------------------------
output_id=".DAJIN_temp/clustering/result_allele_id_${suffix}".txt
output_result=".DAJIN_temp/clustering/result_allele_mutinfo_${suffix}".txt

# ============================================================================
# Report allele mutation info
# ============================================================================

cat "${plot_mutation}" |
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
cat "$allepe_percentage" |
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

# true > .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
# cat "${allepe_percentage}" |
# while read -r input
# do
#     before=$(echo "$input" | cut -d " " -f 1)
#     after=$(echo "$input" | cut -d " " -f 2)
#     #
#     cat "${hdbscan_id}" |
#     awk -v bf="${before}" -v af="${after}" \
#     '$2==bf {$2=af; print}' \
#     >> .DAJIN_temp/clustering/temp/tmp_id_"${suffix}"
# done

# cat .DAJIN_temp/clustering/temp/tmp_id_"${suffix}" |
#     sed "s/ /\t/g" |
#     sort |
# cat - > "${output_id}"


before=$(cat "${allepe_percentage}" | cut -d " " -f 1 | xargs echo)
after=$(cat "${allepe_percentage}" | cut -d " " -f 2 | xargs echo)

cat "${hdbscan_id}" |
    awk -v bf="${before}" -v af="${after}" \
    'BEGIN{
        split(bf,bf_," ")
        split(af,af_," ")}
    {for(i in bf_){
        if($2==bf_[i]) $2=af_[i]
        }
    }1' |
    sed "s/ /\t/g" |
    sort |
cat - > "${output_id}"

# ============================================================================
# 変異情報の同定
# Variant call
# ============================================================================
ref=".DAJIN_temp/fasta/wt.fa"
cat .DAJIN_temp/fasta_ont/"${barcode}".fa |
    minimap2 -ax map-ont "${ref}" - --cs=long 2>/dev/null |
    sort |
cat - > tmp_sam

set $(cat $output_result |
    awk -v cl="${cluster}" \
    '$3==cl {
    gsub("[0-9]*","", $5)
    type=type$5"_"
    site=site$6"_"}
    END{print type, site}')
mutation_type=$(echo "$1" | sed "s/_/ /g") 
mutation_site=$(echo "$2" | sed "s/_/ /g")
echo $mutation_type
echo $mutation_site

cat "${output_id}" |
    grep "${cluster}$" |
    cut -f 1 |
    sort -u |
    join tmp_sam - |
    awk '{print $4, $(NF-1)}' |
    sed "s/cs:Z://g" |
    awk '{seq=""
        padding=$1-1
        for(i=1; i<=padding; i++) seq=seq"-"
        print seq""$2}' |
    # -------------------------------------
    # 条件分岐
    # -------------------------------------
    awk -v type="${mutation_type}" -v site="${mutation_site}" \
    'BEGIN{
        split(type, type_, " ")
        split(site, site_, " ")
        }
    {original_seq = $0
    for(i in type_){
        $0 = original_seq
        if(type_[i] == "I"){
            gsub("*[a-z]", " ", $0)
            gsub("+", " +", $0)
            gsub("[-|=]", " ", $0)
            len=0
            for(j=1; j<=NF;j++) {
                if(len >= site_[i]-1) {print type_[i], len, $j; break}
                else { if($j !~ /+/) len+=length($j) }}
            }
        else {
            gsub("*[a-z]", "=", $0)
            gsub("+[a-z]*", "", $0)
            gsub("[-=]", "", $0)
            print type_[i], site_[i], substr($0, site_[i]+1, 1)
            }
        }
    }' |
    sort |
    uniq -c |
    sed "s/+//g" |
    awk '{if(max[$2] < $1) {max[$2] = $1; out[$2]=$3" "$4}}
        END{for(key in out) print key, out[key]}' |
cat - > tmp_mutation_


label=$(cat $output_result | awk -v cl="${cluster}" '$3==cl {print $1"_"$2"_#"$3}' | head -n 1)
mutation_type=$(cut -d " " -f 1 tmp_mutation_ | xargs echo)
mutation_site=$(cut -d " " -f 2 tmp_mutation_ | xargs echo)
mutation_nuc=$(cut -d " " -f 3 tmp_mutation_ | xargs echo)

cat .DAJIN_temp/fasta/wt.fa |
sed 1d |
awk -F "" -v type="${mutation_type}" -v site="${mutation_site}" -v nuc="${mutation_nuc}" \
    'BEGIN{
        split(type, type_, " ")
        split(site, site_, " ")
        split(nuc, nuc_, " ")
    }
    {for(i in type_){
        if(type_[i] == "S"){
            $(site_[i]+1) = "<span_class=\"Sub\">" nuc_[i] "</span>"
            }
        else if(type_[i] == "I"){
            $(site_[i]+1) = $site_[i] "<span_class=\"Ins\">" nuc_[i] "</span>"
            }
        else if(type_[i] == "D"){
            $(site_[i]+1) = "<span_class=\"Del\">" tolower(nuc_[i]) "</span>"
            }
    }}1' |
sed -e "s/ //g" -e "s/_/ /g"|
sed -e "1i >${label}" |
cat - > test_htmlbody.html



# ============================================================================
# HTMLの作製
# Output HTML file
# ============================================================================
cat << EOF > "${label}".html
<!DOCTYPE html>
<html>
<head>
<style>
p {
    font-family:"Courier New", Courier, monospace;
    color: #585858; 
    width: 50%;
    word-wrap: break-word;
}
.Ins {
    color: white; 
    background-color: #e2041b;
    font-weight: bold;
    font-size: 1.0em;
}
.Del {
    color: white; 
    background-color: #2ca9e1;
    font-weight: bold;
    font-size: 1.0em;
}
.Sub {
    color: white; 
    background-color: #007b43;
    font-weight: bold;
    font-size: 1.0em;
}

</style>
</head>
<body>
<p>
EOF

cat test_htmlbody.html >> "${label}".html

cat << EOF >> "${label}".html
</p>
<hr>
<p>
<span class="Ins">Insertion</span> <span class="Del">Deletion</span>
 <span class="Sub">Substitution</span>
</p>

</body>
</html>
EOF

