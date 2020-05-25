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
# barcode="barcode12"
# alleletype="wt"
# alleletype_original="normal"
# cluster=2
# suffix="${barcode}"_"${alleletype}"
# echo $suffix
# [ "$alleletype" = "normal" ] && alleletype="wt"
# [ "$alleletype" = "abnormal" ] && alleletype="wt"

barcode="${1}"
alleletype="${2}"
alleletype_original="${3}"
cluster="${4}"

suffix="${barcode}"_"${alleletype_original}"
# [ "$alleletype" = "normal" ] && alleletype="wt"
# [ "$alleletype" = "abnormal" ] && alleletype="wt"

mkdir -p ".DAJIN_temp/clustering/temp/" # 念のため

# ----------------------------------------
# Temporal Output
# ----------------------------------------
control_score=".DAJIN_temp/clustering/temp/control_score_${alleletype}"
consensus_mutation=".DAJIN_temp/clustering/temp/consensus_${suffix}"

tmp_mutation=".DAJIN_temp/clustering/temp/tmp_mutation_${suffix}"
# ----------------------------------------------------------
# Output results
# ----------------------------------------------------------
output_id=".DAJIN_temp/clustering/result_allele_id_${suffix}".txt
# output_result=".DAJIN_temp/clustering/result_allele_mutinfo_${suffix}".txt
suffix="${barcode}"_"${alleletype_original}"_"${cluster}"



# ============================================================================
# 変異情報のコンセンサスを得る
# ============================================================================
# true > "${consensus_mutation}"
# #
# cat "${allele_percentage}" |
# cut -d " " -f 1 |
# sort -u |
# while read -r cluster
# do
cat "${output_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 3 |
    # ----------------------------------------
    # 行を「リード指向」から「塩基部位指向」に変換する
    # ----------------------------------------
    awk -F "" \
    '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } }
    END {    
        for(j=1; j<=NF; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){ str=str""a[i,j] }
            print str }
    }' |
    # head -n 740 | tail -n 5 #! -----------------------------
    # ------------------------------------------
    # 各塩基部位において最多の変異をレポートする
    # ------------------------------------------
    awk -F "" '{sequence=$0
        sum[1]=gsub("=","=",sequence)
        sum[2]=gsub("M","M",sequence)
        sum[3]=gsub(/[1-9]|[a-z]/,"@", sequence)
        sum[4]=gsub("D","D",sequence)
        sum[5]=gsub("S","S",sequence)
        max=sum[1]; num=1
        for(i=2; i<=5;i++){if(max<sum[i]){max=sum[i]; num=i}}
        # ------------------------------------------
        # Insertion数をレポートする
        # ------------------------------------------
        max=0; ins_num=0
        if(num==3) {
            for(i=1; i<=NF; i++) { if($i ~ /[0-9]|[a-z]/) array[$i]++ }
            for(key in array){if(max<array[key]) {max=array[key]; ins_num=key}} 
        }
        
        print num, NR, ins_num, "@", (sum[1]+sum[2])/NF,sum[3]/NF,sum[4]/NF,sum[5]/NF
        }' |
    #
    paste - "${control_score}" |
    # head -n 740 | tail -n 5 | #! -----------------------------
    # ------------------------------------------
    # 各塩基部位にたいして「Mの頻度、Iの頻度、Dの頻度、Sの頻度、Iの個数」を表示する
    # ------------------------------------------
    # head test |
    awk 'function abs(v) {return v < 0 ? -v : v}
        $NF==1 {
            I=abs($6-$(NF-3))
            D=abs($7-$(NF-2))
            S=abs($8-$(NF-1))
            M=abs(1-I-D-S)
            print $2, M, I, D, S, $3
        }
        # Sequence errors are annontated as Match
        $NF==2 {
            print $2, 1, 0, 0, 0, 0
        }' |
    # ------------------------------------------
    # 各塩基部位にたいして「最大頻度の変異と挿入塩基数」を表示する
    # ------------------------------------------
    awk '{max=0; num=0
        for(i=2; i<=5;i++){if(max<$i){max=$i; num=i}}
        print num, $1, $NF}' |
    awk -v cl="${cluster}" \
        '{if($1==3) print cl, $2, "I", $NF
        else if($1==4) print cl, $2, "D", $NF
        else if($1==5) print cl, $2, "S", $NF}
        END{print "EOF"}' |
    # ----------------------------------------------------------
    # 挿入塩基数が10以上の場合に数値情報に逆変換する
    # ----------------------------------------------------------
    awk '{
        if($4 ~ /[a-z]/) {
            for (i=10; i<=36; i++) {
                num=i+87
                ins=sprintf("%c", num)
                if($4==ins) $4=i
            }
        }
        print $0}' |
    sed "s/35/>35/g" |    
    # ------------------------------------------
    # 変異がないアリルにはintacと表示する
    # ------------------------------------------
    awk -v cl="${cluster}" \
        '{if(NR==1 && $0=="EOF") print cl, 0, "intact",0
        else print $0}' |
    grep -v "EOF" |
cat - > "${consensus_mutation}"
# done

# # ============================================================================
# # Report allele mutation info
# # ============================================================================

# cat "${consensus_mutation}" |
#     awk '{num=1
#         cluster[$1]=cluster[$1]$3
#         if($3!="M"){
#             cl[NR] = $1
#             loc[NR] = $2
#             mut[NR] = $3
#             ins[NR] = $4
#         }}
#     END{
#         # もし変異がなければintactと表示する
#         for(j in cluster) {
#             if (cluster[j] !~/[I|D|S]/) {print j, 0,"intact"}
#         }
#         for(j in loc){print cl[j], loc[j], mut[j], ins[j]}
#         }' |
#     # ----------------------------------------------------------
#     # 挿入塩基数が10以上の場合に数値情報に逆変換する
#     # ----------------------------------------------------------
#     awk '{
#         if($4 ~ /[a-z]/) {
#             for (i=10; i<=36; i++) {
#                 num=i+87
#                 ins=sprintf("%c", num)
#                 if($4==ins) $4=i
#             }
#         }
#         print $0}' |
#     sed "s/35/>35/g" |
#     sort -t " " -k 1,1n -k 2,2n |
# cat - > ".DAJIN_temp/clustering/temp/tmp_${suffix}"

# # ----------------------------------------------------------
# # barcode, alleletype, クラスター番号, アレル頻度、変異、変異部位、ソート番号を出力する
# # ----------------------------------------------------------
# cat "$allele_percentage" |
#     sort |
#     join - ".DAJIN_temp/clustering/temp/tmp_${suffix}" |
#     sort -t " " -k 1,1n  -k 3,3n |
#     awk '{print $1,$2,$4,$5, NR}' |
#     sed "s/^/${barcode} ${alleletype_original} /g" |
# cat - > "${output_result}"


# ============================================================================
# 変異情報の同定
# Variant call
# ============================================================================
set $(cat "${consensus_mutation}" |
    awk -v cl="${cluster}" \
    '$1==cl {
    type=type$3"_"
    site=site$2"_"}
    END{print type, site}')
mutation_type=$(echo "$1" | sed "s/_/ /g") 
mutation_site=$(echo "$2" | sed "s/_/ /g")
# echo $mutation_type
# echo $mutation_site

ref=".DAJIN_temp/fasta/wt.fa"
cat .DAJIN_temp/fasta_ont/"${barcode}".fa |
    minimap2 -ax map-ont "${ref}" - --cs=long 2>/dev/null |
    sort |
cat - > .DAJIN_temp/clustering/temp/tmp_sam_"${suffix}"

cat "${output_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 1 |
    sort -u |
    join .DAJIN_temp/clustering/temp/tmp_sam_"${suffix}" - |
    awk '{print $4, $(NF-1)}' |
    sed "s/cs:Z://g" |
    awk '{seq=""
        padding=$1-1
        for(i=1; i<=padding; i++) seq=seq"-"
        print seq""$2}' |
    # -------------------------------------
    # 変異塩基の入手
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
            print type_[i], site_[i], substr($0, site_[i], 1)
            }
        }
    }' |
    sort |
    uniq -c |
    sed "s/+//g" |
    awk '{if(max[$2] < $1) {max[$2] = $1; out[$2]=$3" "$4}}
        END{for(key in out) print key, out[key]}' |
cat - > "${tmp_mutation}"

# cat tmp_mutation_
# minimap2 -ax map-ont .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/target.fa --cs |
# ============================================================================
# FASTAファイルの作製
# Output FASTA file
# ============================================================================


# ============================================================================
# HTMLの作製
# Output HTML file
# ============================================================================
mkdir -p .DAJIN_temp/clustering/html

if [ "$(grep -c intact $tmp_mutation)" -eq 1 ]; then
    html_filename=$(echo "${barcode}_allele${cluster}_intact_${alleletype}" |
        sed "s:^:.DAJIN_temp/clustering/html/:g")
else
    html_filename=$(echo "${barcode}_allele${cluster}_indel_${alleletype_original}" |
        sed "s:^:.DAJIN_temp/clustering/html/:g")
fi
mutation_type=$(cut -d " " -f 1 "${tmp_mutation}" | xargs echo)
mutation_site=$(cut -d " " -f 2 "${tmp_mutation}" | xargs echo)
mutation_nuc=$(cut -d " " -f 3 "${tmp_mutation}" | xargs echo)

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
    sed -e "1i >${html_filename}" |
cat - > .DAJIN_temp/clustering/temp/tmp_html_"${suffix}".html

cat << EOF > "${html_filename}".html
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

cat .DAJIN_temp/clustering/temp/tmp_html_"${suffix}".html >> "${html_filename}".html

cat << EOF >> "${html_filename}".html
</p>
<hr>
<p>
<span class="Ins">Insertion</span> <span class="Del">Deletion</span>
 <span class="Sub">Substitution</span>
</p>

</body>
</html>
EOF

echo "$html_filename"
