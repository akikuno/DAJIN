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
# Input arguments
# ----------------------------------------
# barcode="barcode02"
# alleletype="normal"
# cluster=1
# percentage=59
# alleleid=3
# in_suffix="${barcode}"_"${alleletype}"
# out_suffix="${barcode}"_"${alleletype}"_"${alleleid}"

# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

barcode="${1}"
alleletype="${2}"
cluster="${3}"
percentage="${4}"
alleleid="${5}"

in_suffix="${barcode}"_"${alleletype}"
out_suffix="${barcode}"_"${alleletype}"_"${alleleid}"

mapping_alleletype="wt"
# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"


# ----------------------------------------------------------
# Input files
# ----------------------------------------------------------
control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"
allele_id=".DAJIN_temp/clustering/result_allele_id_${in_suffix}".txt

# ----------------------------------------------------------
# Output files
# ----------------------------------------------------------
# temporal --------------------------------
tmp_allele_id=".DAJIN_temp/clustering/temp/allele_id_${out_suffix}"

# results --------------------------------
consensus_mutation=".DAJIN_temp/clustering/temp/consensus_${out_suffix}"
mutation_info=".DAJIN_temp/clustering/temp/mutation_info_${out_suffix}"

# ============================================================================
# 変異情報のコンセンサスを得る
# ============================================================================
#?====================================================================================

cat "${allele_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 3 |
    sed "s/=/M/g" |
    awk -F "" 'BEGIN{OFS=","}{$1=$1}1' |
cat - > "${tmp_allele_id}"

cat "${allele_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 3 |
    # sed "s/z.*/ /g" |
    # awk '{print length($1)}' |
    # sort | uniq -c
    awk '{print substr($0,738,1)}' |
    sort | uniq -c

Rscript DAJIN/src/clustering_variantcall.R "${tmp_allele_id}" "${control_score}" "${cluster}"

if [ -s ".DAJIN_temp/clustering/temp/mutation_${out_suffix}" ]; then
    cat ".DAJIN_temp/clustering/temp/mutation_${out_suffix}" |
        sed "s/^/${cluster} /g" |
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
        sed "s/35$/>35/g" |    
    cat - > "${consensus_mutation}"
else
    echo "${cluster} 0 intact 0" > "${consensus_mutation}"
fi

# done
#?====================================================================================

# start=$(cut -f 2 .DAJIN_temp/data/gggenome_location)

# len=$(samtools view DAJIN_results/BAM/${barcode}_allele${alleleid}.bam |
#     awk '$2==0 || $2==16' |
#     awk '{print length($10)}' |
#     sort | uniq -c |
#     awk '{if(max < $1) {max=$1; len=$2}} END{print len}')

# samtools view DAJIN_results/BAM/${barcode}_allele${alleleid}.bam |
#     awk '$2==0 || $2==16' |
#     awk -v start="${start}" -v len="${len}" '$4==start && length($10) == len' |
#     head -n 1 |
#     awk '{print ">"$1"\n"$10}' |
# cat > test_ref.fa


# cat "${allele_id}" |
#     awk -v cl="${cluster}" '$2==cl' |
#     cut -f 1 |
#     sort |
# cat > tmp_allele_id

# cat .DAJIN_temp/fasta_ont/"${barcode}".fa |
#     awk '$1~/[>@]/ {gsub("@","",$1); printf $1"\t"; next}1' |
#     sort |
#     join - tmp_allele_id |
#     awk '{print ">"$1"\n"$2}' |
# cat > test_que.fa 
# minimap2 -ax map-ont test_ref.fa test_que.fa --cs=long 2>/dev/null |
# cat > test.sam



# samtools view DAJIN_results/BAM/${barcode}_allele${alleleid}.bam |
#     awk '$2==0 || $2==16' |
#     awk -v start="${start}" -v len="${len}" '$4==start && length($10) == len' |
#     awk '{print ">"$1"\n"$10}' |

# minimap2 -ax map-ont .DAJIN_temp/fasta/wt.fa - --cs=long 2>/dev/null |
# awk '{print $(NF-1)}' |
# grep "cs" |
# awk '{cstag=$0
#     gsub("cs:Z:=","",$0)
#     gsub("=", " ", $0)
#     gsub(/[ACGT]/, "M", $0)
#     gsub(/\*[acgt][acgt]/, " S", $0)
#     gsub(/\+[acgt]*/,  " I ", $0)
#     gsub("-",  " ", $0)
#     for(i=1; i<=NF; i++) if($i !~ /[MSI+]/ ){
#         len="%" length($i) "s"
#         D=sprintf(len,""); gsub(/ /," D ",D)
#         $i=D
#         }
#     gsub(" ","", $0)
#     }1' |

# cat > test_MIDS

# cat test_MIDS |
#     awk '{print substr($0,721,1)}' |
#     sort | uniq -c



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
cat - > .DAJIN_temp/clustering/temp/tmp_sam_"${out_suffix}"

cat "${allele_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 1 |
    sort -u |
    join .DAJIN_temp/clustering/temp/tmp_sam_"${out_suffix}" - |
    awk '$2==0 || $2==16' |
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
cat - > "${mutation_info}"

# ============================================================================
# コンセンサス配列の作製
# ============================================================================
mkdir -p .DAJIN_temp/clustering/consensus

mutation_type=$(cut -d " " -f 1 "${mutation_info}" | xargs echo)
mutation_site=$(cut -d " " -f 2 "${mutation_info}" | xargs echo)
mutation_nuc=$(cut -d " " -f 3 "${mutation_info}" | xargs echo)

# -------------------------------
# FASTA file
# -------------------------------

cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    awk -F "" -v type="${mutation_type}" -v site="${mutation_site}" -v nuc="${mutation_nuc}" \
        'BEGIN{
            split(type, type_, " ")
            split(site, site_, " ")
            split(nuc, nuc_, " ")
        }
        {for(i in type_){
            nuc_[i] = toupper(nuc_[i])
            if(type_[i] == "S"){
                $(site_[i]) = nuc_[i]
                }
            else if(type_[i] == "I"){
                $(site_[i]) = $site_[i]""nuc_[i]
                }
            else if(type_[i] == "D"){
                $(site_[i]) = nuc_[i]
                }
        }}1' |
    sed -e "s/ //g" -e "s/_/ /g"|
cat - > .DAJIN_temp/clustering/temp/"${out_suffix}".fa

# -------------------------------
# 出力ファイル名をフォーマット
# -------------------------------
output_filename="${barcode}_allele${alleleid}"

diff_wt=$(cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    diff - .DAJIN_temp/clustering/temp/${out_suffix}.fa |
    wc -l)
diff_target=$(cat .DAJIN_temp/fasta/target.fa |
    sed 1d |
    diff - .DAJIN_temp/clustering/temp/${out_suffix}.fa |
    wc -l)

# if [ "$(awk '$1=="intact"' ${mutation_info} | wc -l)" -eq 0 ]; then
#     include_target=$(
#         cat .DAJIN_temp/data/mutation_points |
#         awk '{print "S",$1+1}' |
#         grep -c - "${mutation_info}")
# else
#     include_target=0
# fi
include_target=1


if [ "${diff_wt}" -eq 0 ]; then
    output_filename="${output_filename}_intact_wt"
elif [ "${diff_target}" -eq 0 ]; then
    output_filename="${output_filename}_intact_target"
elif [ "${include_target}" -gt 0 ]; then
    output_filename="${output_filename}_mutation_target"
elif [ "$(grep -c intact $mutation_info)" -eq 1 ]; then
    output_filename="${output_filename}_intact_${alleletype}"
else
    output_filename="${output_filename}_mutation"
fi

# [ "$(cat .DAJIN_temp/fasta/target.fa | sed 1d |
# diff - .DAJIN_temp/clustering/temp/${out_suffix}.fa |
# wc -l)" -eq 0 ] && output_filename="${output_filename}_intact_target"

# [ "$(cat .DAJIN_temp/fasta/wt.fa | sed 1d |
# diff - .DAJIN_temp/clustering/temp/${out_suffix}.fa |
# wc -l)" -eq 0 ] && output_filename="${output_filename}_intact_wt"

cat .DAJIN_temp/clustering/temp/"${out_suffix}".fa |
    sed -e "1i >${output_filename}_${percentage}%" |
cat - > .DAJIN_temp/clustering/consensus/"${output_filename}".fa

# -------------------------------
# HTML file
# -------------------------------

cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    awk -F "" -v type="${mutation_type}" -v site="${mutation_site}" -v nuc="${mutation_nuc}" \
        'BEGIN{
            split(type, type_, " ")
            split(site, site_, " ")
            split(nuc, nuc_, " ")
        }
        {for(i in type_){
            nuc_[i] = toupper(nuc_[i])
            if(type_[i] == "S"){
                $(site_[i]) = "<span_class=\"Sub\">" nuc_[i] "</span>"
                }
            else if(type_[i] == "I"){
                $(site_[i]) = $site_[i] "<span_class=\"Ins\">" nuc_[i] "</span>"
                }
            else if(type_[i] == "D"){
                $(site_[i]) = "<span_class=\"Del\">" tolower(nuc_[i]) "</span>"
                }
        }}1' |
    sed -e "s/ //g" -e "s/_/ /g"|
    sed -e "1i >${output_filename}_${percentage}%" |
cat - > .DAJIN_temp/clustering/temp/tmp_html_"${out_suffix}".html

cat << EOF > .DAJIN_temp/clustering/consensus/"${output_filename}".html
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

cat .DAJIN_temp/clustering/temp/tmp_html_"${out_suffix}".html |
cat - >> .DAJIN_temp/clustering/consensus/"${output_filename}".html

cat << EOF >> .DAJIN_temp/clustering/consensus/"${output_filename}".html
</p>
<hr>
<p>
<span class="Ins">Insertion</span> <span class="Del">Deletion</span>
 <span class="Sub">Substitution</span>
</p>

</body>
</html>
EOF

ls -l .DAJIN_temp/clustering/consensus/"${output_filename}".fa




# cat "${allele_id}" |
#     awk -v cl="${cluster}" '$2==cl' |
#     cut -f 3 |
#     # ----------------------------------------
#     # 行を「リード指向」から「塩基部位指向」に変換する
#     # ----------------------------------------
#     awk -F "" \
#     '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } }
#     END {    
#         for(j=1; j<=NF; j++) {
#             str=a[1,j]
#             for(i=2; i<=NR; i++){ str=str""a[i,j] }
#             print str }
#     }' |
#     # head -n 740 | tail -n 5 #! -----------------------------
#     # ------------------------------------------
#     # 各塩基部位において最多の変異をレポートする
#     # ------------------------------------------
#     awk -F "" '{sequence=$0
#         sum[1]=gsub("=","=",sequence)
#         sum[2]=gsub("M","M",sequence)
#         sum[3]=gsub(/[1-9]|[a-z]/,"@", sequence)
#         sum[4]=gsub("D","D",sequence)
#         sum[5]=gsub("S","S",sequence)
#         max=sum[1]; num=1
#         for(i=2; i<=5;i++){if(max<sum[i]){max=sum[i]; num=i}}
#         # ------------------------------------------
#         # Insertion数をレポートする
#         # ------------------------------------------
#         max=0; ins_num=0
#         if(num==3) {
#             for(i=1; i<=NF; i++) { if($i ~ /[0-9]|[a-z]/) array[$i]++ }
#             for(key in array){if(max<array[key]) {max=array[key]; ins_num=key}} 
#         }
        
#         print num, NR, ins_num, "@", (sum[1]+sum[2])/NF,sum[3]/NF,sum[4]/NF,sum[5]/NF
#         }' |
#     #
#     paste - "${control_score}" |
#     # head -n 740 | tail -n 5 | #! -----------------------------
#     # ------------------------------------------
#     # 各塩基部位にたいして「Mの頻度、Iの頻度、Dの頻度、Sの頻度、Iの個数」を表示する
#     # ------------------------------------------
#     # head test |
#     awk 'function abs(v) {return v < 0 ? -v : v}
#         $NF==1 {
#             I=abs($6-$(NF-3))
#             D=abs($7-$(NF-2))
#             S=abs($8-$(NF-1))
#             M=abs(1-I-D-S)
#             print $2, M, I, D, S, $3
#         }
#         # Sequence errors are annontated as Match
#         $NF==2 {
#             print $2, 1, 0, 0, 0, 0
#         }' |
#     # ------------------------------------------
#     # 各塩基部位にたいして「最大頻度の変異と挿入塩基数」を表示する
#     # ------------------------------------------
#     awk '{max=0; num=0
#         for(i=2; i<=5;i++){if(max<$i){max=$i; num=i}}
#         print num, $1, $NF}' |
#     awk -v cl="${cluster}" \
#         '{if($1==3) print cl, $2, "I", $NF
#         else if($1==4) print cl, $2, "D", $NF
#         else if($1==5) print cl, $2, "S", $NF}
#         END{print "EOF"}' |
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
#     # ------------------------------------------
#     # 変異がないアリルにはintacと表示する
#     # ------------------------------------------
#     awk -v cl="${cluster}" \
#         '{if(NR==1 && $0=="EOF") print cl, 0, "intact",0
#         else print $0}' |
#     grep -v "EOF" |
# cat - > "${consensus_mutation}"
# # done
