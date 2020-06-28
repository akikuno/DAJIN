#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
export PATH="$(command -p getconf PATH 2>/dev/null)${PATH+:}${PATH-}"
case $PATH in :*) PATH=${PATH#?};; esac
export UNIX_STD=2003  # to make HP-UX conform to POSIX

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST Auguments
#===========================================================
# barcode="barcode12"
# alleletype="wt"
# cluster=1
# percentage=72
# alleleid=1

# in_suffix="${barcode}"_"${alleletype}"
# out_suffix="${barcode}"_"${alleletype}"_"${alleleid}"
# mapping_alleletype="${alleletype}"
# [ "$alleletype" = "normal" ] && mapping_alleletype="wt"
# [ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

#===========================================================
#? Auguments
#===========================================================

barcode="${1}"
alleletype="${2}"
cluster="${3}"
percentage="${4}"
alleleid="${5}"

in_suffix="${barcode}"_"${alleletype}"
out_suffix="${barcode}"_"${alleletype}"_"${alleleid}"
mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"


#===========================================================
#? Input
#===========================================================
control_score=".DAJIN_temp/clustering/temp/control_score_${mapping_alleletype}"
allele_id=".DAJIN_temp/clustering/readid_cl_mids_${in_suffix}"

#===========================================================
#? Output
#===========================================================
mkdir -p .DAJIN_temp/consensus/temp
# .DAJIN_temp/consensus/"${output_filename}".fa
# .DAJIN_temp/consensus/"${output_filename}".html

#===========================================================
#? Temporal
#===========================================================
tmp_allele_id=".DAJIN_temp/consensus/temp/allele_id_${out_suffix}"
consensus_mutation=".DAJIN_temp/consensus/temp/consensus_${out_suffix}"
mutation_info=".DAJIN_temp/consensus/temp/mutation_info_${out_suffix}"
tmp_html=.DAJIN_temp/consensus/temp/tmp_html_"${out_suffix}".html

################################################################################
#! Get consensus sequence
################################################################################

#===========================================================
#? Detecting mutation sites
#===========================================================

cat "${allele_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 3 |
    sed "s/=/M/g" |
    awk -F "" 'BEGIN{OFS=","}{$1=$1}1' |
cat > "${tmp_allele_id}"

# cut -d "," -f 738-740 $tmp_allele_id | sort | uniq -c

Rscript DAJIN/src/consensus.R "${tmp_allele_id}" "${control_score}" "${cluster}"

#===========================================================
#? クラスタ番号、塩基番号、変異の種類、Insertion数の4つをレポートする
#===========================================================

if [ -s ".DAJIN_temp/consensus/temp/mutation_${out_suffix}" ]; then
cat ".DAJIN_temp/consensus/temp/mutation_${out_suffix}" |
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
cat > "${consensus_mutation}"

else
    echo "${cluster} 0 intact 0" > "${consensus_mutation}"
fi


################################################################################
#! 変異塩基の同定
# Variant call
################################################################################

if [ "$(grep -c intact ${consensus_mutation})" -eq 0 ]; then
set $(cat "${consensus_mutation}" |
    awk -v cl="${cluster}" '$1==cl {
    type=type$3"_"
    site=site$2"_"
    size=size$4"_"
    }
    END{print type, site, size}')
mutation_type=$(echo "$1" | sed "s/_/ /g") 
mutation_site=$(echo "$2" | sed "s/_/ /g")
insertion_size=$(echo "$3" | sed "s/_/ /g")
# echo $mutation_type
# echo $mutation_site
# echo $insertion_size

cat .DAJIN_temp/fasta_ont/"${barcode}".fa |
    minimap2 -ax map-ont \
        ".DAJIN_temp/fasta/${mapping_alleletype}.fa" - \
        --cs=long 2>/dev/null |
    sort |
cat > .DAJIN_temp/consensus/temp/tmp_sam_"${out_suffix}"

cat "${allele_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 1 |
    sort -u |
    join .DAJIN_temp/consensus/temp/tmp_sam_"${out_suffix}" - |
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
    awk -v type="${mutation_type}" \
        -v site="${mutation_site}" \
        -v size="${insertion_size}" \
    'BEGIN{
        split(type, type_, " ")
        split(site, site_, " ")
        split(size, size_, " ")
        }
    {original_seq = $0
    for(i in type_){
        $0 = original_seq
        if(type_[i] == "I"){
            gsub(/\*[a-z]/, " ", $0)
            gsub(/\+/, " +", $0)
            gsub(/[-|=]/, " ", $0)                
            len=0
            for(j = 1; j<=NF; j++){
                if(len >= site_[i]-1 && length($j) == size_[i] + 1) {print type_[i], len, $j; break}
                else { if($j !~ /\+/) len+=length($j) }
                }
            }
        else {
            gsub(/\*[a-z]/, "=", $0)
            gsub(/\+[a-z]*/, "", $0)
            gsub(/[-=]/, "", $0)
            print type_[i], site_[i], substr($0, site_[i], 1)
            }
        }
    }' |
    sort |
    uniq -c |
    sed "s/+//g" |
    awk '{if(max[$2] < $1) {max[$2] = $1; out[$2]=$3" "$4}}
        END{for(key in out) print key, out[key]}' |
cat > "${mutation_info}"
else
    echo "intact 0 0" > "${mutation_info}"
fi

# cat "${mutation_info}" #<<<<<

################################################################################
#! コンセンサス配列の作製
################################################################################

mutation_type=$(cut -d " " -f 1 "${mutation_info}" | xargs echo)
mutation_site=$(cut -d " " -f 2 "${mutation_info}" | xargs echo)
mutation_nuc=$(cut -d " " -f 3 "${mutation_info}" | xargs echo)

#===========================================================
#? FASTA file
#===========================================================

cat .DAJIN_temp/fasta/${mapping_alleletype}.fa |
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
                $(site_[i]) = ""
                }
        }}1' |
    sed -e "s/ //g" -e "s/_/ /g"|
cat > .DAJIN_temp/consensus/temp/"${out_suffix}"

#===========================================================
#? 変異情報をもとにOutput file nameをフォーマットする
#===========================================================

output_filename="${barcode}_allele${alleleid}"

diff_target=$(
    cat .DAJIN_temp/fasta/target.fa |
        sed 1d |
        diff - .DAJIN_temp/consensus/temp/${out_suffix} |
    wc -l
    )

diff_wt=$(
    cat .DAJIN_temp/fasta/wt.fa |
        sed 1d |
        diff - .DAJIN_temp/consensus/temp/${out_suffix} |
    wc -l
    )

if [ "${diff_target}" -eq 0 ]; then
    output_filename="${output_filename}_intact_target"
elif [ "${diff_wt}" -eq 0 ]; then
    output_filename="${output_filename}_intact_wt"
elif [ "$(grep -c intact $mutation_info)" -eq 1 ]; then
    output_filename="${output_filename}_intact_${alleletype}"
else
    output_filename="${output_filename}_mutation_${alleletype}"
fi

cat .DAJIN_temp/consensus/temp/"${out_suffix}" |
    fold |
    sed -e "1i >${output_filename}_${percentage}%" |
cat > .DAJIN_temp/consensus/"${output_filename}".fa

#===========================================================
#? HTML file
#===========================================================

cat ".DAJIN_temp/fasta/${mapping_alleletype}.fa" |
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
                $(site_[i]) = "<span_class=\"Del\">" nuc_[i] "</span>"
                }
        }}1' |
    sed -e "s/ //g" -e "s/_/ /g"|
    sed -e "1i >${output_filename}_${percentage}%" |
cat > "${tmp_html}"

cat << EOF > .DAJIN_temp/consensus/"${output_filename}".html
<!DOCTYPE html>
<html>
<head>
<style>
p {
    font-family: Consolas, monaco, "Courier New", Courier, monospace;
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
    color: black; 
    background-color: #66FFFF;
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

cat "${tmp_html}" >> .DAJIN_temp/consensus/"${output_filename}".html

cat << EOF >> .DAJIN_temp/consensus/"${output_filename}".html
</p>
<hr>
<p>
<span class="Ins">Insertion</span> <span class="Del">Deletion</span>
 <span class="Sub">Substitution</span>
</p>

</body>
</html>
EOF

# ls -l .DAJIN_temp/consensus/"${output_filename}".fa
rm "${tmp_allele_id}" "${consensus_mutation}" "${mutation_info}" "${tmp_html}"

# exit 0