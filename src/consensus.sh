#!/bin/sh

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
export UNIX_STD=2003  # to make HP-UX conform to POSIX

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? Auguments
#===========================================================

barcode="${1}"
alleletype="${2}"
cluster="${3}"
percentage="${4}"
alleleid="${5}"

#===========================================================
#? Input
#===========================================================

in_suffix="${barcode}"_"${alleletype}"
out_suffix="${barcode}"_"${alleletype}"_"${alleleid}"
mapping_alleletype="${alleletype}"
[ "$alleletype" = "normal" ] && mapping_alleletype="wt"
[ "$alleletype" = "abnormal" ] && mapping_alleletype="wt"

control_score=".DAJIN_temp/clustering/temp/possible_true_mut_${in_suffix}"
allele_id=".DAJIN_temp/clustering/allele_per/readid_cl_mids_${in_suffix}"

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
mutation_id_loc_type_insnum=".DAJIN_temp/consensus/temp/consensus_${out_suffix}"
mutation_type_site_nuc=".DAJIN_temp/consensus/temp/mutation_type_site_nuc_${out_suffix}"
tmp_html=.DAJIN_temp/consensus/temp/tmp_html_"${out_suffix}"

target_mutation_type=$(
  ref=".DAJIN_temp/fasta/wt.fa"
  que=".DAJIN_temp/fasta/target.fa"
  minimap2 -ax splice "$ref"  "$que"  --cs 2>/dev/null |
  awk '!/@/ {
  cstag=$(NF-1)
  if(cstag ~ "\~") print "D"
  else if(cstag ~ "\+") print "I"
  else if(cstag ~ "\*") print "S"
  }' 2>/dev/null
)

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

if [ -s "${control_score}" ]; then
  Rscript DAJIN/src/consensus.R "${tmp_allele_id}" "${control_score}" 2>/dev/null
else
  true > ".DAJIN_temp/consensus/temp/mutation_${out_suffix}"
fi
#===========================================================
#? Report (1) Cluster ID, (2) Base loc (3) Mutation type (4) Ins num
#===========================================================

if [ -s ".DAJIN_temp/consensus/temp/mutation_${out_suffix}" ]; then
  cat ".DAJIN_temp/consensus/temp/mutation_${out_suffix}" |
  sed "s/^/${cluster} /g" |
  #----------------------------------------------------------
  #* Converte Insertion character to number
  #----------------------------------------------------------
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
  cat > "${mutation_id_loc_type_insnum}"
else
  echo "${cluster} 0 intact 0" > "${mutation_id_loc_type_insnum}"
fi

################################################################################
#! Variant call
################################################################################

if [ "$(grep -c intact ${mutation_id_loc_type_insnum})" -eq 0 ]; then

  set $(
  cat "${mutation_id_loc_type_insnum}" |
  awk -v cl="${cluster}" '$1==cl {
    type=type$3"_"
    site=site$2"_"
    size=size$4"_"
    }
  END{print type, site, size}'
  )

  mutation_type=$(echo "$1" | sed "s/_/ /g")
  mutation_site=$(echo "$2" | sed "s/_/ /g")
  insertion_size=$(echo "$3" | sed "s/_/ /g")

  cat "${allele_id}" |
    awk -v cl="${cluster}" '$2==cl' |
    cut -f 1 |
    sort -u |
    join .DAJIN_temp/consensus/sam/"${barcode}"_"${mapping_alleletype}".sam - |
    awk '$2==0 || $2==16' |
    awk '{print $4, $(NF-1)}' |
    sed "s/cs:Z://g" |
    awk '{seq=""
      padding=$1-1
      for(i=1; i<=padding; i++) seq=seq"-"
      print seq""$2}' |
    #-------------------------------------
    #* Obtain mutation nucreotides
    #-------------------------------------
    awk -v type="${mutation_type}" \
      -v site="${mutation_site}" \
      -v size="${insertion_size}" \
    'BEGIN{
      split(type, type_, " ")
      split(site, site_, " ")
      split(size, size_, " ")
      }
    {
    original_seq = $0
    for(i in type_){
      $0 = original_seq
      if(type_[i] == "I"){
        gsub(/\*[a-z]/, " ", $0)
        gsub(/\+/, " +", $0)
        gsub(/[-|=]/, " ", $0)
        len=0
        for(j = 1; j<=NF; j++){
          if(len >= site_[i]-1 && length($j) == size_[i] + 1) {
            print type_[i], len, $j, i
            break
          } else {
            if($j !~ /\+/) len+=length($j)
          }
        }
      } else {
        $0 = original_seq
        gsub(/\*[a-z]/, "=", $0)
        gsub(/\+[a-z]*/, "", $0)
        gsub(/[-=]/, "", $0)
        print type_[i], site_[i], substr($0, site_[i], 1), i
      }
    }
  }' |
    sort |
    uniq -c |
    grep -v "[ACGT]" |
    sed "s/+//g" |
    awk 'NF==5{
      key=$2"_"$5
      if(max[key] < $1) {max[key] = $1; out[key]=$2" "$3" "$4}}
      END{for(i in out) print out[i]
      }' |
    sort -t " " -k2,2n |
  cat > "${mutation_type_site_nuc}"
else
  echo "intact 0 0" > "${mutation_type_site_nuc}"
fi

################################################################################
#! Report consensus sequence
################################################################################

mutation_type=$(cut -d " " -f 1 "${mutation_type_site_nuc}" | xargs echo)
mutation_site=$(cut -d " " -f 2 "${mutation_type_site_nuc}" | xargs echo)
mutation_nuc=$(cut -d " " -f 3 "${mutation_type_site_nuc}" | xargs echo)

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
  sed -e "s/ //g" -e "s/_/ /g" |
cat > .DAJIN_temp/consensus/temp/"${out_suffix}"

#===========================================================
#? Format output file name
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

#-------------------------------------
#* Annotate "target" when including targetted point mutation
#-------------------------------------

if [ "_${target_mutation_type}" = "_S" ]; then
  mutation_point=$(cut -d " " -f 1 .DAJIN_temp/data/mutation_points)

  contaion_target=$(
    cat "${mutation_type_site_nuc}" |
    awk -v mut="${target_mutation_type}" '$1==mut' |
    cut -d " " -f 2 |
    awk '{print $0-1}' |
    grep -c ${mutation_point} -
    )
else
  contaion_target=0
fi

if [ "${diff_target}" -eq 0 ]; then
  output_filename="${output_filename}_intact_target"
elif [ "${diff_wt}" -eq 0 ]; then
  output_filename="${output_filename}_intact_wt"
elif [ "${contaion_target:-0}" -ne 0 ]; then
  output_filename="${output_filename}_mutation_target"
elif [ "$(grep -c intact ${mutation_type_site_nuc})" -eq 1 ]; then
  output_filename="${output_filename}_intact_${alleletype}"
else
  output_filename="${output_filename}_mutation_${alleletype}"
fi

# echo "$output_filename"

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
  background-color: #FF4B00;
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
  background-color: #006E54; /* #03af7a */
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
<span class="Ins">Insertion</span> <span class="Del">Deletion</span> <span class="Sub">Substitution</span>
</p>

</body>
</html>
EOF

################################################################################
#! Move directory
################################################################################

mkdir -p .DAJIN_temp/consensus/FASTA .DAJIN_temp/consensus/HTML
mv .DAJIN_temp/consensus/"${output_filename}".fa .DAJIN_temp/consensus/FASTA
mv .DAJIN_temp/consensus/"${output_filename}".html .DAJIN_temp/consensus/HTML

exit 0