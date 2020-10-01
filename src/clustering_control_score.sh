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
#? TEST Auguments
#===========================================================

# control=barcode41
# threads=12

#===========================================================
#? Auguments
#===========================================================

control="${1}"
threads="${2}"

#===========================================================
#? Input
#===========================================================

alleletype="wt"

mutation_type=$(
    minimap2 -ax splice \
        .DAJIN_temp/fasta/wt.fa \
        .DAJIN_temp/fasta/target.fa \
        --cs 2>/dev/null |
    grep -v "^@" |
    awk '{
        cstag=$(NF-1)
        if(cstag ~ "~") print "D"
        else if(cstag ~ "\+") print "I"
        else if(cstag ~ "\*") print "S"
        }' 2>/dev/null
)

#===========================================================
#? Output
#===========================================================
mkdir -p ".DAJIN_temp/clustering/temp/"
control_score=".DAJIN_temp/clustering/temp/control_score_${alleletype}"

#===========================================================
#? Temporal
#===========================================================

MIDS_ref=".DAJIN_temp/clustering/temp/MIDS_${control}_${alleletype}"
tmp_MIDS=".DAJIN_temp/clustering/temp/tmp_MIDS_${control}_${alleletype}"
tmp_mask=".DAJIN_temp/clustering/temp/tmp_mask_${control}_${alleletype}"
tmp_control=".DAJIN_temp/clustering/temp/tmp_control_${control}_${alleletype}"
tmp_strecher=".DAJIN_temp/clustering/temp/tmp_strecher_${control}_${alleletype}"
tmp_strecher_wt=".DAJIN_temp/clustering/temp/tmp_strecher_wt_${control}_${alleletype}"
tmp_strecher_label=".DAJIN_temp/clustering/temp/tmp_strecher_label${control}_${alleletype}"

################################################################################
#! Define sequence error
################################################################################

#===========================================================
#? MIDS conversion
#===========================================================

./DAJIN/src/mids_clustering.sh "${control}" "${alleletype}" > "${MIDS_ref}"

#===========================================================
#? Select WT allele
#===========================================================

cat "${MIDS_ref}" |
    sort -k 1,1 |
    join - .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    awk -v alelle="$alleletype" '$NF==alelle' |
    cut -d " " -f 2 |
cat > "${tmp_MIDS}"

#===========================================================
#? Mask repeat sequences
#===========================================================

cat .DAJIN_temp/fasta/wt.fa |
    sed 1d |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    uniq -c |
    awk '$1>=6{$2=tolower($2)}{for(i=1;i<=$1;i++) print $2}' |
cat > "${tmp_mask}"

#===========================================================
#? Define sequence error
#===========================================================

nr=$(cat "${tmp_MIDS}" | wc -l)

cat "${tmp_MIDS}" |
    #----------------------------------------
    #* Transpose matrix
    #----------------------------------------
    awk -v nr="${nr}" -F "" \
    '{for(i=1; i<=NF; i++){
            row[i]=row[i] $i
            if(NR==nr) print row[i]
        }
    }' |
    #----------------------------------------
    #* Define sequence error when control sample has more than 7% mutations
    #----------------------------------------
    awk -F "" '{
        per=7
        INS=gsub(/[1-9]|[a-z]/,"@",$0)
        DEL=gsub("D","D",$0)
        SUB=gsub("S","S",$0)
        if(INS > NF*per/100 || DEL > NF*per/100 || SUB > NF*per/100)
            num=2
        else
            num=1

        print num
    }' |
    paste - "${tmp_mask}" |
    #----------------------------------------
    #* Define sequence error at repeat-masked loci
    #----------------------------------------
    awk '$2 ~ /[acgt]/{$1=2}{print $1}' |
cat > "${control_score}"

################################################################################
#! Generate scores of other possible alleles
################################################################################

cat ".DAJIN_temp/fasta/${alleletype}.fa" |
    sed 1d |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    paste - "${control_score}" |
    awk '{print NR"_"$1, $2}' |
    sort -t " " -k 1,1 |
cat > "${tmp_control}"


find .DAJIN_temp/fasta/ -type f |
    grep -v wt.fa |
    grep -v fasta.fa |
    sed "s:.*/::g" |
    sed "s/.fa.*$//g" |
while read -r label; do

    stretcher \
        -asequence .DAJIN_temp/fasta/wt.fa \
        -bsequence .DAJIN_temp/fasta/"${label}".fa \
        -aformat markx2 \
        -outfile "${tmp_strecher}" 2>/dev/null

    cat "${tmp_strecher}" |
        awk 'NF==2' |
        awk '$0 ~ /[ACGT.-]/' |
        awk 'NR%2==1 {printf $2}' |
        awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    cat > "${tmp_strecher_wt}"

    cat "${tmp_strecher}" |
        awk 'NF==2' |
        awk '$0 ~ /[ACGT.-]/' |
        awk 'NR%2==0 {printf $2}' |
        awk -F "" '{for(i=1;i<=NF;i++) print $i}' |
    cat > "${tmp_strecher_label}"

    splice_num="$(
        minimap2 -ax splice .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/${label}.fa 2>/dev/null |
        grep -c -v "^@"
        )"

    if [ "${splice_num}" -eq 3 ]; then
    # inversion
    set $(
        paste "${tmp_strecher_wt}" "${tmp_strecher_label}" |
        awk '$1=="-" {print "-"; next}
            {refrow++; print $1, $NF}' |
        awk '$2=="-" {$1="-"}1' |
        awk '{printf $1}' |
        sed "s/-\([ACTG].*\)-/-\n\1\n/g" | # flox-ki
        sed "s/--*$//g" | # flox-ki
        awk 'NR==2{gsub(".", "-", $0)} {printf $0}' | # flox-ki
        awk '{mutlen=gsub("-", " ", $0)
            print length($1)+1,length($1)+mutlen+2}'
        )

    cat "${tmp_control}" |
        sort -k 1,1n |
        awk '{print $NF}' |
        sed -e "$1i @" |
        sed -e "$2i @" |
        tr -d "\n" |
        sed "s/@/\n/g" |
        awk -F "" 'NR==2{
            seq=""
            for(i=NF; i>0; i--) seq=seq""$i
            print seq; next}1' |
        tr -d "\n" |
    awk -F "" '{for(i=1;i<=NF;i++) print $i}'
    else
    # deletion/knockin/point mutation
    paste "${tmp_strecher_wt}" "${tmp_strecher_label}" |
        awk '$1=="-" {print $0, NR; next}
            {refrow++; querow=NR; print refrow"_"$0, querow}' |
        sort |
        join -a 1 - "${tmp_control}" |
        awk '$2!="-"' | # deletion
        awk '$2~/[ACGT]/ {$4=1}1' | # point mutation
        awk '$1=="-" {$4=100}1' | # knockin
        sort -t " " -k 3,3n |
        awk '{print $NF}'
    fi |
    cat > ".DAJIN_temp/clustering/temp/control_score_${label}"

done

[ _"${mutation_type}" = "_S" ] &&
mv ".DAJIN_temp/clustering/temp/control_score_target" "${control_score}"

################################################################################
#! remove temporal files
################################################################################

rm "${MIDS_ref}"
rm "${tmp_MIDS}"
rm "${tmp_control}"
rm "${tmp_strecher}"
rm "${tmp_strecher_wt}"
rm "${tmp_strecher_label}"

exit 0