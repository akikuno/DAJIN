#!/bin/sh

################################################################################
# Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C

################################################################################
# I/O naming
################################################################################

#===========================================================
# Auguments
#===========================================================

design=${1}
input_dir=${2}
grna=${3}

################################################################################
# Format FASTA file
################################################################################

cat "${design}" |
    tr -d "\r" |
    awk '{if($1~"^>"){print "\n"$0}
        else {printf $0}}
        END {print ""}' |
    grep -v "^$" |
    cat >.DAJIN_temp/fasta/fasta.fa

design_LF=".DAJIN_temp/fasta/fasta.fa"

#===========================================================
# Separate multiple-FASTA into FASTA files
#===========================================================

cat "${design_LF}" |
    awk 'NR%2==1 {printf $1" "; next}1' |
    awk '{id=$1
        gsub(">","",id)
        output=".DAJIN_temp/fasta/"id".fa"
        print $1"\n"toupper($2) > output
    }'

################################################################################
# Reverse complement if the mutation sites are closer
# to right flanking than left flanking
################################################################################

wt_len=$(awk '!/[>|@]/ {print length}' .DAJIN_temp/fasta/wt.fa)

convert_revcomp=$(
    minimap2 -ax splice .DAJIN_temp/fasta/wt.fa .DAJIN_temp/fasta/target.fa --cs 2>/dev/null |
        awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
        sed -e "s/cs:Z:://g" -e "s/:/ /g" -e "s/~/ /g" |
        tr -d "\~\*\-\+atgc" |
        awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
        awk -v wt_len="$wt_len" '{if(wt_len-$2>$1) print 0; else print 1}'
)

if [ "$convert_revcomp" -eq 1 ]; then
    cat "${design_LF}" |
        ./DAJIN/src/revcomp.sh - |
        cat >.DAJIN_temp/fasta/fasta_revcomp.fa
    design_LF=".DAJIN_temp/fasta/fasta_revcomp.fa"
    grna=$(echo "$grna" | ./DAJIN/src/revcomp.sh -)
fi

#===========================================================
# Separate multiple-FASTA into FASTA files
#===========================================================

cat ${design_LF} |
    awk 'NR%2==1 {printf $1" "; next}1' |
    awk '{id=$1
        gsub(">","",id)
        output=".DAJIN_temp/fasta_conv/"id".fa"
        print $1"\n"toupper($2) > output
    }'

################################################################################
# In the case of Point mutation:
# Generate randome insertion and deletion at gRNA sites as abnormal alleles
################################################################################

#===========================================================
# Define mutation type
# Insertion = I; Deletion = D; Substitution = S
#===========================================================

target_mutation_type=$(cat .DAJIN_temp/target_mutation_type)

#===========================================================
# Generate randome insertion and deletion at gRNA sites
#===========================================================

if [ "_${target_mutation_type}" = "_S" ]; then
    grna_len=$(awk -v grna="$grna" 'BEGIN{print length(grna)}')
    grna_firsthalf=$(awk -v grna="$grna" 'BEGIN{print substr(grna, 1, int(length(grna)/2))}')
    grna_secondhalf=$(awk -v grna="$grna" 'BEGIN{print substr(grna, int(length(grna)/2)+1, length(grna))}')
    # Randome sequence
    ins_seq=$(
        seq_length="$grna_len" &&
            od -A n -t u4 -N $(($seq_length * 100)) /dev/urandom |
            tr -d "\n" |
                sed 's/[^0-9]//g' |
                sed "s/[4-9]//g" |
                sed -e "s/0/A/g" -e "s/1/G/g" -e "s/2/C/g" -e "s/3/T/g" |
                awk -v seq_length="$seq_length" '{print substr($0, 1, seq_length)}'
    )
    # insertion
    cat .DAJIN_temp/fasta_conv/wt.fa |
        sed "s/$grna/$grna_firsthalf,$grna_secondhalf/g" |
        sed "s/,/$ins_seq/g" |
        sed "s/>wt/>wt_ins/g" |
        cat >.DAJIN_temp/fasta_conv/wt_ins.fa
    # deletion
    cat .DAJIN_temp/fasta_conv/wt.fa |
        sed "s/$grna//g" |
        sed "s/>wt/>wt_del/g" |
        cat >.DAJIN_temp/fasta_conv/wt_del.fa
fi

################################################################################
# Format ONT reads into FASTA file
################################################################################

find ${input_dir}/* -type f |
    grep -e ".fq" -e ".fastq" |
    while read -r input; do
        output=$(echo ${input%.f*}.fa | sed "s;${input_dir};.DAJIN_temp/fasta_ont;")
        # Check wheather the files are binary:
        if [ "$(file ${input} | grep -c compressed)" -eq 1 ]; then
            gzip -dc "${input}"
        else
            cat "${input}"
        fi |
            awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' |
            sed "s/^@/>/" |
            cat >"${output}"
    done

exit 0
