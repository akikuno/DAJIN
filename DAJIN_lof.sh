#!/bin/bash


################################################################################
#! Initialize shell environment
################################################################################

umask 0022
export LC_ALL=C
export UNIX_STD=2003  # to make HP-UX conform to POSIX

################################################################################
#! Argument parse
################################################################################

design=prdm14_cas9.fa
control=barcode26
threads=12
genome=mm10
input_dir=fastq

################################################################################
#! Formatting environments
################################################################################

#===========================================================
#? Make temporal directory
#===========================================================
rm -rf ".DAJIN_temp" 2>/dev/null || true
dirs="fasta fasta_conv fasta_ont NanoSim bam data"

echo "${dirs}" |
    sed "s:^:.DAJIN_temp/:g" |
    sed "s: : .DAJIN_temp/:g" |
xargs mkdir -p

#===========================================================
#? Format FASTA file
#===========================================================

cat "${design}" |
    tr -d "\r" |
    awk '{if($1~"^>"){print "\n"$0}
    else {printf $0}}
    END {print ""}' |
    grep -v "^$" |
cat > .DAJIN_temp/fasta/fasta.fa

design_LF=".DAJIN_temp/fasta/fasta.fa"

#---------------------------------------
#* Separate multiple-FASTA into FASTA files
#---------------------------------------

cat ${design_LF} |
    sed "s/^/@/g" |
    tr -d "\n" |
    sed -e "s/@>/\n>/g" \
        -e "s/@/ /g" \
        -e "s/$/\n/g" |
    grep -v "^$" |
    awk '{id=$1
        gsub(">","",id)
        output=".DAJIN_temp/fasta/"id".fa"
        print $1"\n"toupper($2) > output
    }'

#===========================================================
#? Reverse complement if the mutation sites are closer
#? to right flanking than left flanking
#===========================================================

wt_seqlen=$(awk '!/[>|@]/ {print length($0)}' .DAJIN_temp/fasta/wt.fa)

convert_revcomp=$(
    minimap2 -ax splice \
    .DAJIN_temp/fasta/wt.fa \
    .DAJIN_temp/fasta/target.fa --cs 2>/dev/null |
    awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
    sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
    tr -d "\~\*\-\+atgc" |
    awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
    awk -v wt_seqlen="$wt_seqlen" \
    '{if(wt_seqlen-$2>$1) print 0; else print 1}'
    )

if [ "$convert_revcomp" -eq 1 ] ; then
    cat "${design_LF}" |
        ./DAJIN/src/revcomp.sh - |
    cat > .DAJIN_temp/fasta/fasta_revcomp.fa
    design_LF=".DAJIN_temp/fasta/fasta_revcomp.fa"
fi

#---------------------------------------
#* Separate multiple-FASTA into FASTA files
#---------------------------------------
cat ${design_LF} |
    sed "s/^/@/g" |
    tr -d "\n" |
    sed -e "s/@>/\n>/g" \
        -e "s/@/ /g" \
        -e "s/$/\n/g" |
    grep -v "^$" |
    awk '{id=$1
        gsub(">","",id)
        output=".DAJIN_temp/fasta_conv/"id".fa"
        print $1"\n"toupper($2) > output
    }'

#===========================================================
#? Define mutation type
#? Insertion = I; Deletion = D; Substitution = S
#===========================================================

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

#---------------------------------------
#* In the case of Point mutation:
#* Generate randome insertion and deletion at gRNA sites as abnormal alleles
#---------------------------------------

if [ "_${mutation_type}" = "_S" ]; then
    grna_len=$(awk -v grna="$grna" 'BEGIN{print length(grna)}')
    grna_firsthalf=$(awk -v grna="$grna" 'BEGIN{print substr(grna, 1, int(length(grna)/2))}')
    grna_secondhalf=$(awk -v grna="$grna" 'BEGIN{print substr(grna, int(length(grna)/2)+1, length(grna))}')
    # Randome sequence
    ins_seq=$(
        seq_length="$grna_len" &&
        od -A n  -t u4 -N $(($seq_length*100)) /dev/urandom |
        tr -d "\n" |
        sed 's/[^0-9]//g' |
        sed "s/[4-9]//g" |
        sed -e "s/0/A/g" -e "s/1/G/g" -e "s/2/C/g" -e "s/3/T/g" |
        awk -v seq_length=$seq_length '{print substr($0, 1, seq_length)}'
    )
    # insertion
    cat .DAJIN_temp/fasta_conv/wt.fa |
        sed "s/$grna/$grna_firsthalf,$grna_secondhalf/g" |
        sed "s/,/$ins_seq/g" |
        sed "s/>wt/>wt_ins/g" |
    cat > .DAJIN_temp/fasta_conv/wt_ins.fa
    # deletion
    cat .DAJIN_temp/fasta_conv/wt.fa |
        sed "s/$grna//g" |
        sed "s/>wt/>wt_del/g" |
    cat > .DAJIN_temp/fasta_conv/wt_del.fa
fi

#===========================================================
#? Format ONT reads into FASTA file
#===========================================================

for input in ${input_dir}/* ; do
    output=$(
        echo "${input}" |
        sed -e "s#.*/#.DAJIN_temp/fasta_ont/#g" \
            -e "s#\.f.*#.fa#g")
    # Check wheather the files are binary:
    if [ "$(file ${input} | grep -c compressed)" -eq 1 ]
    then
        gzip -dc "${input}"
    else
        cat "${input}"
    fi |
    awk '{if((4+NR)%4==1 || (4+NR)%4==2) print $0}' |
    sed "s/^@/>/g" |
    cat > "${output}"
done

################################################################################
#! NanoSim (v2.5.0)
################################################################################

conda activate DAJIN_nanosim

cat << EOF
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NanoSim read simulation
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

EOF

#===========================================================
#? NanoSim
#===========================================================

./DAJIN/utils/NanoSim/src/read_analysis.py genome \
    -i ".DAJIN_temp/fasta_ont/${control}.fa" \
    -rg .DAJIN_temp/fasta_conv/wt.fa \
    -t ${threads:-1} \
    -o .DAJIN_temp/NanoSim/training

wt_seqlen=$(awk '!/[>|@]/ {print length($0)}' .DAJIN_temp/fasta/wt.fa)

for input in .DAJIN_temp/fasta_conv/*; do
    printf "${input} is now simulating...\n"
    output=$(
        echo "$input" |
        sed -e "s#fasta_conv/#fasta_ont/#g" \
            -e "s/.fasta$//g" -e "s/.fa$//g"
        )
    ## For deletion allele
    input_seqlength=$(
        cat "${input}" |
        awk '!/[>|@]/ {print length($0)-100}'
        )
    if [ "$input_seqlength" -lt "$wt_seqlen" ]; then
        len=${input_seqlength}
    else
        len=${wt_seqlen}
    fi
    ##
    ./DAJIN/utils/NanoSim/src/simulator.py genome \
        -dna_type linear \
        -c .DAJIN_temp/NanoSim/training \
        -rg "${input}" \
        -n 10000 \
        -t "${threads:-1}" \
        -min "${len}" \
        -o "${output}_simulated"
    ##
    rm .DAJIN_temp/fasta_ont/*_error_* .DAJIN_temp/fasta_ont/*_unaligned_* 2>/dev/null || true
done

rm -rf DAJIN/utils/NanoSim/src/__pycache__

printf 'Success!!\nSimulation is finished\n'

################################################################################
#! MIDS conversion
################################################################################

conda activate DAJIN

cat << EOF
++++++++++++++++++++++++++++++++++++++++++
Converting ACGT into MIDS format
++++++++++++++++++++++++++++++++++++++++++
EOF

reference=".DAJIN_temp/fasta_conv/wt.fa"
query=".DAJIN_temp/fasta_conv/target.fa"

#===========================================================
#? Get mutation loci
#===========================================================

cat "${reference}" |
    minimap2 -ax splice - "${query}" --cs 2>/dev/null |
    awk '{for(i=1; i<=NF;i++) if($i ~ /cs:Z/) print $i}' |
    sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
    tr -d "\~\*\-\+atgc" |
    awk '{$NF=0; for(i=1;i<=NF;i++) sum+=$i} END{print $1,sum}' |
cat > .DAJIN_temp/data/mutation_points

#===========================================================
#? MIDS conversion
#===========================================================

find .DAJIN_temp/fasta_ont -type f | sort |
    grep -e "${ont_cont}" -e sim |
    awk '{print "./DAJIN/src/mids_classification.sh", $0, "wt", "&"}' |
    awk -v th=${threads:-1} '{
        if (NR%th==0) gsub("&","&\nwait",$0)
        print}
        END{print "wait"}' |
sh - 2>/dev/null

[ "_${mutation_type}" = "_S" ] && rm .DAJIN_temp/data/MIDS_target*

printf "MIDS conversion was finished...\n"

#===========================================================
#? ACGT conversion
#===========================================================

find .DAJIN_temp/fasta_ont -type f | sort |
    grep -e "${ont_cont}" -e sim |
    awk '{print "./DAJIN/src/acgt_classification.sh", $0, "wt", "&"}' |
    awk -v th=${threads:-1} '{
        if (NR%th==0) gsub("&","&\nwait",$0)
        print}
        END{print "wait"}' |
sh - 2>/dev/null

[ "_${mutation_type}" = "_S" ] && rm .DAJIN_temp/data/MIDS_target*

printf "ACGT conversion was finished...\n"

#===========================================================
#? Prepare data
#===========================================================

#---------------------------------------
#* MIDS
#---------------------------------------

cat .DAJIN_temp/data/MIDS_* |
    grep "_sim" |
    grep -v "^ab" |
    sed -e "s/_aligned_reads//g" |
cat > ".DAJIN_temp/data/DAJIN_MIDS_train.txt"

head -n 1000 .DAJIN_temp/data/MIDS_* |
    grep -e negacon -e del |
    grep -v DAJIN_temp |
    grep -v "^$" |
    grep "^ab" |
    sed -e "s/_aligned_reads//g" |
cat > ".DAJIN_temp/data/DAJIN_MIDS_test.txt"

#---------------------------------------
#* ACGT
#---------------------------------------

cat .DAJIN_temp/data/ACGT_* |
    grep "_sim" |
    grep -v "^ab" |
    sed -e "s/_aligned_reads//g" |
cat > ".DAJIN_temp/data/DAJIN_ACGT_train.txt"

head -n 1000 .DAJIN_temp/data/ACGT_* |
    grep -e negacon -e del |
    grep -v DAJIN_temp |
    grep -v "^$" |
    grep "^ab" |
    sed -e "s/_aligned_reads//g" |
cat > ".DAJIN_temp/data/DAJIN_ACGT_test.txt"

################################################################################
#! Report best metric for LOF
################################################################################

#---------------------------------------
#* MIDS
#---------------------------------------
true > results_lof.csv

iter=20
cat << EOF |
python ./DAJIN/src/ml_metrics.py \
    ".DAJIN_temp/data/DAJIN_MIDS_train.txt" \
    ".DAJIN_temp/data/DAJIN_MIDS_test.txt" \
    "${threads}"
EOF
awk -v iter="${iter}" '{while(i<iter){
        i++
        print $0,i
        }}' |
sh -

cat results_lof.csv |
    grep -v barcodeID |
    awk -F "," 'BEGIN{OFS=","}
    $2=="normal" && $3=="1000" {
        print $0
        $2="abnormal"; $3="0"
        print $0; next}
    $2=="abnormal" && $3=="1000"{
        print $0
        $2="normal"; $3="0"
        print $0; next}1' |
cat > results_lof_add.csv

