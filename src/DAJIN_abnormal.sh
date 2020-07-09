#!/bin/bash

################################################################################
#! Initialize shell environment
################################################################################

set -u
umask 0022
export LC_ALL=C
export UNIX_STD=2003  # to make HP-UX conform to POSIX

error_exit() {
    echo "$@" 1>&2
    exit 1
}

################################################################################
#! I/O settings
################################################################################
#===========================================================
#? Input
#===========================================================

design=DAJIN/misc/abnormal_detection/cables2_anomaldetection.fa
ont_cont=barcode21
ont_dir=fastq

#===========================================================
#? Define threads
#===========================================================

expr "${threads}" + 1 >/dev/null 2>&1
if [ $? -lt 2 ]; then
    :
else
    unset threads
    # Linux and similar...
    threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/1.5+0.5)}')
    # FreeBSD and similar...
    [ -z "$threads" ] && threads=$(getconf NPROCESSORS_ONLN | awk '{print int($0/1.5)+0.5}')
    # Solaris and similar...
    [ -z "$threads" ] && threads=$(ksh93 -c 'getconf NPROCESSORS_ONLN' | awk '{print int($0/1.5+0.5)}')
    # Give up...
    [ -z "$threads" ] && threads=1
fi

################################################################################
#! Setting Conda environment
################################################################################

#===========================================================
#? DAJIN_nanosim
#===========================================================
set +u
type conda > /dev/null 2>&1 || error_exit 'Command "conda" not found'

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [ "$(conda info -e | grep -c DAJIN_nanosim)" -eq 0 ]; then
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda create -y -n DAJIN_nanosim python=3.6
    conda install -y -n DAJIN_nanosim --file ./DAJIN/utils/NanoSim/requirements.txt
    conda install -y -n DAJIN_nanosim minimap2
fi

#===========================================================
#? DAJIN
#===========================================================

if [ "$(conda info -e | cut -d " " -f 1 | grep -c DAJIN$)" -eq 0 ]; then
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda create -y -n DAJIN python=3.7 \
        anaconda nodejs wget \
        tensorflow tensorflow-gpu \
        samtools minimap2 \
        r-essentials r-base
fi
set -u

#===========================================================
#? Required software
#===========================================================
set +u
conda activate DAJIN
set -u

type gzip > /dev/null 2>&1 || error_exit 'Command "gzip" not found'
type wget > /dev/null 2>&1 || error_exit 'Command "wget" not found'
type python > /dev/null 2>&1 || error_exit 'Command "python" not found'
type samtools > /dev/null 2>&1 || error_exit 'Command "samtools" not found'
type minimap2 > /dev/null 2>&1 || error_exit 'Command "minimap2" not found'

python -c "import tensorflow as tf" > /dev/null 2>&1 ||
error_exit '"Tensorflow" not found'

#===========================================================
#? For WSL (Windows Subsystem for Linux)
#===========================================================

uname -a |
grep Microsoft 1>/dev/null 2>/dev/null &&
alias python="python.exe"

################################################################################
#! Formatting environments
################################################################################

#===========================================================
#? Make temporal directory
#===========================================================
rm -rf ".DAJIN_temp" 2>/dev/null || true
dirs="fasta fasta_conv fasta_ont NanoSim bam igvjs data clustering/temp seqlogo/temp"

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
    awk '{print $(NF-1)}' |
    sed -e "s/cs:Z:://g" -e "s/:/\t/g" -e "s/~/\t/g" |
    tr -d "\~\*\-\+atgc" |
    awk '{$NF=0
        for(i=1;i<=NF;i++) sum+=$i
        }END{print $1,sum}' |
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

#---------------------------------------
#* Define mutation type
#---------------------------------------

mutation_type=$(
    minimap2 -ax map-ont \
        .DAJIN_temp/fasta/wt.fa \
        .DAJIN_temp/fasta/target.fa \
        --cs 2>/dev/null |
    grep -v "^@" |
    awk '{
        cstag=$(NF-1)
        if(cstag ~ "-") print "D"
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


for input in ${ont_dir}/* ; do
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
set +u
conda activate DAJIN_nanosim
set -u

cat << EOF
++++++++++++++++++++++++++++++++++++++++++
NanoSim read simulation
++++++++++++++++++++++++++++++++++++++++++
EOF

#===========================================================
#? NanoSim
#===========================================================

printf "Read analysis...\n"
./DAJIN/utils/NanoSim/src/read_analysis.py genome \
    -i ".DAJIN_temp/fasta_ont/${ont_cont}.fa" \
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
set +u
conda activate DAJIN
set -u

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
    awk '{print "./DAJIN/src/acgt_classification.sh", $0, "wt", "&"}' |
    awk -v th=${threads:-1} '{
        if (NR%th==0) gsub("&","&\nwait",$0)
        print}
        END{print "wait"}' |
sh - 2>/dev/null

[ "_${mutation_type}" = "_S" ] && rm .DAJIN_temp/data/MIDS_target*

printf "ACGT conversion was finished...\n"

################################################################################
#! Prediction
################################################################################

cat << EOF
++++++++++++++++++++++++++++++++++++++++++
Allele prediction
++++++++++++++++++++++++++++++++++++++++++
EOF


#===========================================================
#? Prepare simulation data
#===========================================================

#---------------------------------------
#* MIDS
#---------------------------------------

cat .DAJIN_temp/data/MIDS_* |
    grep "_sim" |
    grep -v "^ab_" |
    sed -e "s/_aligned_reads//g" |
cat > ".DAJIN_temp/data/DAJIN_MIDS_sim.txt"

cat .DAJIN_temp/data/MIDS_"${ont_cont}"_wt |
    grep -v "IIIIIIIIII" |
    grep -v "DDDDDDDDDD" |
    grep -v "SSSSSSSSSS" |
    head -n 10000 |
    sed "s/${ont_cont}$/wt_simulated/g" |
cat >> ".DAJIN_temp/data/DAJIN_MIDS_sim.txt"

#---------------------------------------
#* ACGT
#---------------------------------------

cat .DAJIN_temp/data/ACGT_* |
    grep "_sim" |
    grep -v "^ab_" |
    sed -e "s/_aligned_reads//g" |
cat > ".DAJIN_temp/data/DAJIN_ACGT_sim.txt"

cat .DAJIN_temp/data/ACGT_"${ont_cont}"_wt |
    grep -v "IIIIIIIIII" |
    grep -v "DDDDDDDDDD" |
    grep -v "SSSSSSSSSS" |
    head -n 10000 |
    sed "s/${ont_cont}$/wt_simulated/g" |
cat >> ".DAJIN_temp/data/DAJIN_ACGT_sim.txt"

#===========================================================
#? Train model
#===========================================================

#---------------------------------------
#* MIDS
#---------------------------------------

python ./DAJIN/src/ml_simulated.py \
    ".DAJIN_temp/data/DAJIN_MIDS_sim.txt" "${threads}"

#---------------------------------------
#* ACGT
#---------------------------------------

#===========================================================
#? Predict labels
#===========================================================

true > ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt"

find .DAJIN_temp/data/MIDS* |
    grep -v sim |
    sort |
while read -r input; do
    barcode=$(echo $input | cut -d "_" -f 3)
    echo "${barcode} is now processing..."

    python ./DAJIN/src/ml_real.py \
        "${input}" \
        "${mutation_type}" "${threads}" ||
    exit 1
done

cat .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    sort |
cat > .DAJIN_temp/tmp_$$
mv .DAJIN_temp/tmp_$$ .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt

printf "Prediction was finished...\n"

# #===========================================================
# #? Report the percentage of alleles in each sample
# #===========================================================

#---------------------------------------
#* Report the percentage of alleles in each sample
#---------------------------------------

cat .DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt |
    cut -f 2,3 |
    sort |
    uniq -c |
    awk '{barcode[$2]+=$1
        read_info[$2]=$1"____"$3" "read_info[$2]}
    END{for(key in barcode) print key,barcode[key], read_info[key]}' |
    awk '{for(i=3;i<=NF; i++) print $1,$2,$i}' |
    sed "s/____/ /g" |
    awk '{print $1, $3/$2*100, $4}' |
cat > ".DAJIN_temp/tmp_prediction_proportion"


exit 0
