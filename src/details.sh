
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
#? Input
#===========================================================

#===========================================================
#? Output
#===========================================================

mkdir -p .DAJIN_temp/details

################################################################################
#! Generate Details.csv
################################################################################

find .DAJIN_temp/consensus/* -type f |
    grep html |
    sed "s:.*/::g" |
    sed "s/.html//g" |
    sed "s/_/ /g" |
    awk '{print $1"_"$2,$3,$4}' |
    sort |
cat > .DAJIN_temp/details/tmp_nameid

cat .DAJIN_temp/clustering/label* |
    awk '{nr[$1]++; print $0, nr[$1]}' |
    awk '{print $1"_allele"$5, $4, $2}' |
    sort |
    join -a 1 - .DAJIN_temp/details/tmp_nameid |
    sed "s/_/ /" |
    awk '$4=="abnormal" {$5="mutation"}1' |
    awk 'BEGIN{OFS=","}
        {gsub("allele","",$2)
        if($6 == "target") $4 = "target"
        if($4 == "abnormal") $6 ="+"; else $6 = "-"
        gsub("intact","-", $5)
        gsub("mutation","+", $5)

        if($4 == "target" && $5 == "-" && $6 == "-") $7 = "+"
        else $7 = "-"
        }1' |
    sed -e "1i Sample, Allele ID, % of reads, Allele type, Indel, Large indel, Design" |
cat > .DAJIN_temp/details/Details.csv

rm .DAJIN_temp/details/tmp_nameid

################################################################################
#! Plot details.csv
################################################################################

Rscript DAJIN/src/details_plot.R

exit 0