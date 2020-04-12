#!/bin/sh

fastq_pass="tyr-mut-01/20200328_0806_MN33661_FAN31330_f6cf212d/fastq_pass"

fastq_split=$(
    find "${fastq_pass}" -type f |
    wc -l |
    awk -v th="${threads:-1}" \
    '{print int($0/th)}')

find "${fastq_pass}" -type f |
split -l "${fastq_split}" - .DAJIN_temp/qcat_dir/tmp_qcat

find .DAJIN_temp/qcat_dir -name "tmp_qcat*" -type f | sed "s#.*/##g" |
xargs -I @ mkdir -p .DAJIN_temp/qcat_dir/@_dir

for file in $(find .DAJIN_temp/qcat_dir -name "tmp_qcat*" -type f); do
    echo "${file} is now processing..." &&
    #
    cat ${file} | xargs -I @ cat @ | 
    qcat -k PBC096 -b "${file}"_dir &&
    #
    echo "${file} is finished..." &
done
wait

find .DAJIN_temp/qcat_dir -name "barcode*" -type f |
awk -F "/" \
'{barcode[$NF]=barcode[$NF]" "$0}
END{for(key in barcode) print "cat "barcode[key]" > .DAJIN_temp/qcat_dir/"key" &"}' |
awk -v th=${threads:-1} '{
    if (NR%th==0) gsub("&","&\nwait",$0)
    print}
    END{print "wait"}' |
sh -
rm -rf .DAJIN_temp/qcat_dir/tmp_*

total=$(wc -l .DAJIN_temp/qcat_dir/*q | grep total | awk '{print $1}')

wc -l .DAJIN_temp/qcat_dir/*q | grep -v total |
awk -v total="${total}" \
'{if ($1/total*100 > 0.001) print $2}' |
awk -F "/" \
'{if ($NF ~ NR) print }' |
sort \
> .DAJIN_temp/qcat_dir/tmp_fasta

find .DAJIN_temp/qcat_dir -name "barcode*" |
sort |
join -v 1 - .DAJIN_temp/qcat_dir/tmp_fasta |
xargs -I @ rm @
rm .DAJIN_temp/qcat_dir/tmp_fasta
