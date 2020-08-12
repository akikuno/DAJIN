#!/bin/sh

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