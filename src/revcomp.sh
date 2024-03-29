#!/bin/sh

# ======================================
# Initialize shell environment
# ======================================
set -u
umask 0022
export LC_ALL=C

# ======================================

cat "$1" |
    # Complement
    awk 'BEGIN{FS=""}
    {if($1 !~ "^>" && $1 !~ "^@"){
        nuc=""; seq=""
        for(i=1; i<=NF; i++){
            $i=toupper($i)
            if($i=="A") nuc="T"
            else if($i=="C") nuc="G"
            else if($i=="G") nuc="C"
            else if($i=="T") nuc="A"
            else nuc=$i
            seq=seq""nuc
        }
        print seq
    } else print $0}' |
    # Reverse
    awk 'BEGIN{FS=""}
    {if($1!=">" && $1!="@"){
        seq=""
        for(i=NF; i>0; i--){
        seq=seq""$i
        }
    print seq
    } else print $0}'

exit 0
