#!/bin/sh
# ======================================
# Initialize shell environment
# ======================================
set -u
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export PATH="$(command -p getconf PATH)${PATH+:}${PATH-}"
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ======================================

cat "$1" |
awk 'BEGIN{FS=""}
{
    if($1 !~ "^>" && $1 !~ "^@"){
        nuc=""; seq=""; 
        for(i=1; i<=NF; i++){
            $i=toupper($i);
            if($i=="A") nuc="T";
            else if($i=="C") nuc="G";
            else if($i=="G") nuc="C";
            else if($i=="T") nuc="A";
            else nuc=$i;
            seq=seq""nuc
        }
        print seq
    }else print $0
}' |
awk 'BEGIN{FS=""}
{
    if($1!=">" && $1!="@"){
        seq=""
        for(i=NF; i>0; i--){
        seq=seq""$i
        }
    print seq
    }else print $0
}'

exit 0