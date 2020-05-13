#!/bin/sh
# ============================================================================
# Initialize shell environment
# ============================================================================

set -eu
umask 0022
export LC_ALL=C
type command >/dev/null 2>&1 && type getconf >/dev/null 2>&1 &&
export UNIX_STD=2003  # to make HP-UX conform to POSIX

# ============================================================================
# Define the functions for printing usage and error message
# ============================================================================
VERSION=1.0

usuage(){
cat <<- USAGE 1>&2
Usage     : ./DAJIN.sh -f [text file](described at "Input")

Example   : ./DAJIN.sh -f DAJIN/example/example.txt

Input     : Input file should be formatted as below:
            # Example
            ------
            design=DAJIN/example/input.txt
            sequence=DAJIN/example/demultiplex
            control=barcode21
            genome=mm10
            grna=CCCTGCGGCCAGCTTTCAGGCAG
            threads=10
            ------
            - desing: a multi-FASTA file contains sequences of each genotype. ">wt" and ">target" must be included. 
            - sequence: a directory contains FASTA or FASTQ files of long-read sequencing
            - control: control barcode ID
            - genome: reference genome. e.g. mm10, hg38
            - grna: gRNA sequence(s). multiple gRNA sequences must be deliminated by comma. e.g. CCCTGCGGCCAGCTTTCAGGCAG,CCCTGCGGCCAGCTTTCAGGCAG
            - threads: optional. default is two-thirds of available CPU threads.
USAGE
}

usuage_and_exit(){
    usuage
    exit "$1"
}

error(){
    echo "$@" 1>&2
    usage_and_exit 1
}

error_exit() {
  ${2+:} false && echo "${0##*/}: $2" 1>&2
  exit $1
}

# ============================================================================
# Parse arguments
# ============================================================================
[ $# -eq 0 ] && usuage_and_exit 1

while [ $# -gt 0 ]
do
    case "$1" in
        --help | --hel | --he | --h | '--?' | -help | -hel | -he | -h | '-?')
            usuage_and_exit 0
            ;;
        --version | --versio | --versi | --vers | --ver | --ve | --v | \
        -version | -versio | -versi | -vers | -ver | -ve | -v )
            echo "DAJIN version: $VERSION"
            exit 0
            ;;
        --file | -f )
            fasta=$(cat "$2" | grep "design" | sed -e "s/ //g" -e "s/.*=//g")
            ont_dir=$(cat "$2" | grep "sequence" | sed -e "s/ //g" -e "s/.*=//g")
            ont_cont=$(cat "$2" | grep "control" | sed -e "s/ //g" -e "s/.*=//g")
            genome=$(cat "$2" | grep "genome" | sed -e "s/ //g" -e "s/.*=//g")
            grna=$(cat "$2" | grep "grna" | sed -e "s/ //g" -e "s/.*=//g")
            threads=$(cat "$2" | grep "threads" | sed -e "s/ //g" -e "s/.*=//g")
            ;;
        -* )
        error "Unrecognized option : $1"
            ;;
        *)
            break
            ;;
    esac
    shift
done

set +e

if [ -z "$fasta" ] || [ -z "$ont_dir" ] || [ -z "$ont_cont" ] || [ -z "$genome" ] || [ -z "$grna" ]
then
    error_exit 1 "Required arguments are not specified: See ./DAJIN.sh -h"
fi

if [ "$(grep -c '>target' ${fasta})" -eq 0 ] || [ "$(grep -c '>wt' ${fasta})" -eq 0 ]
then
    error_exit 1 "FASTA requires including \">target\" and \">wt\". See ./DAJIN.sh -h"
fi

[ -d "$ont_dir" ] ||
error_exit 1 "No such directory: See ./DAJIN.sh -h"

set -e

# ============================================================================
# Allocate threads
# ============================================================================

set +u
# Linux and similar...
[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/1.5+0.5)}')
# FreeBSD and similar...
[ -z "$threads" ] && threads=$(getconf NPROCESSORS_ONLN | awk '{print int($0/1.5)+0.5}')
# Solaris and similar...
[ -z "$threads" ] && threads=$(ksh93 -c 'getconf NPROCESSORS_ONLN' | awk '{print int($0/1.5+0.5)}')
# Give up...
[ -z "$threads" ] && threads=1
set -u


# ============================================================================
# Required software
# ============================================================================
set +e

type python 1>/dev/null 2>/dev/null || error_exit 1 'Command "python" not found'
type samtools 1>/dev/null 2>/dev/null || error_exit 1 'Command "samtools" not found'
type minimap2 1>/dev/null 2>/dev/null || error_exit 1 'Command "minimap2" not found'
type gzip 1>/dev/null 2>/dev/null || error_exit 1 'Command "gzip" not found'

python -c "import tensorflow as tf" \
1>/dev/null 2>/dev/null ||  error_exit 1 '"Tensorflow" not found'
set -e

# ============================================================================
# For WSL (Windows Subsystem for Linux)
# ============================================================================

uname -a | 
grep Microsoft 1>/dev/null 2>/dev/null &&
alias python="python.exe"