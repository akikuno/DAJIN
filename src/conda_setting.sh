#!/bin/bash

set +u

error_exit() {
    echo "$@" 1>&2
    exit 1
}

type conda > /dev/null 2>&1 || error_exit 'Command "conda" not found'

#===========================================================
#? DAJIN_nanosim
#===========================================================

CONDA_BASE=$(conda info --base)
. "${CONDA_BASE}/etc/profile.d/conda.sh"

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
        anaconda wget \
        tensorflow tensorflow-gpu \
        samtools minimap2 \
        r-essentials r-base r-dbscan
fi

#===========================================================
#? Required software
#===========================================================

conda activate DAJIN

type gzip > /dev/null 2>&1 || error_exit 'Command "gzip" not found'
type wget > /dev/null 2>&1 || error_exit 'Command "wget" not found'
type python > /dev/null 2>&1 || error_exit 'Command "python" not found'
type samtools > /dev/null 2>&1 || error_exit 'Command "samtools" not found'
type minimap2 > /dev/null 2>&1 || error_exit 'Command "minimap2" not found'

python -c "import tensorflow as tf" > /dev/null 2>&1 ||
error_exit '"Tensorflow" not found'