#!/bin/bash

set +u

error_exit() {
    echo "$@" 1>&2
    exit 1
}

conda --version > /dev/null || error_exit 'Command "conda" not found'

#===========================================================
#? DAJIN_nanosim
#===========================================================

CONDA_BASE=$(conda info --base)
. "${CONDA_BASE}/etc/profile.d/conda.sh"

if [ "$(conda info -e | grep -c DAJIN_nanosim)" -eq 0 ]; then
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda update -y conda
    conda create -y -n DAJIN_nanosim python=3.6
    conda install -y -n DAJIN_nanosim --file ./DAJIN/utils/NanoSim/requirements.txt
    conda install -y -n DAJIN_nanosim minimap2
fi

conda activate DAJIN_nanosim

python ./DAJIN/utils/NanoSim/src/simulator.py --version >/dev/null 2>&1 ||
error_exit '"NanoSim" installation has failed'

rm -rf DAJIN/utils/NanoSim/src/__pycache__

#===========================================================
#? DAJIN
#===========================================================

if [ "$(conda info -e | cut -d " " -f 1 | grep -c DAJIN$)" -eq 0 ]; then
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda update -y conda
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

gzip --version > /dev/null || error_exit 'Command "gzip" installation has failed'
wget --version > /dev/null || error_exit 'Command "wget" installation has failed'
python --version > /dev/null || error_exit 'Command "python" installation has failed'
R --version > /dev/null || error_exit 'Command "Rscript" installation has failed'
minimap2 --version > /dev/null || error_exit 'Command "minimap2" installation has failed'

python -c "import tensorflow as tf" > /dev/null ||
error_exit '"Tensorflow" not found'

if samtools --version 2>&1 | grep libcrypto >/dev/null; then
    CONDA_ENV=$(conda info -e | awk '$2=="*"{print $NF}')
    (cd "${CONDA_ENV}"/lib/ && ln -s libcrypto.so.1.1 libcrypto.so.1.0.0)
fi
samtools --version > /dev/null || error_exit 'Command "samtools" installation has failed'