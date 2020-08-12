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
conda config --add channels defaults 2>/dev/null
conda config --add channels bioconda 2>/dev/null
conda config --add channels conda-forge 2>/dev/null

if [ "$(conda info -e | grep -c DAJIN_nanosim)" -eq 0 ]; then
    echo Create "DAJIN_nanosim" environment...
    conda update -y conda 2>/dev/null
    conda create -y -n DAJIN_nanosim python=3.6 2>/dev/null
    conda install -y -n DAJIN_nanosim --file ./DAJIN/utils/NanoSim/requirements.txt 2>/dev/null
    conda install -y -n DAJIN_nanosim minimap2 2>/dev/null
fi

conda activate DAJIN_nanosim

python ./DAJIN/utils/NanoSim/src/simulator.py --version >/dev/null 2>&1 ||
error_exit '"NanoSim" installation has failed'
minimap2 --version > /dev/null || error_exit 'Command "minimap2" installation has failed'

rm -rf DAJIN/utils/NanoSim/src/__pycache__

conda deactivate

#===========================================================
#? DAJIN
#===========================================================

if [ "$(conda info -e | cut -d " " -f 1 | grep -c DAJIN$)" -eq 0 ]; then
    echo Create "DAJIN" environment...
    conda update -y conda >/dev/null 2>&1
    conda create -y -n DAJIN python=3.7 \
        numpy pandas scikit-learn scipy \
        tensorflow tensorflow-gpu \
        wget samtools minimap2 \
        r-essentials r-base r-dbscan >/dev/null 2>&1
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