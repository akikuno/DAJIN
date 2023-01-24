#!/bin/bash

set +u

error_exit() {
    echo "$@" 1>&2
    exit 1
}

###############################################################################
# Update conda
###############################################################################

echo "Update Conda..." >&2

CONDA_BASE=$(conda info --base)
. "${CONDA_BASE}/etc/profile.d/conda.sh"

conda --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "conda" not found'

conda config --add channels defaults >>.DAJIN_temp/log.txt 2>&1
conda config --add channels bioconda >>.DAJIN_temp/log.txt 2>&1
conda config --add channels conda-forge >>.DAJIN_temp/log.txt 2>&1

conda list | grep -q mamba || conda install -y -c conda-forge mamba >>.DAJIN_temp/log.txt 2>&1
mamba update -y -n base conda mamba >>.DAJIN_temp/log.txt 2>&1

###############################################################################
# Setup DAJIN nanosim
###############################################################################

if ! conda info -e | grep -q DAJIN_nanosim; then
    echo Create "DAJIN_nanosim" environment... >&2
    conda create -y -n DAJIN_nanosim python=3.8 >>.DAJIN_temp/log.txt 2>&1
    mamba install -y -n DAJIN_nanosim nanosim samtools==1.10 minimap2 >>.DAJIN_temp/log.txt 2>&1
fi

conda activate DAJIN_nanosim

if ! python ./DAJIN/utils/NanoSim/src/simulator.py --version >>.DAJIN_temp/log.txt 2>&1; then
    CONDA_ENV=$(conda info -e | awk '$2=="*"{print $NF}')
    (
        cd "${CONDA_ENV}"/lib/
        libcrypto=$(ls -l libcrypto.so | awk '{print $NF}')
        ln -sf "$libcrypto" libcrypto.so.1.0.0 >>.DAJIN_temp/log.txt 2>&1
    )
fi

python ./DAJIN/utils/NanoSim/src/simulator.py --version >>.DAJIN_temp/log.txt 2>&1 ||
    error_exit '"NanoSim" installation has failed'
minimap2 --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "minimap2" installation has failed'

rm -rf DAJIN/utils/NanoSim/src/__pycache__

conda deactivate

###############################################################################
# Setup DAJIN
###############################################################################

if ! conda info -e | cut -d " " -f 1 | grep -q "^DAJIN$"; then
    echo Create "DAJIN" environment... >&2
    mamba create -y -n DAJIN >>.DAJIN_temp/log.txt 2>&1
    # tensorflow setting (CPU, GPU w/ RTX, and GPU w/ GTX)
    conda activate DAJIN
    echo "install packages -------------------" >>.DAJIN_temp/log.txt
    if ! type nvidia-smi >>.DAJIN_temp/log.txt 2>&1; then
        mamba install -y -n DAJIN -c anaconda tensorflow \
            numpy pandas scikit-learn joblib hdbscan wget emboss minimap2 \
            r-essentials r-base r-reticulate >>.DAJIN_temp/log.txt 2>&1
    elif nvidia-smi -q | grep -q RTX; then
        mamba install -y -n DAJIN -c conda-forge tensorflow=2.10.0=gpu_py39h039f4ff_0 \
            numpy pandas scikit-learn joblib hdbscan wget emboss minimap2 \
            r-essentials r-base r-tidyverse r-reticulate >>.DAJIN_temp/log.txt 2>&1
    else
        mamba install -y -n DAJIN -c anaconda tensorflow tensorflow-gpu \
            numpy pandas scikit-learn joblib hdbscan wget emboss minimap2 \
            r-essentials r-base r-reticulate >>.DAJIN_temp/log.txt 2>&1
    fi
    echo "install samtools -------------------" >>.DAJIN_temp/log.txt
    mamba install -y -c conda-forge -n DAJIN samtools >>.DAJIN_temp/log.txt 2>&1
fi

conda activate DAJIN

# Install R packages
if ! Rscript -e "installed.packages()" >/dev/null 2>&1 | grep -q pacman; then
    Rscript -e 'install.packages("pacman", repos="https://cloud.r-project.org/")' >>.DAJIN_temp/log.txt 2>&1
fi
Rscript -e 'pacman::p_load("RColorBrewer", "vroom", "furrr", "tidyfast")' >>.DAJIN_temp/log.txt 2>&1

###############################################################################
# Check prerequisites
###############################################################################

gzip --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "gzip" installation has failed'
wget --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "wget" installation has failed'
stretcher --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "stretcher" installation has failed'
python --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "python" installation has failed'
R --version >>.DAJIN_temp/log.txt 2>&1 || error_exit '"R" installation has failed'
minimap2 --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "minimap2" installation has failed'

python -c "import tensorflow as tf" >>.DAJIN_temp/log.txt 2>&1 ||
    error_exit '"Tensorflow" not found'

tf_ver="$(conda list -n DAJIN | awk '$1~/tensorflow/ && $2>1.99')"
[ -z "$tf_ver" ] && error_exit '"Tensorflow 2.x" not found'

if ! samtools --version >>.DAJIN_temp/log.txt 2>&1; then
    CONDA_ENV=$(conda info -e | awk '$2=="*"{print $NF}')
    (
        cd "${CONDA_ENV}"/lib/
        libcrypto=$(ls -l libcrypto.so | awk '{print $NF}')
        ln -sf "$libcrypto" libcrypto.so.1.0.0 >>.DAJIN_temp/log.txt 2>&1
    )
fi
samtools --version >>.DAJIN_temp/log.txt 2>&1 || error_exit 'Command "samtools" installation has failed'
