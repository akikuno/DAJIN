#!/bin/bash

set +u

error_exit() {
    echo "$@" 1>&2
    exit 1
}

CONDA_BASE=$(conda info --base)
. "${CONDA_BASE}/etc/profile.d/conda.sh"

conda --version >/dev/null 2>&1 || error_exit 'Command "conda" not found'

conda update -y -n base conda >/dev/null 2>&1
conda config --add channels defaults 2>/dev/null
conda config --add channels bioconda 2>/dev/null
conda config --add channels conda-forge 2>/dev/null

conda list | grep -q mamba || conda install -y -c conda-forge mamba >/dev/null 2>&1

###############################################################################
# Setup DAJIN nanosim
###############################################################################

if [ "$(conda info -e | grep -c DAJIN_nanosim)" -eq 0 ]; then
    echo Create "DAJIN_nanosim" environment... >&2
    conda create -y -n DAJIN_nanosim python=3.6 >/dev/null 2>&1
    mamba install -y -n DAJIN_nanosim --file ./DAJIN/utils/NanoSim/requirements.txt >/dev/null 2>&1
    mamba install -y -n DAJIN_nanosim minimap2 >/dev/null 2>&1
fi

conda activate DAJIN_nanosim

if ! python ./DAJIN/utils/NanoSim/src/simulator.py --version >/dev/null 2>&1; then
    CONDA_ENV=$(conda info -e | awk '$2=="*"{print $NF}')
    (cd "${CONDA_ENV}"/lib/ && ln -s libcrypto.so.1.1 libcrypto.so.1.0.0 >/dev/null 2>&1)
fi

python ./DAJIN/utils/NanoSim/src/simulator.py --version >/dev/null 2>&1 ||
    error_exit '"NanoSim" installation has failed'
minimap2 --version >/dev/null 2>&1 || error_exit 'Command "minimap2" installation has failed'

rm -rf DAJIN/utils/NanoSim/src/__pycache__

conda deactivate

###############################################################################
# Setup DAJIN
###############################################################################

if [ "$(conda info -e | cut -d " " -f 1 | grep -c DAJIN$)" -eq 0 ]; then
    echo Create "DAJIN" environment... >&2
    conda create -y -n DAJIN python=3.8 >/dev/null 2>&1
    mamba install -y -n DAJIN \
        numpy pandas scikit-learn joblib hdbscan \
        wget emboss samtools minimap2 >/dev/null 2>&1
    mamba install -y -n DAJIN -c conda-forge r-essentials r-base r-reticulate r-vroom r-furrr >/dev/null 2>&1
    # tensorflow setting (CPU, GPU w/ RTX, and GPU w/ GTX)
    conda activate DAJIN
    if ! type nvidia-smi >/dev/null 2>&1; then
        mamba install -y -n DAJIN -c anaconda tensorflow >/dev/null 2>&1
    elif nvidia-smi -q | grep -q RTX; then
        mamba install -y -n DAJIN cudnn=8.0 cudatoolkit=11.0 >/dev/null 2>&1
        pip install -U tensorflow-gpu==2.4.0 >/dev/null 2>&1
    else
        mamba install -y -n DAJIN -c anaconda tensorflow tensorflow-gpu >/dev/null 2>&1
    fi
fi

conda activate DAJIN

Rscript -e 'install.packages(c("pacman", "tidyfast"), repos="https://cloud.r-project.org/")' >/dev/null 2>&1

###############################################################################
# Check prerequisites
###############################################################################

gzip --version >/dev/null 2>&1 || error_exit 'Command "gzip" installation has failed'
wget --version >/dev/null 2>&1 || error_exit 'Command "wget" installation has failed'
stretcher --version >/dev/null 2>&1 || error_exit 'Command "stretcher" installation has failed'
python --version >/dev/null 2>&1 || error_exit 'Command "python" installation has failed'
R --version >/dev/null 2>&1 || error_exit '"R" installation has failed'
minimap2 --version >/dev/null 2>&1 || error_exit 'Command "minimap2" installation has failed'

python -c "import tensorflow as tf" >/dev/null 2>&1 ||
    error_exit '"Tensorflow" not found'

tf_ver="$(conda list -n DAJIN | awk '$1~/tensorflow/ && $2>1.99')"
[ -z "$tf_ver" ] && error_exit '"Tensorflow 2.x" not found'

if ! samtools --version >/dev/null 2>&1; then
    CONDA_ENV=$(conda info -e | awk '$2=="*"{print $NF}')
    (cd "${CONDA_ENV}"/lib/ && ln -s libcrypto.so.1.1 libcrypto.so.1.0.0 >/dev/null 2>&1)
fi
samtools --version >/dev/null 2>&1 || error_exit 'Command "samtools" installation has failed'
