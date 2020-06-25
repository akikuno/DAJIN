#!/bin/bash


conda info -e
conda activate DAJIN_py36
conda info -e
conda deactivate
conda info -e

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [ "$(conda info -e | grep -c DAJIN_nanosim)" -eq 0 ]; then
conda create -y -n DAJIN_nanosim python=3.6
conda install -y -n DAJIN_nanosim -f DAJIN/utils/NanoSim/requirements.txt
conda install -y -n DAJIN_nanosim minimap2
fi

conda activate DAJIN_nanosim

conda activate DAJIN