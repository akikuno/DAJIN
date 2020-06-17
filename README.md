
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">  
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

[日本語はこちら](https://github.com/akikuno/DAJIN/blob/master/misc/README_JP.md)
# DAJIN
A simple, rapid, scalable whole-allelic profile of genome editing aminals using ONT MinION

# Installation
## Linux
The latest [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installation is highly recommended.  

### 1. Setup channels and Install required packages

```
conda config --add channels defaults &&
conda config --add channels bioconda &&
conda config --add channels conda-forge &&
conda update -y -n base conda &&
conda create -y -n DAJIN python=3.6 \
  anaconda git nodejs \
  tensorflow tensorflow-gpu \
  samtools minimap2 \
  r-essentials r-base
```
### 2. Activate the environment
```
conda activate DAJIN
```
### 3. Clone DAJIN repository
```
git clone https://github.com/akikuno/DAJIN.git
```
You need only `2. Activate the environment` from the second time on.

---
## Windows10
Windows user needs WSL and settings by Anaconda Prompt.  
See [this page](https://github.com/akikuno/DAJIN/blob/master/misc/WindowsOS_Setting.md).  

---
## macOS
macOS is not recommended because [Nvidia CUDA will not support it](https://docs.nvidia.com/cuda/cuda-installation-guide-mac-os-x/index.html).  
You can use DAJIN with CPU whereas long computational time.  

# Usuage
```
Usage     : DAJIN.sh -f [text file] (described at "Input")

Input     : Input file should be formatted as below:
            --------------------------------
            design=DAJIN/example/design.txt
            input_dir=DAJIN/example/demultiplex
            control=barcode01
            output_dir=Cables2
            genome=mm10
            grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
            threads=10
            --------------------------------
            - desing: a multi-FASTA file contains sequences of each genotype. ">wt" and ">target" must be included. 
            - input_dir: a directory contains FASTA or FASTQ files of long-read sequencing
            - control: control barcode ID
            - output_dir: output directory name. optional. default is DAJIN_results
            - genome: reference genome. e.g. mm10, hg38
            - grna: gRNA sequence(s). multiple gRNA sequences must be deliminated by comma.
            - threads: optional. default is two-thirds of available CPU threads.
```

# Example

```
./DAJIN/DAJIN.sh -f DAJIN/example/example.txt
```

# Output
## Whole-allelic profile
`results` directory contains a figure of whole-allelic profile.  
This is an example result of three mice.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/sequence_MIDS_prediction_result.png" width="50%">  

Barcode14 and 19 are a founder mice, whose target allele is flox. Barcode21 is a wild-type mice as a control.   
This result shows ~80% of reads from Barcode14 are labeled as "target" (flox), and indicates Barcode14 is the desired mouse that has homozygous floxed allele.

## Alignment viewing using IGV.js
When you want to see the alignment of reads, you can type the following command.  
```
npx live-server results/igvjs/
```
The browser will pop-up the following page:  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/igvjs_localhost.png" width="50%">  

Click `igvjs.html` and you can see the alignment views:  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/igvjs_alignment.png" width="50%">

The barcode14 has two purple sites, where **insertion** occurs.

# Contact
Akihiro Kuno  
akuno@md.tsukuba.ac.jp

# Reference
