
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">  
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

[日本語はこちら](https://github.com/akikuno/DAJIN/blob/master/misc/README_JP.md)

A simple, rapid, scalable whole-allelic profile of genome editing organisms using long-read sequencer

## SETUP

### Linux and Windows10 ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10))

#### 1. Install `git` and `conda`

- [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

or

```sh
# to instal git
sudo apt install git
# to install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

#### 2. Clone DAJIN repository

```sh
git clone https://github.com/akikuno/DAJIN.git
```


<!-- ---
### macOS

macOS is not recommended because [Nvidia CUDA will not support it](https://docs.nvidia.com/cuda/cuda-installation-guide-mac-os-x/index.html).  
You can use DAJIN with CPU though long computational time.   -->

## USAGE

```sh
./DAJIN/DAJIN.sh -i [text file] (described at "Input")
```

### Example usage

```sh
./DAJIN/DAJIN.sh -i DAJIN/example/design.txt
```

### Input file

Input file should be formatted as below:

```
design=DAJIN/example/design.txt
input_dir=DAJIN/example/fastq
control=barcode01
genome=mm10
grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
output_dir=DAJIN_example
threads=10
```

- **desing**: a multi-FASTA file contains sequences of each genotype. ">wt" and ">target" must be included.
- **input_dir**: a directory contains FASTA or FASTQ files of long-read sequencing
- **control**: control barcode ID
- **genome**: reference genome. e.g. mm10, hg38
- **grna**: gRNA sequence(s) including PAM. multiple gRNA sequences must be delimitated by comma.
- **output_dir** (optional): output directory name. default dirname is `DAJIN_results`
- **threads** (optional; integer): default is two-thirds of available CPU threads.

### Output files

DAJIN outputs two files and two folders: `Details.csv`, `Details.pdf`, `BAM`, `Consensus`.

#### Details.csv

`Details.csv` contains allele information.


|Sample   |Allele ID|% of reads|Allele type  |Indel|Large indel|Design|
|---------|---------|----------|-------------|-----|-----------|------|
|barcode01|1        |6.2       |abnormal     |+    |+          |-     |
|barcode01|2        |93.8      |wt           |-    |-          |-     |
|barcode02|1        |100.0     |target       |-    |-          |+     |
|barcode03|1        |38.5      |abnormal     |+    |+          |-     |
|barcode03|2        |13.4      |abnormal     |+    |+          |-     |
|barcode03|3        |7.7       |flox_deletion|+    |-          |-     |
|barcode03|4        |40.4      |flox_deletion|-    |-          |-     |

#### Details.pdf

The output directory contains a figure of whole-allelic profile.  
This is an example result of three samples.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/Details.png" width="75%">  

The barcode01 is a wild-type mice as a control, whereas the barcode02 and barcode03 are founder mice. The target allele is a flox.

This result shows all of Nanopore reads from barcode02 are labeled as "intact target" (flox), and indicates it is the desired sample that has homozygous flox allele.

#### Consensus

The `Conseusus` folder includes FASTA and HTML files which display consensus sequence in each allele.

> in preparation

#### BAM

The `BAM` folder includes BAM files from all and each allele.

> in preparation

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

- Dr. Seiya Mizuno
- Dr. Sinya Ayabe
- Mr. Yoshihisa Ikeda
- Mr. Kotaro Sakamoto
- Ms. Sayaka Suzuki
