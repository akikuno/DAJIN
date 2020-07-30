
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">  
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

[日本語はこちら](https://github.com/akikuno/DAJIN/blob/master/misc/README_JP.md)

A simple, rapid, scalable whole-allelic profile of genome editing aminals using long-read sequencer

## SETUP

### Linux

#### 1. Install `conda`
The [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installation is required.

#### 2. Clone DAJIN repository

```
git clone https://github.com/akikuno/DAJIN.git
```

---
### Windows10
Windows user needs WSL and settings by Anaconda Prompt.  
See [this page](https://github.com/akikuno/DAJIN/blob/master/misc/WindowsOS_Setting.md).  

---
### macOS
macOS is not recommended because [Nvidia CUDA will not support it](https://docs.nvidia.com/cuda/cuda-installation-guide-mac-os-x/index.html).  
You can use DAJIN with CPU though long computational time.  

## USAGE

```sh
./DAJIN/DAJIN.sh -f [text file] (described at "Input")
```

### Example

```
./DAJIN/DAJIN.sh -f DAJIN/example/desing.txt
```

### Input file

Input file should be formatted like below:

```
design=DAJIN/example/design.txt
input_dir=DAJIN/example/demultiplex
control=barcode01
genome=mm10
grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
output_dir=Cables2
threads=10
```

- **desing**: a multi-FASTA file contains sequences of each genotype. ">wt" and ">target" must be included.
- **input_dir**: a directory contains FASTA or FASTQ files of long-read sequencing
- **control**: control barcode ID
- **genome**: reference genome. e.g. mm10, hg38
- **grna**: gRNA sequence(s). multiple gRNA sequences must be deliminated by comma.
- **output_dir** (optional): output directory name. default dirname is `DAJIN_results`
- **threads** (optional): default is two-thirds of available CPU threads.

### Output files

#### Details.csv

`Details.csv` contains sample information and 

|Sample  |Allele ID|% of reads|Allele type  |Indel|Large indel|Design|
|--------|---------|----------|-------------|-----|-----------|------|
|sample01|1        |4.2       |abnormal     |+    |+          |-     |
|sample01|2        |95.8      |wt           |-    |-          |-     |
|sample02|1        |13.2      |abnormal     |+    |+          |-     |
|sample02|2        |12        |abnormal     |+    |+          |-     |
|sample02|3        |74.8      |target       |-    |-          |+     |
|sample03|1        |55.2      |abnormal     |+    |+          |-     |
|sample03|2        |44.8      |flox_deletion|-    |-          |-     |

#### Details.pdf

The output directory contains a figure of whole-allelic profile.  
This is an example result of three samples.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/Details.png" width="75%">  

The sample01 is a wild-type mice as a control, whereas the sample02 and sample03 are founder mice. T target allele is flox.

This result shows ~75% of reads from sample02 are labeled as "target" (flox), and indicates it can be the desired mouse that has homozygous floxed allele.


#### Consensus

The `Conseusus` folder includes FASTA and HTML files which display conseusus sequence in each allele.


#### BAM


## Contact
- Akihiro Kuno akuno@md.tsukuba.ac.jp

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

- Dr. Seiya Mizuno
- Dr. Ayabe Sinya
- Mr. Yoshihisa Ikeda
- Mr. Kotaro Sakamoto
- Ms. Sayaka Suzuki
