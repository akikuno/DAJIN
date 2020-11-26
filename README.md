
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

[日本語はこちら](https://github.com/akikuno/DAJIN/blob/master/misc/README_JP.md)

DAJIN is a genotyping software with simple, rapid and scalable whole-allelic profiling of genome editing organisms using long-read sequencer.  

Here are the DAJIN's features:

- DAJIN automatically identify and classify a diversity of mutations including point mutations, large deletions, inversions, and knock-in
- DAJIN can treat ~100 samples within a day

## Initial setup

We highly recommend Linux OS and NVIDIA GPU to reduce computation time.  
If you have a Windows PC with NVIDIA GPU, please follow the instruction.  
https://docs.nvidia.com/cuda/wsl-user-guide/index.html

FYI: We confirmed DAJIN's operation on [these environments](https://github.com/akikuno/DAJIN/blob/master/misc/TESTED_SYSTEMS.md).

### 1. Install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

```bash
# Instal git
sudo apt install git
# Install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```
### 2. Clone this repository

```sh
git clone https://github.com/akikuno/DAJIN.git
```

## Recommended directory tree

We recommend the following directory tree (described in next section).

```
├── DAJIN
├── input.txt
├── design.fasta
├── fastq
│   ├── barcode01.fastq
│   ├── barcode02.fastq
│   ├── barcode03.fastq
│   ├── ......

```
*You can arbitrary rename `input.txt`, `design.fasta` and `fastq`.

### 3. `input.txt`


`input.txt` should be formatted as below:

```
design=DAJIN/example/example.fa
input_dir=DAJIN/example/fastq
control=barcode01
genome=mm10
grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
output_dir=DAJIN_example
threads=10
```

- **design**: PATH to a multi-FASTA file including sequences of each genotype. **">wt" and ">target" must be included.**
- **input_dir**: PATH to a directory containing demultiplexed FASTQ files
- **control**: control barcode ID
- **genome**: reference genome. e.g. mm10, hg38
- **grna**: gRNA sequence(s) including PAM. multiple gRNA sequences must be delimitated by comma.
- **output_dir** (optional): output directory name. Default is `DAJIN_results`
- **threads** (optional: integer): Default is `2/3` of available CPU threads.
- **filter**  (optional: `on` or `off`): set filter to remove very minor allele (less than 3%). Default is `on`.

`design`, `input_dir`, `control`, `genome`,`grna` are required, but there are in no particular order.

### 2. `design.fasta`

`design.fasta` is a multi-FASTA file as following:

```
>wt
AAAAAA
>target
AATAAA
```
### 3. `fastq` directory

Currently DAJIN accepts [qcat](https://github.com/nanoporetech/qcat)'s demultiplex.  
We plan to update to accepts the output of guppy basecaller.

## USAGE

```sh
./DAJIN/DAJIN.sh -i design.txt
```

### Example usage

```sh
./DAJIN/DAJIN.sh -i DAJIN/example/design.txt
```
You can conduct DAJIN in example small dataset.

### Output files

DAJIN outputs two files and two folders: `Details.csv`, `Details.pdf`, `BAM`, `Consensus`.

#### Details.csv

`Details.csv` contains allele information.


| Sample    |  Allele ID |  % of reads |  Allele type  |  Indel |  Large indel |  Design |
|-----------|------------|-------------|---------------|--------|--------------|---------|
| barcode01 | 1          | 100         | wt            | -      | -            | -       |
| barcode02 | 1          | 11.8        | abnormal      | +      | +            | -       |
| barcode02 | 2          | 88.2        | target        | -      | -            | +       |
| barcode03 | 1          | 9.9         | abnormal      | +      | +            | -       |
| barcode03 | 2          | 38.5        | abnormal      | +      | +            | -       |
| barcode03 | 3          | 51.6        | flox_deletion | -      | -            | -       |

#### Details.pdf

The output directory contains a figure of whole-allelic profile.  
This is an example result of three samples.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/Details.png" width="75%">  

The barcode01 is a wild-type mice as a control, whereas the barcode02 and barcode03 are genome-edited founder mice with a flox knock-in design.

This result shows the most of Nanopore reads of barcode02 are labeled as "intact target" (flox), and indicates the barcode02 is a candidate of the desired homozygous mice.

#### Consensus

The `Conseusus` folder includes FASTA and HTML files which display consensus sequence in each allele.

Here is <a href="https://htmlpreview.github.io/?https://github.com/akikuno/DAJIN/blob/master/misc/images/tyr_c140cg.html" target= _blank rel= noopener> an example of DAJIN consensus sequence</a> using the point mutation.

#### BAM

The `BAM` folder includes BAM files from all and each allele.

The `BAM` files can be visualized by [IGV](http://software.broadinstitute.org/software/igv/).

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

- Dr. Seiya Mizuno
- Dr. Sinya Ayabe
- Mr. Yoshihisa Ikeda
- Mr. Kotaro Sakamoto
- Ms. Sayaka Suzuki
