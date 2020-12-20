
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

[日本語はこちら](https://github.com/akikuno/DAJIN/blob/master/misc/README_JP.md)

DAJIN is a genotyping software with simple, rapid, and scalable whole-allelic profiling of genome editing organisms using a long-read sequencer.  

:dart: DAJIN automatically identify and classify a diversity of mutations, including a point mutation, knock-out, knock-in, and inversion  
:dart: DAJIN uses the nanopore long-read sequencer to capture larger genomic locus (~10 kb) than conventional genotyping methods such as a short-read NGS and Sanger sequencing  
:dart: DAJIN detects various mutations from single-nucleotide variants to large indels (~10 kb)  
:dart: DAJIN can treat ~100 samples within a day  

## Initial Setup

We recommend **Linux OS with NVIDIA GPU** to reduce computation time.  
If you use a Windows PC with NVIDIA GPU, please follow [this instruction](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).  

> We confirmed DAJIN's operation on [these environments](https://github.com/akikuno/DAJIN/blob/master/misc/TESTED_SYSTEMS.md).

### 1. Install conda

```bash
# Install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```
### 2. Clone this repository

```sh
git clone https://github.com/akikuno/DAJIN.git
```
### 3. Prepare design files and fastq directory

We recommend the following directory tree.

```
├── DAJIN
├── design.txt
├── design.fasta
├── fastq
│   ├── barcode01.fastq
│   ├── barcode02.fastq
│   ├── barcode03.fastq
│   ├── ......

```
> You can rename `design.txt`, `design.fasta` and `fastq`.

Descriptions of the files/directory are as following:  
#### 1. `design.txt`


`design.txt` is formatted as below:

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
- **grna**: gRNA sequence(s) *including PAM*. multiple gRNA sequences must be delimitated by comma.
- **output_dir** (optional): output directory name. Default is `DAJIN_results`
- **threads** (optional: integer): Default is `2/3` of available CPU threads.
- **filter**  (optional: `on` or `off`): set filter to remove very minor allele (less than 3%). Default is `on`.

> `design`, `input_dir`, `control`, `genome`,`grna` are required, but there are in no particular dorder.

#### 2. `design.fasta`

`design.fasta` is a multi-FASTA file, which contain a WT and target sequence, as well as byproducts.

[This](https://github.com/akikuno/DAJIN/blob/master/example/example.fa) is an example file of flox design.  
In the case of flox design, 6 allele types (WT, Target, Left LoxP, Right LoxP, flox deletion, and Inversion) can be produced.  
Besides, DAJIN annotates an allele that are different from these allele types as `abnormal` allele.  

#### 3. `fastq` directory

Currently DAJIN accepts [qcat](https://github.com/nanoporetech/qcat)'s demultiplex.  
We plan to update to accepts the output of `guppy` basecaller.

## Usage

```sh
./DAJIN/DAJIN.sh -i design.txt
```

### Example usage

```sh
./DAJIN/DAJIN.sh -i DAJIN/example/design.txt
```
You can conduct DAJIN in example small dataset.

## Results

DAJIN outputs two files and two folders: `Details.csv`, `Details.pdf`, `Consensus`, `BAM`.  

### Details.csv

`Details.csv` contains allele information.
The allele with target mutation is labeled **+** in Design column.

| Sample    |  Allele ID |  % of reads |  Allele type  |  Indel |  Large indel |  Design |
|-----------|------------|-------------|---------------|--------|--------------|---------|
| barcode01 | 1          | 100         | wt            | -      | -            | -       |
| barcode02 | 1          | 11.8        | abnormal      | +      | +            | -       |
| barcode02 | 2          | 88.2        | target        | -      | -            | +       |
| barcode03 | 1          | 9.9         | abnormal      | +      | +            | -       |
| barcode03 | 2          | 38.5        | abnormal      | +      | +            | -       |
| barcode03 | 3          | 51.6        | flox_deletion | -      | -            | -       |

### Details.pdf

The output directory contains a figure of whole-allelic profile.  
This is an example result of three samples.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/Details.png" width="75%">  

The barcode01 is a wild-type mice as a control, whereas the barcode02 and barcode03 are genome-edited founder mice with a flox knock-in design.

This result shows the most of Nanopore reads of barcode02 are labeled as "intact target" (flox), and indicates the barcode02 is a candidate of the desired homozygous mice.

### Consensus

The `Conseusus` folder includes FASTA and HTML files which display consensus sequence in each allele.

Here is <a href="https://htmlpreview.github.io/?https://github.com/akikuno/DAJIN/blob/master/misc/images/tyr_c140cg.html" target= _blank rel= noopener> an example of DAJIN consensus sequence</a> using the point mutation.

### BAM

The `BAM` folder includes BAM files from all and each allele.

The `BAM` files can be visualized by [IGV](http://software.broadinstitute.org/software/igv/).

## License

This project is under the MIT License - see the [LICENSE](https://github.com/akikuno/DAJIN/blob/master/LICENSE) file for details

## Acknowledgments

- Dr. Seiya Mizuno at University of Tsukuba
- Dr. Sinya Ayabe at Riken BioResource Research Center
- Mr. Yoshihisa Ikeda at University of Tsukuba
- Mr. Kotaro Sakamoto at University of Tsukuba
- Ms. Sayaka Suzuki at University of Tsukuba
