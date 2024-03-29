<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">
</p>

****

> [!CAUTION]
> ### DAJIN is deprecated. Please use [DAJIN2](https://github.com/akikuno/DAJIN2), which has become the successor to DAJIN.

[日本語はこちら](https://github.com/akikuno/DAJIN/blob/master/misc/README_JP.md)

DAJIN is a genotyping software for genome-edited organisms using Nanopore sequencer.  
DAJIN is named after 一網**打尽** (Ichimou **DAJIN** or Yī Wǎng **Dǎ jìn**), meaning catching all in one net.  

## Features

- Capturing mutations from SNV to SV, such as a point mutation, knock-out, knock-in, and inversion.  
- Automatic allele clustering and annotation.
- Genotyping ~100 samples with different genome-editing designs at a single run.  

## Setup

We recommend **Linux OS with NVIDIA GPU** to reduce computation time.  
If you use a Windows PC with NVIDIA GPU, please follow [this instruction](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).  

> We confirmed DAJIN's operation on [these environments](https://github.com/akikuno/DAJIN/blob/master/misc/TESTED_SYSTEMS.md).

### 1. Install conda

```bash
# Install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local/
```

### 2. Clone this repository

```sh
git clone https://github.com/akikuno/DAJIN.git
```
### 3. Prepare design files and fastq directory

We recommend the following directory tree.
You can rename `design.txt`, `design.fasta` and `fastq`.

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


## Prepare inputs

### design.txt

`design.txt` is formatted as below:

```
design=DAJIN/example/example.fa
input_dir=DAJIN/example/fastq
control=barcode01
grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
genome=mm10
output_dir=DAJIN_example
threads=10
```

- **design**: PATH to a FASTA file including sequences of each genotype. **">wt" and ">target" must be included.**
- **input_dir**: PATH to a directory containing demultiplexed FASTQ files
- **control**: control barcode ID
- **grna**: gRNA sequence(s) **with PAM**. multiple gRNA sequences must be delimitated by comma.
- **genome**(optional): reference genome. e.g. hg38,mm10 (**not GRCh38,GRCm38**)
- **output_dir** (optional): output directory name. Default is `DAJIN_results`
- **threads** (optional: integer): Default is `2/3` of available CPU threads.
- **filter**  (optional: `on` or `off`): set filter to remove very minor allele (less than 3%). Default is `on`.

> `design`, `input_dir`, `control`, `genome`,`grna` are required, but there are in no particular order.

### design.fasta

`design.fasta` is a multi-FASTA file, which contains a WT and target sequence, as well as byproducts.

[This](https://github.com/akikuno/DAJIN/blob/master/example/example.fa) is an example file of flox design.  
In flox design, 6 allele types (WT, Target, Left LoxP, Right LoxP, flox deletion, and Inversion) may be produced.  
DAJIN defines `SV (structural variants)` allele that is different from these alleles.  


## Usage

```sh
./DAJIN/DAJIN -i design.txt
```

### Example usage

You can test DAJIN by an example small dataset.

```sh
./DAJIN/DAJIN -i DAJIN/example/design.txt
```

## Outputs

DAJIN outputs two files and two folders: `Details.csv`, `Details.pdf`, `Consensus`, `BAM`.  

### Details.csv

`Details.csv` contains allele information.
The allele with target mutation is labeled **+** in the Design column.

| Sample    | Allele ID | % of reads | Allele type   | Indel | Large indel | Design |
| --------- | --------- | ---------- | ------------- | ----- | ----------- | ------ |
| barcode01 | 1         | 100        | wt            | -     | -           | -      |
| barcode02 | 1         | 11.8       | SV            | +     | +           | -      |
| barcode02 | 2         | 88.2       | target        | -     | -           | +      |
| barcode03 | 1         | 9.9        | SV            | +     | +           | -      |
| barcode03 | 2         | 38.5       | SV            | +     | +           | -      |
| barcode03 | 3         | 51.6       | flox_deletion | -     | -           | -      |

### Details.pdf

The output directory contains a figure of whole-allelic profile.  
This is an example result of three samples.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/Details.png" width="75%">  

The barcode01 is a wild-type mouse as a control, whereas the barcode02 and barcode03 are genome-edited founder mice with a flox knock-in design.

This result shows that most Nanopore reads of barcode02 are labeled as "intact target" (flox), and indicates that barcode02 is a candidate for the desired homozygous mice.

### Consensus

The `Conseusus` folder includes FASTA and HTML files that display the consensus sequence in each allele.

Here is <a href="https://htmlpreview.github.io/?https://github.com/akikuno/DAJIN/blob/master/misc/images/tyr_c140cg.html" target= _blank rel= noopener> an example of DAJIN's consensus sequence</a> using the point mutation.

### BAM

The `BAM` folder includes BAM files from all and each allele.

The `BAM` files can be visualized by [IGV](http://software.broadinstitute.org/software/igv/).

## License

This project is under the MIT License - see the [LICENSE](https://github.com/akikuno/DAJIN/blob/master/LICENSE) file for details

## Citation

[PLOS BIOLOGY](https://doi.org/10.1371/journal.pbio.3001507)

```
@article{Kuno_2022,
	title={DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes},
	volume={20},
	url={https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507},
	DOI={10.1371/journal.pbio.3001507},
	number={1},
	journal={PLOS Biology},
	author={Kuno, Akihiro and Ikeda, Yoshihisa and Ayabe, Shinya and Kato, Kanako and Sakamoto, Kotaro and Suzuki, Sayaka R. and Morimoto, Kento and Wakimoto, Arata and Mikami, Natsuki and Ishida, Miyuki and et al.},
	year={2022},
	month={Jan},
	pages={e3001507}
}
```


‌
