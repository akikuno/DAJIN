# miniogenotype
A simple, rapid, scalable whole-allelic profile of genome editing aminals using ONT MinION

# Installation
## Linux
The latest [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installation is highly recommended.  

```
conda create -y -n miniogenotype python=3.6 anaconda git tensorflow-gpu keras tqdm nodejs
conda install -y -n miniogenotype -c bioconda nanosim samtools htslib fasta3 clustalo weblogo
conda activate miniogenotype
git clone https://github.com/akikuno/miniogenotype.git
```
## Windows10
WSLとAnaconda Promptの設定が必要です。  
[こちらのページ](https://github.com/akikuno/miniogenotype/blob/master/misc/WindowsOS_Setting_JP.md)をご覧ください。  

# Usuage
```
allele_profiler.sh \
  -i <FASTA file> (required) \
  -ont_dir <directory> (required)  \
  -ont_ref <FASTA|FASTQ file> (required) \
  -genome <reference genome name> (required) \
  -o <string> (optional) \
  -t <integer> (optional)

Options :   
-i            Multi-FASTA file. It must include ">target" and ">wt"
-ont_dir      Directory containing ONT demultiplexed reads 
-ont_ref      Reference (= wild-type) reads from ONT MinION
-genome       Reference genome name (e.g. hg38, mm10)
              See https://gggenome.dbcls.jp/en/help.html#db_list 
-o            Output file name (default = sequence_MIDS)
-t            Number of threads to NanoSim and BAM compression
              (default =  half of total available thread in CPU)
-h            show the help message
```

# Example

```
./miniogenotype/allele_profiler.sh \
  -i miniogenotype/example/cables2_flox.fa \
  -ont miniogenotype/example/demultiplex \
  -ont_ref miniogenotype/example/demultiplex/barcode21.fastq.gz \
  -genome mm10 \
  -o test \
  -t 8
```

# Output
## Whole-allelic profile
`results` directory contains a figure of whole-allelic profile.  
This is an example result of three mice.  

<img src="https://github.com/akikuno/miniogenotype/blob/images/sequence_MIDS_prediction_result.png" width="50%">  

Barcode14 and 19 are a founder mice, whose target allele is flox. Barcode21 is a wild-type mice as a control.   
This result shows ~80% of reads from Barcode14 are labeled as "target" (flox), and indicates Barcode14 is the desired mouse that has homozygous floxed allele.

## Alignment viewing using IGV.js
When you want to see the alignment of reads, you can type the following command.  
```
npx live-server results/igvjs/
```
The browser will pop-up the following page:  

<img src="https://github.com/akikuno/miniogenotype/blob/images/igvjs_localhost.png" width="50%">  

Click `igvjs.html` and you can see the alignment views:  

<img src="https://github.com/akikuno/miniogenotype/blob/images/igvjs_alignment.png" width="50%">

The barcode14 has two purple sites, where **insertion** occurs.

# Contact
Akihiro Kuno  
akuno@md.tsukuba.ac.jp

# Reference
