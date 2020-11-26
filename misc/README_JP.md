
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

DAJINはゲノム編集生物の遺伝型判別器です. 以下の特徴があります.  

- 多くのゲノム編集デザインに対応しています. 点変異, ノックアウト, ノックイン, 逆位の解析が可能です.
- Nanoporeロングリードシークエンサーによって従来法に比べ広範囲(~10 kb)のゲノム編集領域を探索できます
- 一塩基変異から数kbにおよぶindelを同定できます
- 100サンプル程度までのサンプルについて結果を高速に（GPU環境ならば1日以内で）レポートします

## 推奨環境

計算時間を短くするためにLinux OSとNvidia GPUが使える環境をおすすめします.  
CPUでも実行可能ですが, サンプル数によってはかなり長時間（週単位）がかかってしまう可能性があります.  

以下の環境で動作確認をしています.  

- Ubuntu 18.04, RTX2080ti
- Linux Mint 20.04, GTX1080ti
- Windows10 WSL2 (Ubuntu 18.04) * CPUでの実行

より詳細は計算環境は[こちら](https://github.com/akikuno/DAJIN/blob/master/misc/TESTED_SYSTEMS.md)です.  

> なお, macOSは未検証です.  

## セットアップ

### 1. [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)と[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)をインストールします

```bash
# Instal git
sudo apt install git
# Install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```
### 2. DAJINをダウンロードします

```
git clone https://github.com/akikuno/DAJIN.git
```

## 推奨のディレクトリ構成について

下記のようなディレクトリ構成を推奨しています.  

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

> `input.txt`, `design.fasta` と `fastq`の名前は自由に変更していただけます

ファイル/ディレクトリの説明は以下の通りです.  
### 1. `input.txt`

`input.txt`は以下のような形式のテキストファイルです.  

```
design=DAJIN/example/design.txt
input_dir=DAJIN/example/demultiplex
control=barcode01
genome=mm10
grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
output_dir=DAJIN_cables2
threads=10
filter=on
```

各項目の内容は以下のとおりです。

- **desing**: 考えられる遺伝型の配列を記載したFASTA形式のテキストファイルです。 ">wt"と">target"の2つは含まれている必要があります。
- **input_dir**: demultiplex済みのFASTA/FASTQファイルを含むディレクトリです。
- **control**: 野生型コントロールのバーコード番号です。
- **genome**: `mm10`, `hg38`等の参照ゲノムです。
- **grna**: gRNA配列です。2つ以上の配列はコンマ（,）で区切ります。
- **output_dir（オプショナル）**: 結果を保存するディレクトリの名前です。デフォルトは`DAJIN_results`です。
- **threads（オプショナル）**: DAJINに使用するCPUスレッド数です。デフォルトでは`3分の2`を使用します。
- **filter（オプショナル**: on/off）: マイナーアレル（Targetアレルが1%以下、その他のアレルが3%以下）を解析から除きます。デフォルトは"on"です。

> `design`, `input_dir`, `control`, `genome`,`grna` は必須項目です.
> 各項目は順不同です。

### 2. `design.fasta`

`design.fasta` はマルチFASTA形式のテキストファイルです. WT（ゲノム編集前）とTarget（ゲノム編集後）の配列は必ず含む必要があります.  ほかにも副産物として起こりうる配列を加えることができます.  

floxノックインの例は[こちら](https://github.com/akikuno/DAJIN/blob/master/example/example.fa)です.  
floxノックインの場合は副産物アレルを含めて6つのアレルタイプが考えられます(WT, Target, Left LoxP, Right LoxP, flox deletion, Inversion).  
また, DAJINは`design.fasta`の配列とは違ったアレルを'異常アレル abnormal'としてレポートします.  

### 3. `fastq` directory

現状, DAJINは[qcat](https://github.com/nanoporetech/qcat)によりDemultiplexされたfastqをあつかいます.  
そのためGuppyの場合には出力のディレクトリ構成が異なるため, 変換する必要があります.  
Guppyの出力を扱えるようにアップデートする予定です.  

## DAJINの実行

```bash
./DAJIN/DAJIN.sh -i design.txt
```
### Example usage

```sh
./DAJIN/DAJIN.sh -i DAJIN/example/design.txt
```
:point_up_2:小さいデータセットでDAJINを試すことができます.  


### 出力ファイル

DAJINは解析後, 2つのファイルと2つのディレクトリを作製します. それぞれ`Details.csv`, `Details.pdf`, `Consensus`, `BAM`です.  

#### Details.csv

`Details.csv` は各個体に含まれるアレルの種類と割合が示されています.  
目的の変異をもつアレルはDesignの列が`+`となっています.  


| Sample    |  Allele ID |  % of reads |  Allele type  |  Indel |  Large indel |  Design |
|-----------|------------|-------------|---------------|--------|--------------|---------|
| barcode01 | 1          | 100         | wt            | -      | -            | -       |
| barcode02 | 1          | 11.8        | abnormal      | +      | +            | -       |
| barcode02 | 2          | 88.2        | target        | -      | -            | +       |
| barcode03 | 1          | 9.9         | abnormal      | +      | +            | -       |
| barcode03 | 2          | 38.5        | abnormal      | +      | +            | -       |
| barcode03 | 3          | 51.6        | flox_deletion | -      | -            | -       |

#### Details.pdf

`Details.pdf`は`Details.csv`を可視化したものです.  
下記のような図になります.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/Details.png" width="75%">  

barcode01は野生型コントロールで, barcode02および03がfloxノックインのマウスゲノムです. この図からbarcode02のゲノムはほぼfloxアレルであることから, 目的のノックインがホモで入った個体であると考えられます.  

#### Consensusフォルダ

`Conseusus`フォルダにはFASTAファイルまたはHTMLファイルがあり, 各アレルごとのコンセンサス配列が記載されています.  

とくにHTMLファイルは変異が色付けされているため一瞥して変異箇所とその種類（挿入・欠失・置換）が理解できます.  
<a href="https://htmlpreview.github.io/?https://github.com/akikuno/DAJIN/blob/master/misc/images/tyr_c140cg.html" target= _blank rel= noopener>こちらの一例</a>は点変異のコンセンサス配列です.  

#### BAMフォルダ

`BAM`フォルダには各サンプルごとのBAMファイルがあり, さらに1サンプルすべてのリードか, 各アレルごとのリードのBAMに分かれて保存されています.  
[IGV](http://software.broadinstitute.org/software/igv/)によって可視化できます.  

## ライセンス

DAJINはMITライセンスです. 詳細は[LICENSE](https://github.com/akikuno/DAJIN/blob/master/LICENSE)をご覧ください.

## 謝辞

- 水野聖哉 先生（筑波大学 実験動物学研究室）
- 綾部 信哉 先生（理研バイオリソース研究センター 実験動物開発室）
- 池田 祥久 さん（筑波大学 実験動物学研究室）
- 坂本 航太郎 さん（筑波大学 ヒューマンバイオロジー学位プログラム）
- 鈴木 沙耶香 さん（筑波大学 ヒューマンバイオロジー学位プログラム）
