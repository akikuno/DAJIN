
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

## 推奨環境

以下の環境で動作確認をしています。
- Ubuntu 18.04
- Linux Mint 20.04
- Windows10 WSL2 (Ubuntu 18.04)

## インストール

### `git`と`conda`のインストール

- [git](https://git-scm.com/book/ja/v2/%E4%BD%BF%E3%81%84%E5%A7%8B%E3%82%81%E3%82%8B-Git%E3%81%AE%E3%82%A4%E3%83%B3%E3%82%B9%E3%83%88%E3%83%BC%E3%83%AB)
- [conda](https://docs.conda.io/en/latest/miniconda.html)

##### DAJINをダウンロードします

```
git clone https://github.com/akikuno/DAJIN.git
```

### Windows10

Windows10の場合は[WSL2](https://docs.microsoft.com/ja-jp/windows/wsl/install-win10)が必要です。

---

## 利用方法

### 入力ファイルの用意

以下のようなテキストファイルを作製します。

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

各項目の情報は以下のとおりです。

- **desing**: 考えられる遺伝型の配列を記載したFASTA形式のテキストファイルです。 ">wt"と">target"の2つは含まれている必要があります。
- **input_dir**: demultiplex済みのFASTA/FASTQファイルを含むディレクトリです。
- **control**: 野生型コントロールのバーコード番号です。
- **genome**: `mm10`, `hg38`等の参照ゲノムです。
- **grna**: gRNA配列です。2つ以上の配列はコンマ（,）で区切ります。
- **output_dir（オプショナル）**: 結果を保存するディレクトリの名前です。デフォルトは`DAJIN_results`です。
- **threads（オプショナル）**: DAJINに使用するCPUスレッド数です。デフォルトでは`3分の2`を使用します。
- **filter（オプショナル**: on/off）: マイナーアレル（Targetアレルが1%以下、その他のアレルが3%以下）を解析から除きます。デフォルトは"on"です。

＊各項目は順不同です。


### DAJINの実行

```bash
./DAJIN/DAJIN.sh -f [入力ファイルのPATH]
```

（以下、工事中です…:construction_worker:）