# Windows PCでのGPU環境構築
Windowsでは、WSL(Windows Subsystem for Linux)の設定とGPUを動かすための準備が必要になります。  

## WSLのインストール
[こちら](https://docs.microsoft.com/ja-jp/windows/wsl/install-win10)のページを参考にWSLをインストールします。  
インストールするOSは、`Ubuntu 18.04 LTS`がお勧めです。  

## Linux版Minicondaのインストール
WSLにLinux版のMinicondaをインストールします。  
WSLを起動し、以下のコマンドを実行します。
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

## Windous版Minicondaのインストール
[Miniconda3 Windows 64-bit](https://docs.conda.io/en/latest/miniconda.html#windows-installers)をインストールします。  

**【重要】インストールの際に表示される"Add Anaconda to my PATH environment variable"と  
"Registar Anaconda as my default Python XX"の２項目にチェックをつけてください。（下図）**  


<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/anaconda-install.png" width="50%">  

## 必要なソフトウェアのインストール
Anaconda/Minicondaをインストールした後、スタートメニューから`Anaconda Prompt`を検索して起動します。  
`Anaconda Prompt`を起動したのち、下記のコマンドを実行してください。  
```
:: Install Python packages (Tensorflow-GPU etc) in Windows environment.
conda update -y -n base conda
conda create -y -n DAJIN_win python=3.6 anaconda tensorflow-gpu keras tqdm
conda activate DAJIN_win
bash
# Install other packages (samtools etc) in Linux environment.
echo 'alias python="python.exe"' >> ~/.bashrc
source ~/.bashrc
conda update -y -n base conda
conda create -y -n DAJIN_wsl python=3.6 git nodejs
conda install -y -n DAJIN_wsl -c bioconda nanosim samtools htslib fasta3 clustalo weblogo
conda activate DAJIN_wsl
```

続けて以下のコマンドを実行し、使用しているGPUの情報が表示されれば完了です。  
```
python -c "from tensorflow.python.client import device_lib; 
print(device_lib.list_local_devices())"
```

２回目以降は`Anaconda Prompt`を立ち上げた後に下記のコマンドを実行します。  
```
conda activate DAJIN_win
bash
conda activate DAJIN_wsl
```

# 参考
Using WSL Linux on Windows 10 for Deep Learning Development.  
http://www.erogol.com/using-windows-wsl-for-deep-learning-development/  
