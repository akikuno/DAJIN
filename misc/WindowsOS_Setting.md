# Environment setting for Windows PC
Windows PC requires to set-up WSL (Windows Subsystem for Linux) and Anaconda Prompt to execute DAJIN.  

## Install WSL
[This page](https://docs.microsoft.com/en-us/windows/wsl/install-win10) guides you to install WSL.  
I recommend you to install `Ubuntu 18.04 LTS` on WSL.  

## Install Miniconda to WLS

Install Miniconda of Linux version to WSL.  
Open WSL and execute the following command.  
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

## Install Miniconda to Windows OS
Get [Miniconda3 Windows 64-bit](https://docs.conda.io/en/latest/miniconda.html#windows-installers) and install Miniconda。  

**[IMPORTANT!!] Check both "Add Anaconda to my PATH environment variable" and  
"Registar Anaconda as my default Python XX".**  


<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/anaconda-install.png" width="50%">  

## Set-up environment
Type `Anaconda Prompt` in the search window and open it.  
Then execute the following commands.  
```
# Install Python packages (Tensorflow-GPU etc) in Windows environment.
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

Lastly, when you see the local GPU device by the following command, that is it!
```
python -c "from tensorflow.python.client import device_lib;
print(device_lib.list_local_devices())"
```

From now on, you can execute the following command in `Anaconda Prompt` to conduct DAJIN.  
```
conda activate DAJIN_win
bash
conda activate DAJIN_wsl
```

# Reference
Using WSL Linux on Windows 10 for Deep Learning Development.  
http://www.erogol.com/using-windows-wsl-for-deep-learning-development/  
