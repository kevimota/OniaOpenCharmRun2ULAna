# Analysis Onia + Open Charm

Repository for Quarkonia + Open Charm Analysis using Run2 UL datasets.

# Instalation instructions

## Setting up the conda environment

Miniconda should be installed as in:

https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

```
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b 
~/miniconda3/bin/conda init
echo 'conda deactivate ; conda deactivate' >> ~/.bashrc
rm Miniconda3-latest-Linux-x86_64.sh
```

After installation, reload the session.

## Setup conda Environmment (only once!)

```
conda deactivate ; conda deactivate 
conda create -y -n OniaOpenCharmRun2ULenv python=3.6 xrootd -c conda-forge

conda activate OniaOpenCharmRun2ULenv

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade --ignore-installed --force-reinstall coffea
python3 -m pip install --upgrade --ignore-installed --force-reinstall pick
python3 -m pip install --upgrade --ignore-installed --force-reinstall tabulate
python3 -m pip install --upgrade --ignore-installed --force-reinstall millify
python3 -m pip install --upgrade --ignore-installed --force-reinstall pyyaml
python3 -m pip install --upgrade --ignore-installed --force-reinstall jupyterlab

conda install -y gcc_linux-64 gxx_linux-64
python3 -m pip install --upgrade --ignore-installed --force-reinstall cppyy
python3 -m pip install --upgrade --ignore installed --force-reinstall nbresuse==0.3.3
conda install -y -c conda-forge nodejs

jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install jupyterlab-topbar-extension jupyterlab-system-monitor

#only if CMSSW needed:

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_12
cd CMSSW_10_6_12/src
cmsenv

#######################

mkdir plots
mkdir output

export PYTHONPATH=$CONDA_PREFIX/lib/python3.6/site-packages/
```

## Set the virtual environment after setup:

```
cd OniaOpenCharmRun2ULAna
. quick_setup.sh
```
