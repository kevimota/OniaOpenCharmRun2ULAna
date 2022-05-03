# Analysis Onia + Open Charm

Repository for Quarkonia + Open Charm Analysis using Run2 UL datasets.

# Instalation instructions

## Setting up the conda environment

Using Miniforge as python package manager:

https://github.com/conda-forge/miniforge

To install in a UNIX-like platform:

```
cd ~
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b 
conda config --set auto_activate_base false
rm Miniforge3-$(uname)-$(uname -m).sh
```

After installation, reload the session.

## Setup conda Environmment (only once!)

```
conda deactivate ; conda deactivate 
conda env create -f environment.yml
conda activate OniaOpenCharmRun2ULenv

#only if CMSSW needed:

#export SCRAM_ARCH=slc7_amd64_gcc700

#cmsrel CMSSW_10_6_12
#cd CMSSW_10_6_12/src
#cmsenv

#######################

mkdir plots
mkdir output
```

## Set the virtual environment after setup:

```
cd OniaOpenCharmRun2ULAna
. quick_setup.sh
```
