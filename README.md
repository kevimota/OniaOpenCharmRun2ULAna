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

# Running the Analysis

Events selection: `nanoAODplus_analyzer.py`, can be used with condor with `nanoAODplus_condor.py`. With the selection done, follow the steps

## Do trigger and more tight selection:

Run the code `nanoAODplus_trigger.py` by:
```
python nanoAODplus_trigger.py -p {path} -y {year} -m
```
path = path to the directory where the data is stored

years = ['2016APV', '2016', '2017', '2018']

The data can be plotted by using `nanoAODplus_plotter.py`
```
python nanoAODplus_plotter.py -y {year}
```

Save the data to ROOT TTrees to be able to use RooFit:
```
python tools/save_ttree.py -p {path_in} -o {path_out} -a
```

Fit to the 2D model:
```
python nanoAODplus_fit.py -y {year} -c {particle}
python nanoAODplus_fit.py -y {year} -c {particle} -p
```
particles = ['Upsilon', 'Dstar', 'UpsilonDstar']

## Calculate efficiency from MC

TODO