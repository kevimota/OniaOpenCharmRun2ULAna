# coding: utf-8

import time
import os, sys

from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

from nanoAODplus_processor.EventSelectorProcessor import EventSelectorProcessor
from data.fileset import filesets
import yaml

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

# for argument parsing
import argparse
parser = argparse.ArgumentParser(description="Onia Open Charm NanoAOD analyzer")
parser.add_argument("--name", help="Analyser name", type=str)
args = parser.parse_args()

config_yaml = yaml.load(open("config/local.yaml", "r"), Loader=yaml.FullLoader)

if config_yaml['executor'] == 'futures_executor': 
    executor = processor.futures_executor

tstart = time.time()

files = {'Charmonium2017MINIAOD': filesets['Charmonium2017MINIAOD'][0:1]}

# creating necessary folders into dor output data
os.system("mkdir -p output/" + args.name)
os.system("rm -rf output/" + args.name + "/*")          

output = processor.run_uproot_job(files,
                                  treename='Events',
                                  processor_instance=EventSelectorProcessor(args.name),
                                  executor=executor,
                                  executor_args={'workers': config_yaml['n_cores'], 'flatten': True},
                                  chunksize=config_yaml['chunksize'],
                                 )

elapsed = time.time() - tstart

print("Time elapsed:", elapsed)