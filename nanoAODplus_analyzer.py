# coding: utf-8

import time
import os, sys, subprocess

from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from coffea.util import save, load
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
parser.add_argument("-n", "--name", help="Analyser name", type=str, required=True)
parser.add_argument("-m","--merge", help="Merge the accumulators that were output from a analyzer", action="store_true")
args = parser.parse_args()

if args.merge:
    from tqdm import tqdm
    if (subprocess.run("find output/ -type d -name '" + args.name + "'", shell=True, stdout=subprocess.PIPE).stdout.decode("utf-8") == ''):
        raise Exception("Folder not found!")
    print("Merging files in output/" + args.name)
    files = subprocess.run("ls -d output/" + args.name + "/*", shell=True, stdout=subprocess.PIPE)
    file_list = files.stdout.decode("utf-8").splitlines()
    acc = load(file_list[0])
    for idx, f in tqdm(enumerate(file_list), desc="Merging", unit=" files", total=len(file_list)):
        if (idx == 0): continue
        acc += load(f)
    print("Saving as output/merged/" + args.name + "_merged.coffea")
    os.system("mkdir -p output/merged")
    save(acc, "output/merged/" + args.name + "_merged.coffea")
            
else:
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