# coding: utf-8

import time
import os, sys, subprocess

import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

from nanoAODplus_processor.EventSelectorProcessor import EventSelectorProcessor
from data.fileset import filesets
import yaml

# for argument parsing
import argparse
parser = argparse.ArgumentParser(description="Onia Open Charm NanoAOD analyzer")
parser.add_argument("-n", "--name", help="Analyzer name", type=str, required=True)
parser.add_argument("-s", "--select", help="Do the evt selection", action="store_true")
parser.add_argument("-m","--merge", help="Merge the accumulators that were output from a analyzer", action="store_true")
parser.add_argument("-p","--plots", help="Create the plots from the merged accumulator", action="store_true")
parser.add_argument("-a","--analyze", help="Do the full analysis chain", action="store_true")
args = parser.parse_args()

if (args.select or args.analyze):
    config_yaml = yaml.load(open("config/local.yaml", "r"), Loader=yaml.FullLoader)

    if config_yaml['executor'] == 'futures_executor': 
        executor = processor.futures_executor

    tstart = time.time()

    files = {'MuOnia2017MINIAOD': filesets['MuOnia2017MINIAOD'][0:5]}

    # creating necessary folders into dir output data
    os.system("mkdir -p output/" + args.name)
    os.system("rm -rf output/" + args.name + "/*")          

    output = processor.run_uproot_job(files,
                                    treename='Events',
                                    processor_instance=EventSelectorProcessor(args.name),
                                    executor=executor,
                                    executor_args={'workers': config_yaml['n_cores'], 'flatten': True},
                                    chunksize=config_yaml['chunksize'],
                                    )

    elapsed = round(time.time() - tstart, 2)
    print(f"Process finished in: {elapsed} s")

if (args.merge or args.analyze):
    from tools.merger import merger
    merger(args.name)

if (args.plots or args.analyze):
    from tools.plotter import plotter
    plotter(args.name)    