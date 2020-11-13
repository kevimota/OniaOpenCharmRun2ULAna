# coding: utf-8

import time
import os, sys, subprocess

import coffea.processor as processor
import numpy as np

from nanoAODplus_processor.HistogramingProcessor import HistogramingProcessor
def plotter(name):
    print("Starting plots creation")

    print("Saving histograms in output/" + name + "/hist")
    os.system("mkdir -p output/" + name + "/hist")

    print("Saving plots in plots/" + name)
    os.system("mkdir -p plots/" + name)
    os.system("rm -rf output/" + name + "/hist/*")
    os.system("rm -rf plots/" + name + "/*")

    tstart = time.time()

    _processor = HistogramingProcessor()
    _dummy_accumulator = _processor.accumulator.identity()
    ds = [{"file": "output/" + name + "/merged/" + name + "_merged.coffea", "analyzer_name": name}]
    processor.iterative_executor(ds, _processor.process, _dummy_accumulator)

    elapsed = round(time.time() - tstart, 2)

    print(f"Finished in: {elapsed} s")