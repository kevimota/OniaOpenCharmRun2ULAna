# coding: utf-8

import time
import os, sys

import coffea.processor as processor
from coffea.nanoevents import BaseSchema

from nanoAODplus_processor.EventSelectorProcessor import EventSelectorProcessor
from tools.merger import merger
from tools.plotter import plotter

years = ['2016', '2017', '2018']

tstart = time.time()

job = sys.argv[1]
number = job.split('.')[-2].split('_')[-1]

with open(job) as f:
    file_list = f.read().splitlines()

name = file_list[0].split('/')[-1].split('_')[0]
year = '2017'
for i in years:
    if i in name:
        year = i

name = name + '_' + number

files = {name: file_list}

print(name, year)

# creating necessary folders into dir output data
os.system("mkdir -p output/" + name)
os.system("rm -rf output/" + name + "/*")          

runner = processor.Runner(
                executor=processor.IterativeExecutor(compression=None),
                schema=BaseSchema,
                skipbadfiles=True, 
                chunksize=10000
            )

output = runner(files,
                treename='Events',
                processor_instance=EventSelectorProcessor(name, year),
                )

elapsed = round(time.time() - tstart, 2)
print(f"Process finished in: {elapsed} s")

merger(name)

plotter(name)

print("Transfering files...")
os.system("cp output/" + name + "/* ../.")
