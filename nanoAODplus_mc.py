# coding: utf-8

import time
import os, sys

import coffea.processor as processor
from coffea.nanoevents import BaseSchema

from nanoAODplus_processor.GenParticleProcessor import GenParticleProcessor
import yaml

from tools.utils import get_files

base_folder = "/Users/kevimota/cernbox/CRAB_UserFiles"

if __name__ == '__main__':
    # for argument parsing
    import argparse
    parser = argparse.ArgumentParser(description="Onia Open Charm MC analyzer")
    parser.add_argument("-n", "--name", help="Analyzer name", type=str, required=True)
    args = parser.parse_args()

    path = f"{base_folder}/{args.name}"
    if not os.path.isdir(path):
        print(f"Folder {path} does not exist")
        sys.exit(1)

    if '2016' in args.name:
        if '2016APV' in args.name: year = '2016APV'
        else: year = '2016'
    elif '2017' in args.name: year = '2017'
    elif '2018' in args.name: year = '2018'
    else:
        print("Year not identified!!")
        sys.exit(1)

    config_yaml = yaml.load(open("config/multicore.yaml", "r"), Loader=yaml.FullLoader)

    tstart = time.time()

    files = {args.name: get_files([path])}
    print(files)

    # creating necessary folders into dir output data
    os.system("mkdir -p output/" + args.name)
    os.system("rm -rf output/" + args.name + "/*")       

    if config_yaml['executor'] == 'futures_executor': 
        runner = processor.Runner(
            executor=processor.FuturesExecutor(compression=None, workers=config_yaml['n_cores']),
            schema=BaseSchema,
            skipbadfiles=True, 
            chunksize=config_yaml['chunksize']
        )

    elif config_yaml['executor'] == 'iterative_executor':
        runner = processor.Runner(
            executor=processor.IterativeExecutor(compression=None),
            schema=BaseSchema,
            skipbadfiles=True, 
            chunksize=config_yaml['chunksize']
        )

    output = runner(
        files,
        treename="Events",
        processor_instance=GenParticleProcessor(args.name, year)
    )

    elapsed = round(time.time() - tstart, 2)
    print(f"Process finished in: {elapsed} s")
