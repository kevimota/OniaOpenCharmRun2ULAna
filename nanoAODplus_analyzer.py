# coding: utf-8

import time
import os

import coffea.processor as processor
from coffea.nanoevents import BaseSchema

from nanoAODplus_processor.EventSelectorProcessor import EventSelectorProcessor
#from data.fileset import filesets
import yaml

if __name__ == '__main__':
    # for argument parsing
    import argparse
    parser = argparse.ArgumentParser(description="Onia Open Charm NanoAOD analyzer")
    parser.add_argument("-n", "--name", help="Analyzer name", type=str, required=True)
    parser.add_argument("-y", "--year", help="Year of the dataset", type=str, required=True)
    parser.add_argument("-s", "--select", help="Do the evt selection", action="store_true")
    parser.add_argument("-m","--merge", help="Merge the accumulators that were output from a analyzer", action="store_true")
    parser.add_argument("-a","--analyze", help="Do the full analysis chain", action="store_true")
    args = parser.parse_args()

    if (args.select or args.analyze):
        config_yaml = yaml.load(open("config/multicore.yaml", "r"), Loader=yaml.FullLoader)

        tstart = time.time()

        files = {'MuOniatestAOD': ['MuOniaRun2018A_AOD_994.root']}

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
            processor_instance=EventSelectorProcessor(args.name, args.year)
        )

        elapsed = round(time.time() - tstart, 2)
        print(f"Process finished in: {elapsed} s")

    if (args.merge or args.analyze):
        from tools.merger import merger
        merger(args.name)
  
