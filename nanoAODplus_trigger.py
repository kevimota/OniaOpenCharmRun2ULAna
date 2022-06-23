import time
import os, sys, re

from concurrent.futures import ProcessPoolExecutor
from tqdm.auto import tqdm

import awkward as ak
from coffea.util import load, save

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

from nanoAODplus_processor.Skimmers import Skimmer, FOM

from tools.utils import *

import yaml

def get_coffea_files(path):
    files = []
    with os.scandir(path) as it:
        for f in it:
            if f.name.endswith('.coffea') and (f.name.find('_hists') == -1) and (f.stat().st_size != 0):
                files.append(f.path)

    files.sort(key=natural_keys)
    return files

def merger(path):
    files = get_coffea_files(path)
    for idx, f in tqdm(enumerate(files), desc="Merging", unit=" files", total=len(files)):
        if (idx == 0): 
            acc = load(f)
        else:
            acc += load(f)
        os.system("rm -rf " + f)

    filename = files[0][:files[0].rfind("_")] + ".coffea"
    #filename = f'output/RunII_trigger_processed/{year}/MuOniaRun{year}.coffea'
    save(acc, filename)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Apply trigger and skim dataset")
    parser.add_argument("-p", "--path", help="Path to coffea files", type=str, required=True)
    parser.add_argument("-y", "--year", help="Year of the dataset", type=str, required=True)
    parser.add_argument("-m", "--merge", help="Merge output", action="store_true")
    parser.add_argument("-f", "--fom", help="create the FOM datasets", action="store_true")
    args = parser.parse_args()

    years = ['2016', '2017', '2018']
    if args.year not in years:
        print("Year not in Run II, exiting...")
        sys.exit()

    config_run = yaml.load(open("config/multicore.yaml", "r"), Loader=yaml.FullLoader)
    config_trigger = yaml.load(open("config/skim_trigger.yaml", "r"), Loader=yaml.FullLoader)
    if args.fom:
        config_fom = yaml.load(open("config/fom.yaml", "r"), Loader=yaml.FullLoader)

    folder = args.path[args.path.rfind('/'):args.path.rfind('_')]
    if args.fom: path = f"output/fom/{args.year}{folder}"
    else: path = f"output/RunII_trigger_processed/{args.year}{folder}"

    os.system(f"mkdir -p {path}")
    os.system(f"rm -rf {path}/*")

    if args.fom:
        for var in config_fom:
            if var == 'path': continue
            for lim in config_fom[var]:
                os.system(f"mkdir -p {path}/{var}/{str(lim).replace('.', 'p')}")

    files = get_coffea_files(args.path)  

    if len(files) == 0:
        print("No files were found in the specified directory. Exiting...")
        sys.exit()
    
    if args.fom: s = FOM(config_trigger, config_fom, args.year)
    else: s = Skimmer(config_trigger, args.year)
    print("Starting processing the files...")
    
    tstart = time.time()

    if config_run['executor'] == 'iterative_executor' or len(files) == 1:
        for f in tqdm(files, total=len(files), unit=" files", desc="Processing"):
            s.process(f)
    else:
        with ProcessPoolExecutor(max_workers=config_run['n_cores']) as executor:
            list(tqdm(executor.map(s.process, files), total=len(files), unit=" files", desc="Processing"))

    elapsed = round(time.time() - tstart, 2)
    print(f"Skimming finished in: {elapsed} s")

    if (args.merge):
        print("Starting merging process...")
        if args.fom: 
            for var in config_fom:
                if var == 'path': continue
                print(f"Merging process: {var}")
                for lim in config_fom[var]:
                    l = str(lim).replace('.', 'p')
                    merger(f'{path}/{var}/{l}')
        else: merger(path)

    print("DONE!")