import time
import os, sys, re

from concurrent.futures import ProcessPoolExecutor
from tqdm.auto import tqdm

import awkward as ak
import numpy as np
from coffea.util import load, save
from coffea import processor

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

from tools.utils import *

import yaml

class Skimmer:
    def __init__(self, config, year):
        if config is None or year is None:
            print("Configuration not complete!!!")
        self.config = config
        self.year = year

    def process(self, f):
        # Load the accumulator
        acc = load(f)
        
        Dimu_acc = acc['Dimu']
        Dstar_acc = acc['Dstar']
        Dstar_D0_acc = acc['Dstar_D0']
        DimuDstar_acc = acc['DimuDstar']
        triggers = acc['triggers']

        # Loading variables as awkward arrays
        Dimu = ak.zip({
                'pt': Dimu_acc['pt'].value,
                'eta': Dimu_acc['eta'].value,
                'phi': Dimu_acc['phi'].value,
                'mass': Dimu_acc['mass'].value,
                'rap': Dimu_acc['rap'].value,
                'dl': Dimu_acc['dl'].value,
                'dlSig': Dimu_acc['dlSig'].value,
                'chi2': Dimu_acc['chi2'].value,
                'cosphi': Dimu_acc['cosphi'].value,
                'is_ups': Dimu_acc['is_ups'].value,}, with_name="PtEtaPhiMLorentzVector")

        Dstar = ak.zip({
                'pt': Dstar_acc['pt'].value,
                'eta': Dstar_acc['eta'].value,
                'phi': Dstar_acc['phi'].value,
                'mass': Dstar_acc['mass'].value,
                'rap': Dstar_acc['rap'].value,
                'charge': Dstar_acc['charge'].value,
                'deltam': Dstar_acc['deltam'].value,
                'deltamr': Dstar_acc['deltamr'].value,
                'D0cosphi': Dstar_D0_acc['D0cosphi'].value,
                'D0dlSig': Dstar_D0_acc['D0dlSig'].value,
                'D0pt': Dstar_D0_acc['D0pt'].value,
                'D0eta': Dstar_D0_acc['D0eta'].value,
                'wrg_chg': Dstar_acc['wrg_chg'].value,}, with_name='PtEtaPhiMCandidate')

        DimuDstar_p4 = build_p4(DimuDstar_acc)
        DimuDstar = ak.zip({
                'pt': DimuDstar_p4.pt,
                'eta': DimuDstar_p4.eta,
                'phi': DimuDstar_p4.phi,
                'mass': DimuDstar_p4.mass,
                'charge': DimuDstar_acc['charge'].value,
                'dimu_mass': DimuDstar_acc['Dimu']['mass'].value,
                'dimu_pt': DimuDstar_acc['Dimu']['pt'].value,
                'dimu_eta': DimuDstar_acc['Dimu']['eta'].value,
                'dimu_phi': DimuDstar_acc['Dimu']['phi'].value,
                'dimu_rap': DimuDstar_acc['Dimu']['rap'].value,
                'dstar_deltam': DimuDstar_acc['Dstar']['deltam'].value,
                'dstar_deltamr': DimuDstar_acc['Dstar']['deltamr'].value,
                'dstar_pt': DimuDstar_acc['Dstar']['pt'].value,
                'dstar_eta': DimuDstar_acc['Dstar']['eta'].value,
                'dstar_phi': DimuDstar_acc['Dstar']['phi'].value,
                'dstar_rap': DimuDstar_acc['Dstar']['rap'].value,
                'dstar_d0_pt': DimuDstar_acc['Dstar']['D0pt'].value,
                'dstar_d0_eta': DimuDstar_acc['Dstar']['D0eta'].value,
                'deltarap': DimuDstar_acc['deltarap'].value,
                'deltapt': DimuDstar_acc['deltapt'].value,
                'deltaeta': DimuDstar_acc['deltaeta'].value,
                'deltaphi': DimuDstar_acc['deltaphi'].value,
                'is_ups': DimuDstar_acc['Dimu']['is_ups'].value,
                'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,}, with_name='PtEtaPhiMCandidate')

        # Unflatten
        Dimu = ak.unflatten(Dimu, Dimu_acc['nDimu'].value)
        Dstar = ak.unflatten(Dstar, Dstar_acc['nDstar'].value)
        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)

        #Apply the Trigger and the cuts in the config
        HLT = triggers[self.config['trigger'][self.year]].value
        Dimu = Dimu[HLT]
        Dstar = Dstar[HLT]
        DimuDstar = DimuDstar[HLT]
        
        Dimu = Dimu[Dimu.pt > self.config['limits']['Upsilon_pt']]
        Dimu = Dimu[np.absolute(Dimu.eta) < self.config['limits']['Upsilon_eta']]
        Dimu = Dimu[Dimu.is_ups]
        DimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_pt']]
        DimuDstar = DimuDstar[np.absolute(DimuDstar.dimu_eta) < self.config['limits']['Upsilon_eta']]
        DimuDstar = DimuDstar[DimuDstar.is_ups]

        Dstar = Dstar[Dstar.pt > self.config['limits']['Dstar_pt']]
        Dstar = Dstar[np.absolute(Dstar.eta) < self.config['limits']['Dstar_eta']]
        Dstar = Dstar[Dstar.D0pt > self.config['limits']['Dstar_D0_pt']]
        Dstar = Dstar[np.absolute(Dstar.D0eta) < self.config['limits']['Dstar_D0_eta']]
        Dstar = Dstar[~Dstar.wrg_chg]
        DimuDstar = DimuDstar[DimuDstar.dstar_pt > self.config['limits']['Dstar_pt']]
        DimuDstar = DimuDstar[np.absolute(DimuDstar.dstar_eta) < self.config['limits']['Dstar_eta']]
        DimuDstar = DimuDstar[DimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
        DimuDstar = DimuDstar[np.absolute(DimuDstar.dstar_d0_eta) < self.config['limits']['Dstar_D0_eta']]
        DimuDstar = DimuDstar[~DimuDstar.wrg_chg]

        filename = f"{f[:f.rfind('/')]}/skim{f[f.rfind('/'):]}"
        
        pDimu_acc = build_acc(Dimu)
        pDstar_acc = build_acc(Dstar)
        pDimuDstar_acc = build_acc(DimuDstar)

        output = processor.dict_accumulator({
            'Dimu': pDimu_acc,
            'Dstar': pDstar_acc,
            'DimuDstar': pDimuDstar_acc
        })

        save(output, filename)

def get_coffea_files(path):
    files = []
    with os.scandir(path) as it:
        for f in it:
            if f.name.endswith('.coffea') and (f.name.find('_hists') == -1) and (f.stat().st_size != 0):
                files.append(f.path)

    files.sort(key=natural_keys)
    return files

def merger(path):
    files = get_coffea_files(f"{path}/skim")
    for idx, f in tqdm(enumerate(files), desc="Merging", unit=" files", total=len(files)):
        if (idx == 0): 
            acc = load(f)
        else:
            acc += load(f)
        os.system("rm -rf " + f)

    filename = files[0][:files[0].rfind("_")] + ".coffea"
    save(acc, filename)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Apply trigger and skim dataset")
    parser.add_argument("-p", "--path", help="Analyzer name", type=str, required=True)
    parser.add_argument("-y", "--year", help="Year of the dataset", type=str, required=True)
    args = parser.parse_args()

    years = ['2016', '2017', '2018']
    if args.year not in years:
        print("Year not in Run II, exiting...")
        sys.exit()

    config_run = yaml.load(open("config/multicore.yaml", "r"), Loader=yaml.FullLoader)
    config_trigger = yaml.load(open("config/skim_trigger.yaml", "r"), Loader=yaml.FullLoader)

    os.system(f"mkdir -p {args.path}/skim")
    os.system(f"rm -rf {args.path}/skim/*")

    files = get_coffea_files(args.path)  

    if len(files) == 0:
        print("No files were found in the specified directory. Exiting...")
        sys.exit()
    
    s = Skimmer(config_trigger, args.year)
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

    print("Starting merging process...")
    merger(args.path)

    print("DONE!")