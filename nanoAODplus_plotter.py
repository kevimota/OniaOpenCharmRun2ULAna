import os, sys

from tqdm.auto import tqdm
from coffea.util import load, save

from nanoAODplus_processor.HistogramingProcessor import HistogrammingProcessor

from tools.utils import get_files, get_lumi, get_trigger
from tools.figure import create_plot1d

import matplotlib.pyplot as plt

save_figures = 'plots/variables'

hist_log = ['pt', 'chi2', 'dl', 'dlSig']

def plotter(hists, year, processed_lumi):
    fig, ax = plt.subplots()

    save_folder = f'{save_figures}/{year}'
    print(f'Creating plots in {save_folder}')
    if not os.path.exists(f'{save_folder}/Dimu'): os.makedirs(f'{save_folder}/Dimu')
    if not os.path.exists(f'{save_folder}/Dstar'): os.makedirs(f'{save_folder}/Dstar')
    if not os.path.exists(f'{save_folder}/DimuDstar'): os.makedirs(f'{save_folder}/DimuDstar')

    for p in hists:
        for h in hists[p]:
            log = False
            for i in hist_log:
                if h.find(i) > -1:
                    log = True
                    break
            if log:
                create_plot1d(hists[p][h], ax=ax, log=True, lumi=processed_lumi)
                fig.savefig(f'{save_folder}/{p}/{h}_log_{year}.png')
                ax.clear()
            create_plot1d(hists[p][h], ax=ax, log=False, lumi=processed_lumi)
            fig.savefig(f'{save_folder}/{p}/{h}_{year}.png')
            ax.clear()
            

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Apply trigger and skim dataset")
    parser.add_argument("-y", "--year", help="Year of the dataset", type=str, required=True)
    args = parser.parse_args()

    years = ['2016APV', '2016', '2017', '2018']
    if args.year not in years:
        print("Year not in Run II, exiting...")
        sys.exit()

    folder = f'output/RunII_trigger_processed_vtxfit/{args.year}'
    folders = []
    for it in os.scandir(folder):
        if not it.is_dir(): continue
        folders.append(it.path)
        
    files = get_files(folders, pattern='.coffea')

    h = HistogrammingProcessor()
    print(f"Creating histogram files for year {args.year}")
    for f in tqdm(files, total=len(files), unit=" files", desc="Processing"):
        h.process(f)

    print('Done!')

    hists = get_files(folders, pattern='.hists')
    hists = [load(h) for h in hists]

    processed_lumi = get_lumi(args.year, get_trigger(args.year))

    for i in range(len(hists)):
        if i == 0:
            hs_Dimu = hists[i]['Dimu']
            hs_Dstar = hists[i]['Dstar']
            hs_DimuDstar = hists[i]['DimuDstar']
        else:
            hs_Dimu = {k: hs_Dimu.get(k, 0) + hists[i]['Dimu'].get(k, 0) for k in set(hs_Dimu)}
            hs_Dstar = {k: hs_Dstar.get(k, 0) + hists[i]['Dstar'].get(k, 0) for k in set(hs_Dstar)}
            hs_DimuDstar = {k: hs_DimuDstar.get(k, 0) + hists[i]['DimuDstar'].get(k, 0) for k in set(hs_DimuDstar)}
    
    hists = {
        'Dimu': hs_Dimu,
        'Dstar': hs_Dstar,
        'DimuDstar': hs_DimuDstar,
    }

    plotter(hists, args.year, processed_lumi)