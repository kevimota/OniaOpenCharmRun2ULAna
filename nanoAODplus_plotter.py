from concurrent.futures import ProcessPoolExecutor
import os
import yaml

from tqdm.auto import tqdm
from coffea.util import load

from nanoAODplus_processor.HistogramingProcessor import HistogrammingProcessor

from tools.utils import get_files, get_lumi, get_trigger
from tools.figure import create_plot1d

import matplotlib.pyplot as plt

hist_log = ['pt', 'chi2', 'dl', 'dlSig']

exclude = [
    'UpsilonPt9To30ToMuMuDstarToD0pi_2017-v2_18.root',
]

def plotter(hists, year, processed_lumi, save_folder, is_data=True):
    fig, ax = plt.subplots()

    print(f'Creating plots in {save_folder}')
    for p in hists:
        folder = f'{save_folder}/{p}'
        if not os.path.exists(folder): os.makedirs(folder)

    for p in hists:
        for h in hists[p]:
            log = False
            for i in hist_log:
                if h.find(i) > -1:
                    log = True
                    break
            if log:
                create_plot1d(hists[p][h], ax=ax, log=True, lumi=processed_lumi, is_data=is_data)
                fig.savefig(f'{save_folder}/{p}/{h}_log_{year}.png')
                ax.clear()
            create_plot1d(hists[p][h], ax=ax, log=False, lumi=processed_lumi, is_data=is_data)
            fig.savefig(f'{save_folder}/{p}/{h}_{year}.png')
            ax.clear()
            
years = ['2016APV', '2016', '2017', '2018']

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Apply trigger and skim dataset")
    parser.add_argument("-y", "--year", help="Year of the dataset", type=str, required=True, choices=years)
    parser.add_argument("-s", "--sel", help="Selection", type=str, required=False, default='sel', choices=["sel", "raw"])
    parser.add_argument("-m", "--is_mc", help="is mc", action="store_true")
    args = parser.parse_args()

    config_run = yaml.load(open("config/multicore.yaml", "r"), Loader=yaml.FullLoader)

    if not args.is_mc:
        if args.sel == 'sel':
            folder = f'output/RunII_trigger_processed_vtxfit/{args.year}'
            folders = []
            for it in os.scandir(folder):
                if not it.is_dir(): continue
                folders.append(it.path)
        elif args.sel == 'raw':
            folder = '/Users/kevimota/cernbox/'
            folders = []
            for it in os.scandir(folder):
                if it.name.find('MuOnia') < 0: continue
                if args.year == '2016APV':
                    if it.name.find('HIPM') < 0: continue
                else:
                    if it.name.find(args.year) < 0: continue
                    if it.name.find('HIPM') > -1: continue
                folders.append(it.path)

        files = get_files(folders, pattern='.coffea')

    else:
        folder = '/Users/kevimota/cernbox/CRAB_UserFiles/'
        folders = []
        for it in os.scandir(folder):
            if it.name.find('Upsilon') < 0: continue
            if it.name.find('Dstar') < 0: continue
            if it.name.find(args.year) < 0: continue
            if args.year == '2016':
                if it.name.find('2016APV') > -1: continue
            folders.append(it.path)

        files = get_files(folders, pattern='.root', exclude=exclude)
        
    h = HistogrammingProcessor(args.sel, args.is_mc, args.year)
    print(f"Creating histogram files for year {args.year}")
    if config_run['executor'] == 'iterative_executor' or len(files) == 1:
        for f in tqdm(files, total=len(files), unit=" files", desc="Processing"):
            h.process(f)
    else:
        with ProcessPoolExecutor(max_workers=len(files) if (len(files) < config_run['n_cores']) else config_run['n_cores']) as executor:
            list(tqdm(executor.map(h.process, files), total=len(files), unit=" files", desc="Processing"))

    print('Done!')

    if not args.is_mc:
        if args.sel == 'sel':
            hists = get_files([f"output/sel_data_hists/{args.year}"], pattern='.hists')
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
            save_folder = f'plots/sel_data_hists/{args.year}'
            plotter(hists, args.year, processed_lumi, save_folder)

        elif args.sel == 'raw':
            hists = get_files([f"output/raw_data_hists/{args.year}"], pattern='.hists')
            hists = [load(h) for h in hists]

            processed_lumi = get_lumi(args.year, get_trigger(args.year))

            for i in range(len(hists)):
                if i == 0:
                    hs_Dimu = hists[i]['Dimu']
                    hs_Dstar = hists[i]['Dstar']
                    hs_Dstar_D0 = hists[i]['Dstar_D0']
                    hs_DimuDstar = hists[i]['DimuDstar']
                else:
                    hs_Dimu = {k: hs_Dimu.get(k, 0) + hists[i]['Dimu'].get(k, 0) for k in set(hs_Dimu)}
                    hs_Dstar = {k: hs_Dstar.get(k, 0) + hists[i]['Dstar'].get(k, 0) for k in set(hs_Dstar)}
                    hs_Dstar_D0 = {k: hs_Dstar_D0.get(k, 0) + hists[i]['Dstar_D0'].get(k, 0) for k in set(hs_Dstar_D0)}
                    hs_DimuDstar = {k: hs_DimuDstar.get(k, 0) + hists[i]['DimuDstar'].get(k, 0) for k in set(hs_DimuDstar)}

            hists = {
                'Dimu': hs_Dimu,
                'Dstar': hs_Dstar,
                'Dstar_D0': hs_Dstar_D0,
                'DimuDstar': hs_DimuDstar,
            }

            save_folder = f'plots/raw_data_hists/{args.year}'
            plotter(hists, args.year, processed_lumi, save_folder)
    else:
        hists = get_files([f"output/mc_hists/{args.year}"], pattern='.hists')
        hists = [load(h) for h in hists]

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
        save_folder = f'plots/mc_hists/{args.year}'
        plotter(hists, args.year, None, save_folder, False)
