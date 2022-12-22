import time
import csv
import yaml
from hist.intervals import ratio_uncertainty, poisson_interval
from uncertainties import ufloat

from coffea.nanoevents import BaseSchema

import awkward as ak
import numpy as np
import pandas as pd
import mplhep as hep
from coffea import processor

from hist import Hist
import hist

from nanoAODplus_processor.EfficiencyProcessor import EfficiencyProcessor
from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

import matplotlib.pyplot as plt
plt.style.use(hep.style.CMS)
from matplotlib.text import Text

from tools.utils import *
from tools.collections import *
from tools.figure import create_plot1d, create_plot2d, acceptance_plot

years = ['2016APV', '2016', '2017', '2018']
ps = [
    '/Users/kevimota/cernbox/CRAB_UserFiles/UpsilonPt9To30ToMuMuDstarToD0pi',
    '/Users/kevimota/cernbox/CRAB_UserFiles/UpsilonPt30To60ToMuMuDstarToD0pi',
    '/Users/kevimota/cernbox/CRAB_UserFiles/UpsilonPt60To120ToMuMuDstarToD0pi',
    '/Users/kevimota/cernbox/CRAB_UserFiles/UpsilonPt120ToMuMuDstarToD0pi',
]

exclude = [
    'UpsilonPt9To30ToMuMuDstarToD0pi_2017-v2_18.root',
]


with open('config/efficiency.yaml', 'r') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

def create_eff_csv(num, den, filename):
    uncertainty_type = 'efficiency'
    if len(num) != len(den):
        raise ValueError("num and den must be the same length")

    with open(filename, 'w') as f:
        writer = csv.writer(f)
        header = [
            'x_low', 'x_high', 'y_low', 'y_high', 'n_num', 'n_den', 'eff',
            'eff_err_up', 'eff_err_down',
        ]
        writer.writerow(header)
        
        for row in num:
            if row == 'global': continue
            (x_low, x_high), (y_low, y_high) = [r.split(';') for r in row.split(':')]
            n_num = num[row]
            n_den = den[row]
            try:
                eff = n_num/n_den
                if n_num > n_den:
                    print("Found numerator larger than denominator while calculating binomial uncertainty, switching to poison uncertainty calculation")
                    err_down = np.abs(poisson_interval(eff, n_num / np.square(n_den)) - eff)
                    err_up = err_down
                else:
                    err_down, err_up = ratio_uncertainty(n_num, n_den, uncertainty_type=uncertainty_type)
            except ZeroDivisionError:
                eff = err_down = err_up = n_den = 0
            
            writer.writerow([x_low, x_high, y_low, y_high, n_num, n_den, eff, err_up, err_down,])


def create_eff_plot2D(file_eff, bins_x, bins_y, name_axes, label_axes, ax, with_labels=False, **kwargs):
    hist_eff = (
        Hist.new
        .Var(bins_x, name=name_axes[0], label=label_axes[0])
        .Var(bins_y, name=name_axes[1], label=label_axes[1])
        .Double()
    )
    df = pd.read_csv(file_eff)

    x_high = df.x_high.to_numpy()
    x_low = df.x_low.to_numpy()
    y_high = df.y_high.to_numpy()
    y_low = df.y_low.to_numpy()

    x_center = x_low + (x_high - x_low)/2
    y_center = y_low + (y_high - y_low)/2
    eff = df.eff.to_numpy()
    eff_err_up = df.eff_err_up.to_numpy()
    eff_err_down = df.eff_err_down.to_numpy()

    for i in range(len(x_center)):
        x = x_center[i]
        y = y_center[i]
        hist_eff[hist.loc(x), hist.loc(y)] = eff[i]

    if with_labels:
        n = [len(i.centers) for i in hist_eff.axes]

        labels = []
        for ra, u, d in zip(eff, eff_err_up, eff_err_down):
            ra, u, d = f'{ra:.2f}', f'{u:.2f}', f'{d:.2f}'
            st = '$'+ra+'_{-'+d+'}^{+'+u+'}$'
            labels.append(st)
        labels = np.array(labels).reshape(*n)

        artists = hep.hist2dplot(hist_eff, labels=labels, ax=ax, **kwargs)
        x = ax.get_children()
        for i0 in x:
            if isinstance(i0, Text):
                i0.set_size(10)
                i0.set_rotation(270)
        
        if ('dimu' in file_eff) or ('asso' in file_eff):
            ticks, labels = plt.xticks()
            for idx, i in enumerate(ticks):
                if i == 20.:
                    labels[idx] = None
            
            ax.set_xticklabels(labels)
            #ax.set_xscale('log')
            
                
    else:
        artists = hep.hist2dplot(hist_eff, ax=ax, **kwargs)
    
    return hist_eff


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Create efficiency plots")
    parser.add_argument("-y", "--year", help="Year to create efficiency", type=str, required=True, choices=years)
    parser.add_argument("-p", "--plot", help="Plot efficiency", action='store_true')
    args = parser.parse_args()

    year = args.year

    paths = [f'{p}_{year}-v2' for p in ps]
    files = [get_files([p], exclude=exclude) for p in paths]

    output = []
    tstart = time.time()
    for f in files:
        data = {"test": f[:]}

        upsilon_cut = float(f[0][f[0].find('Pt')+2:f[0].find('To')])
        print(f"Treating file with Pt > {upsilon_cut}")

        output.append(
            processor.run_uproot_job(
                data, 
                treename='Events',
                processor_instance=EfficiencyProcessor(
                    dimu_cut=upsilon_cut,
                    year=year,
                    config=config,
                ),
                executor=processor.futures_executor,
                executor_args={"schema": BaseSchema, 'workers': 8, 'skipbadfiles': True},
                chunksize=360000,
            )
                                    
        )

    print(f"Process finished in: {time.time() - tstart:.2f} s")

    for i, out in enumerate(output):
        if i == 0:
            N_gen_dimu          = out['N_gen_dimu']
            N_reco_dimu         = out['N_reco_dimu']
            N_cuts_dimu         = out['N_cuts_dimu']
            N_trigger_dimu      = out['N_trigger_dimu']
            N_gen_dstar         = out['N_gen_dstar']
            N_reco_dstar        = out['N_reco_dstar']
            N_cuts_dstar        = out['N_cuts_dstar']
            N_num_asso          = out['N_num_asso']
            N_den_asso          = out['N_den_asso']
        else:
            N_gen_dimu          += out['N_gen_dimu']
            N_reco_dimu         += out['N_reco_dimu']
            N_cuts_dimu         += out['N_cuts_dimu']
            N_trigger_dimu      += out['N_trigger_dimu']
            N_gen_dstar         += out['N_gen_dstar']
            N_reco_dstar        += out['N_reco_dstar']
            N_cuts_dstar        += out['N_cuts_dstar']
            N_num_asso          += out['N_num_asso']
            N_den_asso          += out['N_den_asso']

    create_eff_csv(N_reco_dimu, N_gen_dimu, f'output/efficiency/accep_dimu_{year}.csv')
    create_eff_csv(N_reco_dstar, N_gen_dstar, f'output/efficiency/accep_dstar_{year}.csv')
    create_eff_csv(N_cuts_dimu, N_reco_dimu, f'output/efficiency/eff_cuts_dimu_{year}.csv')
    create_eff_csv(N_cuts_dstar, N_reco_dstar, f'output/efficiency/eff_cuts_dstar_{year}.csv')
    create_eff_csv(N_trigger_dimu, N_cuts_dimu, f'output/efficiency/eff_trigger_dimu_{year}.csv')
    create_eff_csv(N_num_asso, N_den_asso, f'output/efficiency/eff_asso_{year}.csv')

    if args.plot:
        for it in os.scandir('output/efficiency'):
            if it.name.find(year) < 0: continue
            if year == '2016' and it.name.find('2016APV') > 0: continue
            print(f'Creating plot for {it.name}')
            fig, ax = plt.subplots()
            if it.name.find('dimu') > -1:
                bins_x = config['bins_pt_dimu']
                bins_y = config['bins_rap_dimu']
                axes_name = ['pt', 'rap']
                axes_label = [r"$p_{T,\mu\mu}$ [GeV]", r"$|y_{\mu\mu}|$"]
            elif it.name.find('dstar') > -1:
                bins_x = config['bins_pt_dstar']
                bins_y = config['bins_rap_dstar']
                axes_name = ['pt', 'rap']
                axes_label = [r"$p_{T,D^*}$ [GeV]", r"$|y_{D^*}|$"]
            elif it.name.find('asso') > -1:
                bins_x = config['bins_pt_dimu']
                bins_y = config['bins_pt_dstar']
                axes_name = ['pt_dimu', 'pt_dstar']
                axes_label = [r"$p_{T,\mu\mu}$ [GeV]", r"$p_{T,D^*}$ [GeV]"]
            else:
                print('Definition not found, skipping...')
                continue
            
            create_eff_plot2D(
                f'{it.path}', 
                bins_x,
                bins_y,
                axes_name,
                axes_label,
                ax,
                with_labels=True, vmin=0, vmax=1,
            )
            savename = it.name.replace('.csv', '.png')
            fig.savefig(f'plots/efficiency/{savename}')
            plt.close()