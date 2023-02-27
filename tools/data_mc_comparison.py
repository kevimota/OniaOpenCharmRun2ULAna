import os
import hist
import mplhep
import numpy as np
import yaml

import matplotlib.pyplot as plt
import figure
from utils import get_lumi, get_trigger

from utils import get_files
from coffea.util import load

years = ['2016APV', '2016', '2017', '2018']
hist_log = ['pt', 'chi2', 'dl', 'dlSig']
rebin = 3j

def plotter(hists, year, processed_lumi, save_folder, is_data=True):
    fig, ax = plt.subplots()

    print(f'Creating plots in {save_folder}')
    if not os.path.exists(save_folder): os.makedirs(save_folder)

    for h in hists['data']:
        data_axis = hists['data'][h].axes[0]
        data_yerrs = np.sqrt(hists['data'][h].values())
        data_xerrs = data_axis.widths * 0.5
        data_centers = data_axis.centers
        data_values = hists['data'][h].values()
        log = False

        h_mc = (
            h
            .replace('d0_', 'D0')
            .replace('_asso_chi2', '_associationchi2')
            .replace('_asso_prob', '_associationProb')
        )

        for i in hist_log:
            if h.find(i) > -1:
                log = True
                break
        if log:
            mplhep.histplot(
                hists['mc'][h_mc],
                histtype="fill",
                edgecolor="black",
                label="DPS MC",
                ax=ax,
                density=True
            )
            ax.errorbar(
                data_centers,
                data_values/np.sum(data_axis.widths*data_values),
                xerr=data_xerrs,
                yerr=data_yerrs/np.sum(data_axis.widths*data_values),
                linestyle="None",
                color="black",
                marker="o",
                label="Data"
            )

            ax.set_yscale('log')

            mplhep.cms.label('Preliminary', loc=0, data=True, lumi=processed_lumi, lumi_format='{0:.2f}')
            ax.legend(fontsize=16)
            down_lim, up_lim = ax.get_ylim()
            up_lim *= 1.1
            ax.set_ylim(down_lim, up_lim)

            fig.savefig(f'{save_folder}/{h}_log_{year}.png')
            ax.clear()

        mplhep.histplot(
            hists['mc'][h_mc],
            histtype="fill",
            edgecolor="black",
            label="DPS MC",
            ax=ax,
            density=True
        )
        ax.errorbar(
            data_centers,
            data_values/np.sum(data_axis.widths*data_values),
            xerr=data_xerrs,
            yerr=data_yerrs/np.sum(data_axis.widths*data_values),
            linestyle="None",
            color="black",
            marker="o",
            label="Data"
        )
        
        mplhep.cms.label('Preliminary', loc=0, data=True, lumi=processed_lumi, lumi_format='{0:.2f}')
        ax.legend(fontsize=16)
        down_lim, up_lim = ax.get_ylim()
        up_lim *= 1.1
        ax.set_ylim(down_lim, up_lim)
        fig.savefig(f'{save_folder}/{h}_{year}.png')
        ax.clear()

for year in years:
    f_hists_mc = get_files([f"output/mc_hists/{year}"], pattern='.hists')
    f_hists_data = get_files([f"output/sel_data_hists/{year}"], pattern='.hists')
    hists_mc = {f[f.rfind('/')+1:]:load(f) for f in f_hists_mc}
    hists_data = {f[f.rfind('/')+1:]:load(f) for f in f_hists_data}

    with open('data/n_signal.yaml') as f:
        n_signal = yaml.load(f, Loader=yaml.FullLoader)

    hs_data_DimuDstar = None
    hs_mc_DimuDstar = None

    for i in hists_data:
        if hs_data_DimuDstar is None:
            hs_data_DimuDstar = hs_data_DimuDstar = {k: hists_data[i]['DimuDstar'].get(k, 0)[::rebin] for k in set(hists_data[i]['DimuDstar'])}
        else:
            hs_data_DimuDstar = {k: hs_data_DimuDstar.get(k, 0) + hists_data[i]['DimuDstar'].get(k, 0)[::rebin] for k in set(hs_data_DimuDstar)}

    """ n_data = {
        '15;30': 0,
        '30;60': 0,
        '60;120': 0,
        '120;150': 0,
    }

    for i in n_data:
        pt_min, pt_max = i.split(';')
        n_data[i] = hs_data_DimuDstar['dimu_pt'][hist.loc(float(pt_min)):hist.loc(float(pt_max))].sum() """

    n_mc = {
        '15;30': 0,
        '30;60': 0,
        '60;120': 0,
        '120;150': 0,
    }

    for i in hists_mc:
        pt_range = ";".join(i[i.rfind('Pt')+2:i.rfind('ToMuMu')].split('To'))
        pt_range = pt_range.replace('9', '15')
        if pt_range == '120':
            pt_range = '120;150'

        n_mc[pt_range] += hists_mc[i]['DimuDstar']['dimu_pt'].sum().value

    for i in hists_mc:
        pt_range = ";".join(i[i.rfind('Pt')+2:i.rfind('ToMuMu')].split('To'))
        pt_range = pt_range.replace('9', '15')
        if pt_range == '120':
            pt_range = '120;150'

        sf = n_signal[year][pt_range]/n_mc[pt_range]
        #sf = n_data[pt_range]/n_mc[pt_range]
        if hs_mc_DimuDstar is None:
            hs_mc_DimuDstar = {k: hists_mc[i]['DimuDstar'].get(k, 0)[::rebin]*sf for k in set(hists_mc[i]['DimuDstar'])}
        else:
            hs_mc_DimuDstar = {k: hs_mc_DimuDstar.get(k, 0) + hists_mc[i]['DimuDstar'].get(k, 0)[::rebin]*sf for k in set(hs_mc_DimuDstar)}

    hists = {
        'data': hs_data_DimuDstar,
        'mc': hs_mc_DimuDstar,
    }
    
    processed_lumi = get_lumi(year, get_trigger(year))

    save_folder = f'plots/comparison/{year}'
    plotter(hists, year, processed_lumi, save_folder, False)

    