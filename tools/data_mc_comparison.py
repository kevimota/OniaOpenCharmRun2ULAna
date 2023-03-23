import os
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

def plotter(hists, year, processed_lumi, save_folder):
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
                [hists['mc_DPS'][h_mc], hists['mc_SPS'][h_mc]],
                histtype="fill",
                edgecolor="black",
                label=["DPS MC", "SPS MC"],
                ax=ax,
                density=True,
                alpha=0.7
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
            [hists['mc_DPS'][h_mc], hists['mc_SPS'][h_mc]],
            histtype="fill",
            edgecolor="black",
            label=["DPS MC", "SPS MC"],
            ax=ax,
            density=True,
            alpha=0.7
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
    f_hists_mc_DPS = get_files([f"output/mc_DPS_hists/{year}"], pattern='.hists')
    f_hists_mc_SPS = get_files([f"output/mc_SPS_hists/{year}"], pattern='.hists')
    f_hists_data = get_files([f"output/sel_data_hists/{year}"], pattern='.hists')
    hists_mc_DPS = {f[f.rfind('/')+1:]:load(f) for f in f_hists_mc_DPS}
    hists_mc_SPS = {f[f.rfind('/')+1:]:load(f) for f in f_hists_mc_SPS}
    hists_data = {f[f.rfind('/')+1:]:load(f) for f in f_hists_data}

    with open('data/n_signal.yaml') as f:
        n_signal = yaml.load(f, Loader=yaml.FullLoader)

    hs_data_DimuDstar = None
    hs_mc_DPS_DimuDstar = None
    hs_mc_SPS_DimuDstar = None

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

    n_mc_DPS = {
        '15;30': 0,
        '30;60': 0,
        '60;120': 0,
        '120;150': 0,
    }

    n_mc_SPS = {
        '15;30': 0,
        '30;60': 0,
        '60;120': 0,
        '120;150': 0,
    }

    for i in hists_mc_DPS:
        if i[i.rfind('/')+1:].startswith('DPS'):
            pt_range = ";".join(i[i.rfind('Filter-')+7:i.rfind(f'_{year}')].split('To'))
        else:
            pt_range = ";".join(i[i.rfind('Pt')+2:i.rfind('ToMuMu')].split('To'))
        pt_range = pt_range.replace('9', '15')
        if pt_range == '120':
            pt_range = '120;150'

        n_mc_DPS[pt_range] += hists_mc_DPS[i]['DimuDstar']['dimu_pt'].sum().value

    for i in hists_mc_SPS:
        pt_range = ";".join(i[i.rfind('Upsilon_')+8:i.rfind('_Dstar')].split('to'))
        pt_range = pt_range.replace('9', '15')
        if pt_range == '120':
            pt_range = '120;150'

        n_mc_SPS[pt_range] += hists_mc_SPS[i]['DimuDstar']['dimu_pt'].sum().value

    for i in hists_mc_DPS:
        if i[i.rfind('/')+1:].startswith('DPS'):
            pt_range = ";".join(i[i.rfind('Filter-')+7:i.rfind(f'_{year}')].split('To'))
        else:
            pt_range = ";".join(i[i.rfind('Pt')+2:i.rfind('ToMuMu')].split('To'))
        pt_range = pt_range.replace('9', '15')
        if pt_range == '120':
            pt_range = '120;150'

        sf = n_signal[year][pt_range]/n_mc_DPS[pt_range]
        #sf = n_data[pt_range]/n_mc_DPS[pt_range]
        if hs_mc_DPS_DimuDstar is None:
            hs_mc_DPS_DimuDstar = {k: hists_mc_DPS[i]['DimuDstar'].get(k, 0)[::rebin]*sf for k in set(hists_mc_DPS[i]['DimuDstar'])}
        else:
            hs_mc_DPS_DimuDstar = {k: hs_mc_DPS_DimuDstar.get(k, 0) + hists_mc_DPS[i]['DimuDstar'].get(k, 0)[::rebin]*sf for k in set(hs_mc_DPS_DimuDstar)}

    for i in hists_mc_SPS:
        pt_range = ";".join(i[i.rfind('Upsilon_')+8:i.rfind('_Dstar')].split('to'))
        pt_range = pt_range.replace('9', '15')
        if pt_range == '120':
            pt_range = '120;150'

        sf = n_signal[year][pt_range]/n_mc_SPS[pt_range]
        #sf = n_data[pt_range]/n_mc_SPS[pt_range]
        if hs_mc_SPS_DimuDstar is None:
            hs_mc_SPS_DimuDstar = {k: hists_mc_SPS[i]['DimuDstar'].get(k, 0)[::rebin]*sf for k in set(hists_mc_SPS[i]['DimuDstar'])}
        else:
            hs_mc_SPS_DimuDstar = {k: hs_mc_SPS_DimuDstar.get(k, 0) + hists_mc_SPS[i]['DimuDstar'].get(k, 0)[::rebin]*sf for k in set(hs_mc_SPS_DimuDstar)}

    hists = {
        'data': hs_data_DimuDstar,
        'mc_DPS': hs_mc_DPS_DimuDstar,
        'mc_SPS': hs_mc_SPS_DimuDstar,
    }
    
    processed_lumi = get_lumi(year, get_trigger(year))

    save_folder = f'plots/comparison/{year}'
    plotter(hists, year, processed_lumi, save_folder)

    