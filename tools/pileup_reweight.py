import uproot
import numpy as np
from hist import Hist

import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep

plt.style.use(mplhep.style.CMS)
    
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
})

pileup_files = {
    'data': {
        '2016': 'data/corrections/pileup_rootfiles/PileupData_2016.root',
        '2017': 'data/corrections/pileup_rootfiles/PileupData_2017.root',
        '2018': 'data/corrections/pileup_rootfiles/PileupData_2018.root',
    },
    'mc': {
        '2016': 'data/corrections/pileup_rootfiles/PileupMC_2016.root',
        '2017': 'data/corrections/pileup_rootfiles/PileupMC_2017.root',
        '2018': 'data/corrections/pileup_rootfiles/PileupMC_2018.root',
    }
}

save_folder_plots = 'plots/pileup_reweight'
save_folder_files = 'data/corrections'

def get_weights(hist_data: Hist, hist_mc: Hist) -> Hist:
    num = hist_data.values()
    den = hist_mc.values()

    scaled_num = num/num.sum()
    scaled_den = den/den.sum()

    weights = np.where(
        (scaled_num > 0) & (scaled_den > 0),
        scaled_num/scaled_den,
        1.0,    
    )

    hist_weights = Hist.new.Regular(100, 0, 100, name='n_vtx', label='Num. of reco vertices').Double()
    hist_weights[...] = weights

    # Save to root
    save_hist = uproot.recreate(f'{save_folder_files}/pileup_reweight_{year}.root')
    save_hist['h_weights'] = hist_weights

    return hist_weights


def create_ratio_plot(hist_data: Hist, hist_mc:Hist, hist_weights: Hist, year: str):
    fig = plt.figure()
    grid = fig.add_gridspec(2, 1, hspace=0, height_ratios=[3, 1])
    main_ax = fig.add_subplot(grid[0])
    subplot_ax = fig.add_subplot(grid[1])
    plt.setp(main_ax.get_xticklabels(), visible=False)
    
    data = hist_data.values()/hist_data.values().sum()
    hist_data_errors = np.sqrt(hist_data.values())/hist_data.values().sum()
    
    mplhep.histplot(
        hist_mc, 
        ax=main_ax, 
        density=True,
    )

    main_ax.errorbar(
        x=hist_data.axes[0].centers,
        y=data,
        yerr=hist_data_errors,
        xerr=hist_data.axes[0].widths,
        fmt="ko", 
        label="Data"
    )

    subplot_ax.errorbar(
        x=hist_weights.axes[0].centers,
        y=hist_weights.values(),
        xerr=hist_weights.axes[0].widths,
        fmt="ko", 
        label="weight"
    )

    subplot_ax.set_ylim(0.8, 1.45)
    subplot_ax.set_xlabel(hist_weights.axes[0].label)
    subplot_ax.set_ylabel('Ratio')

    pathlib.Path(save_folder_plots).mkdir(parents=True, exist_ok=True)

    fig.savefig(f'{save_folder_plots}/ratio_dataMC_{year}.png')

if __name__ == '__main__':
    for year in pileup_files['data']:
        print(f'Running for year {year}')

        file_data = uproot.open(pileup_files['data'][year])
        file_mc = uproot.open(pileup_files['mc'][year])
        hist_data = file_data["pileup"].to_hist()
        hist_mc = file_mc["input_Event/N_TrueInteractions"].to_hist()

        hist_weights = get_weights(hist_data, hist_mc)
        create_ratio_plot(hist_data, hist_mc, hist_weights, year)