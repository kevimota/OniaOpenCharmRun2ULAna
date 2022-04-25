import matplotlib.pyplot as plt
import mplhep
import numpy as np
from matplotlib.text import Text
from hist.intervals import ratio_uncertainty

plt.style.use(mplhep.style.CMS)
    
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
})

def create_plot1d(hist1d, log=False, ax=None, lumi=None):
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

    if ax == None:
        ax = plt.gca()

    hist1d.plot1d(ax=ax, fill=True, ec=(0,0,0,0.5))
    
    if not lumi == None:
        lumi = plt.text(1., 1., lumi + r" fb$^{-1}$ (13 TeV)",
                        fontsize=18,
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        transform=ax.transAxes
                       )
    
    if log:
        ax.set_yscale('log')
        ax.set_ylim(1, None)
    else:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,3), useMathText=True)
    
    centers = hist1d.axes.centers[0]
    values = hist1d.values()
    
    # compute mean and std:
    mean = np.sum(values*centers)/np.sum(values)
    std = np.sqrt(np.sum(values*((centers - mean)**2))/np.sum(values))
    
    annotation = TextArea(f"Total: {np.sum(values):.2e}" \
                    + "\n" + f"Mean: {mean:.2e}" \
                    + "\n" + f"Std: {std:.2e}", textprops=dict(size=14))
    
    at = AnchoredOffsetbox('upper right', child=annotation)
    at.patch.set_facecolor('None')
    ax.add_artist(at)
    
    return ax

def create_plot2d(hist2d, ax):
    if ax == None:
        ax = plt.gca()
    
    hist2d.plot2d(ax=ax)
    
    return ax

def acceptance_plot(hist_reco, hist_gen, ax=None, uncertainty_type ='efficiency', with_labels=True, with_unc=True):
    if ax == None:
        ax = plt.gca()
    
    ratio = hist_reco / hist_gen.values()

    n = [len(i.centers) for i in hist_reco.axes]

    if with_labels:
        try:
            err_down, err_up = ratio_uncertainty(hist_reco.values(), hist_gen.values(), uncertainty_type)
        except ValueError:
            print("Found numerator larger than denominator while calculating binomial uncertainty, switching to poison uncertainty calculation")
            err_down, err_up = ratio_uncertainty(hist_reco.values(), hist_gen.values())
        labels = []
        for ra, u, d in zip(ratio.values().ravel(), err_up.ravel(), err_down.ravel()):
            ra, u, d = f'{ra:.2f}', f'{u:.2f}', f'{d:.2f}'
            st = '$'+ra+'_{-'+d+'}^{+'+u+'}$'
            labels.append(st)
        labels = np.array(labels).reshape(*n)

        mplhep.hist2dplot(ratio, labels=labels, ax=ax)

        x = ax.get_children()
        for i0 in x:
            if isinstance(i0, Text):
                i0.set_size(10)

        return ax
    else:
        mplhep.hist2dplot(ratio, ax=ax)
        return ax