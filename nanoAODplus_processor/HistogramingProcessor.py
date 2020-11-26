import coffea.processor as processor
import boost_histogram as bh

import numpy as np
from coffea.util import save, load

import matplotlib
matplotlib.use('Agg')

def create_plot1d(hist, save_name, log=False):
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.style.CMS)
    # plot
    ax = plt.gca()


    plt.errorbar(hist.axes[0].centers,
             hist.view(),
             np.sqrt(hist.view()),
             fmt='.',
             color='blue',)

    hep.histplot(hist, ax=ax, color='blue')

    if log:
        ax.set_yscale('log')
    else:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,3), useMathText=True)

    ax.set_xlabel(hist.axes[0].metadata, loc='right')
    ax.set_ylabel("Counts", loc='top')

    # compute mean and std:
    mean = (hist.view() * hist.axes[0].centers).sum()/hist.sum()
    std = np.sqrt((hist.view()*((hist.axes[0].centers - mean)**2)).sum()/hist.sum())

    annotation = f"Total {hist.sum()}" \
                    + "\n" + f"Mean: {round(mean,2)}" \
                    + "\n" + f"Std: {round(std,2)}"
    
    ax.annotate(annotation, xy=(0.85, 0.85), xycoords='axes fraction', fontsize = "small",
                    ha='center', annotation_clip=False, bbox=dict(boxstyle='round', fc='None'))

    ax.set_xlim(hist.axes[0].edges[0], hist.axes[0].edges[-1] + hist.axes[0].widths[-1])
    
    fig = ax.get_figure()
    fig.savefig(save_name)
    ax.clear()
    fig.clear()

def create_plot2d(hist, save_name):
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.style.CMS)
    # plot
    ax = plt.gca()

    hep.hist2dplot(hist, ax=ax)
    ax.set_xlabel(hist.axes[0].metadata, loc='right')
    ax.set_ylabel(hist.axes[1].metadata, loc='top')

    fig = ax.get_figure()
    fig.savefig(save_name)
    ax.clear()
    fig.clear()

class HistogramingProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            'foo': processor.defaultdict_accumulator(int)
        })
     
    @property
    def accumulator(self):
        return self._accumulator
     
    def process(self, ds):
        output = self.accumulator.identity()
        acc = load(ds["file"])
        
        ############ Histogram definition
        # Muons
        hist_muon_lead = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,\mu}$ [GeV]"),
                                      bh.axis.Regular(60, -2.5, 2.5, metadata=r"$\eta_{\mu}$"),
                                      bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{\mu}$"),)

        hist_muon_trail = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,\mu}$ [GeV]"),
                                       bh.axis.Regular(60, -2.5, 2.5, metadata=r"$\eta_{\mu}$"),
                                       bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{\mu}$"),)

        #Dimu
        hist_dimu = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,\mu^+\mu^-}$ [GeV]"),
                                 bh.axis.Regular(80, -2.5, 2.5, metadata=r"$\eta_{\mu^+\mu^-}$"),
                                 bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{\mu^+\mu^-}$"),)

        hist_dimu_mass = bh.Histogram(bh.axis.Regular(100, 8.6, 11, metadata=r"$m_{\mu^+\mu^-}$ [GeV]"))

        # D0
        hist_D0 = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,D^0}$ [GeV]"),
                               bh.axis.Regular(80, -2.5, 2.5, metadata=r"$\eta_{D^0}$"),
                               bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{D^0}$"),)

        hist_D0_mass = bh.Histogram(bh.axis.Regular(100, 1.7, 2.0, metadata=r"$m_{D^0}$ [GeV]"))

        hist_D0_eta_mass = bh.Histogram(bh.axis.Regular(80, -2.5, 2.5, metadata=r"$\eta_{D^0}$"),
                                        bh.axis.Regular(100, 1.7, 2.0, metadata=r"$m_{D^0}$ [GeV]"))

        hist_D0_trk = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,D^0 trks}$ [GeV]"),
                                   bh.axis.Regular(80, -2.5, 2.5, metadata=r"$\eta_{D^0 trks}$"),
                                   bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{D^0 trks}$"),)

        # Dstar
        hist_Dstar = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,D*}$ [GeV]"),
                                  bh.axis.Regular(60, -2.5, 2.5, metadata=r"$\eta_{D*}$"),
                                  bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{D*}$"),)

        hist_Dstar_K = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,D* K}$ [GeV]"),
                                    bh.axis.Regular(60, -2.5, 2.5, metadata=r"$\eta_{D* K}$"),
                                    bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{D* K}$"),)

        hist_Dstar_pi = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,D* \pi}$ [GeV]"),
                                     bh.axis.Regular(60, -2.5, 2.5, metadata=r"$\eta_{D* \pi}$"),
                                     bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{D* \pi}$"),)
                                  
        hist_Dstar_pis = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,\pi_s}$ [GeV]"),
                                      bh.axis.Regular(60, -2.5, 2.5, metadata=r"$\eta_{\pi_s}$"),
                                      bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{\pi_s}$"),)

        hist_Dstar_mass = bh.Histogram(bh.axis.Regular(100, 1.8, 2.2, metadata=r"$m_{D*}$ [GeV]"))
        hist_Dstar_mass_refit = bh.Histogram(bh.axis.Regular(100, 1.8, 2.2, metadata=r"$m_{D* refit}$ [GeV]"))
        hist_Dstar_deltamr = bh.Histogram(bh.axis.Regular(50, 0.138, 0.162, metadata=r"$\Delta m_{refit}$ [GeV]"))
        hist_Dstar_deltam = bh.Histogram(bh.axis.Regular(50, 0.138, 0.162, metadata=r"$\Delta m$ [GeV]"))

        # Filling histograms
        hist_muon_lead.fill(acc["Muon_lead"]["__fast_pt"].value,
                            acc["Muon_lead"]["__fast_eta"].value, 
                            acc["Muon_lead"]["__fast_phi"].value)

        hist_muon_trail.fill(acc["Muon_trail"]["__fast_pt"].value,
                             acc["Muon_trail"]["__fast_eta"].value, 
                             acc["Muon_trail"]["__fast_phi"].value)

        hist_dimu.fill(acc["Dimu"]["__fast_pt"].value,
                       acc["Dimu"]["__fast_eta"].value, 
                       acc["Dimu"]["__fast_phi"].value)

        hist_dimu_mass.fill(acc["Dimu"]["__fast_mass"].value) 

        hist_D0.fill(acc["D0"]["__fast_pt"].value,
                     acc["D0"]["__fast_eta"].value, 
                     acc["D0"]["__fast_phi"].value)

        hist_D0_mass.fill(acc["D0"]["__fast_mass"].value)

        hist_D0_eta_mass.fill(acc["D0"]["__fast_eta"].value,
                              acc["D0"]["__fast_mass"].value)

        hist_D0_trk.fill(acc["D0_trk"]["t1_pt"].value,
                         acc["D0_trk"]["t1_eta"].value, 
                         acc["D0_trk"]["t1_phi"].value)

        hist_D0_trk.fill(acc["D0_trk"]["t2_pt"].value,
                         acc["D0_trk"]["t2_eta"].value, 
                         acc["D0_trk"]["t2_phi"].value)

        hist_Dstar.fill(acc["Dstar"]["__fast_pt"].value,
                        acc["Dstar"]["__fast_eta"].value, 
                        acc["Dstar"]["__fast_phi"].value)

        hist_Dstar_mass.fill(acc["Dstar"]["__fast_mass"].value)

        hist_Dstar_mass_refit.fill(acc["Dstar"]["deltamr"].value + acc["Dstar_D0"]["D0_mass"].value)

        hist_Dstar_K.fill(acc["Dstar_trk"]["K_pt"].value,
                          acc["Dstar_trk"]["K_eta"].value, 
                          acc["Dstar_trk"]["K_phi"].value)

        hist_Dstar_pi.fill(acc["Dstar_trk"]["pi_pt"].value,
                           acc["Dstar_trk"]["pi_eta"].value, 
                           acc["Dstar_trk"]["pi_phi"].value)

        hist_Dstar_pis.fill(acc["Dstar_trk"]["pis_pt"].value,
                            acc["Dstar_trk"]["pis_eta"].value, 
                            acc["Dstar_trk"]["pis_phi"].value)


        hist_Dstar_deltamr.fill(acc["Dstar"]["deltamr"].value)
        hist_Dstar_deltam.fill(acc["Dstar"]["deltam"].value)

        # Saving histograms
        save(hist_muon_lead, "output/" + ds['analyzer_name'] + "/hist/hist_Muon_lead.hist")
        save(hist_muon_trail, "output/" + ds['analyzer_name'] + "/hist/hist_Muon_trail.hist")
        save(hist_dimu, "output/" + ds['analyzer_name'] + "/hist/hist_Dimu.hist")
        save(hist_dimu_mass, "output/" + ds['analyzer_name'] + "/hist/hist_Dimu_mass.hist")
        save(hist_D0, "output/" + ds['analyzer_name'] + "/hist/hist_D0.hist")
        save(hist_D0_mass, "output/" + ds['analyzer_name'] + "/hist/hist_D0_mass.hist")
        save(hist_Dstar, "output/" + ds['analyzer_name'] + "/hist/hist_Dstar.hist")
        save(hist_Dstar_mass, "output/" + ds['analyzer_name'] + "/hist/hist_Dstar_mass.hist")
        save(hist_Dstar_mass_refit, "output/" + ds['analyzer_name'] + "/hist/hist_Dstar_mass_refit.hist")
        save(hist_Dstar_deltamr, "output/" + ds['analyzer_name'] + "/hist/hist_Dstar_deltamr.hist")
        save(hist_Dstar_deltam, "output/" + ds['analyzer_name'] + "/hist/hist_Dstar_deltam.hist")

        # Creating plots 1D
        plots_path = "plots/" + ds['analyzer_name'] + "/" 
        create_plot1d(hist_muon_lead[:, sum, sum], plots_path + "Muon_lead_pt.png", log=True)
        create_plot1d(hist_muon_lead[sum, :, sum], plots_path + "Muon_lead_eta.png")
        create_plot1d(hist_muon_lead[sum, sum, :], plots_path + "Muon_lead_phi.png")

        create_plot1d(hist_muon_trail[:, sum, sum], plots_path + "Muon_trail_pt.png", log=True)
        create_plot1d(hist_muon_trail[sum, :, sum], plots_path + "Muon_trail_eta.png")
        create_plot1d(hist_muon_trail[sum, sum, :], plots_path + "Muon_trail_phi.png")

        create_plot1d(hist_dimu[:, sum, sum], plots_path + "Dimu_pt.png", log=True)
        create_plot1d(hist_dimu[sum, :, sum], plots_path + "Dimu_eta.png")
        create_plot1d(hist_dimu[sum, sum, :], plots_path + "Dimu_phi.png")
        create_plot1d(hist_dimu_mass, plots_path + "Dimu_mass.png")

        create_plot1d(hist_D0[:, sum, sum], plots_path + "D0_pt.png", log=True)
        create_plot1d(hist_D0[sum, :, sum], plots_path + "D0_eta.png")
        create_plot1d(hist_D0[sum, sum, :], plots_path + "D0_phi.png")
        create_plot1d(hist_D0_mass, plots_path + "D0_mass.png")
        
        create_plot1d(hist_D0_trk[:, sum, sum], plots_path + "D0_trk_pt.png", log=True)
        create_plot1d(hist_D0_trk[sum, :, sum], plots_path + "D0_trk_eta.png")
        create_plot1d(hist_D0_trk[sum, sum, :], plots_path + "D0_trk_phi.png")

        create_plot1d(hist_Dstar[:, sum, sum], plots_path + "Dstar_pt.png", log=True)
        create_plot1d(hist_Dstar[sum, :, sum], plots_path + "Dstar_eta.png")
        create_plot1d(hist_Dstar[sum, sum, :], plots_path + "Dstar_phi.png")
        create_plot1d(hist_Dstar_mass, plots_path + "Dstar_mass.png")
        create_plot1d(hist_Dstar_mass_refit, plots_path + "Dstar_mass_refit.png")
        create_plot1d(hist_Dstar_deltamr, plots_path + "Dstar_deltamr.png")
        create_plot1d(hist_Dstar_deltam, plots_path + "Dstar_deltam.png")

        create_plot1d(hist_Dstar_K[:, sum, sum], plots_path + "Dstar_K_pt.png", log=True)
        create_plot1d(hist_Dstar_K[sum, :, sum], plots_path + "Dstar_K_eta.png")
        create_plot1d(hist_Dstar_K[sum, sum, :], plots_path + "Dstar_K_phi.png")

        create_plot1d(hist_Dstar_pi[:, sum, sum], plots_path + "Dstar_pi_pt.png", log=True)
        create_plot1d(hist_Dstar_pi[sum, :, sum], plots_path + "Dstar_pi_eta.png")
        create_plot1d(hist_Dstar_pi[sum, sum, :], plots_path + "Dstar_pi_phi.png")

        create_plot1d(hist_Dstar_pis[:, sum, sum], plots_path + "Dstar_pis_pt.png", log=True)
        create_plot1d(hist_Dstar_pis[sum, :, sum], plots_path + "Dstar_pis_eta.png")
        create_plot1d(hist_Dstar_pis[sum, sum, :], plots_path + "Dstar_pis_phi.png")
        
        # Creating plots 2D
        create_plot2d(hist_muon_lead[:,sum,:], plots_path + "Muon_lead_ptXphi")
        create_plot2d(hist_muon_trail[:,sum,:], plots_path + "Muon_trail_ptXphi")

        create_plot2d(hist_D0[:,:,sum], plots_path + "D0_ptXeta.png")
        create_plot2d(hist_D0[sum,:,:], plots_path + "D0_etaXphi.png")
        create_plot2d(hist_D0_eta_mass, plots_path + "D0_etaXmass.png")

        # return dummy accumulator
        return output

    def postprocess(self, accumulator):
        return accumulator      
