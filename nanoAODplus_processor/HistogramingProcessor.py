import coffea.processor as processor
import boost_histogram as bh

import numpy as np
from coffea.util import save, load

import matplotlib
matplotlib.use('Agg')

def create_plot1d(hist, save_name):
   import matplotlib.pyplot as plt
   import mplhep as hep
   plt.style.use(hep.style.CMS)

   # plot 
   ax = plt.gca()
   hep.histplot(hist, ax=ax)

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

   ax.ticklabel_format(axis='y', style='sci', scilimits=(0,3), useMathText=True)
   ax.set_xlim(hist.axes[0].edges[0], hist.axes[0].edges[-1] + hist.axes[0].widths[-1])
   fig = ax.get_figure()
   fig.savefig(save_name)
   ax.clear()

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
      
      # Histogram definition
      hist_muon_lead = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,\mu\mu}$ [GeV]"),
                                    bh.axis.Regular(60, -3.0, 3.0, metadata=r"$\eta_{\mu}$"),
                                    bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{\mu}$"),)

      hist_muon_trail = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,\mu}$ [GeV]"),
                                       bh.axis.Regular(60, -3.0, 3.0, metadata=r"$\eta_{\mu}$"),
                                       bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{\mu}$"),)

      hist_dimu = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,\mu^+\mu^-}$ [GeV]"),
                                 bh.axis.Regular(80, -4.0, 4.0, metadata=r"$\eta_{\mu^+\mu^-}$"),
                                 bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{\mu^+\mu^-}$"),)

      hist_dimu_mass = bh.Histogram(bh.axis.Regular(100, 8, 12, metadata=r"$m_{\mu^+\mu^-}$ [GeV]"))

      hist_D0 = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,D^0}$ [GeV]"),
                              bh.axis.Regular(80, -4.0, 4.0, metadata=r"$\eta_{D^0}$"),
                              bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{D^0}$"),)

      hist_D0_mass = bh.Histogram(bh.axis.Regular(100, 1.8, 2.2, metadata=r"$m_{D^0}$ [GeV]"))

      hist_Dstar = bh.Histogram(bh.axis.Regular(100, 0, 50, metadata=r"$p_{T,D*}$ [GeV]"),
                                 bh.axis.Regular(60, -3.0, 3.0, metadata=r"$\eta_{D*}$"),
                                 bh.axis.Regular(70, -3.5, 3.5, metadata=r"$\phi_{D*}$"),)

      hist_Dstar_mass = bh.Histogram(bh.axis.Regular(100, 1.8, 2.2, metadata=r"$m_{D*}$ [GeV]"))

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

      hist_Dstar.fill(acc["Dstar"]["__fast_pt"].value,
                        acc["Dstar"]["__fast_eta"].value, 
                        acc["Dstar"]["__fast_phi"].value)

      hist_Dstar_mass.fill(acc["Dstar"]["__fast_mass"].value)

      # Saving histograms
      save(hist_muon_lead, "output/" + ds['analyzer_name'] + "/hist/hist_muon_lead.hist")
      save(hist_muon_trail, "output/" + ds['analyzer_name'] + "/hist/hist_muon_trail.hist")
      save(hist_dimu, "output/" + ds['analyzer_name'] + "/hist/hist_dimu.hist")
      save(hist_dimu_mass, "output/" + ds['analyzer_name'] + "/hist/hist_dimu_mass.hist")
      save(hist_D0, "output/" + ds['analyzer_name'] + "/hist/hist_D0.hist")
      save(hist_D0_mass, "output/" + ds['analyzer_name'] + "/hist/hist_D0_mass.hist")
      save(hist_Dstar, "output/" + ds['analyzer_name'] + "/hist/hist_Dstar.hist")
      save(hist_Dstar_mass, "output/" + ds['analyzer_name'] + "/hist/hist_Dstar_mass.hist")

      # Creating plots
      plots_path = "plots/" + ds['analyzer_name'] + "/" 
      create_plot1d(hist_muon_lead[:, sum, sum], plots_path + "muon_lead_pt.png")
      create_plot1d(hist_muon_lead[sum, :, sum], plots_path + "muon_lead_eta.png")
      create_plot1d(hist_muon_lead[sum, sum, :], plots_path + "muon_lead_phi.png")

      create_plot1d(hist_muon_trail[:, sum, sum], plots_path + "muon_trail_pt.png")
      create_plot1d(hist_muon_trail[sum, :, sum], plots_path + "muon_trail_eta.png")
      create_plot1d(hist_muon_trail[sum, sum, :], plots_path + "muon_trail_phi.png")

      create_plot1d(hist_dimu[:, sum, sum], plots_path + "dimu_pt.png")
      create_plot1d(hist_dimu[sum, :, sum], plots_path + "dimu_eta.png")
      create_plot1d(hist_dimu[sum, sum, :], plots_path + "dimu_phi.png")
      create_plot1d(hist_dimu_mass, plots_path + "dimu_mass.png")

      create_plot1d(hist_D0[:, sum, sum], plots_path + "D0_pt.png")
      create_plot1d(hist_D0[sum, :, sum], plots_path + "D0_eta.png")
      create_plot1d(hist_D0[sum, sum, :], plots_path + "D0_phi.png")
      create_plot1d(hist_D0_mass, plots_path + "D0_mass.png")

      create_plot1d(hist_Dstar[:, sum, sum], plots_path + "Dstar_pt.png")
      create_plot1d(hist_Dstar[sum, :, sum], plots_path + "Dstar_eta.png")
      create_plot1d(hist_Dstar[sum, sum, :], plots_path + "Dstar_phi.png")
      create_plot1d(hist_Dstar_mass, plots_path + "Dstar_mass.png")

      # return dummy accumulator
      return output

   def postprocess(self, accumulator):
      return accumulator      
