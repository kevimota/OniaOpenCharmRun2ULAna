from hist import Hist

from coffea.nanoevents.methods import candidate

import awkward as ak
ak.behavior.update(candidate.behavior)

import numpy as np
from coffea.util import load, save

class HistogrammingProcessor:
    def process(self, f) -> None:
        acc = load(f)

        Dimu_acc = acc['Dimu']
        Dstar_acc = acc['Dstar']
        DimuDstar_acc = acc['DimuDstar']

        Dimu_hists = {
            'pt': Hist.new.Regular(100, 0, 300, name='pt', label=r"$p_{T,\mu\mu}$ [GeV]").Double(),
            'eta': Hist.new.Regular(60, -3, 3, name="eta", label=r"$\eta_{\mu\mu}$").Double(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{\mu\mu}$").Double(),
            'mass': Hist.new.Regular(100, 8.6, 11, name="mass", label=r"$m_{\mu\mu}$ [GeV]").Double(),
            'rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$y_{\mu\mu}$").Double(),
            'dl': Hist.new.Regular(100, -1, 1, name="dl", label=r"decay length $\mu\mu$").Double(),
            'dlSig': Hist.new.Regular(100, -20, 20, name="dlSig", label=r"decay length significance $\mu\mu$").Double(),
            'chi2': Hist.new.Regular(100, 0, 10, name="chi2", label=r"$\mu\mu$ vtx fit $\chi^2$").Double(),
            'cosphi': Hist.new.Regular(60, -1, 1, name="cosphi", label=r"cos of pointing angle $\mu\mu$").Double(),
        }

        Dstar_hists = {
            'pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^*}$ [GeV]").Double(),
            'eta': Hist.new.Regular(60, -3, 3, name="eta", label=r"$\eta_{D^*}$").Double(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^*}$").Double(),
            'deltamr': Hist.new.Regular(50, 0.138, 0.162, name="deltam", label=r"$\Delta m_{D^*}$ [GeV]").Double(),
            'rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$y_{D^*}$").Double(),
            'D0cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0 of D^*$").Double(),
            'D0dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0 of D^*$").Double(),
            'D0dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0 of D^*$").Double(),
            'D0pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^0}$ [GeV]").Double(),
            'D0eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^0}$").Double(),
            'D0phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^0}$").Double(),
            'D0mass': Hist.new.Regular(100, 1.82, 1.9, name="mass", label=r"$m_{D^0}$ [GeV]").Double(),
        }
            
        DimuDstar_hists = {
            'pt': Hist.new.Regular(100, 0, 120, name='pt', label=r"$p_{T,\Upsilon D^*}$ [GeV]").Double(),
            'eta': Hist.new.Regular(60, -6, 6, name="eta", label=r"$\eta_{\Upsilon D^*}$").Double(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{\Upsilon D^*}$").Double(),
            'mass': Hist.new.Regular(100, 8.6, 70, name="mass", label=r"$m_{\Upsilon D^*}$ [GeV]").Double(),
            'rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$\eta_{\Upsilon D^*}$").Double(),
            'dimu_mass': Hist.new.Regular(100, 8.6, 11, name="mass", label=r"$m_{\mu\mu}$ [GeV]").Double(),
            'dimu_pt': Hist.new.Regular(100, 0, 300, name='pt', label=r"$p_{T,\mu\mu}$ [GeV]").Double(),
            'dimu_eta': Hist.new.Regular(60, -3, 3, name="eta", label=r"$\eta_{\mu\mu}$").Double(),
            'dimu_phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{\mu\mu}$").Double(),
            'dimu_rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$y_{\mu\mu}$").Double(),
            'dimu_dl': Hist.new.Regular(100, -1, 1, name="dl", label=r"decay length $\mu\mu$").Double(),
            'dimu_dlSig': Hist.new.Regular(100, -20, 20, name="dlSig", label=r"decay length significance $\mu\mu$").Double(),
            'dimu_chi2': Hist.new.Regular(100, 0, 10, name="chi2", label=r"$\mu\mu$ vtx fit $\chi^2$").Double(),
            'dimu_cosphi': Hist.new.Regular(60, -1, 1, name="cosphi", label=r"cos of pointing angle $\mu\mu$").Double(),
            'dstar_deltamr': Hist.new.Regular(50, 0.138, 0.162, name="deltam", label=r"$\Delta m_{D^*}$ [GeV]").Double(),
            'dstar_pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^*}$ [GeV]").Double(),
            'dstar_eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^*}$").Double(),
            'dstar_phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^*}$").Double(),
            'dstar_rap': Hist.new.Regular(60, -4, 4, name="rap", label=r"$y_{D^*}$").Double(),
            'dstar_d0_pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^0}$ [GeV]").Double(),
            'dstar_d0_eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^0}$").Double(),
            'dstar_d0_phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^0}$").Double(),
            'dstar_d0_mass': Hist.new.Regular(100, 1.82, 1.9, name="mass", label=r"$m_{D^0}$ [GeV]").Double(),
            'dstar_d0_cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0 of D^*$").Double(),
            'dstar_d0_dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0 of D^*$").Double(),
            'dstar_d0_dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0 of D^*$").Double(),
            'dstar_asso_chi2': Hist.new.Regular(100, 0, 5, name="assochi2", label=r"Association $\chi^2$").Double(),
            'dstar_asso_prob': Hist.new.Regular(100, 0, 1, name="assoprob", label=r"Association probability").Double(),
            'deltarap': Hist.new.Regular(60, -5, 5, name="deltarap", label=r"$\Delta y_{\Upsilon D^*}$").Double(),
            'deltapt': Hist.new.Regular(100, -30, 80, name='deltapt', label=r"$\Delta p_{T,\Upsilon D^*}$ [GeV]").Double(),
            'deltaeta': Hist.new.Regular(60, -6, 6, name="deltaeta", label=r"$\Delta \eta_{\Upsilon D^*}$").Double(),
            'deltaphi': Hist.new.Regular(60, 0, np.pi, name="deltaphi", label=r"$\Delta \phi_{\Upsilon D^*}$").Double(),
        }

        [Dimu_hists[h].fill(Dimu_acc[h].value) for h in Dimu_hists]
        [Dstar_hists[h].fill(Dstar_acc[h].value) for h in Dstar_hists]
        [DimuDstar_hists[h].fill(DimuDstar_acc[h].value) for h in DimuDstar_hists]

        hists = {
            'Dimu': Dimu_hists,
            'Dstar': Dstar_hists,
            'DimuDstar': DimuDstar_hists,
        }

        filename = f.replace('coffea', 'hists')
        save(hists, filename)
