from hist import Hist
import uproot
import re
import yaml

from coffea.lookup_tools import extractor

from coffea.nanoevents.methods import candidate

import awkward as ak
ak.behavior.update(candidate.behavior)

import numpy as np
from coffea.util import load, save
from tools.utils import *
from tools.collections import *

branch_filter = re.compile('(Dimu_|Muon_|D0_|Dstar_|PVtx_|HLT_Dimuon)')
pileup_file = '/Users/kevimota/CERN/OniaOpenCharmRun2ULAna/data/corrections/pile_up_reweight_{year}.root'
muon_reco_file = '/Users/kevimota/CERN/OniaOpenCharmRun2ULAna/data/corrections/Efficiency_muon_generalTracks_Run{year}_UL_trackerMuon.json'
muon_id_file = '/Users/kevimota/CERN/OniaOpenCharmRun2ULAna/data/corrections/Efficiency_muon_trackerMuon_Run{year}_UL_ID.json'

trgs = {
    "2016APV": "HLT_Dimuon13_Upsilon",
    "2016":    "HLT_Dimuon13_Upsilon",
    "2017":    "HLT_Dimuon24_Upsilon_noCorrL1",
    "2018":    "HLT_Dimuon24_Upsilon_noCorrL1",
} 


def get_weight(evaluator, Muon, PVtx):
    pileup_weight = evaluator['pileup'](ak.num(PVtx))
    mu0_reco_weight = evaluator['NUM_TrackerMuons_DEN_genTracks/abseta_pt_value'](np.absolute(Muon.slot0.eta), Muon.slot0.pt)
    mu1_reco_weight = evaluator['NUM_TrackerMuons_DEN_genTracks/abseta_pt_value'](np.absolute(Muon.slot1.eta), Muon.slot1.pt)
    mu0_id_weight = evaluator['NUM_SoftID_DEN_TrackerMuons/abseta_pt_value'](np.absolute(Muon.slot0.eta), Muon.slot0.pt)
    mu1_id_weight = evaluator['NUM_SoftID_DEN_TrackerMuons/abseta_pt_value'](np.absolute(Muon.slot1.eta), Muon.slot1.pt)
    weight = pileup_weight*mu0_reco_weight*mu1_reco_weight*mu0_id_weight*mu1_id_weight
    
    return weight

class HistogrammingProcessor:
    def __init__(self, selection, is_mc, year) -> None:
        self.selection = selection
        self.year = year
        self.is_mc = is_mc

        with open('config/efficiency.yaml', 'r') as f:
            self.config = yaml.load(f, Loader=yaml.FullLoader)

    def process(self, f) -> None:
        if not self.is_mc:
            if self.selection == 'sel':
                self.process_data_sel(f)
            elif self.selection == 'raw':
                self.process_data_raw(f)
        else:
            self.process_mc(f)

    def process_data_raw(self, f) -> None:
        acc = load(f)

        Dimu_acc = acc['Dimu']
        Muon_lead_acc = acc['Muon_lead']
        Muon_trail_acc = acc['Muon_trail']
        Dstar_acc = acc['Dstar']
        Dstar_D0_acc = acc['Dstar_D0']
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
        }

        Dstar_D0_hists = {
            'D0cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0$ of $D^*$").Double(),
            'D0dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0$ of $D^*$").Double(),
            'D0dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0$ of $D^*$").Double(),
            'D0pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^0}$ [GeV]").Double(),
            'D0eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^0}$").Double(),
            'D0phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^0}$").Double(),
            'D0mass': Hist.new.Regular(100, 1.82, 1.9, name="mass", label=r"$m_{D^0}$ [GeV]").Double(),
        }

        DimuDstar_p4 = build_p4(DimuDstar_acc)

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
            'dstar_D0pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^0}$ [GeV]").Double(),
            'dstar_D0eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^0}$").Double(),
            'dstar_D0phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^0}$").Double(),
            'dstar_D0mass': Hist.new.Regular(100, 1.82, 1.9, name="mass", label=r"$m_{D^0}$ [GeV]").Double(),
            'dstar_D0cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0$ of $D^*$").Double(),
            'dstar_D0dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0$ of $D^*$").Double(),
            'dstar_D0dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0$ of $D^*$").Double(),
            'dstar_associationchi2': Hist.new.Regular(100, 0, 5, name="assochi2", label=r"Association $\chi^2$").Double(),
            'dstar_associationProb': Hist.new.Regular(100, 0, 1, name="assoprob", label=r"Association probability").Double(),
            'deltarap': Hist.new.Regular(60, -5, 5, name="deltarap", label=r"$\Delta y_{\Upsilon D^*}$").Double(),
            'deltapt': Hist.new.Regular(100, -30, 80, name='deltapt', label=r"$\Delta p_{T,\Upsilon D^*}$ [GeV]").Double(),
            'deltaeta': Hist.new.Regular(60, -6, 6, name="deltaeta", label=r"$\Delta \eta_{\Upsilon D^*}$").Double(),
            'deltaphi': Hist.new.Regular(60, 0, np.pi, name="deltaphi", label=r"$\Delta \phi_{\Upsilon D^*}$").Double(),
        }

        [Dimu_hists[h].fill(Dimu_acc[h].value) for h in Dimu_hists]
        [Dstar_hists[h].fill(Dstar_acc[h].value) for h in Dstar_hists]
        [Dstar_D0_hists[h].fill(Dstar_D0_acc[h].value) for h in Dstar_D0_hists]

        DimuDstar_hists['pt'].fill(DimuDstar_p4.pt)
        DimuDstar_hists['eta'].fill(DimuDstar_p4.eta)
        DimuDstar_hists['phi'].fill(DimuDstar_p4.phi)
        DimuDstar_hists['mass'].fill(DimuDstar_p4.mass)
        DimuDstar_hists['rap'].fill(DimuDstar_p4.rap)
        DimuDstar_hists['deltapt'].fill(DimuDstar_acc['deltapt'].value)
        DimuDstar_hists['deltarap'].fill(DimuDstar_acc['deltarap'].value)
        DimuDstar_hists['deltaeta'].fill(DimuDstar_acc['deltaeta'].value)
        DimuDstar_hists['deltaphi'].fill(DimuDstar_acc['deltaphi'].value)
        
        [DimuDstar_hists[h].fill(DimuDstar_acc['Dimu'][h[h.find('_')+1:]].value) for h in DimuDstar_hists if 'dimu' in h]
        [DimuDstar_hists[h].fill(DimuDstar_acc['Dstar'][h[h.find('_')+1:]].value) for h in DimuDstar_hists if 'dstar' in h]

        hists = {
            'Dimu': Dimu_hists,
            'Dstar': Dstar_hists,
            'Dstar_D0': Dstar_D0_hists,
            'DimuDstar': DimuDstar_hists,
        }

        filename = f"output/raw_data_hists/{self.year}{f[f.rfind('/'):]}".replace('coffea', 'hists')
        save(hists, filename)

    def process_data_sel(self, f) -> None:
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
            'D0cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0$ of $D^*$").Double(),
            'D0dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0$ of $D^*$").Double(),
            'D0dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0$ of $D^*$").Double(),
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
            'dstar_d0_cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0$ of $D^*$").Double(),
            'dstar_d0_dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0$ of $D^*$").Double(),
            'dstar_d0_dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0$ of $D^*$").Double(),
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

        filename = f"output/sel_data_hists/{self.year}{f[f.rfind('/'):]}".replace('coffea', 'hists')
        save(hists, filename)

    def process_mc(self, f) -> None:
        events = uproot.open(f)['Events'].arrays(filter_name=branch_filter.match)
        pt_range = f[f.rfind('Pt')+2:f.rfind('ToMuMuDstar')]
        if 'To' in pt_range:
            min_value, max_value = pt_range.split('To')
            dimu_pt_min = float(min_value)
            dimu_pt_max = float(max_value)
        else:
            dimu_pt_min = float(pt_range)
            dimu_pt_max = 9999.

        if len(events) == 0:
            return

        Dimu = ak.zip({**get_vars_dict(events, dimu_cols)}, with_name="PtEtaPhiMCandidate")
        Muon_all = ak.zip({**get_vars_dict(events, muon_cols)}, with_name="PtEtaPhiMCandidate")
        D0 = ak.zip({'mass': events.D0_mass12, **get_vars_dict(events, d0_cols)}, with_name="PtEtaPhiMCandidate")
        Dstar = ak.zip({'mass': (events.Dstar_D0mass + events.Dstar_deltamr),
                        'charge': events.Dstar_pischg,
                        **get_vars_dict(events, dstar_cols)}, 
                        with_name="PtEtaPhiMCandidate")
        PVtx = ak.zip({**get_vars_dict(events, pvtx_cols)})
        HLT = ak.zip({**get_hlt(events, [trgs[self.year]])})
        trigger = HLT[trgs[self.year]]

        # Corrections for Muon and Pileup
        ext = extractor()
        ext.add_weight_sets([f"pileup h_weights {pileup_file.format(year=self.year)}"])
        ext.add_weight_sets([f"* * {muon_reco_file.format(year=self.year)}"])
        ext.add_weight_sets([f"* * {muon_id_file.format(year=self.year)}"])
        ext.finalize()
        evaluator = ext.make_evaluator()

        ############### Phase Space
        Dimu = ak.mask(Dimu, Dimu.charge == 0)
        Dimu = ak.mask(Dimu, (Dimu.mass > self.config['dimu_mass_low']) & (Dimu.mass < self.config['dimu_mass_high']))
        Muon = ak.zip({'0': Muon_all[Dimu.t1muIdx], '1': Muon_all[Dimu.t2muIdx]})
        leading_mu = (Muon.slot0.pt > Muon.slot1.pt)
        Muon = ak.zip({'0': ak.where(leading_mu, Muon.slot0, Muon.slot1), '1': ak.where(~leading_mu, Muon.slot0, Muon.slot1)})

        muon_pt_cut = (Muon.slot0.pt > self.config['muon_pt_min']) & (Muon.slot1.pt > self.config['muon_pt_min'])
        muon_eta_cut = (np.absolute(Muon.slot0.eta) < self.config['muon_eta']) & (np.absolute(Muon.slot0.eta) < self.config['muon_eta'])
        muon_sim_cut = (Muon.slot0.simIdx > -1) & (Muon.slot1.simIdx > -1)
        dimu_pt_cut = (Dimu.pt > max(self.config['dimu_pt_min'], dimu_pt_min)) & (Dimu.pt < min(self.config['dimu_pt_max'], dimu_pt_max))
        dimu_eta_cut = np.absolute(Dimu.rap) < self.config['dimu_rap']
        Muon = ak.mask(Muon, muon_eta_cut & muon_pt_cut & muon_sim_cut & dimu_pt_cut & dimu_eta_cut)
        Dimu = ak.mask(Dimu, muon_eta_cut & muon_pt_cut & muon_sim_cut & dimu_pt_cut & dimu_eta_cut)
        arg_sort = ak.argsort(Dimu.pt, axis=1, ascending=False)
        Dimu = Dimu[arg_sort]
        Muon = Muon[arg_sort]

        Dstar = Dstar[~Dstar.hasMuon]
        Dstar = Dstar[Dstar.Kchg != Dstar.pichg]
        Dstar = Dstar[(Dstar.pt > self.config['dstar_pt_min']) & (Dstar.pt < self.config['dstar_pt_max'])]
        Dstar = Dstar[np.absolute(Dstar.rap) < self.config['dstar_rap']]
        arg_sort = ak.argsort(Dstar.pt, axis=1, ascending=False)
        Dstar = Dstar[arg_sort]

        ############### Cuts
        soft_id = (Muon.slot0.softId > 0) & (Muon.slot1.softId > 0)
        Dimu = ak.mask(Dimu, soft_id)
        Muon = ak.mask(Muon, soft_id)

        Dstar = Dstar[(Dstar.Kpt > self.config['dstar_track_pt_cut']) & (Dstar.pipt > self.config['dstar_track_pt_cut'])]
        Dstar = Dstar[(Dstar.Kchindof < 2.5) & (Dstar.pichindof < 2.5)]
        Dstar = Dstar[(Dstar.KnValid > 4) & (Dstar.pinValid > 4) & (Dstar.KnPix > 1) & (Dstar.pinPix > 1)]
        Dstar = Dstar[(Dstar.Kdxy < 0.5) & (Dstar.pidxy < 0.5)]
        K_theta = 2 * np.arctan(np.exp(-Dstar.Keta))
        pi_theta = 2 * np.arctan(np.exp(-Dstar.pieta))
        Dstar = Dstar[(Dstar.Kdz < 0.5/np.sin(K_theta)) & (Dstar.pidz < 0.5/np.sin(pi_theta))]

        # pis cuts
        Dstar = Dstar[Dstar.pisptr > 0.3]
        Dstar = Dstar[Dstar.pischir < 3]
        Dstar = Dstar[Dstar.pisnValid > 2]

        # D0 of Dstar cuts
        Dstar = Dstar[Dstar.D0cosphi > self.config['dstar_d0_cosphi']]
        Dstar = Dstar[(Dstar.D0mass < D0_PDG_MASS + self.config['dstar_d0_mass']) & (Dstar.D0mass > D0_PDG_MASS - self.config['dstar_d0_mass'])]
        Dstar = Dstar[Dstar.D0pt > self.config['dstar_d0_pt']]
        Dstar = Dstar[Dstar.D0dlSig > self.config['dstar_d0_dlSig']]

        ############### Trigger
        Dimu = Dimu[trigger]
        Muon = Muon[trigger]
        Dstar = Dstar[trigger]
        PVtx = PVtx[trigger]

        ############### Association
        DimuDstar = association(Dimu, Dstar)
        none_cut = remove_none(DimuDstar.slot0.pt)
        DimuDstar = DimuDstar[none_cut]

        Dstar = Dstar[Dstar.associationIdx > -1]
        MuonDstar = ak.zip({'0': Muon[Dstar.associationIdx], '1': Dstar})
        MuonDstar = MuonDstar[none_cut]

        arg_sort = ak.argsort(DimuDstar['cand'].pt, axis=1, ascending=False)
        DimuDstar = DimuDstar[arg_sort]
        MuonDstar = MuonDstar[arg_sort]

        Dimu_hists = {
            'pt': Hist.new.Regular(100, 0, 300, name='pt', label=r"$p_{T,\mu\mu}$ [GeV]").Weight(),
            'eta': Hist.new.Regular(60, -3, 3, name="eta", label=r"$\eta_{\mu\mu}$").Weight(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{\mu\mu}$").Weight(),
            'mass': Hist.new.Regular(100, 8.6, 11, name="mass", label=r"$m_{\mu\mu}$ [GeV]").Weight(),
            'rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$y_{\mu\mu}$").Weight(),
            'dl': Hist.new.Regular(100, -1, 1, name="dl", label=r"decay length $\mu\mu$").Weight(),
            'dlSig': Hist.new.Regular(100, -20, 20, name="dlSig", label=r"decay length significance $\mu\mu$").Weight(),
            'chi2': Hist.new.Regular(100, 0, 10, name="chi2", label=r"$\mu\mu$ vtx fit $\chi^2$").Weight(),
            'cosphi': Hist.new.Regular(60, -1, 1, name="cosphi", label=r"cos of pointing angle $\mu\mu$").Weight(),
        }

        Dstar_hists = {
            'pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^*}$ [GeV]").Weight(),
            'eta': Hist.new.Regular(60, -3, 3, name="eta", label=r"$\eta_{D^*}$").Weight(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^*}$").Weight(),
            'deltamr': Hist.new.Regular(50, 0.138, 0.162, name="deltam", label=r"$\Delta m_{D^*}$ [GeV]").Weight(),
            'rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$y_{D^*}$").Weight(),
            'D0cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0$ of $D^*$").Weight(),
            'D0dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0$ of $D^*$").Weight(),
            'D0dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0$ of $D^*$").Weight(),
            'D0pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^0}$ [GeV]").Weight(),
            'D0eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^0}$").Weight(),
            'D0phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^0}$").Weight(),
            'D0mass': Hist.new.Regular(100, 1.82, 1.9, name="mass", label=r"$m_{D^0}$ [GeV]").Weight(),
        }
            
        DimuDstar_hists = {
            'pt': Hist.new.Regular(100, 0, 120, name='pt', label=r"$p_{T,\Upsilon D^*}$ [GeV]").Weight(),
            'eta': Hist.new.Regular(60, -6, 6, name="eta", label=r"$\eta_{\Upsilon D^*}$").Weight(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{\Upsilon D^*}$").Weight(),
            'mass': Hist.new.Regular(100, 8.6, 70, name="mass", label=r"$m_{\Upsilon D^*}$ [GeV]").Weight(),
            'rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$\eta_{\Upsilon D^*}$").Weight(),
            'dimu_mass': Hist.new.Regular(100, 8.6, 11, name="mass", label=r"$m_{\mu\mu}$ [GeV]").Weight(),
            'dimu_pt': Hist.new.Regular(100, 0, 300, name='pt', label=r"$p_{T,\mu\mu}$ [GeV]").Weight(),
            'dimu_eta': Hist.new.Regular(60, -3, 3, name="eta", label=r"$\eta_{\mu\mu}$").Weight(),
            'dimu_phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{\mu\mu}$").Weight(),
            'dimu_rap': Hist.new.Regular(60, -2.5, 2.5, name="rap", label=r"$y_{\mu\mu}$").Weight(),
            'dimu_dl': Hist.new.Regular(100, -1, 1, name="dl", label=r"decay length $\mu\mu$").Weight(),
            'dimu_dlSig': Hist.new.Regular(100, -20, 20, name="dlSig", label=r"decay length significance $\mu\mu$").Weight(),
            'dimu_chi2': Hist.new.Regular(100, 0, 10, name="chi2", label=r"$\mu\mu$ vtx fit $\chi^2$").Weight(),
            'dimu_cosphi': Hist.new.Regular(60, -1, 1, name="cosphi", label=r"cos of pointing angle $\mu\mu$").Weight(),
            'dstar_deltamr': Hist.new.Regular(50, 0.138, 0.162, name="deltam", label=r"$\Delta m_{D^*}$ [GeV]").Weight(),
            'dstar_pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^*}$ [GeV]").Weight(),
            'dstar_eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^*}$").Weight(),
            'dstar_phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^*}$").Weight(),
            'dstar_rap': Hist.new.Regular(60, -4, 4, name="rap", label=r"$y_{D^*}$").Weight(),
            'dstar_D0pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,D^0}$ [GeV]").Weight(),
            'dstar_D0eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{D^0}$").Weight(),
            'dstar_D0phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{D^0}$").Weight(),
            'dstar_D0mass': Hist.new.Regular(100, 1.82, 1.9, name="mass", label=r"$m_{D^0}$ [GeV]").Weight(),
            'dstar_D0cosphi': Hist.new.Regular(60, 0.95, 1, name="cosphi", label=r"cos of pointing angle $D^0$ of $D^*$").Weight(),
            'dstar_D0dl': Hist.new.Regular(100, 0, 2, name="dl", label=r"decay length $D^0$ of $D^*$").Weight(),
            'dstar_D0dlSig': Hist.new.Regular(100, 0, 20, name="dlSig", label=r"decay length significance $D^0$ of $D^*$").Weight(),
            'dstar_associationchi2': Hist.new.Regular(100, 0, 5, name="assochi2", label=r"Association $\chi^2$").Weight(),
            'dstar_associationProb': Hist.new.Regular(100, 0, 1, name="assoprob", label=r"Association probability").Weight(),
            'deltarap': Hist.new.Regular(60, -5, 5, name="deltarap", label=r"$\Delta y_{\Upsilon D^*}$").Weight(),
            'deltapt': Hist.new.Regular(100, -30, 80, name='deltapt', label=r"$\Delta p_{T,\Upsilon D^*}$ [GeV]").Weight(),
            'deltaeta': Hist.new.Regular(60, -6, 6, name="deltaeta", label=r"$\Delta \eta_{\Upsilon D^*}$").Weight(),
            'deltaphi': Hist.new.Regular(60, 0, np.pi, name="deltaphi", label=r"$\Delta \phi_{\Upsilon D^*}$").Weight(),
        }

        none_cut = remove_none(Dimu.pt)
        Dimu = Dimu[none_cut]
        Muon = Muon[none_cut]
        #Dstar = Dstar

        dimudstar_weight = get_weight(evaluator, MuonDstar.slot0, PVtx)
        dimu_weight = get_weight(evaluator, Muon, PVtx)
        dstar_weight = np.repeat(evaluator['pileup'](ak.num(PVtx)), ak.num(Dstar))

        """ for i0 in range(len(Dimu)):
            if len(Dimu[i0]) == 0: continue
            print(i0)
            for i1 in range(len(Dimu[i0])):
                print(f'nPVtx = {len(PVtx[i0])}, Muon0 = ({Muon[i0][i1].slot0.pt:.2f}, {np.absolute(Muon[i0][i1].slot0.eta):.2f}), Muon1 = ({Muon[i0][i1].slot1.pt:.2f}, {np.absolute(Muon[i0][i1].slot1.eta):.2f}), dimu_weight = {dimu_weight[i0][i1]}')
                print(Dimu.pt[i0][i1]) """

        [Dimu_hists[h].fill(ak.flatten(Dimu[h]), weight=ak.flatten(dimu_weight)) for h in Dimu_hists]
        [Dstar_hists[h].fill(ak.flatten(Dstar[h]), weight=dstar_weight) for h in Dstar_hists]

        DimuDstar_hists['pt'].fill(ak.flatten(DimuDstar.cand.pt), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['eta'].fill(ak.flatten(DimuDstar.cand.eta), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['phi'].fill(ak.flatten(DimuDstar.cand.phi), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['mass'].fill(ak.flatten(DimuDstar.cand.mass), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['rap'].fill(ak.flatten(DimuDstar.rap), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['deltapt'].fill(ak.flatten(DimuDstar.deltapt), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['deltarap'].fill(ak.flatten(DimuDstar.deltarap), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['deltaeta'].fill(ak.flatten(DimuDstar.deltaeta), weight=ak.flatten(dimudstar_weight))
        DimuDstar_hists['deltaphi'].fill(ak.flatten(DimuDstar.deltaphi), weight=ak.flatten(dimudstar_weight))
        
        [DimuDstar_hists[h].fill(ak.flatten(DimuDstar.slot0[h[h.find('_')+1:]]), weight=ak.flatten(dimudstar_weight)) for h in DimuDstar_hists if 'dimu' in h]
        [DimuDstar_hists[h].fill(ak.flatten(DimuDstar.slot1[h[h.find('_')+1:]]), weight=ak.flatten(dimudstar_weight)) for h in DimuDstar_hists if 'dstar' in h]

        filename = f"output/mc_hists/{self.year}{f[f.rfind('/'):]}".replace('root', 'hists')

        hists = {
            'Dimu': Dimu_hists,
            'Dstar': Dstar_hists,
            'DimuDstar': DimuDstar_hists,
        }

        save(hists, filename)