import awkward as ak
import numpy as np
import uproot
from coffea.util import load, save
from coffea import processor

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

from tools.utils import *

class Skimmer:
    def __init__(self, config, year):
        if config is None or year is None:
            print("Configuration not complete!!!")
        self.config = config
        self.year = year

    def process(self, f):
        # Load the accumulator
        acc = load(f)
        
        Dimu_acc = acc['Dimu']
        Dstar_acc = acc['Dstar']
        Dstar_D0_acc = acc['Dstar_D0']
        Dstar_trk_acc = acc['Dstar_trk']
        DimuDstar_acc = acc['DimuDstar']
        triggers = acc['triggers']

        # Loading variables as awkward arrays
        Dimu = ak.zip({
                'pt': Dimu_acc['pt'].value,
                'eta': Dimu_acc['eta'].value,
                'phi': Dimu_acc['phi'].value,
                'mass': Dimu_acc['mass'].value,
                'rap': Dimu_acc['rap'].value,
                'dl': Dimu_acc['dl'].value,
                'dlSig': Dimu_acc['dlSig'].value,
                'chi2': Dimu_acc['chi2'].value,
                'cosphi': Dimu_acc['cosphi'].value,
                'is_ups': Dimu_acc['is_ups'].value,}, with_name="PtEtaPhiMLorentzVector")

        Dstar = ak.zip({
                'pt': Dstar_acc['pt'].value,
                'eta': Dstar_acc['eta'].value,
                'phi': Dstar_acc['phi'].value,
                'mass': Dstar_acc['mass'].value,
                'rap': Dstar_acc['rap'].value,
                'charge': Dstar_acc['charge'].value,
                'deltam': Dstar_acc['deltam'].value,
                'deltamr': Dstar_acc['deltamr'].value,
                'D0cosphi': Dstar_D0_acc['D0cosphi'].value,
                'D0dl': Dstar_D0_acc['D0dl'].value,
                'D0dlSig': Dstar_D0_acc['D0dlSig'].value,
                'D0pt': Dstar_D0_acc['D0pt'].value,
                'D0eta': Dstar_D0_acc['D0eta'].value,
                'D0phi': Dstar_D0_acc['D0phi'].value,
                'D0mass': Dstar_D0_acc['D0mass'].value,
                'Kpt':  Dstar_trk_acc['Kpt'].value,
                'Keta':  Dstar_trk_acc['Keta'].value,
                'Kdxy': Dstar_trk_acc['Kdxy'].value,
                'Kdz': Dstar_trk_acc['Kdz'].value,
                'pipt': Dstar_trk_acc['pipt'].value,
                'pieta': Dstar_trk_acc['pieta'].value,
                'pidxy': Dstar_trk_acc['pidxy'].value,
                'pidz': Dstar_trk_acc['pidz'].value,
                'pispt': Dstar_trk_acc['pispt'].value,
                'piseta': Dstar_trk_acc['piseta'].value,
                'pisdxy': Dstar_trk_acc['pisdxy'].value,
                'pisdz': Dstar_trk_acc['pisdz'].value,
                'wrg_chg': Dstar_acc['wrg_chg'].value,}, with_name='PtEtaPhiMCandidate')

        DimuDstar_p4 = build_p4(DimuDstar_acc)
        DimuDstar = ak.zip({
                'pt': DimuDstar_p4.pt,
                'eta': DimuDstar_p4.eta,
                'phi': DimuDstar_p4.phi,
                'mass': DimuDstar_p4.mass,
                'rap': DimuDstar_p4.rap,
                'charge': DimuDstar_acc['charge'].value,
                'dimu_mass': DimuDstar_acc['Dimu']['mass'].value,
                'dimu_pt': DimuDstar_acc['Dimu']['pt'].value,
                'dimu_eta': DimuDstar_acc['Dimu']['eta'].value,
                'dimu_phi': DimuDstar_acc['Dimu']['phi'].value,
                'dimu_rap': DimuDstar_acc['Dimu']['rap'].value,
                'dimu_dl': DimuDstar_acc['Dimu']['dl'].value,
                'dimu_dlSig': DimuDstar_acc['Dimu']['dlSig'].value,
                'dimu_chi2': DimuDstar_acc['Dimu']['chi2'].value,
                'dimu_cosphi': DimuDstar_acc['Dimu']['cosphi'].value,
                'dstar_deltam': DimuDstar_acc['Dstar']['deltam'].value,
                'dstar_deltamr': DimuDstar_acc['Dstar']['deltamr'].value,
                'dstar_pt': DimuDstar_acc['Dstar']['pt'].value,
                'dstar_eta': DimuDstar_acc['Dstar']['eta'].value,
                'dstar_phi': DimuDstar_acc['Dstar']['phi'].value,
                'dstar_rap': DimuDstar_acc['Dstar']['rap'].value,
                'dstar_d0_pt': DimuDstar_acc['Dstar']['D0pt'].value,
                'dstar_d0_eta': DimuDstar_acc['Dstar']['D0eta'].value,
                'dstar_d0_phi': DimuDstar_acc['Dstar']['D0phi'].value,
                'dstar_d0_mass': DimuDstar_acc['Dstar']['D0mass'].value,
                'dstar_d0_cosphi': DimuDstar_acc['Dstar']['D0cosphi'].value,
                'dstar_d0_dl': DimuDstar_acc['Dstar']['D0dl'].value,
                'dstar_d0_dlSig': DimuDstar_acc['Dstar']['D0dlSig'].value,
                'dstar_k_pt': DimuDstar_acc['Dstar']['Kpt'].value,
                'dstar_k_eta': DimuDstar_acc['Dstar']['Keta'].value,
                'dstar_k_dxy': DimuDstar_acc['Dstar']['Kdxy'].value,
                'dstar_k_dz': DimuDstar_acc['Dstar']['Kdz'].value,
                'dstar_pi_pt': DimuDstar_acc['Dstar']['pipt'].value,
                'dstar_pi_eta': DimuDstar_acc['Dstar']['pieta'].value,
                'dstar_pi_dxy': DimuDstar_acc['Dstar']['pidxy'].value,
                'dstar_pi_dz': DimuDstar_acc['Dstar']['pidz'].value,
                'dstar_pis_pt': DimuDstar_acc['Dstar']['pisptr'].value,
                'dstar_pis_eta': DimuDstar_acc['Dstar']['pisetar'].value,
                'dstar_pis_dxy': DimuDstar_acc['Dstar']['pisdxy'].value,
                'dstar_pis_dz': DimuDstar_acc['Dstar']['pisdz'].value,
                'dstar_asso_chi2': DimuDstar_acc['Dstar']['associationchi2'].value,
                'dstar_asso_prob': DimuDstar_acc['Dstar']['associationProb'].value,
                'deltarap': DimuDstar_acc['deltarap'].value,
                'deltapt': DimuDstar_acc['deltapt'].value,
                'deltaeta': DimuDstar_acc['deltaeta'].value,
                'deltaphi': DimuDstar_acc['deltaphi'].value,
                'is_ups': DimuDstar_acc['Dimu']['is_ups'].value,
                'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,
                'dimu_vtxIdx': DimuDstar_acc['Dimu']['vtxIdx'].value,
                'dstar_vtxIdx': DimuDstar_acc['Dstar']['vtxIdx'].value,}, with_name='PtEtaPhiMCandidate')

        # Unflatten
        Dimu = ak.unflatten(Dimu, Dimu_acc['nDimu'].value)
        Dstar = ak.unflatten(Dstar, Dstar_acc['nDstar'].value)
        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)

        #Apply the Trigger and the cuts in the config
        HLT = triggers[self.config['trigger'][self.year]].value
        Dimu = Dimu[HLT]
        Dstar = Dstar[HLT]
        DimuDstar = DimuDstar[HLT]
        
        Dimu = Dimu[Dimu.pt > self.config['limits']['Upsilon_min_pt']]
        Dimu = Dimu[Dimu.pt < self.config['limits']['Upsilon_max_pt']]
        Dimu = Dimu[np.absolute(Dimu.rap) < self.config['limits']['Upsilon_rap']]
        Dimu = Dimu[Dimu.is_ups]
        DimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_min_pt']]
        DimuDstar = DimuDstar[DimuDstar.dimu_pt < self.config['limits']['Upsilon_max_pt']]
        DimuDstar = DimuDstar[np.absolute(DimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
        DimuDstar = DimuDstar[DimuDstar.is_ups]

        Dstar = Dstar[Dstar.pt > self.config['limits']['Dstar_min_pt']]
        Dstar = Dstar[Dstar.pt < self.config['limits']['Dstar_max_pt']]
        Dstar = Dstar[np.absolute(Dstar.rap) < self.config['limits']['Dstar_rap']]
        Dstar = Dstar[Dstar.D0pt > self.config['limits']['Dstar_D0_pt']]
        Dstar = Dstar[Dstar.D0cosphi > self.config['limits']['Dstar_D0_cosphi']]
        Dstar = Dstar[(Dstar.D0mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (Dstar.D0mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
        Dstar = Dstar[Dstar.D0dlSig > self.config['limits']['Dstar_D0_dlSig']]
        Dstar = Dstar[Dstar.Kpt > self.config['limits']['Dstar_track_pt']]
        Dstar = Dstar[Dstar.pipt > self.config['limits']['Dstar_track_pt']]
        Dstar = Dstar[Dstar.Kdxy < 0.5]
        Dstar.Ktheta = 2 * np.arctan(np.exp(-Dstar.Keta))
        Dstar = Dstar[Dstar.Kdz < 0.5/np.sin(Dstar.Ktheta)]
        Dstar = Dstar[Dstar.pidxy < 0.5]
        Dstar.pitheta = 2 * np.arctan(np.exp(-Dstar.pieta))
        Dstar = Dstar[Dstar.pidz < 0.5/np.sin(Dstar.pitheta)]
        #Dstar.pistheta = 2 * np.arctan(np.exp(-Dstar.piseta))
        Dstar = Dstar[~Dstar.wrg_chg]
        DimuDstar = DimuDstar[DimuDstar.dstar_pt > self.config['limits']['Dstar_min_pt']]
        DimuDstar = DimuDstar[DimuDstar.dstar_pt < self.config['limits']['Dstar_max_pt']]
        DimuDstar = DimuDstar[np.absolute(DimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
        #DimuDstar = DimuDstar[np.absolute(DimuDstar.dstar_eta) < self.config['limits']['Dstar_rap']]
        DimuDstar = DimuDstar[DimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
        DimuDstar = DimuDstar[DimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
        DimuDstar = DimuDstar[(DimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (DimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
        DimuDstar = DimuDstar[DimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
        DimuDstar = DimuDstar[DimuDstar.dstar_k_pt > self.config['limits']['Dstar_track_pt']]
        DimuDstar = DimuDstar[DimuDstar.dstar_pi_pt > self.config['limits']['Dstar_track_pt']]
        DimuDstar = DimuDstar[DimuDstar.dstar_k_dxy < 0.5]
        DimuDstar.dstar_k_theta = 2 * np.arctan(np.exp(-DimuDstar.dstar_k_eta))
        DimuDstar = DimuDstar[DimuDstar.dstar_k_dz < 0.5/np.sin(DimuDstar.dstar_k_theta)]
        #DimuDstar = DimuDstar[DimuDstar.dstar_k_dz < 1]
        DimuDstar = DimuDstar[DimuDstar.dstar_pi_dxy < 0.5]
        DimuDstar.dstar_pi_theta = 2 * np.arctan(np.exp(-DimuDstar.dstar_pi_eta))
        DimuDstar = DimuDstar[DimuDstar.dstar_pi_dz < 0.5/np.sin(DimuDstar.dstar_pi_theta)]
        #DimuDstar = DimuDstar[DimuDstar.dstar_pi_dz < 1]
        DimuDstar = DimuDstar[DimuDstar.dstar_asso_prob > self.config['limits']['Dstar_association_prob']]

        DimuDstar = DimuDstar[DimuDstar.mass > 18]
        DimuDstar = DimuDstar[~DimuDstar.wrg_chg]

        Dstar_D0 = ak.zip({
                'pt': DimuDstar.dstar_d0_pt,
                'eta': DimuDstar.dstar_d0_eta,
                'phi': DimuDstar.dstar_d0_phi,
                'mass': DimuDstar.dstar_d0_mass,
            }, with_name='PtEtaPhiMLorentzVector'
        )
        DimuDstar['dstar_d0_rap'] = np.log((Dstar_D0.energy + Dstar_D0.z)/np.sqrt(Dstar_D0.mass2 + Dstar_D0.pt2))

        folder = f[f.rfind('/'):f.rfind('_')]
        filename = f"output/RunII_trigger_processed_vtxfit/{self.year}{folder}{f[f.rfind('/'):]}"
        
        file_eff = uproot.open(f'output/efficiency/efficiencies_{self.year}.root')
        hists_eff = {h[:h.find(';')]:file_eff[h].to_hist() for h in file_eff}

        eff, eff_err_up, eff_err_down, eff_asso, eff_asso_err_up, eff_asso_err_down, wgt = get_evt_eff(hists_eff, DimuDstar)
        eff = ak.unflatten(eff, ak.num(DimuDstar))
        eff_err_up = ak.unflatten(eff_err_up, ak.num(DimuDstar))
        eff_err_down = ak.unflatten(eff_err_down, ak.num(DimuDstar))
        eff_asso = ak.unflatten(eff_asso, ak.num(DimuDstar))
        eff_asso_err_up = ak.unflatten(eff_asso_err_up, ak.num(DimuDstar))
        eff_asso_err_down = ak.unflatten(eff_asso_err_down, ak.num(DimuDstar))
        wgt = ak.unflatten(wgt, ak.num(DimuDstar))
        DimuDstar['eff'] = eff
        DimuDstar['eff_err_up'] = eff_err_up
        DimuDstar['eff_err_down'] = eff_err_down
        DimuDstar['eff_asso'] = eff_asso
        DimuDstar['eff_asso_err_up'] = eff_asso_err_up
        DimuDstar['eff_asso_err_down'] = eff_asso_err_down
        DimuDstar['wgt'] = wgt

        pDimu_acc = build_acc(Dimu)
        pDstar_acc = build_acc(Dstar)
        pDimuDstar_acc = build_acc(DimuDstar)

        output = processor.dict_accumulator({
            'Dimu': pDimu_acc,
            'Dstar': pDstar_acc,
            'DimuDstar': pDimuDstar_acc
        })
        save(output, filename)

class FOM:
    def __init__(self, config, fom, year):
        if config is None or year is None:
            print("Configuration not complete!!!")
        self.config = config
        self.fom = fom
        self.year = year

    def process(self, f):
        # Load the accumulator
        acc = load(f)

        DimuDstar_acc = acc['DimuDstar']
        triggers = acc['triggers']

        DimuDstar_p4 = build_p4(DimuDstar_acc)
        DimuDstar = ak.zip({
                'pt': DimuDstar_p4.pt,
                'eta': DimuDstar_p4.eta,
                'phi': DimuDstar_p4.phi,
                'mass': DimuDstar_p4.mass,
                'charge': DimuDstar_acc['charge'].value,
                'dimu_mass': DimuDstar_acc['Dimu']['mass'].value,
                'dimu_pt': DimuDstar_acc['Dimu']['pt'].value,
                'dimu_eta': DimuDstar_acc['Dimu']['eta'].value,
                'dimu_phi': DimuDstar_acc['Dimu']['phi'].value,
                'dimu_rap': DimuDstar_acc['Dimu']['rap'].value,
                'dstar_deltam': DimuDstar_acc['Dstar']['deltam'].value,
                'dstar_deltamr': DimuDstar_acc['Dstar']['deltamr'].value,
                'dstar_pt': DimuDstar_acc['Dstar']['pt'].value,
                'dstar_eta': DimuDstar_acc['Dstar']['eta'].value,
                'dstar_phi': DimuDstar_acc['Dstar']['phi'].value,
                'dstar_rap': DimuDstar_acc['Dstar']['rap'].value,
                'dstar_d0_pt': DimuDstar_acc['Dstar']['D0pt'].value,
                'dstar_d0_eta': DimuDstar_acc['Dstar']['D0eta'].value,
                'dstar_d0_mass': DimuDstar_acc['Dstar']['D0mass'].value,
                'dstar_d0_cosphi': DimuDstar_acc['Dstar']['D0cosphi'].value,
                'dstar_d0_dlSig': DimuDstar_acc['Dstar']['D0dlSig'].value,
                'dstar_k_pt': DimuDstar_acc['Dstar']['Kpt'].value,
                'dstar_k_eta': DimuDstar_acc['Dstar']['Keta'].value,
                'dstar_k_dxy': DimuDstar_acc['Dstar']['Kdxy'].value,
                'dstar_k_dz': DimuDstar_acc['Dstar']['Kdz'].value,
                'dstar_pi_pt': DimuDstar_acc['Dstar']['pipt'].value,
                'dstar_pi_eta': DimuDstar_acc['Dstar']['pieta'].value,
                'dstar_pi_dxy': DimuDstar_acc['Dstar']['pidxy'].value,
                'dstar_pi_dz': DimuDstar_acc['Dstar']['pidz'].value,
                'dstar_pis_pt': DimuDstar_acc['Dstar']['pisptr'].value,
                'dstar_pis_eta': DimuDstar_acc['Dstar']['pisetar'].value,
                'dstar_pis_dxy': DimuDstar_acc['Dstar']['pisdxy'].value,
                'dstar_pis_dz': DimuDstar_acc['Dstar']['pisdz'].value,
                'dstar_asso_chi2': DimuDstar_acc['Dstar']['associationchi2'].value,
                'dstar_asso_prob': DimuDstar_acc['Dstar']['associationProb'].value,
                'deltarap': DimuDstar_acc['deltarap'].value,
                'deltapt': DimuDstar_acc['deltapt'].value,
                'deltaeta': DimuDstar_acc['deltaeta'].value,
                'deltaphi': DimuDstar_acc['deltaphi'].value,
                'is_ups': DimuDstar_acc['Dimu']['is_ups'].value,
                'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,
                'dimu_vtxIdx': DimuDstar_acc['Dimu']['vtxIdx'].value,
                'dstar_vtxIdx': DimuDstar_acc['Dstar']['vtxIdx'].value,}, with_name='PtEtaPhiMCandidate')

        # Unflatten
        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)

        #Apply the Trigger and the cuts in the config
        HLT = triggers[self.config['trigger'][self.year]].value
        DimuDstar = DimuDstar[HLT]
        
        for var in self.fom:
            if var == 'D0pt': 
                pDimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dimu_pt < self.config['limits']['Upsilon_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
                pDimuDstar = pDimuDstar[pDimuDstar.is_ups]

                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt > self.config['limits']['Dstar_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt < self.config['limits']['Dstar_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
                pDimuDstar = pDimuDstar[(pDimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (pDimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
                pDimuDstar = pDimuDstar[~pDimuDstar.wrg_chg]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_asso_prob > self.config['limits']['Dstar_association_prob']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dxy < 0.5]
                pDimuDstar.dstar_k_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_k_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dz < 0.5/np.sin(pDimuDstar.dstar_k_theta)]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dxy < 0.5]
                pDimuDstar.dstar_pi_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_pi_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dz < 0.5/np.sin(pDimuDstar.dstar_pi_theta)]
                for i in self.fom[var]:
                    pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_pt > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom_vtxfit/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(pDimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })

                    save(output, filename)

            elif var == 'cosphi': 
                pDimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dimu_pt < self.config['limits']['Upsilon_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
                pDimuDstar = pDimuDstar[pDimuDstar.is_ups]

                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt > self.config['limits']['Dstar_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt < self.config['limits']['Dstar_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
                pDimuDstar = pDimuDstar[(pDimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (pDimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
                pDimuDstar = pDimuDstar[~pDimuDstar.wrg_chg]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_asso_prob > self.config['limits']['Dstar_association_prob']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dxy < 0.5]
                pDimuDstar.dstar_k_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_k_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dz < 0.5/np.sin(pDimuDstar.dstar_k_theta)]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dxy < 0.5]
                pDimuDstar.dstar_pi_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_pi_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dz < 0.5/np.sin(pDimuDstar.dstar_pi_theta)]
                for i in self.fom[var]:
                    pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_cosphi > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom_vtxfit/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(pDimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })

                    save(output, filename)

            elif var == 'dlSig': 
                pDimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dimu_pt < self.config['limits']['Upsilon_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
                pDimuDstar = pDimuDstar[pDimuDstar.is_ups]

                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt > self.config['limits']['Dstar_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt < self.config['limits']['Dstar_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
                pDimuDstar = pDimuDstar[(pDimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (pDimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
                pDimuDstar = pDimuDstar[~pDimuDstar.wrg_chg]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_asso_prob > self.config['limits']['Dstar_association_prob']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dxy < 0.5]
                pDimuDstar.dstar_k_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_k_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dz < 0.5/np.sin(pDimuDstar.dstar_k_theta)]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dxy < 0.5]
                pDimuDstar.dstar_pi_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_pi_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dz < 0.5/np.sin(pDimuDstar.dstar_pi_theta)]
                for i in self.fom[var]:
                    pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_dlSig > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom_vtxfit/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(pDimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })

                    save(output, filename)

            elif var == 'association_prob': 
                pDimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dimu_pt < self.config['limits']['Upsilon_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
                pDimuDstar = pDimuDstar[pDimuDstar.is_ups]

                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt > self.config['limits']['Dstar_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt < self.config['limits']['Dstar_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
                pDimuDstar = pDimuDstar[(pDimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (pDimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
                pDimuDstar = pDimuDstar[~pDimuDstar.wrg_chg]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_pt > self.config['limits']['Dstar_track_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dxy < 0.5]
                pDimuDstar.dstar_k_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_k_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dz < 0.5/np.sin(pDimuDstar.dstar_k_theta)]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dxy < 0.5]
                pDimuDstar.dstar_pi_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_pi_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dz < 0.5/np.sin(pDimuDstar.dstar_pi_theta)]
                for i in self.fom[var]:
                    pDimuDstar = pDimuDstar[pDimuDstar.dstar_asso_prob > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom_vtxfit/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(pDimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })
                    save(output, filename)

            elif var == 'track_pt':
                pDimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dimu_pt < self.config['limits']['Upsilon_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
                pDimuDstar = pDimuDstar[pDimuDstar.is_ups]

                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt > self.config['limits']['Dstar_min_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pt < self.config['limits']['Dstar_max_pt']]
                pDimuDstar = pDimuDstar[np.absolute(pDimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
                pDimuDstar = pDimuDstar[(pDimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (pDimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
                pDimuDstar = pDimuDstar[~pDimuDstar.wrg_chg]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_asso_prob > self.config['limits']['Dstar_association_prob']]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dxy < 0.5]
                pDimuDstar.dstar_k_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_k_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_dz < 0.5/np.sin(pDimuDstar.dstar_k_theta)]
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dxy < 0.5]
                pDimuDstar.dstar_pi_theta = 2 * np.arctan(np.exp(-pDimuDstar.dstar_pi_eta))
                pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_dz < 0.5/np.sin(pDimuDstar.dstar_pi_theta)]
                for i in self.fom[var]:
                    pDimuDstar = pDimuDstar[pDimuDstar.dstar_k_pt > i]
                    pDimuDstar = pDimuDstar[pDimuDstar.dstar_pi_pt > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom_vtxfit/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(pDimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })
                    save(output, filename)
                
