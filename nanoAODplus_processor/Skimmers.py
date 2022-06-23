import awkward as ak
import numpy as np
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
                'D0dlSig': Dstar_D0_acc['D0dlSig'].value,
                'D0pt': Dstar_D0_acc['D0pt'].value,
                'D0eta': Dstar_D0_acc['D0eta'].value,
                'D0mass': Dstar_D0_acc['D0mass'].value,
                'wrg_chg': Dstar_acc['wrg_chg'].value,}, with_name='PtEtaPhiMCandidate')

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
                'deltarap': DimuDstar_acc['deltarap'].value,
                'deltapt': DimuDstar_acc['deltapt'].value,
                'deltaeta': DimuDstar_acc['deltaeta'].value,
                'deltaphi': DimuDstar_acc['deltaphi'].value,
                'is_ups': DimuDstar_acc['Dimu']['is_ups'].value,
                'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,
                'dimu_vtxIdx': DimuDstar_acc['Dimu']['vtxIdx'].value,
                'dstar_vtxIdx': DimuDstar_acc['Dstar']['vtxIdx'].value,
                'd': DimuDstar_acc['d'].value,
                'derr': DimuDstar_acc['derr'].value,
                'dSig': DimuDstar_acc['dSig'].value,}, with_name='PtEtaPhiMCandidate')

        # Unflatten
        Dimu = ak.unflatten(Dimu, Dimu_acc['nDimu'].value)
        Dstar = ak.unflatten(Dstar, Dstar_acc['nDstar'].value)
        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)

        #Apply the Trigger and the cuts in the config
        HLT = triggers[self.config['trigger'][self.year]].value
        Dimu = Dimu[HLT]
        Dstar = Dstar[HLT]
        DimuDstar = DimuDstar[HLT]
        
        Dimu = Dimu[Dimu.pt > self.config['limits']['Upsilon_pt']]
        Dimu = Dimu[np.absolute(Dimu.rap) < self.config['limits']['Upsilon_rap']]
        Dimu = Dimu[Dimu.is_ups]
        DimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_pt']]
        DimuDstar = DimuDstar[np.absolute(DimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
        DimuDstar = DimuDstar[DimuDstar.is_ups]

        Dstar = Dstar[Dstar.pt > self.config['limits']['Dstar_pt']]
        Dstar = Dstar[np.absolute(Dstar.rap) < self.config['limits']['Dstar_rap']]
        Dstar = Dstar[Dstar.D0pt > self.config['limits']['Dstar_D0_pt']]
        Dstar = Dstar[Dstar.D0cosphi > self.config['limits']['Dstar_D0_cosphi']]
        Dstar = Dstar[(Dstar.D0mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (Dstar.D0mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
        Dstar = Dstar[Dstar.D0dlSig > self.config['limits']['Dstar_D0_dlSig']]
        Dstar = Dstar[~Dstar.wrg_chg]
        DimuDstar = DimuDstar[DimuDstar.dstar_pt > self.config['limits']['Dstar_pt']]
        DimuDstar = DimuDstar[np.absolute(DimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
        DimuDstar = DimuDstar[DimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
        DimuDstar = DimuDstar[DimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
        DimuDstar = DimuDstar[(DimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (DimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
        DimuDstar = DimuDstar[DimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
        DimuDstar = DimuDstar[DimuDstar.dimu_vtxIdx == DimuDstar.dstar_vtxIdx]

        DimuDstar = DimuDstar[~DimuDstar.wrg_chg]

        folder = f[f.rfind('/'):f.rfind('_')]
        filename = f"output/RunII_trigger_processed/{self.year}{folder}{f[f.rfind('/'):]}"
        
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
                'deltarap': DimuDstar_acc['deltarap'].value,
                'deltapt': DimuDstar_acc['deltapt'].value,
                'deltaeta': DimuDstar_acc['deltaeta'].value,
                'deltaphi': DimuDstar_acc['deltaphi'].value,
                'is_ups': DimuDstar_acc['Dimu']['is_ups'].value,
                'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,
                'dimu_vtxIdx': DimuDstar_acc['Dimu']['vtxIdx'].value,
                'dstar_vtxIdx': DimuDstar_acc['Dstar']['vtxIdx'].value,
                'd': DimuDstar_acc['d'].value,
                'derr': DimuDstar_acc['derr'].value,
                'dSig': DimuDstar_acc['dSig'].value,}, with_name='PtEtaPhiMCandidate')

        # Unflatten
        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)

        #Apply the Trigger and the cuts in the config
        HLT = triggers[self.config['trigger'][self.year]].value
        DimuDstar = DimuDstar[HLT]
        
        for var in self.fom:
            DimuDstar = DimuDstar[DimuDstar.dimu_pt > self.config['limits']['Upsilon_pt']]
            DimuDstar = DimuDstar[np.absolute(DimuDstar.dimu_rap) < self.config['limits']['Upsilon_rap']]
            DimuDstar = DimuDstar[DimuDstar.is_ups]

            DimuDstar = DimuDstar[DimuDstar.dstar_pt > self.config['limits']['Dstar_pt']]
            DimuDstar = DimuDstar[np.absolute(DimuDstar.dstar_rap) < self.config['limits']['Dstar_rap']]
            DimuDstar = DimuDstar[(DimuDstar.dstar_d0_mass < D0_PDG_MASS + self.config['limits']['Dstar_D0_mass']) & (DimuDstar.dstar_d0_mass > D0_PDG_MASS - self.config['limits']['Dstar_D0_mass'])]
            DimuDstar = DimuDstar[~DimuDstar.wrg_chg]
            if var == 'D0pt': 
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
                DimuDstar = DimuDstar[DimuDstar.dimu_vtxIdx == DimuDstar.dstar_vtxIdx]
                for i in self.fom[var]:
                    DimuDstar = DimuDstar[DimuDstar.dstar_d0_pt > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(DimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })

                    save(output, filename)

            elif var == 'cosphi': 
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
                DimuDstar = DimuDstar[DimuDstar.dimu_vtxIdx == DimuDstar.dstar_vtxIdx]
                for i in self.fom[var]:
                    DimuDstar = DimuDstar[DimuDstar.dstar_d0_cosphi > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(DimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })

                    save(output, filename)

            elif var == 'dlSig': 
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
                DimuDstar = DimuDstar[DimuDstar.dimu_vtxIdx == DimuDstar.dstar_vtxIdx]
                for i in self.fom[var]:
                    DimuDstar = DimuDstar[DimuDstar.dstar_d0_dlSig > i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(DimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })

                    save(output, filename)

            elif var == 'dSig': 
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_pt > self.config['limits']['Dstar_D0_pt']]
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_cosphi > self.config['limits']['Dstar_D0_cosphi']]
                DimuDstar = DimuDstar[DimuDstar.dstar_d0_dlSig > self.config['limits']['Dstar_D0_dlSig']]
                for i in self.fom[var]:
                    DimuDstar = DimuDstar[DimuDstar.dSig < i]
                    folder = f[f.rfind('/'):f.rfind('_')]
                    filename = f"output/fom/{self.year}{folder}/{var}/{str(i).replace('.', 'p')}{f[f.rfind('/'):f.rfind('_')]}{var}{f[f.rfind('_'):]}"
                    
                    pDimuDstar_acc = build_acc(DimuDstar)

                    output = processor.dict_accumulator({
                        'DimuDstar': pDimuDstar_acc
                    })
                    save(output, filename)
                
