import awkward as ak
import numpy as np
import coffea.processor as processor
from coffea.util import save

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

import random, yaml

from tools.collections import *
from tools.utils import *

D0_PDG_MASS = 1.864

class GenParticleProcessor(processor.ProcessorABC):
    def __init__(self, analyzer_name, year):
        self.analyzer_name = analyzer_name
        self.year = year

        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        
        # test if there are any events in the file
        if len(events) == 0:
            return output
        
        ############### Get All the interesting candidates from NTuples
        Dimu = ak.zip({**get_vars_dict(events, dimu_cols)}, with_name="PtEtaPhiMCandidate")
        Muon = ak.zip({**get_vars_dict(events, muon_cols)}, with_name="PtEtaPhiMCandidate")
        D0 = ak.zip({'mass': events.D0_mass12, **get_vars_dict(events, d0_cols)}, with_name="PtEtaPhiMCandidate")
        Dstar = ak.zip({'mass': (events.Dstar_D0mass + events.Dstar_deltamr),
                        'charge': events.Dstar_pischg,
                        **get_vars_dict(events, dstar_cols)}, 
                        with_name="PtEtaPhiMCandidate")
        HLT = ak.zip({**get_hlt(events, hlt_cols[self.year])})
        PVtx = ak.zip({**get_vars_dict(events, pvtx_cols)})
        GenPart = ak.zip({**get_vars_dict(events, gen_part_cols)}, with_name="PtEtaPhiMCandidate")

        ############### Load config files
        with open("config/skim_trigger.yaml", 'r') as f:
            skim_trigger = yaml.load(f, Loader=yaml.FullLoader)
            trigger = skim_trigger['trigger']
            cuts = skim_trigger['limits']

        ############### GenParticles processing
        # Gen_Muon
        GenPart_Muon = GenPart[(np.absolute(GenPart.pdgId) == 13) & (GenPart.genPartIdxMother > -1)]
        gen_ups = (GenPart_Muon.parpdgId == 553) | (GenPart_Muon.parpdgId == 100553) | (GenPart_Muon.parpdgId == 200553)
        GenPart_Muon = GenPart_Muon[gen_ups]
        GenPart_Muon = GenPart_Muon[GenPart_Muon.pt > 3]
        GenPart_Muon = ak.combinations(GenPart_Muon, 2)
        GenPart_Muon = GenPart_Muon[GenPart_Muon.slot0.genPartIdxMother == GenPart_Muon.slot1.genPartIdxMother]

        # Gen_Upsilon
        GenPart_Ups = GenPart[(GenPart_Muon.slot0.genPartIdxMother)]
        GenPart_Ups = GenPart_Ups[GenPart_Ups.pt > cuts['Upsilon_pt']]
        GenPart_Ups = GenPart_Ups[np.absolute(GenPart_Ups.eta) < cuts['Upsilon_eta']] # change to eta cut

        leading_mu = (GenPart_Muon.slot0.pt > GenPart_Muon.slot1.pt)
        GenPart_Muon_lead = ak.where(leading_mu, GenPart_Muon.slot0, GenPart_Muon.slot1)
        GenPart_Muon_trail = ak.where(~leading_mu, GenPart_Muon.slot0, GenPart_Muon.slot1)

        # Gen_Ds
        GenPart_D0 = GenPart[(np.absolute(GenPart.pdgId) == 421) & (np.absolute(GenPart.parpdgId) == 413)]
        GenPart_Dstar = GenPart[GenPart_D0.genPartIdxMother]
        GenPart_Dstar = GenPart_Dstar[GenPart_Dstar.pt > cuts['Dstar_pt']]
        GenPart_Dstar = GenPart_Dstar[np.absolute(GenPart_Dstar.eta) < cuts['Dstar_eta']]
        
        ############### Dimu cuts charge = 0, mass cuts and chi2...
        Dimu = Dimu[Dimu.charge == 0]
        Dimu = Dimu[((Dimu.mass > 8.5) & (Dimu.mass < 11.5))]

        ############### Get the Muons from Dimu, for cuts in their params
        Muon = ak.zip({'0': Muon[Dimu.t1muIdx], '1': Muon[Dimu.t2muIdx]})

        # SoftId and Global Muon cuts
        soft_id = (Muon.slot0.softId > 0) & (Muon.slot1.softId > 0)
        Dimu = Dimu[soft_id]
        Muon = Muon[soft_id]

        # pt and eta cuts
        muon_pt_cut = (Muon.slot0.pt > 3) & (Muon.slot1.pt > 3)
        Dimu = Dimu[muon_pt_cut]
        Muon = Muon[muon_pt_cut]

        muon_eta_cut = (np.absolute(Muon.slot0.eta) <= 2.4) & (np.absolute(Muon.slot1.eta) <= 2.4)
        Dimu = Dimu[muon_eta_cut]
        Muon = Muon[muon_eta_cut]

        Dimu['is_ups'] = (Dimu.mass > 8.5) & (Dimu.mass < 11.5)
        Dimu['is_jpsi'] = (Dimu.mass > 2.95) & (Dimu.mass < 3.25)

        ############### Cuts for D0
        D0 = D0[~D0.hasMuon]
        D0 = D0[(D0.t1pt > 0.8) & (D0.t2pt > 0.8)]
        D0 = D0[(D0.t1chindof < 2.5) & (D0.t2chindof < 2.5)]
        D0 = D0[(D0.t1nValid > 4) & (D0.t2nValid > 4) & (D0.t1nPix > 1) & (D0.t2nPix > 1)]
        D0 = D0[(D0.t1dxy < 0.1) & (D0.t2dxy < 0.1)]
        D0 = D0[(D0.t1dz < 1.) & (D0.t2dz < 1.)]

        # D0 cosphi
        D0 = D0[D0.cosphi > 0.99]

        # D0 dl Significance
        D0 = D0[D0.dlSig > 5.]

        # D0 pt
        D0 = D0[D0.pt > 3.]

        ############### Cuts for Dstar

        # trks cuts
        Dstar = Dstar[~Dstar.hasMuon]
        Dstar = Dstar[(Dstar.Kpt > 0.5) & (Dstar.pipt > 0.5)]
        Dstar = Dstar[(Dstar.Kchindof < 2.5) & (Dstar.pichindof < 2.5)]
        Dstar = Dstar[(Dstar.KnValid > 4) & (Dstar.pinValid > 4) & (Dstar.KnPix > 1) & (Dstar.pinPix > 1)]
        Dstar = Dstar[(Dstar.Kdxy < 0.1) & (Dstar.pidxy < 0.1)]
        Dstar = Dstar[(Dstar.Kdz < 1) & (Dstar.pidz < 1)]

        # pis cuts
        Dstar = Dstar[Dstar.pispt > 0.3]
        Dstar = Dstar[Dstar.pischindof < 3]
        Dstar = Dstar[Dstar.pisnValid > 2]

        # D0 of Dstar cuts
        Dstar = Dstar[Dstar.D0cosphi > 0.99]
        Dstar = Dstar[(Dstar.D0mass < D0_PDG_MASS + 0.025) & (Dstar.D0mass > D0_PDG_MASS - 0.025)]
        Dstar = Dstar[Dstar.D0pt > cuts['Dstar_D0_pt']]
        Dstar = Dstar[Dstar.D0eta < 2.5]
        Dstar = Dstar[Dstar.D0dlSig > 3]

        Dstar['wrg_chg'] = (Dstar.Kchg == Dstar.pichg)
        if cuts['right_charge_only']:
            Dstar = Dstar[~Dstar.wrg_chg]

        ############### Dimu + OpenCharm associations

        DimuDstar = association(Dimu, Dstar)

        ############### Leading and Trailing muon separation
        leading_mu = (Muon.slot0.pt > Muon.slot1.pt)
        Muon_lead = ak.where(leading_mu, Muon.slot0, Muon.slot1)
        Muon_trail = ak.where(~leading_mu, Muon.slot0, Muon.slot1)

        ############### Create the accumulators to save output
        muon_lead_acc = processor.dict_accumulator({})
        for var in Muon_lead.fields:
            muon_lead_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Muon_lead[var])))
        muon_lead_acc["nMuon"] = processor.column_accumulator(ak.to_numpy(ak.num(Muon_lead)))
        output["Muon_lead"] = muon_lead_acc

        muon_trail_acc = processor.dict_accumulator({})
        for var in Muon_trail.fields:
            muon_trail_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Muon_trail[var])))
        muon_trail_acc["nMuon"] = processor.column_accumulator(ak.to_numpy(ak.num(Muon_trail)))
        output["Muon_trail"] = muon_trail_acc

        dimu_acc = processor.dict_accumulator({})
        for var in Dimu.fields:
            if (var.startswith('t')): continue
            dimu_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dimu[var])))
        dimu_acc["nDimu"] = processor.column_accumulator(ak.to_numpy(ak.num(Dimu)))
        output["Dimu"] = dimu_acc

        D0_acc = processor.dict_accumulator({})
        D0_trk_acc = processor.dict_accumulator({})
        for var in D0.fields:
            if (var.startswith('t')):
                D0_trk_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(D0[var])))
            else:
                D0_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(D0[var])))
        D0_acc["nD0"] = processor.column_accumulator(ak.to_numpy(ak.num(D0)))
        output["D0"] = D0_acc
        output["D0_trk"] = D0_trk_acc

        Dstar_acc = processor.dict_accumulator({})
        Dstar_D0_acc = processor.dict_accumulator({})
        Dstar_trk_acc = processor.dict_accumulator({})
        for var in Dstar.fields:
            if var.startswith('D0'):
                Dstar_D0_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dstar[var])))
            elif (var.startswith('K') or var.startswith('pi')):
                Dstar_trk_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dstar[var])))
            else:
                Dstar_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dstar[var])))
        Dstar_acc["nDstar"] = processor.column_accumulator(ak.to_numpy(ak.num(Dstar)))
        output["Dstar"] = Dstar_acc
        output["Dstar_D0"] = Dstar_D0_acc
        output["Dstar_trk"] = Dstar_trk_acc

        DimuDstar_acc = processor.dict_accumulator({})
        DimuDstar_acc['Dimu'] = processor.dict_accumulator({})
        DimuDstar_acc['Dstar'] = processor.dict_accumulator({})
        for var in DimuDstar.fields:
            if (var == '0') or (var =='1'):
                continue
            elif var == 'cand':
                for i0 in DimuDstar[var].fields:
                    DimuDstar_acc[i0] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar[var][i0])))
            else:
                DimuDstar_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar[var])))

        for var in DimuDstar.slot0.fields:
            DimuDstar_acc['Dimu'][var] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar.slot0[var])))

        for var in DimuDstar.slot1.fields:
            DimuDstar_acc['Dstar'][var] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar.slot1[var])))
        DimuDstar_acc['nDimuDstar'] = processor.column_accumulator(ak.to_numpy(ak.num(DimuDstar)))
        output['DimuDstar'] = DimuDstar_acc

        file_hash = str(random.getrandbits(128)) + str(len(events))
        save(output, "output/" + self.analyzer_name + "/" + self.analyzer_name + "_" + file_hash + ".coffea")

        # return dummy accumulator
        return processor.dict_accumulator({
                'cutflow': output['cutflow']
        })

    def postprocess(self, accumulator):
        return accumulator