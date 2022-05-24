import awkward as ak
import numpy as np
import coffea.processor as processor
from coffea.util import save

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

import random

from tools.collections import *
from tools.utils import *

D0_PDG_MASS = 1.864

class EventSelectorProcessor(processor.ProcessorABC):
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
        
        # test if there is any events in the file
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
        PVtx = ak.zip({**get_vars_dict(events, pvtx_cols)})
        HLT = ak.zip({**get_hlt(events, hlt_cols[self.year])})

        output['cutflow']['Number of events']  += len(events)
        output['cutflow']['Number of Dimu']    += ak.sum(ak.num(Dimu))
        output['cutflow']['Number of D0']      += ak.sum(ak.num(D0))
        output['cutflow']['Number of Dstar']   += ak.sum(ak.num(Dstar))

        ############### Dimu cuts charge = 0, mass cuts and chi2...
        Dimu = Dimu[Dimu.charge == 0]
        output['cutflow']['Dimu 0 charge'] += ak.sum(ak.num(Dimu))

        Dimu = Dimu[((Dimu.mass > 8.5) & (Dimu.mass < 11.5)) | ((Dimu.mass > 2.95) & (Dimu.mass < 3.25))]
        output['cutflow']['Quarkonia mass'] += ak.sum(ak.num(Dimu))

        ############### Get the Muons from Dimu, for cuts in their params
        Muon = ak.zip({'0': Muon[Dimu.t1muIdx], '1': Muon[Dimu.t2muIdx]})

        # SoftId and Global Muon cuts
        soft_id = (Muon.slot0.softId > 0) & (Muon.slot1.softId > 0)
        Dimu = Dimu[soft_id]
        Muon = Muon[soft_id]
        output['cutflow']['Dimu muon softId'] += ak.sum(ak.num(Dimu))

        """ global_muon = (Muon.slot0.isGlobal > 0) & (Muon.slot1.isGlobal > 0)
        Dimu = Dimu[global_muon]
        Muon = Muon[global_muon]
        output['cutflow']['Dimu muon global'] += ak.sum(ak.num(Dimu)) """

        # pt and eta cuts
        muon_pt_cut = (Muon.slot0.pt > 3) & (Muon.slot1.pt > 3)
        Dimu = Dimu[muon_pt_cut]
        Muon = Muon[muon_pt_cut]
        output['cutflow']['Dimu muon pt cut'] += ak.sum(ak.num(Dimu))

        muon_eta_cut = (np.absolute(Muon.slot0.eta) < 2.4) & (np.absolute(Muon.slot1.eta) < 2.4)
        Dimu = Dimu[muon_eta_cut]
        Muon = Muon[muon_eta_cut]
        output['cutflow']['Dimu muon eta cut'] += ak.sum(ak.num(Dimu))

        Dimu['is_ups'] = (Dimu.mass > 8.5) & (Dimu.mass < 11.5)
        Dimu['is_jpsi'] = (Dimu.mass > 2.95) & (Dimu.mass < 3.25)

        ############### Cuts for Dstar

        # trks cuts
        Dstar = Dstar[~Dstar.hasMuon]
        output['cutflow']['Dstar trk muon cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[(Dstar.Kpt > 0.5) & (Dstar.pipt > 0.5)]
        output['cutflow']['Dstar trk pt cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[(Dstar.Kchindof < 2.5) & (Dstar.pichindof < 2.5)]
        output['cutflow']['Dstar trk pt cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[(Dstar.KnValid > 4) & (Dstar.pinValid > 4) & (Dstar.KnPix > 1) & (Dstar.pinPix > 1)]
        output['cutflow']['Dstar trk hits cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[(Dstar.Kdxy < 0.1) & (Dstar.pidxy < 0.1)]
        output['cutflow']['Dstar trk pt cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[(Dstar.Kdz < 1) & (Dstar.pidz < 1)]
        output['cutflow']['Dstar trk pt cut'] += ak.sum(ak.num(Dstar))

        # pis cuts
        Dstar = Dstar[Dstar.pispt > 0.3]
        output['cutflow']['Dstar pis pt cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[Dstar.pischindof < 3]
        output['cutflow']['Dstar pis chi2 cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[Dstar.pisnValid > 2]
        output['cutflow']['Dstar pis hits cut'] += ak.sum(ak.num(Dstar))

        # D0 of Dstar cuts
        Dstar = Dstar[Dstar.D0cosphi > 0.99]
        output['cutflow']['Dstar D0 cosphi cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[(Dstar.D0mass < D0_PDG_MASS + 0.025) & (Dstar.D0mass > D0_PDG_MASS - 0.025)]
        output['cutflow']['Dstar D0 mass cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[Dstar.D0pt > 3]
        output['cutflow']['Dstar D0 pt cut'] += ak.sum(ak.num(Dstar))

        Dstar = Dstar[Dstar.D0dlSig > 3]
        output['cutflow']['Dstar D0 dlSig cut'] += ak.sum(ak.num(Dstar))

        Dstar['wrg_chg'] = (Dstar.Kchg == Dstar.pichg)

        ############### Dimu + OpenCharm associations

        DimuDstar = association(Dimu, Dstar)
        calc = ak.zip({
            'x': DimuDstar.slot0.x - DimuDstar.slot1.D0x,
            'y': DimuDstar.slot0.y - DimuDstar.slot1.D0y,
            'z': DimuDstar.slot0.z - DimuDstar.slot1.D0z,
            'Dimu_Covxx': DimuDstar.slot0.Covxx,
            'Dimu_Covyx': DimuDstar.slot0.Covyx,
            'Dimu_Covzx': DimuDstar.slot0.Covzx,
            'Dimu_Covyy': DimuDstar.slot0.Covyy,
            'Dimu_Covzy': DimuDstar.slot0.Covzy,
            'Dimu_Covzz': DimuDstar.slot0.Covzz,
            'Dstar_Covxx': D0[DimuDstar.slot1.D0recIdx].Covxx,
            'Dstar_Covyx': D0[DimuDstar.slot1.D0recIdx].Covyx,
            'Dstar_Covzx': D0[DimuDstar.slot1.D0recIdx].Covzx,
            'Dstar_Covyy': D0[DimuDstar.slot1.D0recIdx].Covyy,
            'Dstar_Covzy': D0[DimuDstar.slot1.D0recIdx].Covzy,
            'Dstar_Covzz': D0[DimuDstar.slot1.D0recIdx].Covzz,
            }, with_name='ThreeVector')

        err = np.sqrt(calc.x*calc.x*(calc.Dimu_Covxx + calc.Dstar_Covxx) + 2*calc.x*calc.y*(calc.Dimu_Covyx + calc.Dstar_Covyx) + 
                calc.y*calc.y*(calc.Dimu_Covyy + calc.Dstar_Covyy) + 2*calc.x*calc.z*(calc.Dimu_Covzx + calc.Dstar_Covzx) +
                calc.z*calc.z*(calc.Dimu_Covzz + calc.Dstar_Covzz) + 2*calc.y*calc.z*(calc.Dimu_Covzy + calc.Dstar_Covzy))

        DimuDstar['d'] = calc.p
        DimuDstar['dSig'] = calc.p/err

        ############### Cuts for D0
        D0 = D0[~D0.hasMuon]
        output['cutflow']['D0 trk muon cut'] += ak.sum(ak.num(D0))

        D0 = D0[(D0.t1pt > 0.8) & (D0.t2pt > 0.8)]
        output['cutflow']['D0 trk pt cut'] += ak.sum(ak.num(D0))

        D0 = D0[(D0.t1chindof < 2.5) & (D0.t2chindof < 2.5)]
        output['cutflow']['D0 trk chi2 cut'] += ak.sum(ak.num(D0))

        D0 = D0[(D0.t1nValid > 4) & (D0.t2nValid > 4) & (D0.t1nPix > 1) & (D0.t2nPix > 1)]
        output['cutflow']['D0 trk hits cut'] += ak.sum(ak.num(D0))

        D0 = D0[(D0.t1dxy < 0.1) & (D0.t2dxy < 0.1)]
        output['cutflow']['D0 trk dxy cut'] += ak.sum(ak.num(D0))

        D0 = D0[(D0.t1dz < 1.) & (D0.t2dz < 1.)]
        output['cutflow']['D0 trk dz cut'] += ak.sum(ak.num(D0))

        # D0 cosphi
        D0 = D0[D0.cosphi > 0.99]
        output['cutflow']['D0 cosphi cut'] += ak.sum(ak.num(D0))

        # D0 dl Significance
        D0 = D0[D0.dlSig > 5.]
        output['cutflow']['D0 dlSig cut'] += ak.sum(ak.num(D0))

        # D0 pt
        D0 = D0[D0.pt > 3.]
        output['cutflow']['D0 pt cut'] += ak.sum(ak.num(D0))

        ############### Final computation of number of objects
        output['cutflow']['Dimu final']    += ak.sum(ak.num(Dimu))
        output['cutflow']['D0 final']      += ak.sum(ak.num(D0))
        output['cutflow']['Dstar final']   += ak.sum(ak.num(Dstar))
        output['cutflow']['Dimu Dstar Associated'] += ak.sum(ak.num(DimuDstar))

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

        evt_info_acc = processor.dict_accumulator({})
        evt_info_acc['event'] = processor.column_accumulator(ak.to_numpy(events.event))
        evt_info_acc['run'] = processor.column_accumulator(ak.to_numpy(events.run))
        evt_info_acc['luminosityBlock'] = processor.column_accumulator(ak.to_numpy(events.luminosityBlock))
        output['event_info'] = evt_info_acc

        PVtx_acc = processor.dict_accumulator({})
        for var in PVtx.fields:
            PVtx_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(PVtx[var])))
        PVtx_acc['nPVtx'] = processor.column_accumulator(ak.to_numpy(ak.num(PVtx)))
        output['PVtx'] = PVtx_acc

        triggers_acc = processor.dict_accumulator({})
        for var in HLT.fields:
            triggers_acc[var] = processor.column_accumulator(ak.to_numpy(HLT[var]))
        output['triggers'] = triggers_acc

        file_hash = str(random.getrandbits(128)) + str(len(events))
        save(output, "output/" + self.analyzer_name + "/" + self.analyzer_name + "_" + file_hash + ".coffea")

        # return dummy accumulator
        return processor.dict_accumulator({
                'cutflow': output['cutflow']
        })

    def postprocess(self, accumulator):
        return accumulator
