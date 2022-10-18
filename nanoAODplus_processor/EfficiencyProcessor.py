import awkward as ak
import numpy as np
import coffea.processor as processor

from hist import Hist

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

from tools.collections import *
from tools.utils import save_kin_hists, D0_PDG_MASS, remove_none

def n_events_bin(bins_x, bins_y, cand_x, cand_y, cand, acc):
    for i in range(len(bins_x)):
        if i == 0: continue
        for j in range(len(bins_y)):
            if j == 0: continue
            bin_x_high = bins_x[i]
            bin_x_low = bins_x[i-1]
            bin_rap_high = bins_y[j]
            bin_rap_low = bins_y[j-1]
            x_cut = (cand_x >= bin_x_low) & (cand_x < bin_x_high)
            y_cut = (np.absolute(cand_y) >= bin_rap_low) & (np.absolute(cand_y) < bin_rap_high)
            cand_cut = cand[x_cut & y_cut]
            acc[f'{bin_x_low};{bin_x_high}:{bin_rap_low};{bin_rap_high}'] += len(cand_cut[ak.num(cand_cut) > 0])

class EfficiencyProcessor(processor.ProcessorABC):
    def __init__(self, config, year, dimu_cut=0):
        self.config = config
        self.year = year
        self.dimu_cut = dimu_cut

        Gen_Dimu_hists = {
            'pt': Hist.new.Regular(100, 0, 300, name='pt', label=r"$p_{T,Gen \mu\mu}$ [GeV]").Double(),
            'eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{Gen \mu\mu}$").Double(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{Gen \mu\mu}$").Double(),
        }

        Reco_Dimu_hists = {
            'pt': Hist.new.Regular(100, 0, 300, name='pt', label=r"$p_{T,Reco \mu\mu}$ [GeV]").Double(),
            'eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{Reco \mu\mu}$").Double(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{Reco \mu\mu}$").Double(),
            'mass': Hist.new.Regular(100, 8.6, 11, name="mass", label=r"$m_{Reco \mu\mu}$ [GeV]").Double(),
        }

        Gen_Dstar_hists = {
            'pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,Gen D^*}$ [GeV]").Double(),
            'eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{Gen D^*}$").Double(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{Gen D^*}$").Double(),
        }

        Reco_Dstar_hists = {
            'pt': Hist.new.Regular(100, 0, 100, name='pt', label=r"$p_{T,Reco D^*}$ [GeV]").Double(),
            'eta': Hist.new.Regular(60, -4, 4, name="eta", label=r"$\eta_{Reco D^*}$").Double(),
            'phi': Hist.new.Regular(60, -np.pi, np.pi, name="phi", label=r"$\phi_{Reco D^*}$").Double(),
            'deltam': Hist.new.Regular(50, 0.138, 0.162, name="deltam", label=r"$\Delta m_{Reco D^*}$ [GeV]").Double(),
        }

        self._accumulator = {
            'N_evt': processor.value_accumulator(int),
            'N_gen_dstar': processor.defaultdict_accumulator(int),
            'N_reco_dstar': processor.defaultdict_accumulator(int),
            'N_gen_dimu': processor.defaultdict_accumulator(int),
            'N_reco_dimu': processor.defaultdict_accumulator(int),
            'N_cuts_dstar': processor.defaultdict_accumulator(int),
            'N_cuts_dimu': processor.defaultdict_accumulator(int),
            'N_trigger_dimu': processor.defaultdict_accumulator(int),
            'N_num_asso': processor.defaultdict_accumulator(int),
            'N_den_asso': processor.defaultdict_accumulator(int),
            'Gen_Dimu': Gen_Dimu_hists,
            'Reco_Dimu': Reco_Dimu_hists,
            'Gen_Dstar': Gen_Dstar_hists,
            'Reco_Dstar': Reco_Dstar_hists,
        }

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator
        
        # test if there are any events in the file
        if len(events) == 0:
            return output
        
        output['N_evt'] += len(events)

        ############### Get All the interesting candidates from NTuples
        Dimu = ak.zip({**get_vars_dict(events, dimu_cols)}, with_name="PtEtaPhiMCandidate")
        Muon = ak.zip({**get_vars_dict(events, muon_cols)}, with_name="PtEtaPhiMCandidate")
        D0 = ak.zip({'mass': events.D0_mass12, **get_vars_dict(events, d0_cols)}, with_name="PtEtaPhiMCandidate")
        Dstar = ak.zip({'mass': (events.Dstar_D0mass + events.Dstar_deltamr),
                        'charge': events.Dstar_pischg,
                        **get_vars_dict(events, dstar_cols)}, 
                        with_name="PtEtaPhiMCandidate")
        GenPart = ak.zip({**get_vars_dict(events, gen_part_cols)}, with_name="PtEtaPhiMCandidate")
        HLT = ak.zip({**get_hlt(events, [self.config['trigger'][self.year]])})

        #GenParticles for Dimu and for Dstar
        GenPart_Muon = GenPart[(np.absolute(GenPart.pdgId) == 13) & (GenPart.genPartIdxMother > -1)]

        ############### Acceptance

        if self.config['particle'] == 'J/psi':
            gen_pdgid = (GenPart_Muon.parpdgId == 443)
        elif self.config['particle'] == 'Y':
            gen_pdgid = (GenPart_Muon.parpdgId == 553) | (GenPart_Muon.parpdgId == 100553) | (GenPart_Muon.parpdgId == 200553)
        
        GenPart_Muon = GenPart_Muon[gen_pdgid]
        GenPart_Muon = GenPart_Muon[GenPart_Muon.pt > self.config['muon_pt_min']]
        GenPart_Muon = GenPart_Muon[np.absolute(GenPart_Muon.eta) < self.config['muon_eta']]
        GenPart_Muon = ak.combinations(GenPart_Muon, 2)
        leading_mu = (GenPart_Muon.slot0.pt > GenPart_Muon.slot1.pt)
        GenPart_Muon = ak.zip({'0': ak.where(leading_mu, GenPart_Muon.slot0, GenPart_Muon.slot1), '1': ak.where(~leading_mu, GenPart_Muon.slot0, GenPart_Muon.slot1)})
        GenPart_Muon = GenPart_Muon[GenPart_Muon.slot0.genPartIdxMother == GenPart_Muon.slot1.genPartIdxMother]
        
        GenPart_Dimu = GenPart[GenPart_Muon.slot0.genPartIdxMother]
        GenPart_Dimu = GenPart_Dimu[(GenPart_Dimu.pt > max(self.config['dimu_pt_min'], self.dimu_cut)) & (GenPart_Dimu.pt < self.config['dimu_pt_max'])]
        GenPart_Dimu['rap'] = np.log((GenPart_Dimu.t+GenPart_Dimu.z)/(GenPart_Dimu.t-GenPart_Dimu.z))/2       
        GenPart_Dimu = GenPart_Dimu[np.absolute(GenPart_Dimu.rap) < self.config['dimu_rap']]        

        GenPart_D0 = GenPart[(np.absolute(GenPart.pdgId) == 421) & (np.absolute(GenPart.parpdgId) == 413)]
        GenPart_Dstar = GenPart[GenPart_D0.genPartIdxMother]
        GenPart_Dstar = GenPart_Dstar[(GenPart_Dstar.pt > self.config['dstar_pt_min']) & (GenPart_Dstar.pt < self.config['dstar_pt_max'])]
        GenPart_Dstar['rap'] = np.log((GenPart_Dstar.t+GenPart_Dstar.z)/(GenPart_Dstar.t-GenPart_Dstar.z))/2
        GenPart_Dstar = GenPart_Dstar[np.absolute(GenPart_Dstar.rap) < self.config['dstar_rap']]

        Dimu = ak.mask(Dimu, Dimu.charge == 0)
        Dimu = ak.mask(Dimu, (Dimu.mass > self.config['dimu_mass_low']) & (Dimu.mass < self.config['dimu_mass_high']))
        Muon = ak.zip({'0': Muon[Dimu.t1muIdx], '1': Muon[Dimu.t2muIdx]})
        leading_mu = (Muon.slot0.pt > Muon.slot1.pt)
        Muon = ak.zip({'0': ak.where(leading_mu, Muon.slot0, Muon.slot1), '1': ak.where(~leading_mu, Muon.slot0, Muon.slot1)})

        muon_pt_cut = (Muon.slot0.pt > self.config['muon_pt_min']) & (Muon.slot1.pt > self.config['muon_pt_min'])
        muon_eta_cut = (np.absolute(Muon.slot0.eta) < self.config['muon_eta']) & (np.absolute(Muon.slot0.eta) < self.config['muon_eta'])
        muon_sim_cut = (Muon.slot0.simIdx > -1) & (Muon.slot1.simIdx > -1)
        dimu_pt_cut = (Dimu.pt > max(self.config['dimu_pt_min'], self.dimu_cut)) & (Dimu.pt < self.config['dimu_pt_max'])
        dimu_eta_cut = np.absolute(Dimu.rap) < self.config['dimu_rap']
        Muon = ak.mask(Muon, muon_eta_cut & muon_pt_cut & muon_sim_cut & dimu_pt_cut & dimu_eta_cut)
        Dimu = ak.mask(Dimu, muon_eta_cut & muon_pt_cut & muon_sim_cut & dimu_pt_cut & dimu_eta_cut)
        #muon_sim_cut_2 = (GenPart[Muon.slot0.simIdx] == GenPart_Muon.slot0) & (GenPart[Muon.slot1.simIdx] == GenPart_Muon.slot1)
        """ Muon0 = ak.cartesian([GenPart[Muon.slot0.simIdx], GenPart_Muon.slot0])
        Muon1 = ak.cartesian([GenPart[Muon.slot1.simIdx], GenPart_Muon.slot1])
        muon_sim_cut = ak.any(Muon0.slot0 == Muon0.slot1, axis=1) & ak.any(Muon1.slot0 == Muon1.slot1, axis=1)
        Muon = ak.mask(Muon, muon_sim_cut)
        Dimu = ak.mask(Dimu, muon_sim_cut)
        print(Dimu.pt) """
        
        #Muon = ak.mask(Muon, muon_sim_cut_2)
        #Dimu = ak.mask(Dimu, muon_sim_cut_2)

        Dstar = Dstar[~Dstar.hasMuon]
        Dstar = Dstar[Dstar.Kchg != Dstar.pichg]
        Dstar = Dstar[(Dstar.pt > self.config['dstar_pt_min']) & (Dstar.pt < self.config['dstar_pt_max'])]
        Dstar = Dstar[np.absolute(Dstar.rap) < self.config['dstar_rap']]
        Dstar_sim = Dstar[Dstar.simIdx > -1]

        Dimu_aux = Dimu[remove_none(Dimu.pt)]

        output['N_gen_dimu']['global'] += len(GenPart_Dimu[ak.num(GenPart_Dimu) > 0])
        output['N_reco_dimu']['global'] += len(Dimu_aux[ak.num(Dimu_aux) > 0])
        output['N_gen_dstar']['global'] += len(GenPart_Dstar[ak.num(GenPart_Dstar) > 0])
        output['N_reco_dstar']['global'] += len(Dstar_sim[ak.num(Dstar_sim) > 0])

        n_events_bin(self.config['bins_pt_dimu'], self.config['bins_rap_dimu'], GenPart_Dimu.pt, GenPart_Dimu.rap, GenPart_Dimu, output['N_gen_dimu'])
        n_events_bin(self.config['bins_pt_dimu'], self.config['bins_rap_dimu'], Dimu_aux.pt, Dimu_aux.rap, Dimu_aux, output['N_reco_dimu'])
        n_events_bin(self.config['bins_pt_dstar'], self.config['bins_rap_dstar'], GenPart_Dstar.pt, GenPart_Dstar.rap, GenPart_Dstar, output['N_gen_dstar'])
        n_events_bin(self.config['bins_pt_dstar'], self.config['bins_rap_dstar'], Dstar_sim.pt, Dstar_sim.rap, Dstar_sim, output['N_reco_dstar'])

        ############### /Acceptance

        ############### Cuts Efficiency
        soft_id = (Muon.slot0.softId > 0) & (Muon.slot1.softId > 0)
        Dimu = ak.mask(Dimu, soft_id)
        Muon = ak.mask(Muon, soft_id)

        Dstar = Dstar[(Dstar.Kpt > self.config['dstar_track_pt_cut']) & (Dstar.pipt > self.config['dstar_track_pt_cut'])]
        Dstar = Dstar[(Dstar.Kchindof < 2.5) & (Dstar.pichindof < 2.5)]
        Dstar = Dstar[(Dstar.KnValid > 4) & (Dstar.pinValid > 4) & (Dstar.KnPix > 1) & (Dstar.pinPix > 1)]
        Dstar = Dstar[(Dstar.Kdxy < 0.1) & (Dstar.pidxy < 0.1)]
        Dstar = Dstar[(Dstar.Kdz < 1) & (Dstar.pidz < 1)]

        # pis cuts
        Dstar = Dstar[Dstar.pispt > 0.3]
        Dstar = Dstar[Dstar.pischindof < 3]
        Dstar = Dstar[Dstar.pisnValid > 2]

        # D0 of Dstar cuts
        Dstar = Dstar[Dstar.D0cosphi > self.config['dstar_d0_cosphi']]
        Dstar = Dstar[(Dstar.D0mass < D0_PDG_MASS + self.config['dstar_d0_mass']) & (Dstar.D0mass > D0_PDG_MASS - self.config['dstar_d0_mass'])]
        Dstar = Dstar[Dstar.D0pt > self.config['dstar_d0_pt']]
        Dstar = Dstar[Dstar.D0dlSig > self.config['dstar_d0_dlSig']]

        Dimu_aux = Dimu[remove_none(Dimu.pt)]

        output['N_cuts_dimu']['global'] += len(Dimu_aux[ak.num(Dimu_aux) > 0])
        output['N_cuts_dstar']['global'] += len(Dstar[ak.num(Dstar) > 0])

        n_events_bin(self.config['bins_pt_dimu'], self.config['bins_rap_dimu'], Dimu_aux.pt, Dimu_aux.rap, Dimu_aux, output['N_cuts_dimu'])
        n_events_bin(self.config['bins_pt_dstar'], self.config['bins_rap_dstar'], Dstar.pt, Dstar.rap, Dstar, output['N_cuts_dstar'])

        ############### /Cuts Efficiency 

        ############### Trigger Efficiency

        trigger = HLT[self.config['trigger'][self.year]]
        Dimu = Dimu[trigger]
        Muon = Muon[trigger]
        Dstar = Dstar[trigger]

        Dimu_aux = Dimu[remove_none(Dimu.pt)]

        output['N_trigger_dimu']['global'] += len(Dimu_aux[ak.num(Dimu_aux) > 0])

        n_events_bin(self.config['bins_pt_dimu'], self.config['bins_rap_dimu'], Dimu_aux.pt, Dimu_aux.rap, Dimu_aux, output['N_trigger_dimu'])

        ############### /Trigger Efficiency 

        ############### Efficiency Association

        Dstar = Dstar[Dstar.associationIdx >= 0]
        DimuDstar = ak.zip({'0': Dimu[Dstar.associationIdx], '1': Dstar})
        DimuDstar = DimuDstar[remove_none(DimuDstar.slot0.pt)]
        
        output['N_den_asso']['global'] += len(DimuDstar[ak.num(DimuDstar) > 0])
        n_events_bin(self.config['bins_pt_dimu'], self.config['bins_pt_dstar'], DimuDstar.slot0.pt, DimuDstar.slot1.pt, DimuDstar, output['N_den_asso'])

        DimuDstar = DimuDstar[DimuDstar.slot1.associationProb > self.config['vertex_probability_cut']]
        output['N_num_asso']['global'] += len(DimuDstar[ak.num(DimuDstar) > 0])
        n_events_bin(self.config['bins_pt_dimu'], self.config['bins_pt_dstar'], DimuDstar.slot0.pt, DimuDstar.slot1.pt, DimuDstar, output['N_num_asso'])

        ############### /Efficiency Association
    
        #save the hists:
        save_kin_hists(output['Gen_Dimu'], GenPart_Dimu, gen=True)
        save_kin_hists(output['Reco_Dimu'], Dimu_aux)
        save_kin_hists(output['Gen_Dstar'], GenPart_Dstar, gen=True, get_deltam=True)
        save_kin_hists(output['Reco_Dstar'], Dstar, get_deltam=True)

        return output

    def postprocess(self, accumulator):
        return accumulator
