from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np
from coffea.util import save
import random

from tools.collections import *

D0_PDG_MASS = 1.864

class EventSelectorProcessor(processor.ProcessorABC):
    def __init__(self, analyzer_name):
        self.analyzer_name = analyzer_name

        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        output = self.accumulator.identity()
        if df.size == 0:
            return processor.dict_accumulator({
                'foo': processor.defaultdict_accumulator(int),
                'cutflow': output['cutflow']
          })

        # Dimu candidates
        if df['nDimu'].size != 0:
            Dimu = JaggedCandidateArray.candidatesfromcounts(
                    df['nDimu'],
                    **get_vars_dict(df, dimu_cols)
                    )
        else:
            Dimu = JaggedCandidateArray.candidatesfromcounts(
                    np.array([]),
                    **get_vars_dict(df, dimu_cols)
                    )

        # Muon candidates
        if df['nMuon'].size != 0:
            Muon = JaggedCandidateArray.candidatesfromcounts(
                    df['nMuon'],
                    **get_vars_dict(df, muon_cols)
                    )
        else:
            Muon = JaggedCandidateArray.candidatesfromcounts(
                    np.array([]),
                    **get_vars_dict(df, muon_cols)
                    )

        # D0 candidates
        if df['nD0'].size != 0:
            D0 = JaggedCandidateArray.candidatesfromcounts(
                    df['nD0'],
                    mass=df['D0_mass12'],
                    **get_vars_dict(df, d0_cols)
                    )
        else:
            D0 = JaggedCandidateArray.candidatesfromcounts(
                    np.array([]),
                    mass=np.array([]),
                    **get_vars_dict(df, d0_cols)
                    )

        # Dstar candidates
        if df['nDstar'].size != 0:
            Dstar = JaggedCandidateArray.candidatesfromcounts(
                    df['nDstar'],
                    mass=df['Dstar_deltam'] + df['DstarD0_mass'],
                    **get_vars_dict(df, dstar_cols)
                    )
        else:
            Dstar = JaggedCandidateArray.candidatesfromcounts(
                    np.array([]),
                    mass=np.array([]),
                    **get_vars_dict(df, dstar_cols)
                    )

        output['cutflow']['all events']  += Muon.size
        output['cutflow']['all Dimu']    += Dimu.counts.sum()
        output['cutflow']['all D0']      += D0.counts.sum()
        output['cutflow']['all Dstar']   += Dstar.counts.sum()

        ############### Cuts
        # Dimu cuts: charge = 0, mass cuts and chi2...
        Dimu = Dimu[Dimu.charge == 0]
        output['cutflow']['Dimu 0 charge'] += Dimu.counts.sum()

        dimu_mass_cut = (Dimu.mass > 8.5) & (Dimu.mass < 11.5)
        Dimu = Dimu[dimu_mass_cut]
        output['cutflow']['Upsilon mass'] += Dimu.counts.sum()

        dimu_chi2_cut = (Dimu.chi2 < 5.)
        Dimu = Dimu[dimu_chi2_cut]
        output['cutflow']['chi2 cut'] = Dimu.counts.sum()

        ############### Get the Muons from Dimu, for cuts in their params
        if df['nDimu'].size != 0:
            mu1Idx = (Dimu.t1_muIdx + Muon.starts).content
            mu2Idx = (Dimu.t2_muIdx + Muon.starts).content
            Muon1 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[mu1Idx])
            Muon2 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[mu2Idx])
        else:
            Muon1 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[np.array([], dtype='int64')])
            Muon2 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[np.array([], dtype='int64')])

        # SoftId and Global Muon cuts
        soft_id = (Muon1.softId > 0) & (Muon2.softId > 0)
        Dimu = Dimu[soft_id]
        Muon1 = Muon1[soft_id]
        Muon2 = Muon2[soft_id]
        output['cutflow']['Dimu muon softId'] += Dimu.counts.sum()

        global_muon = (Muon1.isGlobal > 0) & (Muon2.isGlobal > 0)
        Dimu = Dimu[global_muon]
        Muon1 = Muon1[global_muon]
        Muon2 = Muon2[global_muon]
        output['cutflow']['Dimu muon global'] += Dimu.counts.sum()

        # pt and eta cuts
        muon_pt_cut = (Muon1.pt > 3) & (Muon2.pt > 3)
        Dimu = Dimu[muon_pt_cut]
        Muon1 = Muon1[muon_pt_cut]
        Muon2 = Muon2[muon_pt_cut]
        output['cutflow']['Dimu muon pt cut'] += Dimu.counts.sum()

        muon_eta_cut = (np.absolute(Muon1.eta) <= 2.4) & (np.absolute(Muon2.eta) <= 2.4)
        Dimu = Dimu[muon_eta_cut]
        Muon1 = Muon1[muon_eta_cut]
        Muon2 = Muon2[muon_eta_cut]
        output['cutflow']['Dimu muon eta cut'] += Dimu.counts.sum()

        ############### Cuts for D0

        # trk cuts
        D0_trk_muon_cut = ~D0.hasMuon
        D0 = D0[D0_trk_muon_cut]
        output['cutflow']['D0 trk muon cut'] += D0.counts.sum()

        D0_trk_pt_cut = (D0.t1_pt > 0.8) & (D0.t2_pt > 0.8)
        D0 = D0[D0_trk_pt_cut]
        output['cutflow']['D0 trk pt cut'] += D0.counts.sum()

        D0_trk_chi2_cut = (D0.t1_chindof < 2.5) & (D0.t2_chindof < 2.5)
        D0 = D0[D0_trk_chi2_cut]
        output['cutflow']['D0 trk chi2 cut'] += D0.counts.sum()

        D0_trk_hits_cut = (D0.t1_nValid > 4) & (D0.t2_nValid > 4) & (D0.t1_nPix > 1) & (D0.t2_nPix > 1)
        D0 = D0[D0_trk_hits_cut]
        output['cutflow']['D0 trk hits cut'] += D0.counts.sum()

        D0_trk_dxy_cut = (D0.t1_dxy < 0.1) & (D0.t2_dxy < 0.1)
        D0 = D0[D0_trk_dxy_cut]
        output['cutflow']['D0 trk dxy cut'] += D0.counts.sum()

        D0_trk_dz_cut = (D0.t1_dz < 1.) & (D0.t2_dz < 1.)
        D0 = D0[D0_trk_dz_cut]
        output['cutflow']['D0 trk dz cut'] += D0.counts.sum()

        # D0 cosphi
        D0_cosphi_cut = (D0.cosphi > 0.99)
        D0 = D0[D0_cosphi_cut]
        output['cutflow']['D0 cosphi cut'] += D0.counts.sum()

        # D0 dl Significance
        D0_dlSig_cut = (D0.dlSig > 5.)
        D0 = D0[D0_dlSig_cut]
        output['cutflow']['D0 dlSig cut'] += D0.counts.sum()

        # D0 pt
        D0_pt_cut = (D0.pt > 3.)
        D0 = D0[D0_pt_cut]
        output['cutflow']['D0 pt cut'] += D0.counts.sum()

        ############### Cuts for Dstar

        # trks cuts
        Dstar_trk_muon_cut = ~Dstar.hasMuon
        Dstar = Dstar[Dstar_trk_muon_cut]
        output['cutflow']['Dstar trk muon cut'] += Dstar.counts.sum()

        Dstar_trk_pt_cut = (Dstar.K_pt > 0.5) & (Dstar.pi_pt > 0.5)
        Dstar = Dstar[Dstar_trk_pt_cut]
        output['cutflow']['Dstar trk pt cut'] += Dstar.counts.sum()

        Dstar_trk_chi2_cut = (Dstar.K_chindof < 2.5) & (Dstar.pi_chindof < 2.5)
        Dstar = Dstar[Dstar_trk_chi2_cut]
        output['cutflow']['Dstar trk pt cut'] += Dstar.counts.sum()

        Dstar_trk_hits_cut = (Dstar.K_nValid > 4) & (Dstar.pi_nValid > 4) & (Dstar.K_nValid > 1) & (Dstar.pi_nValid > 1)
        Dstar = Dstar[Dstar_trk_hits_cut]
        output['cutflow']['Dstar trk hits cut'] += Dstar.counts.sum()

        Dstar_trk_dxy_cut = (Dstar.K_dxy < 0.1) & (Dstar.pi_dxy < 0.1)
        Dstar = Dstar[Dstar_trk_dxy_cut]
        output['cutflow']['Dstar trk pt cut'] += Dstar.counts.sum()

        Dstar_trk_dz_cut = (Dstar.K_dz < 1) & (Dstar.pi_dz < 1)
        Dstar = Dstar[Dstar_trk_dz_cut]
        output['cutflow']['Dstar trk pt cut'] += Dstar.counts.sum()

        # pis cuts
        Dstar_pis_pt_cut = (Dstar.pis_pt > 0.3)
        Dstar = Dstar[Dstar_pis_pt_cut]
        output['cutflow']['Dstar pis pt cut'] += Dstar.counts.sum()

        Dstar_pis_chi2_cut = (Dstar.pis_chindof < 3)
        Dstar = Dstar[Dstar_pis_chi2_cut]
        output['cutflow']['Dstar pis chi2 cut'] += Dstar.counts.sum()

        Dstar_pis_hits_cut = (Dstar.pis_nValid > 2)
        Dstar = Dstar[Dstar_pis_hits_cut]
        output['cutflow']['Dstar pis hits cut'] += Dstar.counts.sum()

        # D0 of Dstar cuts
        DstarD0_cosphi_cut = (Dstar.D0_cosphi > 0.99)
        Dstar = Dstar[DstarD0_cosphi_cut]
        output['cutflow']['Dstar D0 cosphi cut'] += Dstar.counts.sum()

        DstarD0_mass_cut = (Dstar.D0_mass < D0_PDG_MASS + 0.025) & (Dstar.D0_mass > D0_PDG_MASS - 0.025)
        Dstar = Dstar[DstarD0_mass_cut]
        output['cutflow']['Dstar D0 mass cut'] += Dstar.counts.sum()

        DstarD0_pt_cut = (Dstar.D0_pt > 3)
        Dstar = Dstar[DstarD0_pt_cut]
        output['cutflow']['Dstar D0 pt cut'] += Dstar.counts.sum()

        DstarD0_dlSig_cut = (Dstar.D0_dlSig > 3)
        Dstar = Dstar[DstarD0_dlSig_cut]
        output['cutflow']['Dstar D0 dlSig cut'] += Dstar.counts.sum()

        """ Dstar_wrong_charge_cut = (Dstar.K_chg != Dstar.pi_chg)
        Dstar = Dstar[Dstar_wrong_charge_cut]
        output['cutflow']['Dstar wrong charge cut'] += Dstar.counts.sum() """

        ############### Upsilon + Dstar association

        ############### Final computation of number of objects
        output['cutflow']['Dimu final']    += Dimu.counts.sum()
        output['cutflow']['D0 final']      += D0.counts.sum()
        output['cutflow']['Dstar final']   += Dstar.counts.sum()

        ############### Leading and Trailing muon separation
        leading_mu = (Muon1.pt.content > Muon2.pt.content)
        Muon_lead = JaggedCandidateArray.candidatesfromoffsets(Dimu.offsets,
                                                               pt=np.where(leading_mu, Muon1.pt.content, Muon2.pt.content),
                                                               eta=np.where(leading_mu, Muon1.eta.content, Muon2.eta.content),
                                                               phi=np.where(leading_mu, Muon1.phi.content, Muon2.phi.content),
                                                               mass=np.where(leading_mu, Muon1.mass.content, Muon2.mass.content),)

        Muon_trail = JaggedCandidateArray.candidatesfromoffsets(Dimu.offsets,
                                                                pt=np.where(~leading_mu, Muon1.pt.content, Muon2.pt.content),
                                                                eta=np.where(~leading_mu, Muon1.eta.content, Muon2.eta.content),
                                                                phi=np.where(~leading_mu, Muon1.phi.content, Muon2.phi.content),
                                                                mass=np.where(~leading_mu, Muon1.mass.content, Muon2.mass.content),)



        ############### Create the accumulators to save output
        muon_lead_acc = processor.dict_accumulator({})
        for var in Muon_lead.columns:
            if var == 'p4': continue
            muon_lead_acc[var] = processor.column_accumulator(np.array(Muon_lead[var].flatten()))
        muon_lead_acc["nMuon"] = processor.column_accumulator(Muon_lead.counts)
        output["Muon_lead"] = muon_lead_acc

        muon_trail_acc = processor.dict_accumulator({})
        for var in Muon_trail.columns:
            if var == 'p4': continue
            muon_trail_acc[var] = processor.column_accumulator(np.array(Muon_trail[var].flatten()))
        muon_trail_acc["nMuon"] = processor.column_accumulator(Muon_trail.counts)
        output["Muon_trail"] = muon_trail_acc

        dimu_acc = processor.dict_accumulator({})
        for var in Dimu.columns:
            if (var == 'p4' or var.startswith('t')): continue
            dimu_acc[var] = processor.column_accumulator(np.array(Dimu[var].flatten()))
        dimu_acc["nDimu"] = processor.column_accumulator(Dimu.counts)
        output["Dimu"] = dimu_acc

        D0_acc = processor.dict_accumulator({})
        D0_trk_acc = processor.dict_accumulator({})
        for var in D0.columns:
            if (var == 'p4'): continue
            elif ( var.startswith('t')):
                D0_trk_acc[var] = processor.column_accumulator(np.array(D0[var].flatten()))
            else:
                D0_acc[var] = processor.column_accumulator(np.array(D0[var].flatten()))
        D0_acc["nD0"] = processor.column_accumulator(D0.counts)
        output["D0"] = D0_acc
        output["D0_trk"] = D0_trk_acc

        Dstar_acc = processor.dict_accumulator({})
        Dstar_D0_acc = processor.dict_accumulator({})
        Dstar_trk_acc = processor.dict_accumulator({})
        for var in Dstar.columns:
            if (var == 'p4' ): continue
            elif var.startswith('D0'):
                Dstar_D0_acc[var] = processor.column_accumulator(np.array(Dstar[var].flatten()))
            elif (var.startswith('K') or var.startswith('pi')):
                Dstar_trk_acc[var] = processor.column_accumulator(np.array(Dstar[var].flatten()))
            else:
                Dstar_acc[var] = processor.column_accumulator(np.array(Dstar[var].flatten()))
        Dstar_acc["nDstar"] = processor.column_accumulator(Dstar.counts)
        output["Dstar"] = Dstar_acc
        output["Dstar_D0"] = Dstar_D0_acc
        output["Dstar_trk"] = Dstar_trk_acc

        file_hash = str(random.getrandbits(128)) + str(df.size)
        save(output, "output/" + self.analyzer_name + "/" + self.analyzer_name + "_" + file_hash + ".coffea")

        # return dummy accumulator
        return processor.dict_accumulator({
                'foo': processor.defaultdict_accumulator(int),
                'cutflow': output['cutflow']
          })


    def postprocess(self, accumulator):
        return accumulator
