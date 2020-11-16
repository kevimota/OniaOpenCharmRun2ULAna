from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np
from coffea.util import save
import random

class EventSelectorProcessor(processor.ProcessorABC):
   def __init__(self, analyzer_name):
      self.analyzer_name = analyzer_name
      
      self._accumulator = processor.dict_accumulator({
         'cutflow': processor.defaultdict_accumulator(int),
         'dataset': processor.value_accumulator(str)
      })
    
   @property
   def accumulator(self):
      return self._accumulator
    
   def process(self, df):
      output = self.accumulator.identity()

      output['dataset'] = df['dataset']

      # Dimu candidates
      if df['nDimu'].size != 0:
         Dimu = JaggedCandidateArray.candidatesfromcounts(
               df['nDimu'],
               pt=df['Dimu_pt'],
               eta=df['Dimu_eta'],
               phi=df['Dimu_phi'],
               mass=df['Dimu_mass'],
               charge=df['Dimu_charge'],
               t1Idx=df['Dimut1_muIdx'],
               t2Idx=df['Dimut2_muIdx'],
               )
      else:
         Dimu = JaggedCandidateArray.candidatesfromcounts(
               np.array([]),
               pt=np.array([]),
               eta=np.array([]),
               phi=np.array([]),
               mass=np.array([]),
               charge=np.array([]),
               t1Idx=np.array([]),
               t2Idx=np.array([]),
               )


      # Muon candidates
      if df['nMuon'].size != 0:
         Muon = JaggedCandidateArray.candidatesfromcounts(
               df['nMuon'],
               pt=df['Muon_pt'],
               eta=df['Muon_eta'],
               phi=df['Muon_phi'],
               mass=df['Muon_mass'],
               charge=df['Muon_charge'],
               isGlobal=df['Muon_isGlobal'],
               softId=df['Muon_softId'],
               )
      else:  
         Muon = JaggedCandidateArray.candidatesfromcounts(
               np.array([]),
               pt=np.array([]),
               eta=np.array([]),
               phi=np.array([]),
               mass=np.array([]),
               charge=np.array([]),
               isGlobal=np.array([]),
               softId=np.array([]),
               )        
      
      # D0 candidates
      if df['nD0'].size != 0:
         D0 = JaggedCandidateArray.candidatesfromcounts(
               df['nD0'],
               pt=df['D0_pt'],
               eta=df['D0_eta'],
               phi=df['D0_phi'],
               mass=df['D0_mass12'],
               mass12=df['D0_mass12'],
               mass21=df['D0_mass21'],
               vtxIdx=df['D0_vtxIdx'],
               cosphi=df['D0_cosphi'],
               dlSig=df['D0_dlSig'],
               t1_pt=df['D0t1_pt'],
               t1_chindof=df['D0t1_chindof'],
               t1_nValid=df['D0t1_nValid'],
               t1_nPix=df['D0t1_nPix'],
               t1_dz=df['D0t1_dz'],
               t1_dxy=df['D0t1_dxy'],
               t2_pt=df['D0t2_pt'],
               t2_chindof=df['D0t2_chindof'],
               t2_nValid=df['D0t2_nValid'],
               t2_nPix=df['D0t2_nPix'],
               t2_dz=df['D0t2_dz'],
               t2_dxy=df['D0t2_dxy'],
               )
      else:  
         D0 = JaggedCandidateArray.candidatesfromcounts(
               np.array([]),
               pt=np.array([]),
               eta=np.array([]),
               phi=np.array([]),
               mass=np.array([]),
               mass12=np.array([]),
               mass21=np.array([]),
               vtxIdx=np.array([]),
               cosphi=np.array([]),
               dlSig=np.array([]),
               t1_chindof=np.array([]),
               t1_nValid=np.array([]),
               t1_nPix=np.array([]),
               t1_dz=np.array([]),
               t1_dxy=np.array([]),
               t2_chindof=np.array([]),
               t2_nValid=np.array([]),
               t2_nPix=np.array([]),
               t2_dz=np.array([]),
               t2_dxy=np.array([]),
               )

      # Dstar candidates
      if df['nDstar'].size != 0:
         Dstar = JaggedCandidateArray.candidatesfromcounts(
               df['nDstar'],
               pt=df['Dstar_pt'],
               eta=df['Dstar_eta'],
               phi=df['Dstar_phi'],
               mass=df['Dstar_deltam'] + df['DstarD0_mass'],
               vtxIdx=df['Dstar_vtxIdx'],
               K_pt=df['DstarK_pt'],
               K_chindof=df['DstarK_chindof'],
               K_nValid=df['DstarK_nValid'],
               K_nPix=df['DstarK_nPix'],
               K_dxy=df['DstarK_dxy'],
               K_dz=df['DstarK_dz'],
               pi_pt=df['Dstarpi_pt'],
               pi_chindof=df['Dstarpi_chindof'],
               pi_nValid=df['Dstarpi_nValid'],
               pi_nPix=df['Dstarpi_nPix'],
               pi_dxy=df['Dstarpi_dxy'],
               pi_dz=df['Dstarpi_dz'],
               pis_pt=df['Dstarpis_pt'],
               pis_chindof=df['Dstarpis_chindof'],
               pis_nValid=df['Dstarpis_nValid'],
               pis_dxy=df['Dstarpis_dxy'],
               pis_dz=df['Dstarpis_dz'],
               D0_cosphi=df['DstarD0_cosphi'],
               D0_pt=df['DstarD0_pt'],
               D0_dlSig=df['DstarD0_dlSig'],
               )
      else:  
         Dstar = JaggedCandidateArray.candidatesfromcounts(
               np.array([]),
               pt=np.array([]),
               eta=np.array([]),
               phi=np.array([]),
               mass=np.array([]),
               vtxIdx=np.array([]),
               K_pt=np.array([]),
               K_chindof=np.array([]),
               K_nValid=np.array([]),
               K_nPix=np.array([]),
               K_dxy=np.array([]),
               K_dz=np.array([]),
               pi_pt=np.array([]),
               pi_chindof=np.array([]),
               pi_nValid=np.array([]),
               pi_nPix=np.array([]),
               pi_dxy=np.array([]),
               pi_dz=np.array([]),
               pis_pt=np.array([]),
               pis_chindof=np.array([]),
               pis_nValid=np.array([]),
               pis_dxy=np.array([]),
               pis_dz=np.array([]),
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

      ############### Get the Muons from Dimu, for cuts in their params

      # For Muon 1:
      mu1Idx = (Dimu.t1Idx + Muon.starts).content
      mu2Idx = (Dimu.t2Idx + Muon.starts).content
      Muon1 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[mu1Idx])
      Muon2 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[mu2Idx])

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

      # Cuts into events with at least 1 Dimu
      ndimu_cut = Dimu.counts > 0
      Dimu = Dimu[ndimu_cut]
      Muon1 = Muon1[ndimu_cut]
      Muon2 = Muon2[ndimu_cut]
      D0 = D0[ndimu_cut]
      Dstar = Dstar[ndimu_cut]

      ############### Cuts for D0

      # trk cuts
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

      D0_trk_dz_cut = (D0.t1_dz < 1) & (D0.t2_dz < 1)
      D0 = D0[D0_trk_dz_cut]
      output['cutflow']['D0 trk dz cut'] += D0.counts.sum()

      # D0 cosphi
      D0_cosphi_cut = (D0.cosphi > 0.99)
      D0 = D0[D0_cosphi_cut]
      output['cutflow']['D0 cosphi cut'] += D0.counts.sum()

      # D0 dl Significance
      D0_dlSig_cut = (D0.dlSig > 5)
      D0 = D0[D0_dlSig_cut]
      output['cutflow']['D0 dlSig cut'] += D0.counts.sum()

      # D0 pt
      D0_pt_cut = (D0.pt > 3)
      D0 = D0[D0_pt_cut]
      output['cutflow']['D0 pt cut'] += D0.counts.sum()

      ############### Cuts for Dstar

      # trks cuts
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

      DstarD0_pt_cut = (Dstar.D0_pt > 3)
      Dstar = Dstar[DstarD0_pt_cut]
      output['cutflow']['Dstar D0 pt cut'] += Dstar.counts.sum()

      DstarD0_dlSig_cut = (Dstar.D0_dlSig > 3)
      Dstar = Dstar[DstarD0_dlSig_cut]
      output['cutflow']['Dstar D0 dlSig cut'] += Dstar.counts.sum() 

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
      columns = ['__fast_pt', '__fast_phi', '__fast_eta', '__fast_mass']
      muon_lead_acc = processor.dict_accumulator({})
      for var in columns:
         muon_lead_acc[var] = processor.column_accumulator(np.array(Muon_lead[var].flatten()))
      muon_lead_acc["nMuon"] = processor.column_accumulator(Muon_lead.counts)
      output["Muon_lead"] = muon_lead_acc

      muon_trail_acc = processor.dict_accumulator({})
      for var in columns:
         muon_trail_acc[var] = processor.column_accumulator(np.array(Muon_trail[var].flatten()))
      muon_trail_acc["nMuon"] = processor.column_accumulator(Muon_trail.counts)
      output["Muon_trail"] = muon_trail_acc
      
      dimu_acc = processor.dict_accumulator({})
      for var in columns:
         dimu_acc[var] = processor.column_accumulator(np.array(Dimu[var].flatten()))
      dimu_acc["nDimu"] = processor.column_accumulator(Dimu.counts)
      output["Dimu"] = dimu_acc
      
      D0_acc = processor.dict_accumulator({})
      for var in columns:
         D0_acc[var] = processor.column_accumulator(np.array(D0[var].flatten()))
      D0_acc["nD0"] = processor.column_accumulator(D0.counts)
      output["D0"] = D0_acc

      Dstar_acc = processor.dict_accumulator({})
      for var in columns:
         Dstar_acc[var] = processor.column_accumulator(np.array(Dstar[var].flatten()))
      Dstar_acc["nDstar"] = processor.column_accumulator(Dstar.counts)
      output["Dstar"] = Dstar_acc

      file_hash = str(random.getrandbits(128)) + str(df.size)
      save(output, "output/" + self.analyzer_name + "/" + self.analyzer_name + "_" + file_hash + ".coffea")

      # return dummy accumulator
      return processor.dict_accumulator({
            'foo': processor.defaultdict_accumulator(int),
        })


   def postprocess(self, accumulator):
      return accumulator
    