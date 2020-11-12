from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray, Table
import numpy as np
from coffea.util import save, load
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
               vtxIdx=df['D0_vtxIdx'],
               x=df['D0_x'],
               y=df['D0_y'],
               z=df['D0_z'],
               )
      else:  
         D0 = JaggedCandidateArray.candidatesfromcounts(
               np.array([]),
               pt=np.array([]),
               eta=np.array([]),
               phi=np.array([]),
               mass=np.array([]),
               vtxIdx=np.array([]),
               x=np.array([]),
               y=np.array([]),
               z=np.array([]),
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
               D0x=df['DstarD0_x'],
               D0y=df['DstarD0_y'],
               D0z=df['DstarD0_z'],
               )
      else:  
         Dstar = JaggedCandidateArray.candidatesfromcounts(
               np.array([]),
               pt=np.array([]),
               eta=np.array([]),
               phi=np.array([]),
               mass=np.array([]),
               vtxIdx=np.array([]),
               D0x=np.array([]),
               D0y=np.array([]),
               D0z=np.array([]),
               )

      output['cutflow']['all events']  += Muon.size
      output['cutflow']['all Dimu']    += Dimu.counts.sum()
      output['cutflow']['all D0']      += D0.counts.sum()
      output['cutflow']['all Dstar']   += Dstar.counts.sum()

      # Dimu cuts: charge = 0...
      Dimu = Dimu[Dimu.charge == 0]
      output['cutflow']['Dimu 0 charge'] += Dimu.counts.sum()

      mass_cut = (Dimu.mass > 8.5) & (Dimu.mass < 11.5)
      Dimu = Dimu[mass_cut]
      output['cutflow']['Upsilon mass'] += Dimu.counts.sum()

      ############### Get the Muons from Dimu, for cuts in their params

      # For Muon 1:
      mu1Idx = (Dimu.t1Idx + Muon.starts).content
      mu2Idx = (Dimu.t2Idx + Muon.starts).content
      Muon1 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[mu1Idx])
      Muon2 = JaggedCandidateArray.fromoffsets(Dimu.offsets, Muon.content[mu2Idx])

      # SoftId and Global Muon cuts
      soft_id1 = (Muon1.softId > 0)
      Dimu = Dimu[soft_id1]
      Muon1 = Muon1[soft_id1]
      Muon2 = Muon2[soft_id1]
      output['cutflow']['Dimu muon1 softId'] += Dimu.counts.sum()

      soft_id2 = (Muon2.softId > 0)
      Dimu = Dimu[soft_id2]
      Muon1 = Muon1[soft_id2]
      Muon2 = Muon2[soft_id2]
      output['cutflow']['Dimu muon2 softId'] += Dimu.counts.sum()

      global_muon1 = (Muon1.isGlobal > 0)
      Dimu = Dimu[global_muon1]
      Muon1 = Muon1[global_muon1]
      Muon2 = Muon2[global_muon1]
      output['cutflow']['Dimu muon1 global'] += Dimu.counts.sum()

      global_muon2 = (Muon2.isGlobal > 0)
      Dimu = Dimu[global_muon2]
      Muon1 = Muon1[global_muon2]
      Muon2 = Muon2[global_muon2]
      output['cutflow']['Dimu muon2 global'] += Dimu.counts.sum()

      # pt and eta cuts
      pt_cut1 = (Muon1.pt > 3)
      Dimu = Dimu[pt_cut1]
      Muon1 = Muon1[pt_cut1]
      Muon2 = Muon2[pt_cut1]
      output['cutflow']['Muon1 pt cut'] += Dimu.counts.sum()

      pt_cut2 = (Muon2.pt > 3)
      Dimu = Dimu[pt_cut2]
      Muon1 = Muon1[pt_cut2]
      Muon2 = Muon2[pt_cut2]
      output['cutflow']['Muon2 pt cut'] += Dimu.counts.sum()

      eta_cut1 = (np.absolute(Muon1.eta) <= 2.4)
      Dimu = Dimu[eta_cut1]
      Muon1 = Muon1[eta_cut1]
      Muon2 = Muon2[eta_cut1]
      output['cutflow']['Muon1 eta cut'] += Dimu.counts.sum()

      eta_cut2 = (np.absolute(Muon2.eta) <= 2.4)
      Dimu = Dimu[eta_cut2]
      Muon1 = Muon1[eta_cut2]
      Muon2 = Muon2[eta_cut2]
      output['cutflow']['Muon2 eta cut'] += Dimu.counts.sum()

      # Cuts into events with at least 1 Dimu
      ndimu_cut = Dimu.counts > 0
      Dimu = Dimu[ndimu_cut]
      Muon1 = Muon1[ndimu_cut]
      Muon2 = Muon2[ndimu_cut]
      D0 = D0[ndimu_cut]
      Dstar = Dstar[ndimu_cut]

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
    