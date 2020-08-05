from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
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
      
      # Muon candidates
      output['dataset'] = df['dataset']
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
               vtxIdx=df['Muon_vtxIdx'],
               pfRelIso04_all=df['Muon_pfRelIso04_all'],
               x=df['Muon_x'],
               y=df['Muon_y'],
               z=df['Muon_z'],
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
               vtxIdx=np.array([]),
               pfRelIso04_all=np.array([]),
               x=np.array([]),
               y=np.array([]),
               z=np.array([]),
               )        
      
      # Dzero candidates
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
      output['cutflow']['all muons']   += Muon.counts.sum()
      output['cutflow']['all D0']      += D0.counts.sum()
      output['cutflow']['all Dstar']   += Dstar.counts.sum()
      
      # global and soft muon
      soft_id = (Muon.softId > 0)
      Muon = Muon[soft_id]
      output['cutflow']['soft muon'] += soft_id.sum().sum()

      global_muon = (Muon.isGlobal > 0)
      Muon = Muon[global_muon]
      output['cutflow']['global muon'] += global_muon.sum().sum()

      #pt and eta cuts
      pt_cut = (Muon.pt > 3)
      Muon = Muon[pt_cut]
      output['cutflow']['pt cut'] += Muon.counts.sum()

      eta_cut = (np.absolute(Muon.eta) <= 2.4)
      Muon = Muon[eta_cut]
      output['cutflow']['eta cut'] += Muon.counts.sum()

      #isolated muon
      iso_muon = (Muon.pfRelIso04_all < 0.4)
      Muon = Muon[iso_muon]
      output['cutflow']['iso muon'] += Muon.counts.sum()

      #valid vtx
      valid_vtx = (Muon.vtxIdx != -1)
      Muon = Muon[valid_vtx]
      output['cutflow']['valid vtx'] += Muon.counts.sum()

      #dimuon
      twomuons = (Muon.counts > 1)
      Muon = Muon[twomuons]
      D0 = D0[twomuons]
      Dstar = Dstar[twomuons]
      output['cutflow']['two muons']         += Muon.counts.sum()
      output['cutflow']['D0 two muons']      += D0.counts.sum()
      output['cutflow']['Dstar two muons']   += Dstar.counts.sum()

      Dimuon = Muon.distincts()
      output['cutflow']['all dimuons'] += Dimuon.counts.sum()

      #opposite charge muons
      opposite_charge = (Dimuon.i0['charge'] * Dimuon.i1['charge'] < 0)
      Dimuon = Dimuon[opposite_charge]
      output['cutflow']['opposite charge'] += Dimuon.counts.sum()

      #same vtx or close in z
      same_vtx = (Dimuon.i0['vtxIdx'] == Dimuon.i1['vtxIdx']) | (np.absolute(Dimuon.i0['z'] - Dimuon.i1['z']) < 0.2)
      Dimuon = Dimuon[same_vtx]
      output['cutflow']['same vtx'] += Dimuon.counts.sum()

      # Only events with at least 1 dimuon
      evtcut = (Dimuon.counts > 0)
      Dimuon = Dimuon[evtcut]
      D0 = D0[evtcut]
      Dstar = Dstar[evtcut]

      output['cutflow']['D0 evt cut'] += D0.counts.sum()
      output['cutflow']['Dstar evt cut'] += Dstar.counts.sum()
      
      # Leading and Trailing muon separation
      leading_mu = (Dimuon.i0.pt.content > Dimuon.i1.pt.content)
      Muon_lead = JaggedCandidateArray.candidatesfromoffsets(Dimuon.offsets, 
                                                       pt=np.where(leading_mu, Dimuon.i0.pt.content, Dimuon.i1.pt.content),
                                                       eta=np.where(leading_mu, Dimuon.i0.eta.content, Dimuon.i1.eta.content),
                                                       phi=np.where(leading_mu, Dimuon.i0.phi.content, Dimuon.i1.phi.content),
                                                       mass=np.where(leading_mu, Dimuon.i0.mass.content, Dimuon.i1.mass.content),)

      Muon_trail = JaggedCandidateArray.candidatesfromoffsets(Dimuon.offsets, 
                                                       pt=np.where(~leading_mu, Dimuon.i0.pt.content, Dimuon.i1.pt.content),
                                                       eta=np.where(~leading_mu, Dimuon.i0.eta.content, Dimuon.i1.eta.content),
                                                       phi=np.where(~leading_mu, Dimuon.i0.phi.content, Dimuon.i1.phi.content),
                                                       mass=np.where(~leading_mu, Dimuon.i0.mass.content, Dimuon.i1.mass.content),)

      #Create the accumulators to save output
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
      
      dimuon_acc = processor.dict_accumulator({})
      for var in columns:
         dimuon_acc[var] = processor.column_accumulator(np.array(Dimuon[var].flatten()))
      dimuon_acc["nDimuon"] = processor.column_accumulator(Dimuon.counts)
      output["Dimuon"] = dimuon_acc
      
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