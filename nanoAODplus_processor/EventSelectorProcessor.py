from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

class EventSelectorProcessor(processor.ProcessorABC):
   def __init__(self):
      dataset_axis = hist.Cat("dataset", "Primary dataset")

      Muon_pt_axis = hist.Bin("pt", r"$p_{T,\mu}$ [GeV]", 3000, 0.25, 300)
      Muon_eta_axis = hist.Bin("eta", r"$\eta_{\mu}$", 60, -3.0, 3.0)
      Muon_phi_axis = hist.Bin("phi", r"$\phi_{\mu}$", 70, -3.5, 3.5)

      Dimuon_mass_axis = hist.Bin("mass", r"$m_{\mu\mu}$ [GeV]", 3600, 0.25, 120)
      Dimuon_pt_axis = hist.Bin("pt", r"$p_{T,\mu\mu}$ [GeV]", 3000, 0.25, 300)
      Dimuon_eta_axis = hist.Bin("eta", r"$\eta_{\mu\mu}$", 100, -5.0, 5.0)
      Dimuon_phi_axis = hist.Bin("phi", r"$\phi_{\mu\mu}$", 70, -3.5, 3.5)

      D0_mass_axis = hist.Bin("mass", r"$m_{D^0}$ [GeV]", 100, 0.25, 3)
      D0_pt_axis = hist.Bin("pt", r"$p_{T,D^0}$ [GeV]", 3000, 0.25, 300)
      D0_eta_axis = hist.Bin("eta", r"$\eta_{D^0}$", 100, -5.0, 5.0)
      D0_phi_axis = hist.Bin("phi", r"$\phi_{D^0}$", 70, -3.5, 3.5)

      Dstar_mass_axis = hist.Bin("mass", r"$m_{D^*}$ [GeV]", 100, 0.25, 3)
      Dstar_pt_axis = hist.Bin("pt", r"$p_{T,D^*}$ [GeV]", 3000, 0.25, 300)
      Dstar_eta_axis = hist.Bin("eta", r"$\eta_{D^*}$", 100, -5.0, 5.0)
      Dstar_phi_axis = hist.Bin("phi", r"$\phi_{D^*}$", 70, -3.5, 3.5)
      
      self._accumulator = processor.dict_accumulator({
         'Muon_pt': hist.Hist("Counts", dataset_axis, Muon_pt_axis),
         'Muon_eta': hist.Hist("Counts", dataset_axis, Muon_eta_axis),
         'Muon_phi': hist.Hist("Counts", dataset_axis, Muon_phi_axis),
         'Dimuon_mass': hist.Hist("Counts", dataset_axis, Dimuon_mass_axis),
         'Dimuon_pt': hist.Hist("Counts", dataset_axis, Dimuon_pt_axis),
         'Dimuon_eta': hist.Hist("Counts", dataset_axis, Dimuon_eta_axis),
         'Dimuon_phi': hist.Hist("Counts", dataset_axis, Dimuon_phi_axis),
         'D0_mass': hist.Hist("Counts", dataset_axis, D0_mass_axis),
         'D0_pt': hist.Hist("Counts", dataset_axis, D0_pt_axis),
         'D0_eta': hist.Hist("Counts", dataset_axis, D0_eta_axis),
         'D0_phi': hist.Hist("Counts", dataset_axis, D0_phi_axis),
         'Dstar_mass': hist.Hist("Counts", dataset_axis, Dstar_mass_axis),
         'Dstar_pt': hist.Hist("Counts", dataset_axis, Dstar_pt_axis),
         'Dstar_eta': hist.Hist("Counts", dataset_axis, Dstar_eta_axis),
         'Dstar_phi': hist.Hist("Counts", dataset_axis, Dstar_phi_axis),
         'cutflow': processor.defaultdict_accumulator(int),
      })
    
   @property
   def accumulator(self):
      return self._accumulator
    
   def process(self, df):
      output = self.accumulator.identity()
      
      # Muon candidates
      dataset = df['dataset']
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
               mass=df['Dstar_mass'],
               vtxIdx=df['Dstar_vtxIdx'],
               x=df['Dstar_x'],
               y=df['Dstar_y'],
               z=df['Dstar_z'],
               )
      else:  
         Dstar = JaggedCandidateArray.candidatesfromcounts(
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
      output['cutflow']['two muons']         += Muons.counts.sum()
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
      
      output['muon_pt'].fill(dataset=dataset, pt=Muon.pt.flatten())
      output['muon_eta'].fill(dataset=dataset, eta=Muon.eta.flatten())
      output['muon_phi'].fill(dataset=dataset, phi=Muon.phi.flatten())

      output['dimu_mass'].fill(dataset=dataset, mass=Dimuon.mass.flatten())
      output['dimu_pt'].fill(dataset=dataset, pt=Dimuon.pt.flatten())
      output['dimu_eta'].fill(dataset=dataset, eta=Dimuon.eta.flatten())
      output['dimu_phi'].fill(dataset=dataset, phi=Dimuon.phi.flatten())

      output['D0_mass'].fill(dataset=dataset, mass=D0.mass.flatten())
      output['D0_pt'].fill(dataset=dataset, pt=D0.pt.flatten())
      output['D0_eta'].fill(dataset=dataset, eta=D0.eta.flatten())
      output['D0_phi'].fill(dataset=dataset, phi=D0.phi.flatten())

      output['Dstar_mass'].fill(dataset=dataset, mass=Dstar.mass.flatten())
      output['Dstar_pt'].fill(dataset=dataset, pt=Dstar.pt.flatten())
      output['Dstar_eta'].fill(dataset=dataset, eta=Dstar.eta.flatten())
      output['Dstar_phi'].fill(dataset=dataset, phi=Dstar.phi.flatten())
      
      return output

   def postprocess(self, accumulator):
      return accumulator