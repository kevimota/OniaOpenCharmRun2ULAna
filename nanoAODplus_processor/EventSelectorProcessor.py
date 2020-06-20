from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

class EventSelectorProcessor(processor.ProcessorABC):
   def __init__(self):
      dataset_axis = hist.Cat("dataset", "Primary dataset")

      muon_pt_axis = hist.Bin("pt", r"$p_{T,\mu}$ [GeV]", 3000, 0.25, 300)
      muon_eta_axis = hist.Bin("eta", r"$\eta_{\mu}$", 60, -3.0, 3.0)
      muon_phi_axis = hist.Bin("phi", r"$\phi_{\mu}$", 70, -3.5, 3.5)

      dimu_mass_axis = hist.Bin("mass", r"$m_{\mu\mu}$ [GeV]", 3600, 0.25, 120)
      dimu_pt_axis = hist.Bin("pt", r"$p_{T,\mu\mu}$ [GeV]", 3000, 0.25, 300)
      dimu_eta_axis = hist.Bin("eta", r"$\eta_{\mu\mu}$", 100, -5.0, 5.0)
      dimu_phi_axis = hist.Bin("phi", r"$\phi_{\mu\mu}$", 70, -3.5, 3.5)

      D0_mass_axis = hist.Bin("mass", r"$m_{D^0}$ [GeV]", 100, 0.25, 3)
      D0_pt_axis = hist.Bin("pt", r"$p_{T,D^0}$ [GeV]", 3000, 0.25, 300)
      D0_eta_axis = hist.Bin("eta", r"$\eta_{D^0}$", 100, -5.0, 5.0)
      D0_phi_axis = hist.Bin("phi", r"$\phi_{D^0}$", 70, -3.5, 3.5)
      
      self._accumulator = processor.dict_accumulator({
         'muon_pt': hist.Hist("Counts", dataset_axis, muon_pt_axis),
         'muon_eta': hist.Hist("Counts", dataset_axis, muon_eta_axis),
         'muon_phi': hist.Hist("Counts", dataset_axis, muon_phi_axis),
         'dimu_mass': hist.Hist("Counts", dataset_axis, dimu_mass_axis),
         'dimu_pt': hist.Hist("Counts", dataset_axis, dimu_pt_axis),
         'dimu_eta': hist.Hist("Counts", dataset_axis, dimu_eta_axis),
         'dimu_phi': hist.Hist("Counts", dataset_axis, dimu_phi_axis),
         'D0_mass': hist.Hist("Counts", dataset_axis, D0_mass_axis),
         'D0_pt': hist.Hist("Counts", dataset_axis, D0_pt_axis),
         'D0_eta': hist.Hist("Counts", dataset_axis, D0_eta_axis),
         'D0_phi': hist.Hist("Counts", dataset_axis, D0_phi_axis),
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
         muons = JaggedCandidateArray.candidatesfromcounts(
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
         muons = JaggedCandidateArray.candidatesfromcounts(
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

      output['cutflow']['all events'] += muons.size
      output['cutflow']['all muons'] += muons.counts.sum()
      output['cutflow']['all D0'] += D0.counts.sum()
      
      # global and soft muon
      soft_id = (muons.softId > 0)
      muons = muons[soft_id]
      output['cutflow']['soft muon'] += soft_id.sum().sum()

      global_muon = (muons.isGlobal > 0)
      muons = muons[global_muon]
      output['cutflow']['global muon'] += global_muon.sum().sum()

      #pt and eta cuts
      pt_cut = (muons.pt > 3)
      muons = muons[pt_cut]
      output['cutflow']['pt cut'] += pt_cut.sum().sum()

      eta_cut = (np.absolute(muons.eta) <= 2.4)
      muons = muons[eta_cut]
      output['cutflow']['eta cut'] += eta_cut.sum().sum()

      #isolated muon
      iso_muon = (muons.pfRelIso04_all < 0.4)
      muons = muons[iso_muon]
      output['cutflow']['iso muon'] += iso_muon.sum().sum()

      #valid vtx
      valid_vtx = (muons.vtxIdx != -1)
      muons = muons[valid_vtx]
      output['cutflow']['valid vtx'] += valid_vtx.sum().sum()

      #dimuon
      twomuons = (muons.counts >= 2)
      output['cutflow']['two muons'] += twomuons.sum()
      
      D0 = D0[twomuons]
      output['cutflow']['D0 two muons'] += D0.counts.sum()
      dimuons = muons[twomuons].distincts()

      opposite_charge = (dimuons.i0['charge'] * dimuons.i1['charge'] < 0)
      dimuons = dimuons[opposite_charge]
      output['cutflow']['opposite charge'] += opposite_charge.any().sum()

      #same vtx or close in z
      same_vtx = (dimuons.i0['vtxIdx'] == dimuons.i1['vtxIdx']) | (np.absolute(dimuons.i0['z'] - dimuons.i1['z']) < 0.2)
      dimuons = dimuons[same_vtx]
      output['cutflow']['same vtx'] += same_vtx.any().sum()
      
      output['muon_pt'].fill(dataset=dataset, pt=muons.pt.flatten())
      output['muon_eta'].fill(dataset=dataset, eta=muons.eta.flatten())
      output['muon_phi'].fill(dataset=dataset, phi=muons.phi.flatten())

      output['dimu_mass'].fill(dataset=dataset,mass=dimuons.mass.flatten())
      output['dimu_pt'].fill(dataset=dataset, pt=dimuons.pt.flatten())
      output['dimu_eta'].fill(dataset=dataset, eta=dimuons.eta.flatten())
      output['dimu_phi'].fill(dataset=dataset, phi=dimuons.phi.flatten())

      output['D0_mass'].fill(dataset=dataset,mass=D0.mass.flatten())
      output['D0_pt'].fill(dataset=dataset, pt=D0.pt.flatten())
      output['D0_eta'].fill(dataset=dataset, eta=D0.eta.flatten())
      output['D0_phi'].fill(dataset=dataset, phi=D0.phi.flatten())
      
      return output

   def postprocess(self, accumulator):
      return accumulator