from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np
from coffea.util import save

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
         'Muon_pt': processor.column_accumulator(np.zeros(shape=(0,))),
         'Muon_eta': processor.column_accumulator(np.zeros(shape=(0,))),
         'Muon_phi': processor.column_accumulator(np.zeros(shape=(0,))),
         'cutflow': processor.defaultdict_accumulator(int),
      })
    
   @property
   def accumulator(self):
      return self._accumulator
    
   def process(self, events):
      dataset = events.metadata['dataset']
      output = self.accumulator.identity()
      
      Muons = events.Muon
      D0 = events.D0
      Dstar = events.Dstar

      output['cutflow']['all events']  += Muons.size
      output['cutflow']['all muons']   += Muons.counts.sum()
      output['cutflow']['all D0']      += D0.counts.sum()
      output['cutflow']['all Dstar']   += Dstar.counts.sum()
      
      # global and soft muon
      soft_id = (Muons.softId > 0)
      Muons = Muons[soft_id]
      output['cutflow']['soft muon'] += Muons.counts.sum()

      global_muon = (Muons.isGlobal > 0)
      Muons = Muons[global_muon]
      output['cutflow']['global muon'] += Muons.counts.sum()

      #pt and eta cuts
      pt_cut = (Muons.pt > 3)
      Muons = Muons[pt_cut]
      output['cutflow']['pt cut'] += Muons.counts.sum()

      eta_cut = (np.absolute(Muons.eta) <= 2.4)
      Muons = Muons[eta_cut]
      output['cutflow']['eta cut'] += Muons.counts.sum()

      #isolated muon
      iso_muon = (Muons.pfRelIso04_all < 0.4)
      Muons = Muons[iso_muon]
      output['cutflow']['iso muon'] += Muons.counts.sum()

      #valid vtx
      valid_vtx = (Muons.vtxIdx != -1)
      Muons = Muons[valid_vtx]
      output['cutflow']['valid vtx'] += Muons.counts.sum()

      #dimuon
      twomuons = (Muons.counts > 1)
      Muons = Muons[twomuons]
      D0 = D0[twomuons]
      Dstar = Dstar[twomuons]
      output['cutflow']['two muons']         += Muons.counts.sum()
      output['cutflow']['D0 two muons']      += D0.counts.sum()
      output['cutflow']['Dstar two muons']   += Dstar.counts.sum()
      Dimuons = Muons.distincts()
      output['cutflow']['all dimuons'] += Dimuons.counts.sum()

      opposite_charge = (Dimuons.i0['charge'] * Dimuons.i1['charge'] < 0)
      Dimuons = Dimuons[opposite_charge]
      output['cutflow']['opposite charge'] += Dimuons.counts.sum()

      #same vtx or close in z
      same_vtx = (Dimuons.i0['vtxIdx'] == Dimuons.i1['vtxIdx']) | (np.absolute(Dimuons.i0['z'] - Dimuons.i1['z']) < 0.2)
      Dimuons = Dimuons[same_vtx]
      output['cutflow']['same vtx'] += Dimuons.counts.sum()

      # Only events with at least 1 dimuon
      evtcut = (Dimuons.counts > 0)
      Dimuons = Dimuons[evtcut]
      D0 = D0[evtcut]
      Dstar = Dstar[evtcut]
      output['cutflow']['D0 evt cut'] += D0.counts.sum()
      output['cutflow']['Dstar evt cut'] += Dstar.counts.sum()

      Dimuons = Dimuons.i0 + Dimuons.i1 

      output['Muon_pt'] += processor.column_accumulator(Muons.pt.flatten())
      output['Muon_eta'] += processor.column_accumulator(Muons.eta.flatten())
      output['Muon_phi'] += processor.column_accumulator(Muons.phi.flatten())
      
      output['muon_pt'].fill(dataset=dataset, pt=Muons.pt.flatten())
      output['muon_eta'].fill(dataset=dataset, eta=Muons.eta.flatten())
      output['muon_phi'].fill(dataset=dataset, phi=Muons.phi.flatten())

      output['dimu_mass'].fill(dataset=dataset,mass=Dimuons.mass.flatten())
      output['dimu_pt'].fill(dataset=dataset, pt=Dimuons.pt.flatten())
      output['dimu_eta'].fill(dataset=dataset, eta=Dimuons.eta.flatten())
      output['dimu_phi'].fill(dataset=dataset, phi=Dimuons.phi.flatten())

      output['D0_mass'].fill(dataset=dataset,mass=D0.mass12.flatten())
      output['D0_pt'].fill(dataset=dataset, pt=D0.pt.flatten())
      output['D0_eta'].fill(dataset=dataset, eta=D0.eta.flatten())
      output['D0_phi'].fill(dataset=dataset, phi=D0.phi.flatten())
      
      return output

   def postprocess(self, accumulator):
      return accumulator