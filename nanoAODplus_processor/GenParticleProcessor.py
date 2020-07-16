from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

class GenParticleProcessor(processor.ProcessorABC):
   def __init__(self):
      dataset_axis = hist.Cat("dataset", "Primary dataset")

      Muon_pt_axis = hist.Bin("pt", r"$p_{T,\mu}$ [GeV]", 3000, 0.25, 300)
      Muon_eta_axis = hist.Bin("eta", r"$\eta_{\mu}$", 60, -3.0, 3.0)
      Muon_phi_axis = hist.Bin("phi", r"$\phi_{\mu}$", 70, -3.5, 3.5)

      Upsilon_mass_axis = hist.Bin("mass", r"$m_{\mu\mu}$ [GeV]", 16, 9, 11)
      Upsilon_pt_axis = hist.Bin("pt", r"$p_{T,\mu\mu}$ [GeV]", 3000, 0.25, 300)
      Upsilon_eta_axis = hist.Bin("eta", r"$\eta_{\mu\mu}$", 100, -5.0, 5.0)
      Upsilon_phi_axis = hist.Bin("phi", r"$\phi_{\mu\mu}$", 70, -3.5, 3.5)

      D0_mass_axis = hist.Bin("mass", r"$m_{D^0}$ [GeV]", 10, 1.2, 2.2)
      D0_pt_axis = hist.Bin("pt", r"$p_{T,D^0}$ [GeV]", 3000, 0.25, 300)
      D0_eta_axis = hist.Bin("eta", r"$\eta_{D^0}$", 100, -5.0, 5.0)
      D0_phi_axis = hist.Bin("phi", r"$\phi_{D^0}$", 70, -3.5, 3.5)
      
      self._accumulator = processor.dict_accumulator({
         'Muon_pt': hist.Hist("Counts", dataset_axis, Muon_pt_axis),
         'Muon_eta': hist.Hist("Counts", dataset_axis, Muon_eta_axis),
         'Muon_phi': hist.Hist("Counts", dataset_axis, Muon_phi_axis),
         'Upsilon_mass': hist.Hist("Counts", dataset_axis, Upsilon_mass_axis),
         'Upsilon_pt': hist.Hist("Counts", dataset_axis, Upsilon_pt_axis),
         'Upsilon_eta': hist.Hist("Counts", dataset_axis, Upsilon_eta_axis),
         'Upsilon_phi': hist.Hist("Counts", dataset_axis, Upsilon_phi_axis),
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
      
      # GenParticles
      dataset = df['dataset']
      if df['nGenPart'].size != 0:
         GenPart = JaggedCandidateArray.candidatesfromcounts(
               df['nGenPart'],
               pt=df['GenPart_pt'],
               eta=df['GenPart_eta'],
               phi=df['GenPart_phi'],
               mass=df['GenPart_mass'],
               charge=df['GenPart_charge'],
               pdgId=df['GenPart_pdgId'],
               vx=df['GenPart_vx'],
               vy=df['GenPart_vy'],
               vz=df['GenPart_vz'],
               mpdgId=df['GenPart_mpdgId'],
               mvx=df['GenPart_mvx'],
               mvy=df['GenPart_mvy'],
               mvz=df['GenPart_mvz'],
               )
      else: 
         GenPart = JaggedCandidateArray.candidatesfromcounts(
               np.array([]),
               pt=np.array([]),
               eta=np.array([]),
               phi=np.array([]),
               mass=np.array([]),
               charge=np.array([]),
               pdgId=np.array([]),
               vx=np.array([]),
               vy=np.array([]),
               vz=np.array([]),
               mpdgId=np.array([]),
               mvx=np.array([]),
               mvy=np.array([]),
               mvz=np.array([]),
               ) 

         GenPart = JaggedCandidateArray.candidatesfromcounts(
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

      muonid = (np.absolute(GenPart.pdgId) == 13)
      Muon = GenPart[muonid]

      upsilonid = (np.absolute(GenPart.pdgId) == 553)
      Upsilon = GenPart[upsilonid]

      d0id = (np.absolute(GenPart.pdgId) == 421)
      D0 = GenPart[d0id]
      
      output['Muon_pt'].fill(dataset=dataset, pt=Muon.pt.flatten())
      output['Muon_eta'].fill(dataset=dataset, eta=Muon.eta.flatten())
      output['Muon_phi'].fill(dataset=dataset, phi=Muon.phi.flatten())

      output['Upsilon_mass'].fill(dataset=dataset, mass=Upsilon.mass.flatten())
      output['Upsilon_pt'].fill(dataset=dataset, pt=Upsilon.pt.flatten())
      output['Upsilon_eta'].fill(dataset=dataset, eta=Upsilon.eta.flatten())
      output['Upsilon_phi'].fill(dataset=dataset, phi=Upsilon.phi.flatten())

      output['D0_mass'].fill(dataset=dataset, mass=D0.mass.flatten())
      output['D0_pt'].fill(dataset=dataset, pt=D0.pt.flatten())
      output['D0_eta'].fill(dataset=dataset, eta=D0.eta.flatten())
      output['D0_phi'].fill(dataset=dataset, phi=D0.phi.flatten())
      
      return output

   def postprocess(self, accumulator):
      return accumulator