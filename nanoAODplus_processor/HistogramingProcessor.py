from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
import boost_histogram as bh
from awkward import JaggedArray
import numpy as np
from coffea.util import save, load

class HistogramingProcessor(processor.ProcessorABC):
   def __init__(self, analyzer_name):
      self.analyzer_name = analyzer_name
      
      self._accumulator = processor.dict_accumulator({
         'foo': processor.defaultdict_accumulator(int),
      })
    
   @property
   def accumulator(self):
      return self._accumulator
    
   def process(self, df):
      output = self.accumulator.identity()
      acc = load("output/" + self.analyzer_name + "/merged/" + self.analyzer_name + "_merged.coffea")
      
      # Muon candidates
      Muon_lead = JaggedCandidateArray.candidatesfromcounts(acc["Muon_lead"]['nMuon'].value,
                                                            pt=acc["Muon_lead"]['__fast_pt'].value,
                                                            eta=acc["Muon_lead"]['__fast_eta'].value,
                                                            phi=acc["Muon_lead"]['__fast_phi'].value,
                                                            mass=acc["Muon_lead"]['__fast_mass'].value,)
      
      Muon_trail = JaggedCandidateArray.candidatesfromcounts(acc["Muon_trail"]['nMuon'].value,
                                                            pt=acc["Muon_trail"]['__fast_pt'].value,
                                                            eta=acc["Muon_trail"]['__fast_eta'].value,
                                                            phi=acc["Muon_trail"]['__fast_phi'].value,
                                                            mass=acc["Muon_trail"]['__fast_mass'].value,)

      Dimuon = JaggedCandidateArray.candidatesfromcounts(acc["Dimuon"]['nDimuon'].value,
                                                            pt=acc["Dimuon"]['__fast_pt'].value,
                                                            eta=acc["Dimuon"]['__fast_eta'].value,
                                                            phi=acc["Dimuon"]['__fast_phi'].value,
                                                            mass=acc["Dimuon"]['__fast_mass'].value,)

      D0 = JaggedCandidateArray.candidatesfromcounts(acc["D0"]['nD0'].value,
                                                            pt=acc["D0"]['__fast_pt'].value,
                                                            eta=acc["D0"]['__fast_eta'].value,
                                                            phi=acc["D0"]['__fast_phi'].value,
                                                            mass=acc["D0"]['__fast_mass'].value,)

      Dstar = JaggedCandidateArray.candidatesfromcounts(acc["Dstar"]['nDstar'].value,
                                                            pt=acc["Dstar"]['__fast_pt'].value,
                                                            eta=acc["Dstar"]['__fast_eta'].value,
                                                            phi=acc["Dstar"]['__fast_phi'].value,
                                                            mass=acc["Dstar"]['__fast_mass'].value,)

      

      # return dummy accumulator
      return output


   def postprocess(self, accumulator):
      return accumulator