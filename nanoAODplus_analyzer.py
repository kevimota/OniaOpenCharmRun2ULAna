# coding: utf-8

import time

from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

from nanoAODplus_handler.AnalyzerProcessor import AnalyzerProcessor
from data.fileset import filesets

import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)

tstart = time.time()

files = {'Charmonium2017MINIAOD': filesets['Charmonium2017MINIAOD'][0:10], 
           'MuOnia2017MINIAOD': filesets['MuOnia2017MINIAOD'][0:10], 
           'DoubleMuon2017AOD': filesets['DoubleMuon2017AOD'][0:100]
          }

output = processor.run_uproot_job(files,
                                  treename='Events',
                                  processor_instance=AnalyzerProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 8, 'flatten': True},
                                  #executor=processor.iterative_executor,
                                  #executor_args={'flatten': True},
                                  chunksize=10000,
                                 )

elapsed = time.time() - tstart

ax_mass = hist.plot1d(output['mass'], overlay='dataset')
ax_mass.set_xlim(0,12)
fig = ax_mass.get_figure()
fig.savefig("plots/mass.png")
fig.clf()

ax_pt = hist.plot1d(output['pt'], overlay='dataset')
ax_pt.set_xlim(0,20)
fig = ax_pt.get_figure()
fig.savefig("plots/pt.png")
fig.clf()

ax_eta = hist.plot1d(output['eta'], overlay='dataset')
fig = ax_eta.get_figure()
fig.savefig("plots/eta.png")
fig.clf()

ax_phi = hist.plot1d(output['phi'], overlay='dataset')
fig = ax_phi.get_figure()
fig.savefig("plots/phi.png")
fig.clf()

print("Events/s:", output['cutflow']['all events']/elapsed, "Time elapsed:", elapsed)
print(output['cutflow'])