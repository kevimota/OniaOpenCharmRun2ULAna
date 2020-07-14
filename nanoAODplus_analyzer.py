# coding: utf-8

import time

from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

from nanoAODplus_processor.EventSelectorProcessor import EventSelectorProcessor
from data.fileset import filesets
import yaml

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

config_yaml = yaml.load(open("config/local.yaml", "r"), Loader=yaml.FullLoader)

if config_yaml['executor'] == 'futures_executor': 
    executor = processor.futures_executor

tstart = time.time()

files = {'Charmonium2017MINIAOD': filesets['Charmonium2017MINIAOD'][0:10], 
           'MuOnia2017MINIAOD': filesets['MuOnia2017MINIAOD'][0:10], 
           'DoubleMuon2017AOD': filesets['DoubleMuon2017AOD'][0:30]
          }

output = processor.run_uproot_job(files,
                                  treename='Events',
                                  processor_instance=EventSelectorProcessor(),
                                  executor=executor,
                                  executor_args={'workers': config_yaml['n_cores'], 'flatten': True},
                                  chunksize=config_yaml['chunksize'],
                                 )

elapsed = time.time() - tstart

ax_muon_pt = hist.plot1d(output['muon_pt'], overlay='dataset')
ax_muon_pt.set_xlim(0,100)
fig = ax_muon_pt.get_figure()
fig.savefig("plots/muon_pt.png")
fig.clf()

ax_muon_eta = hist.plot1d(output['muon_eta'], overlay='dataset')
fig = ax_muon_eta.get_figure()
fig.savefig("plots/muon_eta.png")
fig.clf()

ax_muon_phi = hist.plot1d(output['muon_phi'], overlay='dataset')
fig = ax_muon_phi.get_figure()
fig.savefig("plots/muon_phi.png")
fig.clf()

ax_dimu_mass = hist.plot1d(output['dimu_mass'], overlay='dataset')
ax_dimu_mass.set_xlim(0,12)
fig = ax_dimu_mass.get_figure()
fig.savefig("plots/dimu_mass.png")
fig.clf()

ax_dimu_pt = hist.plot1d(output['dimu_pt'], overlay='dataset')
ax_dimu_pt.set_xlim(0,20)
fig = ax_dimu_pt.get_figure()
fig.savefig("plots/dimu_pt.png")
fig.clf()

ax_dimu_eta = hist.plot1d(output['dimu_eta'], overlay='dataset')
fig = ax_dimu_eta.get_figure()
fig.savefig("plots/dimu_eta.png")
fig.clf()

ax_dimu_phi = hist.plot1d(output['dimu_phi'], overlay='dataset')
fig = ax_dimu_phi.get_figure()
fig.savefig("plots/dimu_phi.png")
fig.clf()

print("Events/s:", output['cutflow']['all events']/elapsed, "Time elapsed:", elapsed)
print(output['cutflow'])