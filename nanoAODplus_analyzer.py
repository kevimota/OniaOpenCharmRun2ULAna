# coding: utf-8

import time

from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

from nanoAODplus_handler.AnalyzerProcessor import AnalyzerProcessor

import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)

tstart = time.time()    
'''fileset = {
    'Data10': [
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_500.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_501.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_502.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_503.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_504.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_505.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_506.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_507.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_508.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_509.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_510.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_511.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_512.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_513.root',
        'root://t2-cms-xrootd01.desy.de:1094//store/user/nujomhar/Data2010/CRAB_UserFiles/Data10MuOniaRunAtrigobj4/200301_205327/0000/Data10_Mu_trigobj4_514.root'

    ]
}'''

path = open("datasets/data2010_0000_path.txt")
files = path.read().splitlines()
files = files[0:80]

fileset = {'Data10': files}

output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=AnalyzerProcessor(),
                                  executor=processor.futures_executor,
                                  #executor_args={'workers': 6, 'flatten': True},
                                  executor_args={'workers': 8, 'flatten': True},
                                  chunksize=1000000000000,
                                 )

elapsed = time.time() - tstart

""" ax_mass = hist.plot1d(output['mass'], overlay='dataset')
fig = ax_mass.get_figure()
fig.savefig("plots/mass.png")

ax_pt = hist.plot1d(output['pt'], overlay='dataset')
ax_pt.set_xlim(0,20)
fig = ax_pt.get_figure()
fig.savefig("plots/pt.png")


ax_eta = hist.plot1d(output['eta'], overlay='dataset')
fig = ax_eta.get_figure()
fig.savefig("plots/eta.png")

ax_phi = hist.plot1d(output['phi'], overlay='dataset')
fig = ax_phi.get_figure()
fig.savefig("plots/phi.png") """

print("Events/s:", output['cutflow']['all events']/elapsed, "Time elapsed:", elapsed)
print(output['cutflow'])



