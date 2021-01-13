from coffea import processor, hist

import awkward1 as ak
from coffea.util import load

'''

def create_plot1d(hist, save_name, log=False):
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.style.CMS)

    ax = plt.gca()


    plt.errorbar(hist.axes[0].centers,
             hist.view(),
             np.sqrt(hist.view()),
             fmt='.',
             color='blue',)

    hep.histplot(hist, ax=ax, color='blue')

    if log:
        ax.set_yscale('log')
    else:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,3), useMathText=True)

    ax.set_xlabel(hist.axes[0].metadata, loc='right')
    ax.set_ylabel("Counts", loc='top')

    # compute mean and std:
    mean = (hist.view() * hist.axes[0].centers).sum()/hist.sum()
    std = np.sqrt((hist.view()*((hist.axes[0].centers - mean)**2)).sum()/hist.sum())

    annotation = f"Total {hist.sum()}" \
                    + "\n" + f"Mean: {round(mean,2)}" \
                    + "\n" + f"Std: {round(std,2)}"
    
    ax.annotate(annotation, xy=(0.85, 0.85), xycoords='axes fraction', fontsize = "small",
                    ha='center', annotation_clip=False, bbox=dict(boxstyle='round', fc='None'))

    ax.set_xlim(hist.axes[0].edges[0], hist.axes[0].edges[-1] + hist.axes[0].widths[-1])
    
    fig = ax.get_figure()
    fig.savefig(save_name)
    ax.clear()
    fig.clear()

def create_plot2d(hist, save_name):
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.style.CMS)
    
    ax = plt.gca()

    hep.hist2dplot(hist, ax=ax)
    ax.set_xlabel(hist.axes[0].metadata, loc='right')
    ax.set_ylabel(hist.axes[1].metadata, loc='top')

    fig = ax.get_figure()
    fig.savefig(save_name)
    ax.clear()
    fig.clear()
'''
def build_p4(acc):
    p4 = ak.zip({'x': acc['x'].value, 
                 'y': acc['y'].value,
                 'z': acc['z'].value,
                 't': acc['t'].value}, with_name="LorentzVector")

    return p4

class HistogramingProcessor(processor.ProcessorABC):
    def __init__(self, analyzer_name):
        self.analyzer_name = analyzer_name
        
        self._accumulator = processor.dict_accumulator({
            'Muon_lead_p': hist.Hist("Events", 
                                   hist.Bin("pt", "$p_{T,\mu}$ [GeV]", 100, 0, 50),
                                   hist.Bin("eta", "$\eta_{\mu}$", 60, -2.5, 2.5),
                                   hist.Bin("phi", "$\phi_{\mu}$", 70, -3.5, 3.5)),
            'Muon_trail_p': hist.Hist("Events", 
                                       hist.Bin("pt", "$p_{T,\mu}$ [GeV]", 100, 0, 50),
                                       hist.Bin("eta", "$\eta_{\mu}$", 60, -2.5, 2.5),
                                       hist.Bin("phi", "$\phi_{\mu}$", 70, -3.5, 3.5)),
            'Upsilon_mass': hist.Hist("Events", hist.Bin("mass", "$m_{\mu^+\mu^-}$ [GeV]", 100, 8.6, 11)),
            'Upsilon_p': hist.Hist("Events", 
                                   hist.Bin("pt", "$p_{T,\mu^+\mu^-}$ [GeV]", 100, 0, 50),
                                   hist.Bin("eta", "$\eta_{\mu^+\mu^-}$", 60, -2.5, 2.5),
                                   hist.Bin("phi", "$\phi_{\mu^+\mu^-}$", 70, -3.5, 3.5)),
            'Jpsi_mass': hist.Hist("Events", hist.Bin("mass", "$m_{\mu^+\mu^-}$ [GeV]", 100, 2.95, 3.25)),
            'Jpsi_p': hist.Hist("Events", 
                                   hist.Bin("pt", "$p_{T,\mu^+\mu^-}$ [GeV]", 100, 0, 50),
                                   hist.Bin("eta", "$\eta_{\mu^+\mu^-}$", 60, -2.5, 2.5),
                                   hist.Bin("phi", "$\phi_{\mu^+\mu^-}$", 70, -3.5, 3.5)),
            'D0_mass12': hist.Hist("Events", hist.Bin("mass", "$m_{D^0, 12}$ [GeV]", 100, 1.7, 2.0)),
            'D0_mass21': hist.Hist("Events", hist.Bin("mass", "$m_{D^0, 21}$ [GeV]", 100, 1.7, 2.0)),
            'D0_p': hist.Hist("Events", 
                              hist.Bin("pt", "$p_{T,D^0}$ [GeV]", 100, 0, 50),
                              hist.Bin("eta", "$\eta_{D^0}$", 80, -2.5, 2.5),
                              hist.Bin("phi", "$\phi_{D^0}$", 70, -3.5, 3.5)),
            'D0_eta_mass': hist.Hist("Events",
                                     hist.Bin("eta", "$\eta_{D^0}$", 80, -2.5, 2.5),
                                     hist.Bin("mass", "$m_{D^0}$ [GeV]", 100, 1.7, 2.0)),
            'D0_trk_p': hist.Hist("Events", 
                                  hist.Bin("pt", "$p_{T,D^0 trks}$ [GeV]", 100, 0, 50),
                                  hist.Bin("eta", "$\eta_{D^0 trks}$", 80, -2.5, 2.5),
                                  hist.Bin("phi", "$\phi_{D^0 trks}$", 70, -3.5, 3.5)),
            'Dstar_p': hist.Hist("Events", 
                                 hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                 hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                 hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
            'Dstar_deltam': hist.Hist("Events", hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
            'Dstar_deltamr': hist.Hist("Events", hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
            'Dstar_pw': hist.Hist("Events", 
                                 hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                 hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                 hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
            'Dstar_deltamw': hist.Hist("Events", hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
            'Dstar_deltamrw': hist.Hist("Events", hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
            'Dstar_K_p': hist.Hist("Events", 
                                   hist.Bin("pt", "$p_{T,D* K}$ [GeV]", 100, 0, 30),
                                   hist.Bin("eta", "$\eta_{D* K}$", 60, -2.5, 2.5),
                                   hist.Bin("phi", "$\phi_{D* K}$", 70, -3.5, 3.5)),
            'Dstar_pi_p': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,D* \pi}$ [GeV]", 100, 0, 30),
                                    hist.Bin("eta", "$\eta_{D* \pi}$", 60, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{D* \pi}$", 70, -3.5, 3.5)),
            'Dstar_pis_p': hist.Hist("Events", 
                                     hist.Bin("pt", "$p_{T,D* \pi_s}$ [GeV]", 100, 0, 20),
                                     hist.Bin("eta", "$\eta_{D* \pi_s}$", 60, -2.5, 2.5),
                                     hist.Bin("phi", "$\phi_{D* \pi_s}$", 70, -3.5, 3.5)),
            'UpsilonDstar': processor.dict_accumulator({
                'Upsilon_mass': hist.Hist("Events", hist.Bin("mass", "$m_{\mu^+\mu^-}$ [GeV]", 100, 8.6, 11)),
                'Upsilon_p': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,\mu^+\mu^-}$ [GeV]", 100, 0, 50),
                                    hist.Bin("eta", "$\eta_{\mu^+\mu^-}$", 60, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{\mu^+\mu^-}$", 70, -3.5, 3.5)),
                'Upsilon_deltarap': hist.Hist("Events", hist.Bin("deltarap", "$\Delta y$", 50, -5, 5)),
                'UpsilonDstar_mass': hist.Hist("Events", hist.Bin("mass", "$m_{\Upsilon D*}$ [GeV]", 1000, 0, 20)),
                'Dstar_p': hist.Hist("Events", 
                                 hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                 hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                 hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
                'Dstar_deltam': hist.Hist("Events", hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
                'Dstar_deltamr': hist.Hist("Events", hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
                'Dstar_pw': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                    hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
                'Dstar_deltamw': hist.Hist("Events", hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
                'Dstar_deltamrw': hist.Hist("Events", hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
            'JpsiDstar': processor.dict_accumulator({
                'Jpsi_mass': hist.Hist("Events", hist.Bin("mass", "$m_{\mu^+\mu^-}$ [GeV]", 100, 2.95, 3.25)), 
                'Jpsi_p': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,\mu^+\mu^-}$ [GeV]", 100, 0, 50),
                                    hist.Bin("eta", "$\eta_{\mu^+\mu^-}$", 60, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{\mu^+\mu^-}$", 70, -3.5, 3.5)),
                'Jpsi_deltarap': hist.Hist("Events", hist.Bin("deltarap", "$\Delta y$", 50, -5, 5)),
                'JpsiDstar_mass': hist.Hist("Events", hist.Bin("mass", "$m_{J/\psi D*}$ [GeV]", 1000, 0, 20)),
                'Dstar_p': hist.Hist("Events", 
                                 hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                 hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                 hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
                'Dstar_deltam': hist.Hist("Events", hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
                'Dstar_deltamr': hist.Hist("Events", hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
                'Dstar_pw': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                    hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
                'Dstar_deltamw': hist.Hist("Events", hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
                'Dstar_deltamrw': hist.Hist("Events", hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
            }),
            }),
        })
     
    @property
    def accumulator(self):
        return self._accumulator
     
    def process(self, file):
        output = self.accumulator.identity()
        acc = load(file)

        Muon_lead_acc = acc['Muon_lead']
        Muon_trail_acc = acc['Muon_trail']
        Dimu_acc = acc['Dimu']
        D0_acc = acc['D0']
        D0_trk_acc = acc['D0_trk']
        Dstar_acc = acc['Dstar']
        Dstar_trk_acc = acc['Dstar_trk']
        DimuDstar_acc = acc['DimuDstar']

        DimuDstar_p4 = build_p4(DimuDstar_acc)

        # Filling histograms
        output['Muon_lead_p'].fill(pt=Muon_lead_acc['pt'].value,
                                   eta=Muon_lead_acc['eta'].value,
                                   phi=Muon_lead_acc['phi'].value)

        output['Muon_trail_p'].fill(pt=Muon_trail_acc['pt'].value,
                                    eta=Muon_trail_acc['eta'].value,
                                    phi=Muon_trail_acc['phi'].value)

        output['Upsilon_mass'].fill(mass=Dimu_acc['mass'].value[Dimu_acc['is_ups'].value])

        output['Upsilon_p'].fill(pt=Dimu_acc['pt'].value[Dimu_acc['is_ups'].value],
                                 eta=Dimu_acc['eta'].value[Dimu_acc['is_ups'].value],
                                 phi=Dimu_acc['phi'].value[Dimu_acc['is_ups'].value])

        output['Jpsi_mass'].fill(mass=Dimu_acc['mass'].value[Dimu_acc['is_jpsi'].value])

        output['Jpsi_p'].fill(pt=Dimu_acc['pt'].value[Dimu_acc['is_jpsi'].value],
                              eta=Dimu_acc['eta'].value[Dimu_acc['is_jpsi'].value],
                              phi=Dimu_acc['phi'].value[Dimu_acc['is_jpsi'].value])

        output['D0_mass12'].fill(mass=D0_acc['mass12'].value)
        output['D0_mass21'].fill(mass=D0_acc['mass21'].value)

        output['D0_p'].fill(pt=D0_acc['pt'].value,
                            eta=D0_acc['eta'].value,
                            phi=D0_acc['phi'].value)

        output['D0_eta_mass'].fill(eta=D0_acc['eta'].value,
                                   mass=D0_acc['mass'].value)

        output['D0_trk_p'].fill(pt=D0_trk_acc['t1_pt'].value,
                                eta=D0_trk_acc['t1_eta'].value,
                                phi=D0_trk_acc['t1_phi'].value)

        output['D0_trk_p'].fill(pt=D0_trk_acc['t2_pt'].value,
                                eta=D0_trk_acc['t2_eta'].value,
                                phi=D0_trk_acc['t2_phi'].value)
        
        output['Dstar_deltamr'].fill(deltamr=Dstar_acc['deltamr'].value)
        output['Dstar_deltam'].fill(deltam=Dstar_acc['deltam'].value)

        output['Dstar_p'].fill(pt=Dstar_acc['pt'].value,
                               eta=Dstar_acc['eta'].value,
                               phi=Dstar_acc['phi'].value)

        output['Dstar_deltamrw'].fill(deltamr=Dstar_acc['deltamr'].value[Dstar_acc['wrg_chg'].value])
        output['Dstar_deltamw'].fill(deltam=Dstar_acc['deltam'].value[Dstar_acc['wrg_chg'].value])

        output['Dstar_pw'].fill(pt=Dstar_acc['pt'].value[Dstar_acc['wrg_chg'].value],
                                eta=Dstar_acc['eta'].value[Dstar_acc['wrg_chg'].value],
                                phi=Dstar_acc['phi'].value[Dstar_acc['wrg_chg'].value])

        output['Dstar_K_p'].fill(pt=Dstar_trk_acc['K_pt'].value,
                                 eta=Dstar_trk_acc['K_eta'].value,
                                 phi=Dstar_trk_acc['K_phi'].value)

        output['Dstar_pi_p'].fill(pt=Dstar_trk_acc['pi_pt'].value,
                                  eta=Dstar_trk_acc['pi_eta'].value,
                                  phi=Dstar_trk_acc['pi_phi'].value)

        output['Dstar_pis_p'].fill(pt=Dstar_trk_acc['pis_pt'].value,
                                   eta=Dstar_trk_acc['pis_eta'].value,
                                   phi=Dstar_trk_acc['pis_phi'].value)

        ############# DimuDstar
        is_ups_dstar = DimuDstar_acc['Dimu']['is_ups'].value
        is_jpsi_dstar = DimuDstar_acc['Dimu']['is_jpsi'].value
        dimudstar_wrg_chg = DimuDstar_acc['Dstar']['wrg_chg'].value

        # Upsilon
        output['UpsilonDstar']['Upsilon_mass'].fill(mass=DimuDstar_acc['Dimu']['mass'].value[is_ups_dstar])

        output['UpsilonDstar']['Upsilon_p'].fill(pt=DimuDstar_acc['Dimu']['pt'].value[is_ups_dstar],
                                                 eta=DimuDstar_acc['Dimu']['eta'].value[is_ups_dstar],
                                                 phi=DimuDstar_acc['Dimu']['phi'].value[is_ups_dstar])

        output['UpsilonDstar']['Dstar_deltamr'].fill(deltamr=DimuDstar_acc['Dstar']['deltamr'].value[is_ups_dstar])
        output['UpsilonDstar']['Dstar_deltam'].fill(deltam=DimuDstar_acc['Dstar']['deltam'].value[is_ups_dstar])

        output['UpsilonDstar']['Dstar_p'].fill(pt=DimuDstar_acc['Dstar']['pt'].value[is_ups_dstar],
                                               eta=DimuDstar_acc['Dstar']['eta'].value[is_ups_dstar],
                                               phi=DimuDstar_acc['Dstar']['phi'].value[is_ups_dstar])

        output['UpsilonDstar']['Dstar_deltamrw'].fill(deltamr=DimuDstar_acc['Dstar']['deltamr'].value[is_ups_dstar & dimudstar_wrg_chg])
        output['UpsilonDstar']['Dstar_deltamw'].fill(deltam=DimuDstar_acc['Dstar']['deltam'].value[is_ups_dstar & dimudstar_wrg_chg])

        output['UpsilonDstar']['Dstar_pw'].fill(pt=DimuDstar_acc['Dstar']['pt'].value[is_ups_dstar & dimudstar_wrg_chg],
                                             eta=DimuDstar_acc['Dstar']['eta'].value[is_ups_dstar & dimudstar_wrg_chg],
                                             phi=DimuDstar_acc['Dstar']['phi'].value[is_ups_dstar & dimudstar_wrg_chg])

        output['UpsilonDstar']['Upsilon_deltarap'].fill(deltarap=DimuDstar_acc['deltarap'].value[is_ups_dstar])

        # Jpsi
        output['JpsiDstar']['Jpsi_mass'].fill(mass=DimuDstar_acc['Dimu']['mass'].value[is_jpsi_dstar])

        output['JpsiDstar']['Jpsi_p'].fill(pt=DimuDstar_acc['Dimu']['pt'].value[is_jpsi_dstar],
                                                 eta=DimuDstar_acc['Dimu']['eta'].value[is_jpsi_dstar],
                                                 phi=DimuDstar_acc['Dimu']['phi'].value[is_jpsi_dstar])

        output['JpsiDstar']['Dstar_deltamr'].fill(deltamr=DimuDstar_acc['Dstar']['deltamr'].value[is_jpsi_dstar])
        output['JpsiDstar']['Dstar_deltam'].fill(deltam=DimuDstar_acc['Dstar']['deltam'].value[is_jpsi_dstar])

        output['JpsiDstar']['Dstar_p'].fill(pt=DimuDstar_acc['Dstar']['pt'].value[is_jpsi_dstar],
                                               eta=DimuDstar_acc['Dstar']['eta'].value[is_jpsi_dstar],
                                               phi=DimuDstar_acc['Dstar']['phi'].value[is_jpsi_dstar])

        output['JpsiDstar']['Dstar_deltamrw'].fill(deltamr=DimuDstar_acc['Dstar']['deltamr'].value[is_jpsi_dstar & dimudstar_wrg_chg])
        output['JpsiDstar']['Dstar_deltamw'].fill(deltam=DimuDstar_acc['Dstar']['deltam'].value[is_jpsi_dstar & dimudstar_wrg_chg])

        output['JpsiDstar']['Dstar_pw'].fill(pt=DimuDstar_acc['Dstar']['pt'].value[is_jpsi_dstar & dimudstar_wrg_chg],
                                             eta=DimuDstar_acc['Dstar']['eta'].value[is_jpsi_dstar & dimudstar_wrg_chg],
                                             phi=DimuDstar_acc['Dstar']['phi'].value[is_jpsi_dstar & dimudstar_wrg_chg])

        output['JpsiDstar']['Jpsi_deltarap'].fill(deltarap=DimuDstar_acc['deltarap'].value[is_jpsi_dstar])

        return output

    def postprocess(self, accumulator):
        return accumulator      
