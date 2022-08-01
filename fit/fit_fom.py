import ROOT
ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

import tools.root_plotting as plot
from tools.figure import create_fom

import os
import numpy as np

colors_hex = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
TColor = ROOT.TColor()

colors = [TColor.GetColor(color) for color in colors_hex]
styles = [1, 1, 2, 2, 2]

def fit_fom(param, path, value, fit_params, alpha_CB, n_CB, year):
    chain = ROOT.TChain("UpsilonDstar")
    for f in path:
        chain.Add(f"{f}/{param}/{value}/UpsilonDstar.root")

    wspace = ROOT.RooWorkspace(f"upsilondstar_fom_{param}{value}_{year}")
    upsilon_mass = ROOT.RooRealVar("dimu_mass", "Mass Upsilon", 8.7, 11.2)
    dstar_deltamr = ROOT.RooRealVar("dstar_deltamr", "Dstar Delta m", 0.141, 0.158)
    data = ROOT.RooDataSet("data", 
                       "Data 2D Upsilon + Dstar", 
                       ROOT.RooArgSet(upsilon_mass, dstar_deltamr), 
                       ROOT.RooFit.Import(chain))

    save_path = path[0][:path[0].rfind('/')] + f'/{param}/{value}'

    # Signal PDFs
    m1S = ROOT.RooRealVar("m1S","PDG mass Upsilon(1S)", fit_params['Upsilon']['M1S'])
    m2S = ROOT.RooRealVar("m2S","PDG mass Upsilon(2S)", fit_params['Upsilon']['M2S'])
    m3S = ROOT.RooRealVar("m3S","PDG mass Upsilon(3S)", fit_params['Upsilon']['M3S'])
    mscale = ROOT.RooRealVar("mscale", "mass scale factor", *fit_params['UpsilonDstar']['mscale'])
    mean1S  = ROOT.RooFormulaVar("mean1S","@0*@1", ROOT.RooArgList(m1S, mscale))
    mean2S  = ROOT.RooFormulaVar("mean2S","@0*@1", ROOT.RooArgList(m2S, mscale))
    mean3S  = ROOT.RooFormulaVar("mean3S","@0*@1", ROOT.RooArgList(m3S, mscale))
    sigma1S  = ROOT.RooRealVar("sigma1S","sigma1S", *fit_params['UpsilonDstar']['sigma1S'])
    sigma2S  = ROOT.RooFormulaVar("sigma2S","@0*@1/@2", ROOT.RooArgList(sigma1S, m2S, m1S))
    sigma3S  = ROOT.RooFormulaVar("sigma3S","@0*@1/@2", ROOT.RooArgList(sigma1S, m3S, m1S))
    upsilon1S_frac = ROOT.RooRealVar("upsilon1S_frac", "Upsilon(1S) fraction", *fit_params['UpsilonDstar']['upsilon1S_frac'])
    upsilon2S_frac = ROOT.RooRealVar("upsilon2S_frac", "Upsilon(2S) fraction", *fit_params['UpsilonDstar']['upsilon2S_frac']) 
    alpha_CB = ROOT.RooRealVar("alpha_CB", "alpha CB Upsilon", alpha_CB)
    n_CB = ROOT.RooRealVar("n_CB", "n CB Upsilon", n_CB)

    signal1S = ROOT.RooCBShape("signal1S", "Upsilon(1S) signal PDF", upsilon_mass, mean1S, sigma1S, alpha_CB, n_CB)
    signal2S = ROOT.RooCBShape("signal2S", "Upsilon(2S) signal PDF", upsilon_mass, mean2S, sigma2S, alpha_CB, n_CB)
    signal3S = ROOT.RooCBShape("signal3S", "Upsilon(3S) signal PDF", upsilon_mass, mean3S, sigma3S, alpha_CB, n_CB)
    upsilon_signal = ROOT.RooAddPdf("upsilon_signal", "Upsilon Signal PDF", ROOT.RooArgList(signal1S, signal2S, signal3S), 
                                    ROOT.RooArgList(upsilon1S_frac, upsilon2S_frac), ROOT.kTRUE)

    # Upsilon Bkg
    bkg_a1  = ROOT.RooRealVar("bkg_a1", "bkg_a1", *fit_params['UpsilonDstar']['bkg_a1'])
    bkg_a2  = ROOT.RooRealVar("bkg_a2", "bkg_a2", *fit_params['UpsilonDstar']['bkg_a2'])

    upsilon_bkg = ROOT.RooChebychev("upsilon_bkg", "Upsilon Background PDF", upsilon_mass, ROOT.RooArgList(bkg_a1, bkg_a2))

    # Dstar Signal Double Gaussian - same mean
    dstar_mean = ROOT.RooRealVar("dstar_mean", "Dstar Gaussian Mean", *fit_params['UpsilonDstar']['dstar_mean'])
    dstar_sigma1 = ROOT.RooRealVar("dstar_sigma1", "Dstar Gaussian 1 Sigma", *fit_params['UpsilonDstar']['dstar_sigma1'])
    dstar_sigma2 = ROOT.RooRealVar("dstar_sigma2", "Dstar Gaussian 2 Sigma", *fit_params['UpsilonDstar']['dstar_sigma2'])
    gauss1_frac = ROOT.RooRealVar("gauss1_frac", "Dstar Gaussian Fraction", *fit_params['UpsilonDstar']['gauss1_frac'])

    g1 = ROOT.RooGaussian("g1", "Dstar Gaussian 1", dstar_deltamr, dstar_mean, dstar_sigma1)
    g2 = ROOT.RooGaussian("g2", "Dstar Gaussian 2", dstar_deltamr, dstar_mean, dstar_sigma2)

    # Dstar Background = Phenomenological Threshold Function 
    p0 = ROOT.RooRealVar("p0","", *fit_params['UpsilonDstar']['p0'])
    p1 = ROOT.RooRealVar('p1',"", *fit_params['UpsilonDstar']['p1'])
    p2 = ROOT.RooRealVar('p2',"", *fit_params['UpsilonDstar']['p2'])

    dstar_signal = ROOT.RooAddPdf("dstar_model", "Dstar Model", ROOT.RooArgList(g1, g2),
                                ROOT.RooArgList(gauss1_frac), ROOT.kTRUE)

    dstar_bkg = ROOT.RooGenericPdf("dstar_bkg","Dstar Background PDF","(1 - exp(-(@0 -0.13957)/@1)) * (@0/0.13957)**@2 + @3 * (@0/0.13957 - 1)",
                                ROOT.RooArgList(dstar_deltamr,p0,p1,p2))

    signal = ROOT.RooProdPdf("signal", "Signal of 2D model", ROOT.RooArgList(upsilon_signal, dstar_signal))
    bkg1 = ROOT.RooProdPdf("bkg1", "Bkg1 of 2D model", ROOT.RooArgList(upsilon_signal, dstar_bkg))
    bkg2 = ROOT.RooProdPdf("bkg2", "Bkg2 of 2D model", ROOT.RooArgList(upsilon_bkg, dstar_signal))
    bkg3 = ROOT.RooProdPdf("bkg3", "Bkg3 of 2D model", ROOT.RooArgList(upsilon_bkg, dstar_bkg))

    signal_frac = ROOT.RooRealVar("signal_frac", "signal fraction", *fit_params['UpsilonDstar']['signal_frac'])
    bkg1_frac = ROOT.RooRealVar("bkg1_frac", "bkg1 fraction", *fit_params['UpsilonDstar']['bkg1_frac'])
    bkg2_frac = ROOT.RooRealVar("bkg2_frac", "bkg2 fraction", *fit_params['UpsilonDstar']['bkg2_frac'])

    model2D = ROOT.RooAddPdf("model2D", "2D Model Upsilon + Dstar", ROOT.RooArgList(signal, bkg1, bkg2, bkg3),
                                ROOT.RooArgList(signal_frac, bkg1_frac, bkg2_frac), ROOT.kTRUE)

    result = model2D.fitTo(data, ROOT.RooFit.BatchMode("cpu"), ROOT.RooFit.Save())

    print("Fit DONE. Saving workspace and params to " + save_path)
    wspace.Import(data)
    wspace.Import(model2D)
    wspace.Import(result)

    os.system('mkdir -p ' + save_path)

    wspace.writeToFile(save_path + "/UpsilonDstar_fit.root")

def plot_results(param, value, path, year, lumi):
    print("Saving plots to: plots/fom/" + year)

    plot.ModTDRStyle(width=800)
    plot.lumi_13TeV = f"{lumi:.2f} " + "fb^{-1}"
    if not os.path.exists(f"{path}/{param}/{value}/UpsilonDstar_fit.root"): 
        print(f"File {path}/{param}/{value}/UpsilonDstar_fit.root does not exist")
        return

    f = ROOT.TFile.Open(f"{path}/{param}/{value}/UpsilonDstar_fit.root")
    wspace = f.Get(f"upsilondstar_fom_{param}{value}_{year}")

    all_vars = []
    upsilon_mass = wspace.var("dimu_mass")
    dstar_deltamr = wspace.var("dstar_deltamr")
    data = wspace.data("data")
    model2D = wspace.pdf("model2D")
    n_params = len(model2D.getParameters(data).selectByAttrib('Constant', False))
    all_vars.append([upsilon_mass, dstar_deltamr, data, model2D])

    f.Close()

    # Canvas for plotting Upsilon projection 
    c1 = ROOT.TCanvas("Upsilon mass")
    frame_upsilon = upsilon_mass.frame()
    frame_upsilon.GetXaxis().SetTitle("M(\mu^+\mu^-) [GeV/c^2]")
    
    # Plot the Data
    data.plotOn(frame_upsilon, Name="Data", DataError="SumW2")

    # Plot the Model components
    model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Signal"), ROOT.RooFit.Components("signal"),
                ROOT.RooFit.LineStyle(styles[1]), ROOT.RooFit.LineColor(colors[1]))
    model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Background 1"), ROOT.RooFit.Components("bkg1"),
                ROOT.RooFit.LineStyle(styles[2]), ROOT.RooFit.LineColor(colors[2]))
    model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Background 2"), ROOT.RooFit.Components("bkg2"),
                ROOT.RooFit.LineStyle(styles[3]), ROOT.RooFit.LineColor(colors[3]))
    model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Background 3"), ROOT.RooFit.Components("bkg3"),
                ROOT.RooFit.LineStyle(styles[4]), ROOT.RooFit.LineColor(colors[4]))
    model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Model2D"), ROOT.RooFit.LineStyle(styles[0]),
                ROOT.RooFit.LineColor(colors[0]))

    chi2_upsilon = frame_upsilon.chiSquare("Model2D", "Data", n_params)

    leg_upsilon = ROOT.TLegend(0.73, 0.69, 0.94, 0.91)
    leg_upsilon.AddEntry(frame_upsilon.findObject("Data"), "Data", "LEP")
    leg_upsilon.AddEntry(frame_upsilon.findObject("Model2D"), "Model Fit", "L")
    leg_upsilon.AddEntry(frame_upsilon.findObject("Signal"), "Signal Fit", "L")
    leg_upsilon.AddEntry(frame_upsilon.findObject("Background 1"), "Background - #Upsilon_{sig} D^{*}_{bkg}", "L")
    leg_upsilon.AddEntry(frame_upsilon.findObject("Background 2"), "Background - #Upsilon_{bkg} D^{*}_{sig}", "L")
    leg_upsilon.AddEntry(frame_upsilon.findObject("Background 3"), "Background - #Upsilon_{bkg} D^{*}_{bkg}", "L")

    frame_upsilon.Draw()
    leg_upsilon.Draw("same")
    c1.Draw()

    plot.CMS_lumi(c1, 4, 1)

    chi2_label_ups = ROOT.TLatex()
    chi2_label_ups.SetNDC()
    chi2_label_ups.SetTextFont(43)
    chi2_label_ups.SetTextSize(15)
    chi2_label_ups.DrawLatex(.80, 0.65,"#chi^{2} = " + f'{chi2_upsilon:.2f}')

    c1.SaveAs(f"plots/fom/{year}/fit2D_upsilon_proj_{param}{value}_{year}.png")
    
    # Canvas for Dstar 
    c2 = ROOT.TCanvas("Dstar deltamr")

    # Frame
    frame_dstar = dstar_deltamr.frame(ROOT.RooFit.Title("Dstar delta m"))
    frame_dstar.GetXaxis().SetTitle("M(k\pi\pi_s) - M(k\pi) [GeV/c^2]")
    
    # Data
    data.plotOn(frame_dstar, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    
    # Plot the Model components
    model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Signal"), ROOT.RooFit.Components("signal"),
                ROOT.RooFit.LineStyle(styles[1]), ROOT.RooFit.LineColor(colors[1]))
    model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Background 1"), ROOT.RooFit.Components("bkg1"),
                ROOT.RooFit.LineStyle(styles[2]), ROOT.RooFit.LineColor(colors[2]))
    model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Background 2"), ROOT.RooFit.Components("bkg2"),
                ROOT.RooFit.LineStyle(styles[3]), ROOT.RooFit.LineColor(colors[3]))
    model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Background 3"), ROOT.RooFit.Components("bkg3"),
                ROOT.RooFit.LineStyle(styles[4]), ROOT.RooFit.LineColor(colors[4]))
    model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Model2D"), ROOT.RooFit.LineStyle(styles[0]),
                ROOT.RooFit.LineColor(colors[0]))

    chi2_dstar = frame_dstar.chiSquare("Model2D", "Data", n_params)

    # Legends
    leg_dstar = ROOT.TLegend(0.73, 0.69, 0.94, 0.91)
    leg_dstar.AddEntry(frame_dstar.findObject("Data"), "Data", "LEP")
    leg_dstar.AddEntry(frame_dstar.findObject("Model2D"), "Model Fit", "L")
    leg_dstar.AddEntry(frame_dstar.findObject("Signal"), "Signal Fit", "L")
    leg_dstar.AddEntry(frame_dstar.findObject("Background 1"), "Background - #Upsilon_{sig} D^{*}_{bkg}", "L")
    leg_dstar.AddEntry(frame_dstar.findObject("Background 2"), "Background - #Upsilon_{bkg} D^{*}_{sig}", "L")
    leg_dstar.AddEntry(frame_dstar.findObject("Background 3"), "Background - #Upsilon_{bkg} D^{*}_{bkg}", "L")

    frame_dstar.Draw()
    leg_dstar.Draw("same")
    c2.Draw()
    plot.CMS_lumi(c2, 4, 1)

    chi2_label_dstar = ROOT.TLatex()
    chi2_label_dstar.SetNDC()
    chi2_label_dstar.SetTextFont(43)
    chi2_label_dstar.SetTextSize(15)
    chi2_label_dstar.DrawLatex(.80, 0.65,"#chi^{2} = " + f'{chi2_dstar:.2f}')

    c2.SaveAs(f"plots/fom/{year}/fit2D_dstar_proj_{param}{value}_{year}.png")

    if (chi2_dstar > 3) or (chi2_upsilon > 3): return -1
    else: return 0


def plot_fom(param, config, path, year, lumi):
    n_sig_frac = np.array([])
    n_sig_frac_error = np.array([])
    n_evt = np.array([])

    for i in config:
        value = str(i).replace('.', 'p')
        f = ROOT.TFile.Open(f"{path}/{param}/{value}/UpsilonDstar_fit.root")
        wspace = f.Get(f"upsilondstar_fom_{param}{value}_{year}")

        data = wspace.data("data")
        model2D = wspace.pdf("model2D")
        signal_frac = wspace.var("signal_frac")
        
        f.Close()

        n_evt = np.append(n_evt, data.sumEntries())
        n_sig_frac = np.append(n_sig_frac, signal_frac.getVal())
        n_sig_frac_error = np.append(n_sig_frac_error, signal_frac.getErrorHi())

    n_sig = n_evt*n_sig_frac
    n_sig_error = n_evt*n_sig_frac_error
    n_bkg = n_evt - n_sig
    fom = n_sig/np.sqrt(n_sig + n_bkg)

    create_fom(fom, config, param, year, lumi=lumi)