import ROOT
ROOT.EnableImplicitMT()

ERAS = ['B', 'C', 'D', 'E', 'F']
#ERAS = ['C']
chain = ROOT.TChain("UpsilonDstar")

for i in ERAS:
     save_path = "/eos/user/k/kmotaama/MuOniaRun2017" + i + "_output/"
     dataset = "MuOniaRun2017" + i
     chain.Add(save_path + dataset + "_UpsilonDstar.root")

wspace = ROOT.RooWorkspace("UpsilonDstar_fit")

colors_hex = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
TColor = ROOT.TColor()

colors = [TColor.GetColor(color) for color in colors_hex]

'''for color in colors_hex:
     colors.append(TColor.GetColor(color))'''

#file = ROOT.TFile.Open(save_path + dataset + "_UpsilonDstar.root")

upsilon_mass = ROOT.RooRealVar("upsilon_mass", "Mass Upsilon", 8.7, 11.2)
dstar_deltam = ROOT.RooRealVar("dstar_deltam", "Dstar Delta m ", 0.14, 0.16)

data = ROOT.RooDataSet("data", 
                       "Data 2D Upsilon + Dstar", 
                       ROOT.RooArgSet(upsilon_mass, dstar_deltam), 
                       ROOT.RooFit.Import(chain))

# Upsilon Signal Crystal Ball

M1S = 9.46   # Upsilon(1S) PDG mass
M2S = 10.02  
M3S = 10.35 

# Signal PDFs
m1S = ROOT.RooRealVar("m1S","PDG mass Upsilon(1S)", M1S)
m2S = ROOT.RooRealVar("m2S","PDG mass Upsilon(2S)", M2S)
m3S = ROOT.RooRealVar("m3S","PDG mass Upsilon(3S)", M3S)
mscale = ROOT.RooRealVar("mscale", "mass scale factor", 1, 0, 1.2)

mean1S  = ROOT.RooFormulaVar("mean1S","@0*@1", ROOT.RooArgList(m1S, mscale))
mean2S  = ROOT.RooFormulaVar("mean2S","@0*@1", ROOT.RooArgList(m2S, mscale))
mean3S  = ROOT.RooFormulaVar("mean3S","@0*@1", ROOT.RooArgList(m3S, mscale))

sigma1S  = ROOT.RooRealVar("sigma1S","sigma1S", 0.08, 0, 1.2)
sigma2S  = ROOT.RooFormulaVar("sigma2S","@0*@1/@2", ROOT.RooArgList(sigma1S, m2S, m1S))
sigma3S  = ROOT.RooFormulaVar("sigma3S","@0*@1/@2", ROOT.RooArgList(sigma1S, m3S, m1S))

upsilon_1s_frac = ROOT.RooRealVar("upsilon_1s_frac", "Upsilon(1S) fraction", 0.3, 0, 1)
upsilon_2s_frac = ROOT.RooRealVar("upsilon_2s_frac", "Upsilon(2S) fraction", 0.15, 0, 1) 
#upsilon_3s_frac = ROOT.RooRealVar("upsilon_3s_frac", "Upsilon(3S) fraction", 0.1, 0, 1)

alpha = ROOT.RooRealVar("alpha", "alpha CB Upsilon", 1.2126e+00)
n = ROOT.RooRealVar("n", "n CB Upsilon", 8.7842e+00)

signal1S = ROOT.RooCBShape("signal1S", "Upsilon(1S) signal PDF", upsilon_mass, mean1S, sigma1S, alpha, n)
signal2S = ROOT.RooCBShape("signal2S", "Upsilon(2S) signal PDF", upsilon_mass, mean2S, sigma2S, alpha, n)
signal3S = ROOT.RooCBShape("signal3S", "Upsilon(3S) signal PDF", upsilon_mass, mean3S, sigma3S, alpha, n)

# Upsilon Bkg
bkg_a1  = ROOT.RooRealVar("bkg_a1", "bkg_a1", -0.06, -1, 1)
bkg_a2  = ROOT.RooRealVar("bkg_a2", "bkg_a2", -0.08, -1, 1)

upsilon_signal = ROOT.RooAddPdf("upsilon_signal", "Upsilon Signal PDF", ROOT.RooArgList(signal1S, signal2S, signal3S), 
                                ROOT.RooArgList(upsilon_1s_frac, upsilon_2s_frac), ROOT.kTRUE)

upsilon_bkg = ROOT.RooChebychev("upsilon_bkg", "Upsilon Background PDF", upsilon_mass, ROOT.RooArgList(bkg_a1, bkg_a2))

'''upsilon_model = ROOT.RooAddPdf("upsilon_model", "Upsilon Model", ROOT.RooArgList(signal1S, signal2S, signal3S, upsilon_bkg),
                               ROOT.RooArgList(upsilon_1s_frac, upsilon_2s_frac, upsilon_3s_frac), ROOT.kTRUE)'''

# Dstar Signal Double Gaussian - same mean
dstar_mean = ROOT.RooRealVar("mean", "Dstar Gaussian Mean", 0.1455, 0.142, 0.158)
dstar_sigma1 = ROOT.RooRealVar("sigma_1", "Dstar Gaussian 1 Sigma", 0.0001, 0.000001, 0.001)
dstar_sigma2 = ROOT.RooRealVar("sigma_2", "Dstar Gaussian 2 Sigma", 0.00001, 0.000001, 0.001)
gauss1_frac = ROOT.RooRealVar("gauss1_frac", "Dstar Gaussian Fraction", 0.3, 0, 1)
#gauss2_frac = ROOT.RooRealVar("gauss2_frac", "Dstar Signal Fraction", 0.3, 0, 1)

g1 = ROOT.RooGaussian("g1", "Dstar Gaussian 1", dstar_deltam, dstar_mean, dstar_sigma1)
g2 = ROOT.RooGaussian("g2", "Dstar Gaussian 2", dstar_deltam, dstar_mean, dstar_sigma2)

# Dstar Background = Phenomenological Threshold Function 
p0 = ROOT.RooRealVar("p0","", 0.005, -1, 1)
p1 = ROOT.RooRealVar('p1',"", -10, -50, 5)
p2 = ROOT.RooRealVar('p2',"", 5, 0, 10)

dstar_signal = ROOT.RooAddPdf("dstar_model", "Dstar Model", ROOT.RooArgList(g1, g2),
                              ROOT.RooArgList(gauss1_frac), ROOT.kTRUE)

dstar_bkg = ROOT.RooGenericPdf("dstar_bkg","Dstar Background PDF","(1 - exp(-(@0 -0.13957)/@1)) * (@0/0.13957)**@2 + @3 * (@0/0.13957 - 1)",
                               ROOT.RooArgList(dstar_deltam,p0,p1,p2))

'''dstar_model = ROOT.RooAddPdf("dstar_model", "Dstar Model", ROOT.RooArgList(g1, g2, dstar_bkg),
                       ROOT.RooArgList(gauss1_frac, gauss2_frac), ROOT.kTRUE)'''

signal = ROOT.RooProdPdf("signal", "Signal of 2D model", ROOT.RooArgList(upsilon_signal, dstar_signal))
bkg1 = ROOT.RooProdPdf("bkg1", "Bkg1 of 2D model", ROOT.RooArgList(upsilon_signal, dstar_bkg))
bkg2 = ROOT.RooProdPdf("bkg2", "Bkg2 of 2D model", ROOT.RooArgList(upsilon_bkg, dstar_signal))
bkg3 = ROOT.RooProdPdf("bkg3", "Bkg3 of 2D model", ROOT.RooArgList(upsilon_bkg, dstar_bkg))

signal_frac = ROOT.RooRealVar("signal_frac", "signal fraction", 0.3, 0, 1)
bkg1_frac = ROOT.RooRealVar("bkg1_frac", "bkg1 fraction", 0.3, 0, 1)
bkg2_frac = ROOT.RooRealVar("bkg2_frac", "bkg2 fraction", 0.3, 0, 1)

'''model2D = ROOT.RooProdPdf("model2D", "2D Model Upsilon + Dstar", ROOT.RooArgList(upsilon_model, dstar_model))'''
model2D = ROOT.RooAddPdf("model2D", "2D Model Upsilon + Dstar", ROOT.RooArgList(signal, bkg1, bkg2, bkg3),
                              ROOT.RooArgList(signal_frac, bkg1_frac, bkg2_frac), ROOT.kTRUE)
model2D.fitTo(data, ROOT.RooFit.Save())
result = model2D.fitTo(data, ROOT.RooFit.Save())

result.floatParsFinal().Print("S")

c1 = ROOT.TCanvas("c1")
frame_upsilon = upsilon_mass.frame(ROOT.RooFit.Title("Dimuon Invariant mass"))
data.plotOn(frame_upsilon, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Upsilon Model"))
frame_upsilon.Draw()
c1.Draw()

c2 = ROOT.TCanvas("c2")
frame_dstar = dstar_deltam.frame(ROOT.RooFit.Title("Dstar delta m"))
data.plotOn(frame_dstar, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Dstar Model"))

frame_dstar.Draw()
c2.Draw()

c3 = ROOT.TCanvas("c3")
ph2 = dstar_deltam.createHistogram("dstar vs upsilon pdf", upsilon_mass)
model2D.fillHistogram(ph2,ROOT.RooArgList(dstar_deltam,upsilon_mass))
ph2.Draw("SURF")
c3.Draw()

c4 = ROOT.TCanvas("c4")
dh2 = dstar_deltam.createHistogram("dstar vs upsilon data", upsilon_mass)
data.fillHistogram(dh2,ROOT.RooArgList(dstar_deltam,upsilon_mass))
dh2.Rebin2D(3,3)
dh2.Draw("LEGO")
c4.Draw()

getattr(wspace, "import")(data)
getattr(wspace, "import")(model2D)
getattr(wspace, "import")(result)

wspace.writeToFile("/eos/user/k/kmotaama/MuOniaRun2017_UpsilonDstar_fits.root")

x = input()
