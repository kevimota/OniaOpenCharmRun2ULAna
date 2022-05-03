import ROOT
ROOT.EnableImplicitMT()

colors_hex = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
colors = []

TColor = ROOT.TColor()

for color in colors_hex:
     colors.append(TColor.GetColor(color))

file = ROOT.TFile.Open("/eos/user/k/kmotaama/MuOniaRun2017_UpsilonDstar_fits.root")
wspace = file.Get("UpsilonDstar_fit")

upsilon_mass = wspace.var("upsilon_mass")
dstar_deltam = wspace.var("dstar_deltam")
data = wspace.data("data")
model2D = wspace.pdf("model2D")
result = wspace.obj("fitresult_model2D_data")

c1 = ROOT.TCanvas("c1")
frame_upsilon = upsilon_mass.frame(ROOT.RooFit.Title("Dimuon Invariant mass"))
data.plotOn(frame_upsilon, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
#model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Signal"), ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(colors[1]))
#model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Bkg1"), ROOT.RooFit.Components("bkg1"), ROOT.RooFit.LineColor(colors[2]))
#model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Bkg2"), ROOT.RooFit.Components("bkg2"), ROOT.RooFit.LineColor(colors[3]))
#model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Bkg3"), ROOT.RooFit.Components("bkg3"), ROOT.RooFit.LineColor(colors[4]))
#model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Upsilon Model"), ROOT.RooFit.LineColor(colors[0]))
model2D.plotOn(frame_upsilon, ROOT.RooFit.Name("Upsilon Model"))

#leg_upsilon = ROOT.TLegend(0.7, 0.7, 0.93, 0.92)
#leg_upsilon.AddEntry(frame_upsilon.findObject("Data"), "Data", "LEP")
#leg_upsilon.AddEntry(frame_upsilon.findObject("Upsilon Model"), "Fit", "L")
#leg_upsilon.AddEntry(frame_upsilon.findObject("Signal"), "Signal Component", "L")
#leg_upsilon.AddEntry(frame_upsilon.findObject("Bkg1"), "Bkg1 Component", "L")
#leg_upsilon.AddEntry(frame_upsilon.findObject("Bkg2"), "Bkg2 Component", "L")
#leg_upsilon.AddEntry(frame_upsilon.findObject("Bkg3"), "Bkg3 Component", "L")

frame_upsilon.Draw()
#leg_upsilon.Draw("same")
c1.Draw()
c1.SaveAs('upsilon_fit.png')

c2 = ROOT.TCanvas("c2")
frame_dstar = dstar_deltam.frame(ROOT.RooFit.Title("Dstar delta m"))
data.plotOn(frame_dstar, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
#model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Signal"), ROOT.RooFit.Components("signal"), ROOT.RooFit.LineColor(colors[1]))
#model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Bkg1"), ROOT.RooFit.Components("bkg1"), ROOT.RooFit.LineColor(colors[2]))
#model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Bkg2"), ROOT.RooFit.Components("bkg2"), ROOT.RooFit.LineColor(colors[3]))
#model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Bkg3"), ROOT.RooFit.Components("bkg3"), ROOT.RooFit.LineColor(colors[4]))
#model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Upsilon Model"), ROOT.RooFit.LineColor(colors[0]))
model2D.plotOn(frame_dstar, ROOT.RooFit.Name("Upsilon Model"))

#leg_dstar = ROOT.TLegend(0.7, 0.7, 0.93, 0.92)
#leg_dstar.AddEntry(frame_dstar.findObject("Data"), "Data", "LEP")
#leg_dstar.AddEntry(frame_dstar.findObject("Upsilon Model"), "Fit", "L")
#leg_dstar.AddEntry(frame_dstar.findObject("Signal"), "Signal Component", "L")
#leg_dstar.AddEntry(frame_dstar.findObject("Bkg1"), "Bkg1 Component", "L")
#leg_dstar.AddEntry(frame_dstar.findObject("Bkg2"), "Bkg2 Component", "L")
#leg_dstar.AddEntry(frame_dstar.findObject("Bkg3"), "Bkg3 Component", "L")

frame_dstar.Draw()
#leg_dstar.Draw("same")
c2.Draw()
c2.SaveAs('dstar_fit.png')

c3 = ROOT.TCanvas("c3")
ph2 = dstar_deltam.createHistogram("dstar vs upsilon pdf", upsilon_mass)
model2D.fillHistogram(ph2,ROOT.RooArgList(dstar_deltam,upsilon_mass))
ph2.Draw("SURF")
c3.Draw()

c4 = ROOT.TCanvas("c4")
dh2 = dstar_deltam.createHistogram("dstar vs upsilon data", upsilon_mass)
data.fillHistogram(dh2,ROOT.RooArgList(dstar_deltam,upsilon_mass))
#dh2.Rebin2D(3,3)
dh2.SetStats(0)
dh2.Draw("COLZ")
c4.Draw()

result.floatParsFinal().Print("S")
value = result.floatParsFinal().find("signal_frac").getVal()
error = result.floatParsFinal().find("signal_frac").getError()

print(f"{data.sumEntries()*value} +/- {data.sumEntries()*error}")

x = input()

