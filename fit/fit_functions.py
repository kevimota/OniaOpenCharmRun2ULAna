import ROOT
ROOT.EnableImplicitMT()

import tools.root_plotting as plot

import yaml, os

colors_hex = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
TColor = ROOT.TColor()

colors = [TColor.GetColor(color) for color in colors_hex]
styles = [1, 1, 2, 2, 2]

def save_fit_params(path, result, wspace, channel):
    fit_params = {}
    for p in result.floatParsFinal():
        name = p.GetName()
        param = wspace.var(name)
        fit_params[name] = {
                'value': param.getVal(),
                'error_hi': param.getErrorHi(),
                'error_lo': param.getErrorLo(),
                }

    data = wspace.data("data")
    n_evt = data.sumEntries()
    n_params = len(result.floatParsFinal())

    fit_params['n_evt'] = n_evt
    fit_params['n_params'] = n_params

    with open(f'{path}/{channel}_fit_params.yaml', 'w') as f:
        yaml.dump(fit_params, f)

def fit_upsilon(config, year):
    chain = ROOT.TChain("Upsilon")
    for f in config['path'][year]:
        chain.Add(f"{f}/Upsilon.root")

    wspace = ROOT.RooWorkspace(f"upsilon_fit_{year}")
    mass = ROOT.RooRealVar("mass", "Mass Upsilon", 8.7, 11.2)
    data = ROOT.RooDataSet("data", "Data Upsilon", ROOT.RooArgSet(mass), ROOT.RooFit.Import(chain))

    # Signal Variables
    m1S = ROOT.RooRealVar("m1S","PDG mass Upsilon(1S)", config['Upsilon']['M1S'])
    m2S = ROOT.RooRealVar("m2S","PDG mass Upsilon(2S)", config['Upsilon']['M2S'])
    m3S = ROOT.RooRealVar("m3S","PDG mass Upsilon(3S)", config['Upsilon']['M3S'])
    mscale = ROOT.RooRealVar("mscale", "mass scale factor", *config['Upsilon']['mscale'])
    mean1S  = ROOT.RooFormulaVar("mean1S","@0*@1", ROOT.RooArgList(m1S, mscale))
    mean2S  = ROOT.RooFormulaVar("mean2S","@0*@1", ROOT.RooArgList(m2S, mscale))
    mean3S  = ROOT.RooFormulaVar("mean3S","@0*@1", ROOT.RooArgList(m3S, mscale))
    sigma1S  = ROOT.RooRealVar("sigma1S","sigma1S", *config['Upsilon']['sigma1S'])
    sigma2S  = ROOT.RooFormulaVar("sigma2S","@0*@1/@2", ROOT.RooArgList(sigma1S, m2S, m1S))
    sigma3S  = ROOT.RooFormulaVar("sigma3S","@0*@1/@2", ROOT.RooArgList(sigma1S, m3S, m1S))
    alpha_CB = ROOT.RooRealVar("alpha_CB", "alpha CB Upsilon", *config['Upsilon']['alpha_CB'])
    n_CB = ROOT.RooRealVar("n_CB", "n CB Upsilon", *config['Upsilon']['n_CB'])

    # Signal pdfs
    signal1S = ROOT.RooCBShape("signal1S", "Upsilon(1S) signal PDF", mass, mean1S, sigma1S, alpha_CB, n_CB)
    signal2S = ROOT.RooCBShape("signal2S", "Upsilon(2S) signal PDF", mass, mean2S, sigma2S, alpha_CB, n_CB)
    signal3S = ROOT.RooCBShape("signal3S", "Upsilon(3S) signal PDF", mass, mean3S, sigma3S, alpha_CB, n_CB)

    # Background variables
    bkg_a1  = ROOT.RooRealVar("bkg_a1", "bkg_a1", *config['Upsilon']['bkg_a1'])
    bkg_a2  = ROOT.RooRealVar("bkg_a2", "bkg_a2", *config['Upsilon']['bkg_a2'])

    # Background pdf
    bkg = ROOT.RooChebychev("bkg", "Upsilon Background PDF", mass, ROOT.RooArgList(bkg_a1, bkg_a2))

    # Model
    upsilon1S_frac = ROOT.RooRealVar("upsilon1S_frac", "Upsilon(1S) fraction", *config['Upsilon']['upsilon1S_frac'])
    upsilon2S_frac = ROOT.RooRealVar("upsilon2S_frac", "Upsilon(2S) fraction", *config['Upsilon']['upsilon2S_frac'])
    upsilon3S_frac = ROOT.RooRealVar("upsilon3S_frac", "Upsilon(2S) fraction", *config['Upsilon']['upsilon3S_frac'])

    model = ROOT.RooAddPdf("model", "Upsilon Model", ROOT.RooArgList(signal1S, signal2S, signal3S, bkg),
                       ROOT.RooArgList(upsilon1S_frac, upsilon2S_frac, upsilon3S_frac), ROOT.kTRUE)

    result = model.fitTo(data, ROOT.RooFit.BatchMode("cpu"), ROOT.RooFit.Save()) #fitting
    
    print("Fit DONE. Saving workspace and params")
    getattr(wspace, "import")(data)
    getattr(wspace, "import")(model)
    getattr(wspace, "import")(result)

    save_path = config['path'][year][0][:config['path'][year][0].rfind('/')]
    wspace.writeToFile(save_path + "/Upsilon_fit.root")
    save_fit_params(save_path, result, wspace, "Upsilon")

def fit_dstar(config, year):
    chain = ROOT.TChain("Dstar")
    for f in config['path'][year]:
        chain.Add(f"{f}/Dstar.root")

    wspace = ROOT.RooWorkspace(f"dstar_fit_{year}")
    deltamr = ROOT.RooRealVar("deltamr", "Dstar Delta m ", 0.141, 0.158)
    data = ROOT.RooDataSet("data", "Data Dstar", ROOT.RooArgSet(deltamr), ROOT.RooFit.Import(chain))

    # Signal Double Gaussian - same mean
    mean = ROOT.RooRealVar("mean", "Dstar Gaussian Mean", *config['Dstar']['mean'])
    sigma1 = ROOT.RooRealVar("sigma1", "Dstar Gaussian 1 Sigma", *config['Dstar']['sigma1'])
    sigma2 = ROOT.RooRealVar("sigma2", "Dstar Gaussian 2 Sigma", *config['Dstar']['sigma2'])

    g1 = ROOT.RooGaussian("g1", "Dstar Gaussian 1", deltamr, mean, sigma1)
    g2 = ROOT.RooGaussian("g2", "Dstar Gaussian 2", deltamr, mean, sigma2)

    # Background = Phenomenological Threshold Function 
    p0 = ROOT.RooRealVar("p0","", *config['Dstar']['p0'])
    p1 = ROOT.RooRealVar('p1',"", *config['Dstar']['p1'])
    p2 = ROOT.RooRealVar('p2',"", *config['Dstar']['p2'])

    bkg = ROOT.RooGenericPdf("bkg","Dstar Background PDF","(1 - exp(-(@0 -0.13957)/@1)) * (@0/0.13957)**@2 + @3 * (@0/0.13957 - 1)",
                            ROOT.RooArgList(deltamr, p0, p1, p2))

    gauss1_frac = ROOT.RooRealVar("gauss1_frac", "Dstar Gaussian Fraction", *config['Dstar']['gauss1_frac'])
    gauss2_frac = ROOT.RooRealVar("gauss2_frac", "Dstar Signal Fraction", *config['Dstar']['gauss2_frac'])

    model = ROOT.RooAddPdf("model", "Dstar Model", ROOT.RooArgList(g1, g2, bkg),
                       ROOT.RooArgList(gauss1_frac, gauss2_frac), ROOT.kTRUE)

    result = model.fitTo(data, ROOT.RooFit.BatchMode("cpu"), ROOT.RooFit.Save()) #fitting

    print("Fit DONE. Saving workspace and params")
    getattr(wspace, "import")(data)
    getattr(wspace, "import")(model)
    getattr(wspace, "import")(result)

    save_path = config['path'][year][0][:config['path'][year][0].rfind('/')]
    wspace.writeToFile(save_path + "/Dstar_fit.root")
    save_fit_params(save_path, result, wspace, "Dstar")

def fit_upsilondstar(config, year):
    chain = ROOT.TChain("UpsilonDstar")
    for f in config['path'][year]:
        chain.Add(f"{f}/UpsilonDstar.root")

    wspace = ROOT.RooWorkspace(f"upsilondstar_fit_{year}")
    upsilon_mass = ROOT.RooRealVar("dimu_mass", "Mass Upsilon", 8.7, 11.2)
    dstar_deltamr = ROOT.RooRealVar("dstar_deltamr", "Dstar Delta m ", 0.141, 0.158)
    data = ROOT.RooDataSet("data", 
                       "Data 2D Upsilon + Dstar", 
                       ROOT.RooArgSet(upsilon_mass, dstar_deltamr), 
                       ROOT.RooFit.Import(chain))

    # Load previous Upsilon results
    save_path = config['path'][year][0][:config['path'][year][0].rfind('/')]
    with open(f"{save_path}/Upsilon_fit_params.yaml") as f:
        upsilon_fit = yaml.load(f, Loader=yaml.FullLoader)

    # Signal PDFs
    m1S = ROOT.RooRealVar("m1S","PDG mass Upsilon(1S)", config['Upsilon']['M1S'])
    m2S = ROOT.RooRealVar("m2S","PDG mass Upsilon(2S)", config['Upsilon']['M2S'])
    m3S = ROOT.RooRealVar("m3S","PDG mass Upsilon(3S)", config['Upsilon']['M3S'])
    mscale = ROOT.RooRealVar("mscale", "mass scale factor", *config['UpsilonDstar']['mscale'])
    mean1S  = ROOT.RooFormulaVar("mean1S","@0*@1", ROOT.RooArgList(m1S, mscale))
    mean2S  = ROOT.RooFormulaVar("mean2S","@0*@1", ROOT.RooArgList(m2S, mscale))
    mean3S  = ROOT.RooFormulaVar("mean3S","@0*@1", ROOT.RooArgList(m3S, mscale))
    sigma1S  = ROOT.RooRealVar("sigma1S","sigma1S", *config['UpsilonDstar']['sigma1S'])
    sigma2S  = ROOT.RooFormulaVar("sigma2S","@0*@1/@2", ROOT.RooArgList(sigma1S, m2S, m1S))
    sigma3S  = ROOT.RooFormulaVar("sigma3S","@0*@1/@2", ROOT.RooArgList(sigma1S, m3S, m1S))
    upsilon1S_frac = ROOT.RooRealVar("upsilon1S_frac", "Upsilon(1S) fraction", *config['UpsilonDstar']['upsilon1S_frac'])
    upsilon2S_frac = ROOT.RooRealVar("upsilon2S_frac", "Upsilon(2S) fraction", *config['UpsilonDstar']['upsilon2S_frac']) 
    alpha_CB = ROOT.RooRealVar("alpha_CB", "alpha CB Upsilon", upsilon_fit['alpha_CB']['value'])
    n_CB = ROOT.RooRealVar("n_CB", "n CB Upsilon", upsilon_fit['n_CB']['value'])
    """ alpha_CB = ROOT.RooRealVar("alpha_CB", "alpha CB Upsilon", *config['Upsilon']['alpha_CB'])
    n_CB = ROOT.RooRealVar("n_CB", "n CB Upsilon", *config['Upsilon']['n_CB']) """

    signal1S = ROOT.RooCBShape("signal1S", "Upsilon(1S) signal PDF", upsilon_mass, mean1S, sigma1S, alpha_CB, n_CB)
    signal2S = ROOT.RooCBShape("signal2S", "Upsilon(2S) signal PDF", upsilon_mass, mean2S, sigma2S, alpha_CB, n_CB)
    signal3S = ROOT.RooCBShape("signal3S", "Upsilon(3S) signal PDF", upsilon_mass, mean3S, sigma3S, alpha_CB, n_CB)
    upsilon_signal = ROOT.RooAddPdf("upsilon_signal", "Upsilon Signal PDF", ROOT.RooArgList(signal1S, signal2S, signal3S), 
                                    ROOT.RooArgList(upsilon1S_frac, upsilon2S_frac), ROOT.kTRUE)

    # Upsilon Bkg
    bkg_a1  = ROOT.RooRealVar("bkg_a1", "bkg_a1", *config['UpsilonDstar']['bkg_a1'])
    bkg_a2  = ROOT.RooRealVar("bkg_a2", "bkg_a2", *config['UpsilonDstar']['bkg_a2'])

    upsilon_bkg = ROOT.RooChebychev("upsilon_bkg", "Upsilon Background PDF", upsilon_mass, ROOT.RooArgList(bkg_a1, bkg_a2))

    # Dstar Signal Double Gaussian - same mean
    dstar_mean = ROOT.RooRealVar("dstar_mean", "Dstar Gaussian Mean", *config['UpsilonDstar']['dstar_mean'])
    dstar_sigma1 = ROOT.RooRealVar("dstar_sigma1", "Dstar Gaussian 1 Sigma", *config['UpsilonDstar']['dstar_sigma1'])
    dstar_sigma2 = ROOT.RooRealVar("dstar_sigma2", "Dstar Gaussian 2 Sigma", *config['UpsilonDstar']['dstar_sigma2'])
    gauss1_frac = ROOT.RooRealVar("gauss1_frac", "Dstar Gaussian Fraction", *config['UpsilonDstar']['gauss1_frac'])

    g1 = ROOT.RooGaussian("g1", "Dstar Gaussian 1", dstar_deltamr, dstar_mean, dstar_sigma1)
    g2 = ROOT.RooGaussian("g2", "Dstar Gaussian 2", dstar_deltamr, dstar_mean, dstar_sigma2)

    # Dstar Background = Phenomenological Threshold Function 
    p0 = ROOT.RooRealVar("p0","", *config['UpsilonDstar']['p0'])
    p1 = ROOT.RooRealVar('p1',"", *config['UpsilonDstar']['p1'])
    p2 = ROOT.RooRealVar('p2',"", *config['UpsilonDstar']['p2'])

    dstar_signal = ROOT.RooAddPdf("dstar_model", "Dstar Model", ROOT.RooArgList(g1, g2),
                                ROOT.RooArgList(gauss1_frac), ROOT.kTRUE)

    dstar_bkg = ROOT.RooGenericPdf("dstar_bkg","Dstar Background PDF","(1 - exp(-(@0 -0.13957)/@1)) * (@0/0.13957)**@2 + @3 * (@0/0.13957 - 1)",
                                ROOT.RooArgList(dstar_deltamr,p0,p1,p2))

    signal = ROOT.RooProdPdf("signal", "Signal of 2D model", ROOT.RooArgList(upsilon_signal, dstar_signal))
    bkg1 = ROOT.RooProdPdf("bkg1", "Bkg1 of 2D model", ROOT.RooArgList(upsilon_signal, dstar_bkg))
    bkg2 = ROOT.RooProdPdf("bkg2", "Bkg2 of 2D model", ROOT.RooArgList(upsilon_bkg, dstar_signal))
    bkg3 = ROOT.RooProdPdf("bkg3", "Bkg3 of 2D model", ROOT.RooArgList(upsilon_bkg, dstar_bkg))

    signal_frac = ROOT.RooRealVar("signal_frac", "signal fraction", *config['UpsilonDstar']['signal_frac'])
    bkg1_frac = ROOT.RooRealVar("bkg1_frac", "bkg1 fraction", *config['UpsilonDstar']['bkg1_frac'])
    bkg2_frac = ROOT.RooRealVar("bkg2_frac", "bkg2 fraction", *config['UpsilonDstar']['bkg2_frac'])

    model2D = ROOT.RooAddPdf("model2D", "2D Model Upsilon + Dstar", ROOT.RooArgList(signal, bkg1, bkg2, bkg3),
                                ROOT.RooArgList(signal_frac, bkg1_frac, bkg2_frac), ROOT.kTRUE)
    result = model2D.fitTo(data, ROOT.RooFit.BatchMode("cpu"), ROOT.RooFit.Save())

    print("Fit DONE. Saving workspace and params")
    getattr(wspace, "import")(data)
    getattr(wspace, "import")(model2D)
    getattr(wspace, "import")(result)

    wspace.writeToFile(save_path + "/UpsilonDstar_fit.root")
    save_fit_params(save_path, result, wspace, "UpsilonDstar")

def plot_fit(config, year, channel):
    # Get Lumi files
    with open("config/lumi.yaml", 'r') as f:
        lumi = yaml.load(f, Loader=yaml.FullLoader)[year]
    with open("config/skim_trigger.yaml", 'r') as f:
        trigger = yaml.load(f, Loader=yaml.FullLoader)['trigger'][year]

    processed_lumi = 0
    for era in lumi:
        processed_lumi += lumi[era][trigger]

    path = config['path'][year][0][:config['path'][year][0].rfind('/')]
    if channel == "Upsilon":
        if not os.path.exists(path + "/Upsilon_fit.root"): 
            print(f"File {path + '/Upsilon_fit.root'} does not exist")

        with open(f'{path}/{channel}_fit_params.yaml', 'r') as f:
            fit_params = yaml.load(f, Loader=yaml.FullLoader)
        
        n_params = fit_params['n_params']

        f = ROOT.TFile.Open(path + "/Upsilon_fit.root")
        wspace = f.Get(f"upsilon_fit_{year}")

        mass = wspace.var("mass")
        data = wspace.data("data")
        model = wspace.pdf("model")
        result = wspace.obj("fitresult_model_data")

        frame = mass.frame(ROOT.RooFit.Title("Dimuon Invariant mass"))
        data.plotOn(frame, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
        model.plotOn(frame, ROOT.RooFit.Name("#Upsilon(1S) fit"), ROOT.RooFit.Components("signal1S"), ROOT.RooFit.LineColor(colors[1]))
        model.plotOn(frame, ROOT.RooFit.Name("#Upsilon(2S) fit"), ROOT.RooFit.Components("signal2S"), ROOT.RooFit.LineColor(colors[2]))
        model.plotOn(frame, ROOT.RooFit.Name("#Upsilon(3S) fit"), ROOT.RooFit.Components("signal3S"), ROOT.RooFit.LineColor(colors[3]))
        model.plotOn(frame, ROOT.RooFit.Name("Combinatorial fit"), ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineColor(colors[4]), ROOT.RooFit.LineStyle(2))
        model.plotOn(frame, ROOT.RooFit.Name("Model"), ROOT.RooFit.LineColor(colors[0]))
        
        n_evt = data.sumEntries()
        n_params = len(result.floatParsFinal())
        chi2 = frame.chiSquare("Model", "Data", n_params)

        leg = ROOT.TLegend(0.7, 0.7, 0.93, 0.92)
        leg.AddEntry(frame.findObject("Data"), "Data", "LEP")
        leg.AddEntry(frame.findObject("Model"), "Fit", "L")
        leg.AddEntry(frame.findObject("#Upsilon(1S) fit"), "#Upsilon(1S) Fit", "L")
        leg.AddEntry(frame.findObject("#Upsilon(2S) fit"), "#Upsilon(2S) Fit", "L")
        leg.AddEntry(frame.findObject("#Upsilon(3S) fit"), "#Upsilon(3S) Fit", "L")
        leg.AddEntry(frame.findObject("Combinatorial fit"), "Combinatorial Fit", "L")

        canvas = ROOT.TCanvas("c1")

        frame.Draw()
        leg.Draw("same")

        canvas.Draw()
        canvas.SaveAs(f"plots/fit/{channel.lower()}_{year}_fit.png")

        # Save to file some interesting params
        fit_params['chi2'] = chi2

        with open(f'{path}/{channel}_fit_params.yaml', 'w') as f:
            yaml.dump(fit_params, f)


    if channel == "Dstar":
        if not os.path.exists(path + "/Dstar_fit.root"): 
            print(f"File {path + '/Dstar_fit.root'} does not exist")

        with open(f'{path}/{channel}_fit_params.yaml', 'r') as f:
            fit_params = yaml.load(f, Loader=yaml.FullLoader)
        
        n_params = fit_params['n_params']
        f = ROOT.TFile.Open(path + "/Dstar_fit.root")
        wspace = f.Get(f"dstar_fit_{year}")

        deltamr = wspace.var("deltamr")
        data = wspace.data("data")
        model = wspace.pdf("model")
        result = wspace.obj("fitresult_model_data")

        frame = deltamr.frame(ROOT.RooFit.Title("D* #Delta m fit"))
        data.plotOn(frame, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
        model.plotOn(frame, ROOT.RooFit.Name("D* fit"), ROOT.RooFit.Components("g*"), ROOT.RooFit.LineColor(colors[1]))
        model.plotOn(frame, ROOT.RooFit.Name("Combinatorial fit"), ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineColor(colors[2]), ROOT.RooFit.LineStyle(2))
        model.plotOn(frame, ROOT.RooFit.Name("Model"), ROOT.RooFit.LineColor(colors[0]))

        n_evt = data.sumEntries()
        n_params = len(result.floatParsFinal())
        chi2 = frame.chiSquare("Model", "Data", n_params)

        leg = ROOT.TLegend(0.7, 0.7, 0.93, 0.92)
        leg.AddEntry(frame.findObject("Data"), "Data", "LEP")
        leg.AddEntry(frame.findObject("Model"), "Fit", "L")
        leg.AddEntry(frame.findObject("D* fit"), "D* fit Fit", "L")
        leg.AddEntry(frame.findObject("Combinatorial fit"), "Combinatorial Fit", "L")

        canvas = ROOT.TCanvas("c1")

        frame.Draw()
        leg.Draw("same")

        canvas.Draw()
        canvas.SaveAs(f"plots/fit/{channel.lower()}_{year}_fit.png")

        # Save to file some interesting params
        fit_params['chi2'] = chi2

        with open(f'{path}/{channel}_fit_params.yaml', 'w') as f:
            yaml.dump(fit_params, f)

    if channel == "UpsilonDstar":
        plots_upsilondstar(path, year, processed_lumi)

def plots_upsilondstar(path, year, lumi):
    plot.ModTDRStyle(width=800)
    plot.lumi_13TeV = f"{lumi:.2f} " + "fb^{-1}"
    if not os.path.exists(path + "/UpsilonDstar_fit.root"): 
        print(f"File {path + '/UpsilonDstar_fit.root'} does not exist")
        return

    with open(f'{path}/UpsilonDstar_fit_params.yaml', 'r') as f:
        fit_params = yaml.load(f, Loader=yaml.FullLoader)
    n_params = fit_params['n_params']

    f = ROOT.TFile.Open(path + "/UpsilonDstar_fit.root")
    wspace = f.Get(f"upsilondstar_fit_{year}")

    upsilon_mass = wspace.var("dimu_mass")
    dstar_deltamr = wspace.var("dstar_deltamr")
    data = wspace.data("data")
    model2D = wspace.pdf("model2D")
    result = wspace.obj("fitresult_model_data")
    
    # Canvas for plotting Upsilon projection 
    c1 = ROOT.TCanvas("Upsilon mass")
    frame_upsilon = upsilon_mass.frame(ROOT.RooFit.Title("Dimuon Invariant mass"))
    frame_upsilon.GetXaxis().SetTitle("M(\mu^+\mu^-) [GeV/c^2]")
    
    # Plot the Data
    data.plotOn(frame_upsilon, ROOT.RooFit.Name("Data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

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

    c1.SaveAs(f"plots/fit/fit2D_upsilon_proj_{year}.png")

    ## Upsilon pull distribution
    cp_upsilon = ROOT.TCanvas("Upsilon Pull", "", 0, 0, 800, 250)

    histpull_upsilon = frame_upsilon.pullHist()

    # New frame to draw pull distribution
    frame_pull_upsilon = upsilon_mass.frame(ROOT.RooFit.Title("Pull Distribution"))
    frame_pull_upsilon.GetXaxis().SetTitle("")
    frame_pull_upsilon.GetXaxis().SetTitleSize(0.15)    
    frame_pull_upsilon.GetXaxis().SetLabelSize(0.12) 
    frame_pull_upsilon.GetYaxis().SetTitleSize(0.15)    
    frame_pull_upsilon.GetYaxis().SetLabelSize(0.12) 

    # Add the distribution to the frame
    frame_pull_upsilon.addPlotable(histpull_upsilon, "P")

    frame_pull_upsilon.Draw()
    line = ROOT.TLine(frame_pull_upsilon.GetXaxis().GetXmax(), 0, frame_pull_upsilon.GetXaxis().GetXmin(), 0)
    line.SetLineStyle(2)
    line.Draw()

    cp_upsilon.SaveAs(f"plots/fit/fit2D_upsilon_pull_{year}.png")
    
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

    c2.SaveAs(f"plots/fit/fit2D_dstar_proj_{year}.png")

    ## Dstar pull distribution

    # Canvas for Dstar pull 
    cpd = ROOT.TCanvas("Dstar Pull", "", 0, 0, 800, 250)

    histpull_dstar = frame_dstar.pullHist()

    # New frame to draw pull distribution
    frame_pull_dstar = dstar_deltamr.frame(ROOT.RooFit.Title("Pull Distribution"))
    frame_pull_dstar.GetXaxis().SetTitle("")
    frame_pull_dstar.GetXaxis().SetTitleSize(0.15)    
    frame_pull_dstar.GetXaxis().SetLabelSize(0.12)    
    frame_pull_dstar.GetYaxis().SetTitleSize(0.15)    
    frame_pull_dstar.GetYaxis().SetLabelSize(0.12)    

    # Add the distribution to the frame
    frame_pull_dstar.addPlotable(histpull_dstar, "P")

    frame_pull_dstar.Draw()
    line = ROOT.TLine(frame_pull_dstar.GetXaxis().GetXmax(), 0, frame_pull_dstar.GetXaxis().GetXmin(), 0)
    line.SetLineStyle(2)
    line.Draw()

    cpd.SaveAs(f"plots/fit/fit2D_dstar_pull_{year}.png")

    #################################### 2D fitting
    plot.ModTDRStyle(width=1000, height=1000)
    c3 = ROOT.TCanvas("2D fit", '')
    c3.Draw()

    ph2 = dstar_deltamr.createHistogram("", upsilon_mass)
    ph2.SetStats(False)
    ph2.SetTitle("2D fit - #Upsilon D*")
    ph2.GetXaxis().SetLabelSize(0.02)
    ph2.GetYaxis().SetLabelSize(0.02)
    ph2.GetZaxis().SetLabelSize(0.02)
    ph2.GetXaxis().SetTitle('')
    ph2.GetYaxis().SetTitle('')
    model2D.fillHistogram(ph2, ROOT.RooArgList(dstar_deltamr, upsilon_mass))
    ph2.Draw("SURF")

    # D* axis
    mumu = ROOT.TLatex()
    mumu.SetNDC()
    mumu.SetTextFont(43)
    mumu.SetTextSize(20)
    mumu.SetTextAlign(31)
    mumu.SetTextAngle(12)
    mumu.DrawLatex(0.8, 0.15, "M(k\pi\pi_s) - M(k\pi) [GeV/c^2]")

    # Jpsi axis
    mumu = ROOT.TLatex()
    mumu.SetNDC()
    mumu.SetTextFont(43)
    mumu.SetTextSize(20)
    mumu.SetTextAlign(31)
    mumu.SetTextAngle(-33)
    mumu.DrawLatex(0.22, 0.17, "M(\mu^+\mu^-) [GeV]")

    # Candidates
    mumu = ROOT.TLatex()
    mumu.SetNDC()
    mumu.SetTextFont(43)
    mumu.SetTextSize(20)
    mumu.SetTextAlign(31)
    mumu.SetTextAngle(90)
    mumu.DrawLatex(0.027, 0.68, "Candidates")

    right = ROOT.TLatex()
    right.SetNDC()
    right.SetTextFont(43)
    right.SetTextSize(30)
    right.SetTextAlign(13)
    right.SetTextAngle(13)
    right.DrawLatex(.10, 0.84,"#bf{CMS} #scale[0.7]{#it{Preliminary}}")

    right.SetTextSize(25)
    right.SetTextAngle(12)
    right.DrawLatex(.55, .95, f"{lumi:.2f}" + " fb^{-1}")

    ROOT.gPad.Update()

    """ PaletteAxis = ph2.GetListOfFunctions().FindObject("palette")
    PaletteAxis.SetX1NDC(1.4)
    PaletteAxis.SetX2NDC(0.95) """
    
    c3.Update()
    c3.SaveAs(f"plots/fit/fit2D_pdf_{year}.png")

    #################################### Mass correlation plot
    with open('config/fit.yaml') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    chain = ROOT.TChain("UpsilonDstar")
    for f in config['path'][year]:
        chain.Add(f"{f}/UpsilonDstar.root")

    hist = ROOT.TH2D('h1', 'Associated Particles; ;; ', 30, 8.7, 11.2, 30, 0.141, 0.16)
    hist.SetTitle("Mass correlation - #Upsilon D*")
    hist.GetXaxis().SetLabelSize(0.02)
    hist.GetYaxis().SetLabelSize(0.02)
    hist.GetZaxis().SetLabelSize(0.02)
    chain.Draw("dstar_deltamr:dimu_mass>>h1")

    # Canvas Definition
    c4 = ROOT.TCanvas("Mass correlation", '')
    c4.Draw()

    hist.SetStats(False)
    hist.Draw("lego2z")

    # Jpsi axis
    mumu = ROOT.TLatex()
    mumu.SetNDC()
    mumu.SetTextFont(43)
    mumu.SetTextSize(25)
    mumu.SetTextAlign(31)
    mumu.SetTextAngle(11)
    mumu.DrawLatex(0.85, 0.14, "M(\mu^+\mu^-) [GeV/c^2]")

    # D* axis
    mumu = ROOT.TLatex()
    mumu.SetNDC()
    mumu.SetTextFont(43)
    mumu.SetTextSize(25)
    mumu.SetTextAlign(31)
    mumu.SetTextAngle(-35)
    mumu.DrawLatex(0.30, 0.10, "M(k\pi\pi_s) - M(k\pi) [GeV/c^2]")

    # Candidates
    mumu = ROOT.TLatex()
    mumu.SetNDC()
    mumu.SetTextFont(43)
    mumu.SetTextSize(25)
    mumu.SetTextAlign(31)
    mumu.SetTextAngle(90)
    mumu.DrawLatex(0.03, 0.65, "Candidates")

    right = ROOT.TLatex()
    right.SetNDC()
    right.SetTextFont(43)
    right.SetTextSize(35)
    right.SetTextAlign(13)
    right.SetTextAngle(13)
    right.DrawLatex(.10, 0.84,"#bf{CMS} #scale[0.7]{#it{Preliminary}}")

    right.SetTextSize(25)
    right.SetTextAngle(13)
    right.DrawLatex(.54, .94, f"{lumi:.2f}" + " fb^{-1}")

    ROOT.gPad.Update()

    PaletteAxis = hist.GetListOfFunctions().FindObject("palette")
    PaletteAxis.SetX1NDC(1.4)
    PaletteAxis.SetX2NDC(0.95)

    c4.Update()
    c4.SaveAs(f"plots/fit/fit2D_mass_corr_{year}.png")

    # Save to file some interesting params
    fit_params['chi2_upsilon'] = chi2_upsilon
    fit_params['chi2_dstar'] = chi2_dstar

    with open(f'{path}/UpsilonDstar_fit_params.yaml', 'w') as f:
        yaml.dump(fit_params, f)