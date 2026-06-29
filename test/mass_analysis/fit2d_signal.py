import ROOT

# Import the Double Crystal Ball PDF 
ROOT.gROOT.ProcessLineSync(".L ./mass_analysis/dCB/RooDoubleCBFast.cc+")

# Suppress the opening of many Canvas's 
ROOT.gROOT.SetBatch(True) 

#input
file_path = "./histos/H/SR_HRho_OldTrigger_BDT_Signal.root"
tree_name = "tree_output"
infile = ROOT.TFile.Open(file_path)
tree = infile.Get(tree_name)

#Vars
bosonMass = ROOT.RooRealVar("bosonMass", "m(H) [GeV]", 100, 170)
mesonMass = ROOT.RooRealVar("mesonMass", "m(#rho) [GeV]", 0.62, 0.92)

# Create RooDataSet (Unbinned)
data = ROOT.RooDataSet("data", "Signal", ROOT.RooArgSet(bosonMass, mesonMass), ROOT.RooFit.Import(tree))
print(f"Number of entries loaded in dataset: {data.numEntries()}")

# Higgs pdf
mean_H     = ROOT.RooRealVar("mean_H", "Higgs Mean", 125.0, 120.0, 130.0)
sigma_H    = ROOT.RooRealVar("sigma_H", "Higgs Sigma", 3.0, 0.5, 8.0)
alphaL_H   = ROOT.RooRealVar("alphaL_H", "Higgs Alpha Left", 1.5, 0.1, 15.0)
nL_H       = ROOT.RooRealVar("nL_H", "Higgs n Left", 2.0, 0.1, 50.0)
alphaR_H   = ROOT.RooRealVar("alphaR_H", "Higgs Alpha Right", 1.5, 0.1, 15.0)
nR_H       = ROOT.RooRealVar("nR_H", "Higgs n Right", 2.0, 0.1, 50.0)

pdf_H = ROOT.RooDoubleCBFast("pdf_H", "DSCB Higgs", bosonMass, mean_H, sigma_H, alphaL_H, nL_H, alphaR_H, nR_H)

# Meson pdf
mean_rho   = ROOT.RooRealVar("mean_rho", "Rho Mean", 0.770, 0.74, 0.80)
sigma_rho  = ROOT.RooRealVar("sigma_rho", "Rho Sigma", 0.05, 0.02, 0.12)
alphaL_rho = ROOT.RooRealVar("alphaL_rho", "Rho Alpha Left", 1.0, 0.05, 20.0)
nL_rho     = ROOT.RooRealVar("nL_rho", "Rho n Left", 1.5, 0.05, 100.0)
alphaR_rho = ROOT.RooRealVar("alphaR_rho", "Rho Alpha Right", 1.0, 0.05, 20.0)
nR_rho     = ROOT.RooRealVar("nR_rho", "Rho n Right", 1.5, 0.05, 100.0)

pdf_rho = ROOT.RooDoubleCBFast("pdf_rho", "DSCB Rho", mesonMass, mean_rho, sigma_rho, alphaL_rho, nL_rho, alphaR_rho, nR_rho)

# 2D PDF
pdf_2D = ROOT.RooProdPdf("pdf_2D", "2D Signal PDF (Higgs x Meson)", ROOT.RooArgList(pdf_H, pdf_rho))

# 2D Unbinned Maximum Likelihood Fit
fit_result = pdf_2D.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(1))

print("\n" + "="*50)
print("FIT RESULTS:")
fit_result.Print()
print("="*50 + "\n")

# Correlation
print("Evaluating correlation between observables...")
h2_data = ROOT.TH2F("h2_data", "Data Correlation", 70, 100, 170, 50, 0.62, 0.92)
data.fillHistogram(h2_data, ROOT.RooArgList(bosonMass, mesonMass))
correlation_data = h2_data.GetCorrelationFactor()
print(f"--> Linear correlation (Pearson) in DATA: {correlation_data:.4f}")

canvas = ROOT.TCanvas("canvas", "Fit Projections & Correlation", 1800, 600)
canvas.Divide(3, 1)

# Higgs Mass proj
canvas.cd(1)
frame_H = bosonMass.frame(ROOT.RooFit.Title("Higgs Mass Projection"))
data.plotOn(frame_H, ROOT.RooFit.Name("data_H"))
pdf_2D.plotOn(frame_H, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("pdf_H"))
chi2_H = frame_H.chiSquare("pdf_H", "data_H", 6)
frame_H.Draw()

leg_H = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
leg_H.SetBorderSize(0)
leg_H.SetFillStyle(0)
leg_H.SetTextSize(0.035)
leg_H.AddEntry("data_H", "Signal", "ep")
leg_H.AddEntry("pdf_H", "DCB fit", "l")
leg_H.AddEntry("", f"#chi^{{2}} / ndf = {chi2_H:.2f}", "")
leg_H.Draw()

# Rho Meson proj
canvas.cd(2)
frame_rho = mesonMass.frame(ROOT.RooFit.Title("Meson Mass Projection"))
data.plotOn(frame_rho, ROOT.RooFit.Name("data_rho"))
pdf_2D.plotOn(frame_rho, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf_rho"))
chi2_rho = frame_rho.chiSquare("pdf_rho", "data_rho", 6)
frame_rho.Draw()

leg_rho = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
leg_rho.SetBorderSize(0)
leg_rho.SetFillStyle(0)
leg_rho.SetTextSize(0.035)
leg_rho.AddEntry("data_rho", "Signal", "ep")
leg_rho.AddEntry("pdf_rho", "DCB Fit", "l")
leg_rho.AddEntry("", f"#chi^{{2}} / ndf = {chi2_rho:.2f}", "")
leg_rho.Draw()

# --- Panel 3: 2D Heatmap (COLZ) ---
canvas.cd(3)
ROOT.gPad.SetRightMargin(0.15)

h2_heatmap = data.createHistogram("h2_heatmap", bosonMass, ROOT.RooFit.Binning(50), ROOT.RooFit.YVar(mesonMass, ROOT.RooFit.Binning(50)))
h2_heatmap.SetTitle("Scatter plot")
h2_heatmap.GetXaxis().SetTitle("m(H) [GeV]")
h2_heatmap.GetYaxis().SetTitle("m(#rho) [GeV]")

ROOT.gStyle.SetPalette(ROOT.kBird) 
h2_heatmap.Draw("COLZ")

leg_corr = ROOT.TLegend(0.18, 0.82, 0.65, 0.88)
leg_corr.SetBorderSize(0)
leg_corr.SetFillColor(ROOT.kWhite)
leg_corr.SetFillStyle(1001)
leg_corr.SetTextSize(0.035)
leg_corr.AddEntry("", f"Pearson Corr. = {correlation_data:.4f}", "")
leg_corr.Draw()

# Save plots
output_path_pdf = "/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/mass_analysis/fit_projections_2D.pdf"
output_path_png = "/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/mass_analysis/fit_projections_2D.png"
canvas.SaveAs(output_path_pdf)
canvas.SaveAs(output_path_png)
print("Plots saved successfully with 2D Heatmap.")