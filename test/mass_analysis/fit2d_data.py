import ROOT

# Suppress the opening of many Canvas's 
ROOT.gROOT.SetBatch(True) 

#input
file_path = "../histos/H/FullDataset_HRho2024_ORTrigger_BDTCat1_Data.root"
tree_name = "tree_output"
infile = ROOT.TFile.Open(file_path)
tree = infile.Get(tree_name)

#Vars
bosonMass = ROOT.RooRealVar("bosonMass", "m(H) [GeV]", 110., 140.)
mesonMass = ROOT.RooRealVar("mesonMass", "m(#rho) [GeV]", 0.55, 0.95)

# Blind region definition
bosonMass.setRange("LowSideband",110.,120.)
bosonMass.setRange("HighSideband",130.,140.)
bosonMass.setRange("full",110.,140.)


# Create RooDataSet (Unbinned)
data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(bosonMass, mesonMass), ROOT.RooFit.Import(tree))
print(f"Number of entries loaded in dataset: {data.numEntries()}")

#Give the blind range - just for plotting ------------------------------------------------------------------------------------------------------------------------
data_blinded = data.reduce("bosonMass < 120. || bosonMass > 130.")
print(f"Entries after blinding: {data_blinded.numEntries()}")

# Higgs pdf - chebychev
#a_H = ROOT.RooRealVar("a_H_chebychev","a_H",-1.,-3.,0.)
#b_H = ROOT.RooRealVar("b_H_chebychev","b_H",0.3,-3.,3.)
#c_H = ROOT.RooRealVar("c_H_chebychev","c_H",-0.01,-1.,1.)
a_H = ROOT.RooRealVar("a_H_chebychev", "a_H", -0.430604)
a_H.setConstant(True)
b_H = ROOT.RooRealVar("b_H_chebychev","b_H",0.0459796)
b_H.setConstant(True)

pdf_H = ROOT.RooChebychev("pdf_H", "Chebychev Higgs", bosonMass, ROOT.RooArgList(a_H, b_H))#, c_H))

'''
#Meson pdf without background
# Meson pdf - voigtian+background

#mass_min = 0.55
#mass_max = 1.
#mean_central ,  mean_min,  mean_max = 0.77, mass_min, mass_max
#sigma_central, sigma_min, sigma_max = 0.022, 0.0001, 0.1
#width_central, width_min, width_max = 0.047, 0.0, 0.2
#mean  = ROOT.RooRealVar("mean","The mean of the gaussian pdf", mean_central, mean_min, mean_max)
#sigma = ROOT.RooRealVar("sigma","The sigma of the gaussian pdf", sigma_central, sigma_min, sigma_max)
#width = ROOT.RooRealVar("width","The width of the gaussian pdf", width_central, width_min, width_max)

#pdf_rho = ROOT.RooVoigtian("signalPDF","Voigtian pdf",mesonMass,mean,width,sigma)
'''
# Meson pdf - voigtian
mass_min = 0.55
mass_max = 0.95
mean_central, mean_min, mean_max = 0.77, 0.74, 0.79
sigma_central, sigma_min, sigma_max = 0.055, 0.01, 0.2
width_central, width_min, width_max = 0.047, 0.0, 0.2

#mean  = ROOT.RooRealVar("mean","The mean of the gaussian pdf", mean_central, mean_min, mean_max)
#sigma = ROOT.RooRealVar("sigma","The sigma of the gaussian pdf", sigma_central, sigma_min, sigma_max)
#width = ROOT.RooRealVar("width","The width of the gaussian pdf", width_central, width_min, width_max)
mean = ROOT.RooRealVar("mean", "The mean of the gaussian pdf", 0.751561)
mean.setConstant(True)
sigma = ROOT.RooRealVar("sigma","The sigma of the gaussian pdf", 0.054935)
sigma.setConstant(True)

signalPDF = ROOT.RooGaussian("signalPDF","Gauss pdf",mesonMass,mean,sigma)
#signalPDF = ROOT.RooVoigtian("signalPDF","Voigtian pdf",mass,mean,width,sigma)

# Construct background pdf
#a1_central, a1_min, a1_max = 0.3,-0.5,0.5
#a2_central, a2_min, a2_max = -0.1,-2.,2.
#a3_central, a3_min, a3_max = -0.02,-2.,2.

#a1_off = ROOT.RooRealVar("a1_off","The a1 of background", a1_central, a1_min, a1_max)
#a2_off = ROOT.RooRealVar("a2_off","The a2 of background", a2_central, a2_min, a2_max)
#a3_off = ROOT.RooRealVar("a3_off","The a3 of background", a3_central, a3_min, a3_max)
a1_off = ROOT.RooRealVar("a1_off","The a1 of background", 0.006057)
a1_off.setConstant(True)
a2_off = ROOT.RooRealVar("a2_off","The a2 of background",-0.023487)
a2_off.setConstant(True)
a3_off = ROOT.RooRealVar("a3_off","The a3 of background", -0.002644)
a3_off.setConstant(True)

backgroundPDF = ROOT.RooChebychev("backgroundPDF","The background PDF",mesonMass,ROOT.RooArgList(a1_off,a2_off,a3_off))

# N events
# ---------------------------------------------------------------
#data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(mesonMass), ROOT.RooFit.Import(tree))
nData = data.numEntries()
#print(f"Number of entries loaded in dataset: {nData}")

nsig = ROOT.RooRealVar("nsig", "signal yield",0.2*nData, 0.,nData)
nbkg = ROOT.RooRealVar("nbkg", "background yield",0.8*nData,0.,nData)

# Create the model
# ---------------------------------------------------------------
pdf_rho = ROOT.RooAddPdf("model", "Signal and Background PDF", ROOT.RooArgList(signalPDF, backgroundPDF), ROOT.RooArgList(nsig, nbkg))



# 2D PDF
pdf_2D = ROOT.RooProdPdf("pdf_2D", "2D Signal PDF (Higgs x Meson)", ROOT.RooArgList(pdf_H, pdf_rho))

mask = ROOT.RooFormulaVar("mask", "(bosonMass < 120 || bosonMass > 130) ? 1.0 : 0.0", ROOT.RooArgList(bosonMass))

pdf_2D_masked = ROOT.RooGenericPdf("pdf_2D_masked", "@0 * @1", ROOT.RooArgList(pdf_2D, mask))


# 2D Unbinned Fit
###attention -- fitting on the whole dataset at the moment. try to implement the range (don't use data_blinded)
fit_result = pdf_2D.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(1)) #ROOT.RooFit.Range("LowSideband,HighSideband")

print("\n" + "="*60)
print("FIT RESULTS:")
fit_result.Print()
print("="*60 + "\n")

# Correlation
print("Evaluating correlation between observables...")
h2_data = ROOT.TH2F("h2_data", "Data Correlation", 60, 105, 165, 50, 0.65, 0.90)
data.fillHistogram(h2_data, ROOT.RooArgList(bosonMass, mesonMass))
correlation_data = h2_data.GetCorrelationFactor()
print(f"--> Linear correlation (Pearson) in DATA: {correlation_data:.4f}")

canvas = ROOT.TCanvas("canvas", "Fit Projections & Correlation", 1800, 600)
canvas.Divide(3, 1)

# Higgs Mass proj
canvas.cd(1)
frame_H = bosonMass.frame(ROOT.RooFit.Title("Higgs Mass Projection"), ROOT.RooFit.Bins(60))
data_blinded.plotOn(frame_H, ROOT.RooFit.Name("data_H"), ROOT.RooFit.CutRange("LowSideband,HighSideband"))
pdf_2D.plotOn(frame_H, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("pdf_H"), ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.NormRange("LowSideband,HighSideband"))
chi2_H = frame_H.chiSquare("pdf_H", "data_H")
frame_H.Draw()

leg_H = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
leg_H.SetBorderSize(0)
leg_H.SetFillStyle(0)
leg_H.SetTextSize(0.035)
leg_H.AddEntry("data_H", "Signal", "ep")
leg_H.AddEntry("pdf_H", "Chebychev fit", "l")
leg_H.AddEntry("", f"#chi^{{2}} / ndf = {chi2_H:.2f}", "")
leg_H.Draw()

# Rho Meson proj
canvas.cd(2)
frame_rho = mesonMass.frame(ROOT.RooFit.Title("Meson Mass Projection"))
data.plotOn(frame_rho, ROOT.RooFit.Name("data_rho"))
pdf_2D.plotOn(frame_rho, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf_rho"))
chi2_rho = frame_rho.chiSquare("pdf_rho", "data_rho")
frame_rho.Draw()

leg_rho = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
leg_rho.SetBorderSize(0)
leg_rho.SetFillStyle(0)
leg_rho.SetTextSize(0.035)
leg_rho.AddEntry("data_rho", "Signal", "ep")
leg_rho.AddEntry("pdf_rho", "Voigtian Fit", "l")
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
output_path_pdf = "/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/mass_analysis/fit_projections_data_2D.pdf"
output_path_png = "/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/mass_analysis/fit_projections_data_2D.png"
canvas.SaveAs(output_path_pdf)
canvas.SaveAs(output_path_png)
print("Plots saved successfully with 2D Heatmap.")

