import ROOT
import cmsstyle
import argparse

# Suppress the opening of many Canvas's 
ROOT.gROOT.SetBatch(True) 

#input
file_path = "../../histos/H/FullDataset_HRho2024_ORTrigger_BDTCat0_Data.root"
tree_name = "tree_output"
infile = ROOT.TFile.Open(file_path)
tree = infile.Get(tree_name)

CHANNEL = "Rho"

#CMS-style plotting ---------------------------------------------------------------
cmsstyle.setCMSStyle()

# set the luminosity and the COM energy to show in the canvases
cmsstyle.SetLumi(108.95, run=None)
cmsstyle.SetEnergy(13.6)
iPeriod = 4
iPos = 11

#Vars
bosonMass = ROOT.RooRealVar("bosonMass", "m(H) [GeV]", 110., 140.)
mesonMass = ROOT.RooRealVar("mesonMass", "m(#rho) [GeV]", 0.55, 0.95)

# Blind region definition
xLowRange  = 110.
xHighRange = 140.
xSidLeft = 120.
xSidRight = 130.
bosonMass.setRange("LowSideband",110.,120.)
bosonMass.setRange("HighSideband",130.,140.)
bosonMass.setRange("full",110.,140.)


# Create RooDataSet (Unbinned)
data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(bosonMass, mesonMass), ROOT.RooFit.Import(tree), ROOT.RooFit.Cut("(bosonMass < 120 || bosonMass > 130)"))
print(f"Number of entries loaded in dataset: {data.numEntries()}")


#dataset to be imported in the workspace, without blinding
observed_data = ROOT.RooDataSet("observed_data","observed_data",ROOT.RooArgSet(bosonMass, mesonMass),ROOT.RooFit.Import(tree))
nEntries = observed_data.numEntries() 

#Give the blind range ------------------------------------------------------------------------------------------------------------------------
data_blinded = data.reduce("bosonMass < 120. || bosonMass > 130.")
#print(f"Number of entries loaded in dataset: {data.numEntries()}")
print(f"Entries in Higgs sidebands: {data_blinded.numEntries()}")

# Higgs pdfs - chebychev

#Initialize Chebychev pdfs
a_bkg = ROOT.RooRealVar("a_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","a_bkg",-0.4961)
a_bkg.setConstant(True)
b_bkg = ROOT.RooRealVar("b_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","b_bkg",0.062)
b_bkg.setConstant(True)
c_bkg = ROOT.RooRealVar("c_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","c_bkg",-0.0103)
c_bkg.setConstant(True)

a1_bkg = ROOT.RooRealVar("a1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","a1_bkg",-0.4947)
a1_bkg.setConstant(True)
b1_bkg = ROOT.RooRealVar("b1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","b1_bkg",0.060)
b1_bkg.setConstant(True)
c1_bkg = ROOT.RooRealVar("c1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","c1_bkg",-0.0114)
c1_bkg.setConstant(True)
d1_bkg = ROOT.RooRealVar("d1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","d1_bkg",0.009)
d1_bkg.setConstant(True)

a2_bkg = ROOT.RooRealVar("a2_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","a2_bkg",-0.328)
a2_bkg.setConstant(True)
b2_bkg = ROOT.RooRealVar("b2_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","b2_bkg",-0.001)
b2_bkg.setConstant(True)
c2_bkg = ROOT.RooRealVar("c2_bkg_chebychev_"+CHANNEL+"_GFcat_bdt0","c2_bkg",-0.01,-1.,1.)
#c2_bkg.setConstant(True)

bkgPDF_chebychev  = ROOT.RooChebychev("chebychev_GFcat_bdt0_bkg","bkgPDF",bosonMass,ROOT.RooArgList(a_bkg, b_bkg,c_bkg))
bkgPDF_chebychev_Nplus1 = ROOT.RooChebychev("chebychev_Nplus1_GFcat_bdt0_bkg","bkgPDF",bosonMass,ROOT.RooArgList(a1_bkg,b1_bkg,c1_bkg,d1_bkg))
bkgPDF_chebychev_Nminus1 = ROOT.RooChebychev("chebychev_Nminus1_GFcat_bdt0_bkg","bkgPDF",bosonMass,ROOT.RooArgList(a2_bkg, b2_bkg))

#Initialize an exponential pdf
e_bkg = ROOT.RooRealVar("e_bkg_exponential_"+CHANNEL+"_GFcat_bdt0","e_bkg",-0.064 )
e_bkg.setConstant(True)
f_bkg = ROOT.RooRealVar("f_bkg_exponential_"+CHANNEL+"_GFcat_bdt0","f_bkg",0.3)
f_bkg.setConstant(True)
bkgPDF_exponential = ROOT.RooExponential("exponential_GFcat_bdt0_bkg","bkgPDF",bosonMass,e_bkg)

bkgPDF_chebyplusexp = ROOT.RooAddPdf("chebyplusexp_GFcat_bdt0_bkg","bkgPDF",ROOT.RooArgSet(bkgPDF_exponential,bkgPDF_chebychev_Nminus1),ROOT.RooArgList(f_bkg))

# Meson pdf - voigtian
mean = ROOT.RooRealVar("mean", "The mean of the gaussian pdf", 0.751)
mean.setConstant(True)
sigma = ROOT.RooRealVar("sigma","The sigma of the gaussian pdf", 0.053)
sigma.setConstant(True)

signalPDF = ROOT.RooGaussian("signalPDF","Gauss pdf",mesonMass,mean,sigma)
#signalPDF = ROOT.RooVoigtian("signalPDF","Voigtian pdf",mass,mean,width,sigma)

# Construct background pdf
a1_off = ROOT.RooRealVar("a1_off","The a1 of background", 0.015)
a1_off.setConstant(True)
a2_off = ROOT.RooRealVar("a2_off","The a2 of background", -0.035)
a2_off.setConstant(True)
a3_off = ROOT.RooRealVar("a3_off","The a3 of background", -0.0009)
a3_off.setConstant(True)

backgroundPDF = ROOT.RooChebychev("backgroundPDF","The background PDF",mesonMass,ROOT.RooArgList(a1_off,a2_off,a3_off))

# N events
# ---------------------------------------------------------------
#data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(mesonMass), ROOT.RooFit.Import(tree))
nData = data.numEntries()
#print(f"Number of entries loaded in dataset: {nData}")

nsig = ROOT.RooRealVar("nsig", "signal yield",0.2*nData,0.,nData)
nbkg = ROOT.RooRealVar("nbkg", "background yield",0.8*nData,0.,nData)

# Create the model for the meson distribution: voigtian signal over polinomial background
# ---------------------------------------------------------------
pdf_rho = ROOT.RooAddPdf("model", "Signal and Background PDF", ROOT.RooArgList(signalPDF, backgroundPDF), ROOT.RooArgList(nsig, nbkg))

# 2D PDF
pdf_2D_cheby = ROOT.RooProdPdf("pdf_2D_cheby", "2D Signal PDF (Higgs x Meson)", ROOT.RooArgList(bkgPDF_chebychev, pdf_rho))
pdf_2D_cheby_Nplus1 = ROOT.RooProdPdf("pdf_2D_cheby_Nplus1", "2D Signal PDF (Higgs x Meson)", ROOT.RooArgList(bkgPDF_chebychev_Nplus1, pdf_rho))
pdf_2D_chebyplusexp = ROOT.RooProdPdf("pdf_2D_chebyplusexp", "2D Signal PDF (Higgs x Meson)", ROOT.RooArgList(bkgPDF_chebyplusexp, pdf_rho))


#Do the fit ------------------------------------------------------------------------------------------------------------------------------
fitResult_chebychev = pdf_2D_cheby.fitTo(data_blinded,ROOT.RooFit.Save(),ROOT.RooFit.PrintLevel(1))
fitResult_chebychev_Nplus1 = pdf_2D_cheby_Nplus1.fitTo(data_blinded,ROOT.RooFit.Save(),ROOT.RooFit.PrintLevel(1))
fitResult_exponential = pdf_2D_chebyplusexp.fitTo(data_blinded,ROOT.RooFit.Save(),ROOT.RooFit.PrintLevel(1))

#print("\n" + "="*60)
#print("FIT RESULTS:")
#fit_result.Print()
#print("="*60 + "\n")

# Correlation
print("Evaluating correlation between observables...")
h2_data = ROOT.TH2F("h2_data", "Data Correlation", 60, 105, 165, 50, 0.65, 0.90)
data.fillHistogram(h2_data, ROOT.RooArgList(bosonMass, mesonMass))
correlation_data = h2_data.GetCorrelationFactor()
print(f"--> Linear correlation (Pearson) in DATA: {correlation_data:.4f}")

#Plot ------------------------------------------------------------------------------------------------------------------------
'''
def draw_parameters(parameters, x1=0.68, y1=0.55, x2=0.94, y2=0.90):

    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillStyle(0)
    box.SetBorderSize(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(0.025)

    for par in parameters:
        if par.isConstant():
            box.AddText(f"{par.GetName()} = {par.getVal():.4g}")
        else:
            box.AddText(
                f"{par.GetName()} = {par.getVal():.4g} #pm {par.getError():.2g}"
            )

    box.Draw()
    return box
'''


canvas_chebychev = ROOT.TCanvas("canvas", "Fit Projections & Correlation", 1800, 600)
canvas_chebychev.Divide(3, 1)

for i in range(1, 4):
    canvas_chebychev.cd(i)
    ROOT.gPad.SetLeftMargin(0.2)

# Higgs Mass proj
canvas_chebychev.cd(1)

#Chebychev frame
xframe_chebychev = bosonMass.frame(30)
data_blinded.plotOn(xframe_chebychev, ROOT.RooFit.CutRange("LowSideband,HighSideband"))
pdf_2D_cheby.plotOn(xframe_chebychev, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.NormRange("LowSideband,HighSideband"), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(1))
xframe_chebychev.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
xframe_chebychev.GetXaxis().SetTitle("m_{M,#gamma} [GeV]")
xframe_chebychev.SetMaximum(1.3*xframe_chebychev.GetMaximum())
xframe_chebychev.GetYaxis().SetTitleOffset(1.8)
xframe_chebychev.Draw()
'''
chi2_H = frame_H.chiSquare("pdf_H", "data_H")
binWidth_H = frame_H.getFitRangeBinW()
frame_H.GetYaxis().SetTitle(f"Events / {binWidth_H:.1f} GeV")
frame_H.GetXaxis().SetTitle("m_{Trk^{+}Trk^{-}#gamma} [GeV]")
frame_H.SetMaximum(frame_H.GetMaximum() * 1.30)
frame_H.Draw()
'''
#Calculate Chi square and parameters 
nParam_cheby = fitResult_chebychev.floatParsFinal().getSize()
print('n param cheby = ',nParam_cheby)
chi2_cheby = xframe_chebychev.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xSidLeft-xLowRange+xHighRange-xSidRight) - nParam_cheby
cut_chi2_cheby = "{:.2f}".format(chi2_cheby*my_ndof) #Crop the chi2 to 2 decimal digits
print('Chi square cheby = ',chi2_cheby)
print('')

leg1 = ROOT.TLegend(0.5,0.22,0.72,0.90) #right positioning
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.AddEntry(cut_chi2_cheby,"#chi^{2}/ndof = " + cut_chi2_cheby + " / " + str(my_ndof),"brNDC")
leg1.Draw()

#box_cheby = draw_parameters([a_bkg, b_bkg, c_bkg],x1=0.42, y1=0.70, x2=0.90, y2=0.90)

# Rho Meson proj
canvas_chebychev.cd(2)
frame_rho_cheby = mesonMass.frame(30)
data_blinded.plotOn(frame_rho_cheby)
pdf_2D_cheby.plotOn(frame_rho_cheby, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(1), ROOT.RooFit.Name("pdf_rho"))
frame_rho_cheby.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
frame_rho_cheby.GetXaxis().SetTitle("m_{M} [GeV]")
frame_rho_cheby.SetMaximum(1.3*frame_rho_cheby.GetMaximum())
frame_rho_cheby.GetYaxis().SetTitleOffset(1.8)
frame_rho_cheby.Draw()

'''
chi2_rho = frame_rho_cheby.chiSquare("pdf_rho", "data_rho")
binWidth_rho = frame_rho_cheby.getFitRangeBinW()
frame_rho_cheby.GetYaxis().SetTitle(f"Events / {binWidth_rho:.3f} GeV")
frame_rho_cheby.GetXaxis().SetTitle("m_{Trk^{+}Trk^{-}} [GeV]")
frame_rho_cheby.SetMaximum(frame_rho_cheby.GetMaximum() * 1.20)
frame_rho_cheby.Draw()
'''

#Calculate Chi square and parameters 
nParam_voig = fitResult_chebychev.floatParsFinal().getSize()
print('n param voig = ',nParam_voig)
chi2_voig = frame_rho_cheby.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof_rho = 30 - nParam_voig
cut_chi2_voig = "{:.2f}".format(chi2_voig*my_ndof_rho) #Crop the chi2 to 2 decimal digits
print('Chi square voig = ',chi2_voig)
print('')
leg_rho = ROOT.TLegend(0.5,0.22,0.72,0.90) #right positioning
leg_rho.SetHeader(" ")
leg_rho.SetNColumns(1)
leg_rho.SetFillColorAlpha(0,0.)
leg_rho.SetBorderSize(0)
leg_rho.SetLineColor(1)
leg_rho.SetLineStyle(1)
leg_rho.SetLineWidth(1)
leg_rho.SetFillStyle(1001)
leg_rho.AddEntry(cut_chi2_voig,"#chi^{2}/ndof = " + cut_chi2_voig + " / " + str(my_ndof),"brNDC")
leg_rho.Draw()

#box_rho = draw_parameters([nsig, nbkg, mean, sigma, a1_off, a2_off, a3_off],x1=0.68, y1=0.55, x2=0.94, y2=0.90)

# --- Panel 3: 2D Heatmap (COLZ) ---
canvas_chebychev.cd(3)
#ROOT.gPad.SetRightMargin(0.)

h2_heatmap = data.createHistogram("h2_heatmap", bosonMass, ROOT.RooFit.Binning(50), ROOT.RooFit.YVar(mesonMass, ROOT.RooFit.Binning(50)))
h2_heatmap.SetTitle("Scatter plot")
h2_heatmap.GetXaxis().SetTitle("m_{M#gamma} [GeV]")
h2_heatmap.GetYaxis().SetTitle("m_{M} [GeV]")

ROOT.gStyle.SetPalette(ROOT.kBird) 
ROOT.gStyle.SetOptStat(0)
h2_heatmap.Draw("COLZ")

leg_corr = ROOT.TLegend(0.18, 0.82, 0.65, 0.88)
leg_corr.SetBorderSize(0)
leg_corr.SetFillColor(ROOT.kWhite)
leg_corr.SetFillStyle(1001)
leg_corr.SetTextSize(0.035)
leg_corr.AddEntry("", f"Pearson Corr. = {correlation_data:.4f}", "")
leg_corr.Draw()

# Save plots
output_path_pdf = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit_projections_data_rho_cheby.pdf"
output_path_png = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit_projections_data_rho_cheby.png"
canvas_chebychev.SaveAs(output_path_pdf)
canvas_chebychev.SaveAs(output_path_png)


#-----------------------------------------------------------------------------------------------


canvas_exponential = ROOT.TCanvas("canvas", "Fit Projections & Correlation", 1800, 600)
canvas_exponential.Divide(3, 1)

for i in range(1, 4):
    canvas_exponential.cd(i)
    ROOT.gPad.SetLeftMargin(0.2)

# Higgs Mass proj
canvas_exponential.cd(1)

#Exponential frame
xframe_exponential = bosonMass.frame(30)
data_blinded.plotOn(xframe_exponential, ROOT.RooFit.CutRange("LowSideband,HighSideband"))
pdf_2D_chebyplusexp.plotOn(xframe_exponential, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.NormRange("LowSideband,HighSideband"), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(1))
xframe_exponential.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
xframe_exponential.GetXaxis().SetTitle("m_{M,#gamma} [GeV]")
xframe_exponential.SetMaximum(1.3*xframe_exponential.GetMaximum())
xframe_exponential.GetYaxis().SetTitleOffset(1.8)
xframe_exponential.Draw()
'''
chi2_H = frame_H.chiSquare("pdf_H", "data_H")
binWidth_H = frame_H.getFitRangeBinW()
frame_H.GetYaxis().SetTitle(f"Events / {binWidth_H:.1f} GeV")
frame_H.GetXaxis().SetTitle("m_{Trk^{+}Trk^{-}#gamma} [GeV]")
frame_H.SetMaximum(frame_H.GetMaximum() * 1.30)
frame_H.Draw()
'''
#Calculate Chi square and parameters 
nParam_exp = fitResult_exponential.floatParsFinal().getSize()
print('n param cheby = ',nParam_exp)
chi2_exp = xframe_exponential.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xSidLeft-xLowRange+xHighRange-xSidRight) - nParam_exp
cut_chi2_exp = "{:.2f}".format(chi2_exp*my_ndof) #Crop the chi2 to 2 decimal digits
print('Chi square exp = ',chi2_exp)
print('')

leg2 = ROOT.TLegend(0.5,0.22,0.72,0.90) #right positioning
leg2.SetHeader(" ")
leg2.SetNColumns(1)
leg2.SetFillColorAlpha(0,0.)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)
leg2.AddEntry(cut_chi2_exp,"#chi^{2}/ndof = " + cut_chi2_exp + " / " + str(my_ndof),"brNDC")
leg2.Draw()

#box_exp = draw_parameters([a2_bkg, b2_bkg, e_bkg, f_bkg],x1=0.42, y1=0.70, x2=0.90, y2=0.90)

# Rho Meson proj
canvas_exponential.cd(2)
frame_rho_exp = mesonMass.frame(30)
data_blinded.plotOn(frame_rho_exp)
pdf_2D_chebyplusexp.plotOn(frame_rho_exp, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(1), ROOT.RooFit.Name("pdf_rho"))
frame_rho_exp.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
frame_rho_exp.GetXaxis().SetTitle("m_{M} [GeV]")
frame_rho_exp.SetMaximum(1.3*frame_rho_exp.GetMaximum())
frame_rho_exp.GetYaxis().SetTitleOffset(1.8)
frame_rho_exp.Draw()

'''
chi2_rho = frame_rho_cheby.chiSquare("pdf_rho", "data_rho")
binWidth_rho = frame_rho_cheby.getFitRangeBinW()
frame_rho_cheby.GetYaxis().SetTitle(f"Events / {binWidth_rho:.3f} GeV")
frame_rho_cheby.GetXaxis().SetTitle("m_{Trk^{+}Trk^{-}} [GeV]")
frame_rho_cheby.SetMaximum(frame_rho_cheby.GetMaximum() * 1.20)
frame_rho_cheby.Draw()
'''

#Calculate Chi square and parameters 
nParam_voig = fitResult_exponential.floatParsFinal().getSize()
print('n param voig = ',nParam_voig)
chi2_voig = frame_rho_exp.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof_rho = 30 - nParam_voig
cut_chi2_voig = "{:.2f}".format(chi2_voig*my_ndof_rho) #Crop the chi2 to 2 decimal digits
print('Chi square voig = ',chi2_voig)
print('')
leg_rho = ROOT.TLegend(0.5,0.22,0.72,0.90) #right positioning
leg_rho.SetHeader(" ")
leg_rho.SetNColumns(1)
leg_rho.SetFillColorAlpha(0,0.)
leg_rho.SetBorderSize(0)
leg_rho.SetLineColor(1)
leg_rho.SetLineStyle(1)
leg_rho.SetLineWidth(1)
leg_rho.SetFillStyle(1001)
leg_rho.AddEntry(cut_chi2_voig,"#chi^{2}/ndof = " + cut_chi2_voig + " / " + str(my_ndof),"brNDC")
leg_rho.Draw()

#box_rho = draw_parameters([nsig, nbkg, mean, sigma, a1_off, a2_off, a3_off],x1=0.68, y1=0.55, x2=0.94, y2=0.90)

# --- Panel 3: 2D Heatmap (COLZ) ---
canvas_exponential.cd(3)
#ROOT.gPad.SetRightMargin(0.)

h2_heatmap = data.createHistogram("h2_heatmap", bosonMass, ROOT.RooFit.Binning(50), ROOT.RooFit.YVar(mesonMass, ROOT.RooFit.Binning(50)))
h2_heatmap.SetTitle("Scatter plot")
h2_heatmap.GetXaxis().SetTitle("m_{M#gamma} [GeV]")
h2_heatmap.GetYaxis().SetTitle("m_{M} [GeV]")

ROOT.gStyle.SetPalette(ROOT.kBird) 
ROOT.gStyle.SetOptStat(0)
h2_heatmap.Draw("COLZ")

leg_corr = ROOT.TLegend(0.18, 0.82, 0.65, 0.88)
leg_corr.SetBorderSize(0)
leg_corr.SetFillColor(ROOT.kWhite)
leg_corr.SetFillStyle(1001)
leg_corr.SetTextSize(0.035)
leg_corr.AddEntry("", f"Pearson Corr. = {correlation_data:.4f}", "")
leg_corr.Draw()

# Save plots
output_path_pdf = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit_projections_data_rho_exp.pdf"
output_path_png = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit_projections_data_rho_exp.png"
canvas_exponential.SaveAs(output_path_pdf)
canvas_exponential.SaveAs(output_path_png)

#----------------------------------------------------------------------

canvas_chebychev_Nplus1 = ROOT.TCanvas("canvas", "Fit Projections & Correlation", 1800, 600)
canvas_chebychev_Nplus1.Divide(3, 1)

for i in range(1, 4):
    canvas_chebychev_Nplus1.cd(i)
    ROOT.gPad.SetLeftMargin(0.2)

# Higgs Mass proj
canvas_chebychev_Nplus1.cd(1)

#Chebychev_Nplus1 frame
xframe_chebychev_Nplus1 = bosonMass.frame(30)
data_blinded.plotOn(xframe_chebychev_Nplus1, ROOT.RooFit.CutRange("LowSideband,HighSideband"))
pdf_2D_cheby_Nplus1.plotOn(xframe_chebychev_Nplus1, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.NormRange("LowSideband,HighSideband"), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(1))
xframe_chebychev_Nplus1.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
xframe_chebychev_Nplus1.GetXaxis().SetTitle("m_{M,#gamma} [GeV]")
xframe_chebychev_Nplus1.SetMaximum(1.3*xframe_chebychev_Nplus1.GetMaximum())
xframe_chebychev_Nplus1.GetYaxis().SetTitleOffset(1.8)
xframe_chebychev_Nplus1.Draw()
'''
chi2_H = frame_H.chiSquare("pdf_H", "data_H")
binWidth_H = frame_H.getFitRangeBinW()
frame_H.GetYaxis().SetTitle(f"Events / {binWidth_H:.1f} GeV")
frame_H.GetXaxis().SetTitle("m_{Trk^{+}Trk^{-}#gamma} [GeV]")
frame_H.SetMaximum(frame_H.GetMaximum() * 1.30)
frame_H.Draw()
'''
#Calculate Chi square and parameters 
nParam_cheby_Nplus1 = fitResult_chebychev_Nplus1.floatParsFinal().getSize()
print('n param cheby_Nplus1 = ',nParam_cheby_Nplus1)
chi2_cheby_Nplus1 = xframe_chebychev_Nplus1.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xSidLeft-xLowRange+xHighRange-xSidRight) - nParam_cheby_Nplus1
cut_chi2_cheby_Nplus1 = "{:.2f}".format(chi2_cheby_Nplus1*my_ndof) #Crop the chi2 to 2 decimal digits
print('Chi square cheby_Nplus1 = ',chi2_cheby_Nplus1)
print('')

leg3 = ROOT.TLegend(0.5,0.22,0.72,0.90) #right positioning
leg3.SetHeader(" ")
leg3.SetNColumns(1)
leg3.SetFillColorAlpha(0,0.)
leg3.SetBorderSize(0)
leg3.SetLineColor(1)
leg3.SetLineStyle(1)
leg3.SetLineWidth(1)
leg3.SetFillStyle(1001)
leg3.AddEntry(cut_chi2_cheby_Nplus1,"#chi^{2}/ndof = " + cut_chi2_cheby_Nplus1 + " / " + str(my_ndof),"brNDC")
leg3.Draw()

#box_cheby_Nplus1 = draw_parameters([a1_bkg, b2_bkg, c1_bkg, d1_bkg],x1=0.42, y1=0.70, x2=0.90, y2=0.90)

# Rho Meson proj
canvas_chebychev_Nplus1.cd(2)
frame_rho_cheby_Nplus1 = mesonMass.frame(30)
data_blinded.plotOn(frame_rho_cheby_Nplus1)
pdf_2D_cheby_Nplus1.plotOn(frame_rho_cheby_Nplus1, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(1), ROOT.RooFit.Name("pdf_rho"))
frame_rho_cheby_Nplus1.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
frame_rho_cheby_Nplus1.GetXaxis().SetTitle("m_{M} [GeV]")
frame_rho_cheby_Nplus1.SetMaximum(1.3*frame_rho_cheby_Nplus1.GetMaximum())
frame_rho_cheby_Nplus1.GetYaxis().SetTitleOffset(1.8)
frame_rho_cheby_Nplus1.Draw()

'''
chi2_rho = frame_rho_cheby.chiSquare("pdf_rho", "data_rho")
binWidth_rho = frame_rho_cheby.getFitRangeBinW()
frame_rho_cheby.GetYaxis().SetTitle(f"Events / {binWidth_rho:.3f} GeV")
frame_rho_cheby.GetXaxis().SetTitle("m_{Trk^{+}Trk^{-}} [GeV]")
frame_rho_cheby.SetMaximum(frame_rho_cheby.GetMaximum() * 1.20)
frame_rho_cheby.Draw()
'''

#Calculate Chi square and parameters 
nParam_voig = fitResult_chebychev_Nplus1.floatParsFinal().getSize()
print('n param voig = ',nParam_voig)
chi2_voig = frame_rho_cheby_Nplus1.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof_rho = 30 - nParam_voig
cut_chi2_voig = "{:.2f}".format(chi2_voig*my_ndof_rho) #Crop the chi2 to 2 decimal digits
print('Chi square voig = ',chi2_voig)
print('')
leg_rho = ROOT.TLegend(0.5,0.22,0.72,0.90) #right positioning
leg_rho.SetHeader(" ")
leg_rho.SetNColumns(1)
leg_rho.SetFillColorAlpha(0,0.)
leg_rho.SetBorderSize(0)
leg_rho.SetLineColor(1)
leg_rho.SetLineStyle(1)
leg_rho.SetLineWidth(1)
leg_rho.SetFillStyle(1001)
leg_rho.AddEntry(cut_chi2_voig,"#chi^{2}/ndof = " + cut_chi2_voig + " / " + str(my_ndof),"brNDC")
leg_rho.Draw()

#box_rho = draw_parameters([nsig, nbkg, mean, sigma, a1_off, a2_off, a3_off],x1=0.68, y1=0.55, x2=0.94, y2=0.90)

# --- Panel 3: 2D Heatmap (COLZ) ---
canvas_chebychev_Nplus1.cd(3)
#ROOT.gPad.SetRightMargin(0.)

h2_heatmap = data.createHistogram("h2_heatmap", bosonMass, ROOT.RooFit.Binning(50), ROOT.RooFit.YVar(mesonMass, ROOT.RooFit.Binning(50)))
h2_heatmap.SetTitle("Scatter plot")
h2_heatmap.GetXaxis().SetTitle("m_{M#gamma} [GeV]")
h2_heatmap.GetYaxis().SetTitle("m_{M} [GeV]")

ROOT.gStyle.SetPalette(ROOT.kBird) 
ROOT.gStyle.SetOptStat(0)
h2_heatmap.Draw("COLZ")

leg_corr = ROOT.TLegend(0.18, 0.82, 0.65, 0.88)
leg_corr.SetBorderSize(0)
leg_corr.SetFillColor(ROOT.kWhite)
leg_corr.SetFillStyle(1001)
leg_corr.SetTextSize(0.035)
leg_corr.AddEntry("", f"Pearson Corr. = {correlation_data:.4f}", "")
leg_corr.Draw()

# Save plots
output_path_pdf = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit_projections_data_rho_cheby_Nplus1.pdf"
output_path_png = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit_projections_data_rho_cheby_Nplus1.png"
canvas_chebychev_Nplus1.SaveAs(output_path_pdf)
canvas_chebychev_Nplus1.SaveAs(output_path_png)


# Multipdf ------------------------------------------------------------------------------------------------------------------------------
a_bkg.setConstant(True);
b_bkg.setConstant(True);
c_bkg.setConstant(True);
a1_bkg.setConstant(True);
b1_bkg.setConstant(True);
c1_bkg.setConstant(True);
d1_bkg.setConstant(True);
a2_bkg.setConstant(True);
b2_bkg.setConstant(True);
c2_bkg.setConstant(True);
e_bkg.setConstant(True);
f_bkg.setConstant(True);


cat = ROOT.RooCategory("pdf_index_GFcat_bdt0","Index of Pdf which is active")
mypdfs = ROOT.RooArgList()
mypdfs.add(pdf_2D_cheby)
mypdfs.add(pdf_2D_cheby_Nplus1)
mypdfs.add(pdf_2D_chebyplusexp)

multipdf = ROOT.RooMultiPdf("multipdf_"+CHANNEL+"_GFcat_bdt0_bkg","All Pdfs",cat,mypdfs)

#create Workspace ------------------------------------------------------------------------------------------------------------------------------
norm     = nEntries #fileInput.Get("h_InvMass_TwoTrk_Photon").Integral() #get the normalization of ggH signal (area under ggH signal)
bkg_norm = ROOT.RooRealVar(multipdf.GetName()+ "_norm", multipdf.GetName()+ "_norm", norm,0.5*norm, 2*norm)

print('************************************** n. events = ',nEntries)
workspace = ROOT.RooWorkspace("workspace_STAT_"+CHANNEL+"_GFcat_bdt0_2024")
getattr(workspace,'import')(cat)
getattr(workspace,'import')(multipdf)
getattr(workspace,'import')(observed_data)
getattr(workspace,'import')(bkg_norm)
print('integral BKG',bkg_norm.getVal())
#workspace.Print()

inputWS = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_bdt0_2024.root","RECREATE") #there's only one ws for both ggH and VBF 
inputWS.cd()
workspace.Write()
#print('-------------------------------------------"
#print('Final print to check the workspace update:"
#workspace.Print()
inputWS.Close()