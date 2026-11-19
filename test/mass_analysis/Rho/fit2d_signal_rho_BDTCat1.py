import ROOT
import cmsstyle
import argparse

# Import the Double Crystal Ball PDF 
#ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

# Suppress the opening of many Canvas's 
ROOT.gROOT.SetBatch(True) 

CHANNEL = "Rho"

#CMS-style plotting ---------------------------------------------------------------
cmsstyle.setCMSStyle()
iPeriod = 4
iPos = 11

'''
#input
file_path = "../../histos/H/SR_HRho_ORTrigger_BDTCat1_Signal.root"
tree_name = "tree_output"
infile = ROOT.TFile.Open(file_path)
tree = infile.Get(tree_name)
'''
#Signal input rootfiles ---------------------------------------------------------------
fileInput = ROOT.TFile("../../histos/H/SR_HRhoPU_ORTrigger_BDTCat1_Signal.root")
tree = fileInput.Get("tree_output")

#Define the observable ---------------------------------------------------------------
xLowRange  = 110.
xHighRange = 140.
bosonMass = ROOT.RooRealVar("bosonMass", "m(H) [GeV]", 110, 140)
mesonMass = ROOT.RooRealVar("mesonMass", "m(#rho) [GeV]", 0.62, 0.92)

# Create RooDataSet (Unbinned)
data = ROOT.RooDataSet("data", "Signal", ROOT.RooArgSet(bosonMass, mesonMass), ROOT.RooFit.Import(tree))
print(f"Number of entries loaded in dataset: {data.numEntries()}")
#nEntries = data.numEntries() 


#Double Crystal Ball definition for Higgs ---------------------------------------------------------------
dCB_pole_H  = ROOT.RooRealVar("dCB_pole_H_"+CHANNEL+"_GFcat_bdt1_ggH", "Higgs Mean", 125.0, 120.0, 130.0)
dCB_width_H = ROOT.RooRealVar("dCB_width_H_"+CHANNEL+"_GFcat_bdt1_ggH", "Higgs Sigma", 3.0, 0.5, 8.0)
dCB_aL_H    = ROOT.RooRealVar("dCB_aL_H_"+CHANNEL+"_GFcat_bdt1_ggH", "Higgs Alpha Left", 1.5, 0.1, 15.0)
dCB_aR_H    = ROOT.RooRealVar("dCB_aR_H_"+CHANNEL+"_GFcat_bdt1_ggH", "Higgs Alpha Right", 1.5, 0.1, 15.0)
dCB_nL_H    = ROOT.RooRealVar("dCB_nL_H_"+CHANNEL+"_GFcat_bdt1_ggH", "Higgs n Left", 2.0, 0.1, 50.0)
dCB_nR_H    = ROOT.RooRealVar("dCB_nR_H_"+CHANNEL+"_GFcat_bdt1_ggH", "Higgs n Right", 2.0, 0.1, 50.0)

signalPDF_H = ROOT.RooDoubleCBFast("signalPDF_H_"+CHANNEL+"_GFcat_bdt1_ggH", "DSCB Higgs", bosonMass, dCB_pole_H, dCB_width_H, dCB_aL_H, dCB_nL_H, dCB_aR_H, dCB_nR_H)

#Double Crystal Ball definition for Rho ---------------------------------------------------------------
dCB_pole_rho  = ROOT.RooRealVar("dCB_pole_M_"+CHANNEL+"_GFcat_bdt1_ggH", "Rho Mean", 0.770, 0.74, 0.80)
dCB_width_rho = ROOT.RooRealVar("dCB_width_M_"+CHANNEL+"_GFcat_bdt1_ggH", "Rho Sigma", 0.05, 0.02, 0.12)
dCB_aL_rho    = ROOT.RooRealVar("dCB_aL_M_"+CHANNEL+"_GFcat_bdt1_ggH", "Rho Alpha Left", 1.0, 0.05, 20.0)
dCB_aR_rho    = ROOT.RooRealVar("dCB_aR_M_"+CHANNEL+"_GFcat_bdt1_ggH", "Rho Alpha Right", 1.0, 0.05, 20.0)
dCB_nL_rho    = ROOT.RooRealVar("dCB_nL_M_"+CHANNEL+"_GFcat_bdt1_ggH", "Rho n Left", 1.5, 0.05, 100.0)
dCB_nR_rho    = ROOT.RooRealVar("dCB_nR_M_"+CHANNEL+"_GFcat_bdt1_ggH", "Rho n Right", 1.5, 0.05, 100.0)

signalPDF_M = ROOT.RooDoubleCBFast("signalPDF_M_"+CHANNEL+"_GFcat_bdt1_ggH", "DSCB Rho", mesonMass, dCB_pole_rho, dCB_width_rho, dCB_aL_rho, dCB_nL_rho, dCB_aR_rho, dCB_nR_rho)

# 2D PDF
signalPDF_2D = ROOT.RooProdPdf("signalPDF_2D_"+CHANNEL+"_GFcat_bdt1_ggH", "2D Signal PDF (Higgs x Meson)", ROOT.RooArgList(signalPDF_H, signalPDF_M))

# 2D Unbinned Maximum Likelihood Fit
fitResult = signalPDF_2D.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(1))

print("\n" + "="*50)
print("FIT RESULTS:")
fitResult.Print()
print("="*50 + "\n")

# Correlation
print("Evaluating correlation between observables...")
h2_data = ROOT.TH2F("h2_data", "Data Correlation", 70, 110, 140, 50, 0.62, 0.92)
data.fillHistogram(h2_data, ROOT.RooArgList(bosonMass, mesonMass))
correlation_data = h2_data.GetCorrelationFactor()
print(f"--> Linear correlation (Pearson) in DATA: {correlation_data:.4f}")

def draw_parameters(parameters, x1=0.55, y1=0.55, x2=0.93, y2=0.90):

    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillStyle(0)
    box.SetBorderSize(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(0.022)

    for par in parameters:
        box.AddText(
            f"{par.GetTitle()} = {par.getVal():.4g} #pm {par.getError():.2g}"
        )

    box.Draw()
    return box


canvas = ROOT.TCanvas("canvas", "Fit Projections & Correlation", 1800, 600)
canvas.Divide(3, 1)

for i in range(1, 4):
    canvas.cd(i)
    ROOT.gPad.SetLeftMargin(0.15)

# Higgs Mass proj
canvas.cd(1)
xframe_H = bosonMass.frame()
data.plotOn(xframe_H)
signalPDF_2D.plotOn(xframe_H, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(1))
xframe_H.SetMaximum(1.2 * xframe_H.GetMaximum())
xframe_H.GetXaxis().SetTitle("m_{M,#gamma} [GeV]")
xframe_H.Draw()

# ChiSquare test ---------------------------------------------------------------------------------------------------------------------------
nParam_H = fitResult.floatParsFinal().getSize()
chi2_H = xframe_H.chiSquare()  # Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xHighRange-xLowRange)*4 - nParam_H
cut_chi2_H = "{:.2f}".format(chi2_H*my_ndof) #Crop the chi2 to 2 decimal digits

print('nParam_H =', nParam_H)
print('cut_chi2_H =', cut_chi2_H)

# Legend ----------------------------------------
leg1 = ROOT.TLegend(0.6, 0.29, 0.87, 0.70)  # right positioning
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0, 0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.AddEntry(cut_chi2_H,"#chi^{2} / ndof = " + cut_chi2_H + " / " + str(my_ndof),"brNDC")
leg1.Draw()       

# Higgs DCB parameters
box_H = draw_parameters(
    [dCB_pole_H, dCB_width_H, dCB_aL_H, dCB_nL_H,
     dCB_aR_H, dCB_nR_H],
    x1=0.52, y1=0.57, x2=0.93, y2=0.90
)

# Rho Meson proj
canvas.cd(2)
xframe_M = mesonMass.frame()
data.plotOn(xframe_M)
signalPDF_2D.plotOn(xframe_M, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(1))
xframe_M.SetMaximum(1.2 * xframe_M.GetMaximum())
xframe_M.GetXaxis().SetTitle("m_{M} [GeV]")
xframe_M.Draw()

# ChiSquare test ---------------------------------------------------------------------------------------------------------------------------
nParam_M = fitResult.floatParsFinal().getSize()
chi2_M = xframe_M.chiSquare()  # Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(30)*4 - nParam_H
cut_chi2_M = "{:.2f}".format(chi2_M*my_ndof) #Crop the chi2 to 2 decimal digits

print('nParam_M =', nParam_M)
print('cut_chi2_M =', cut_chi2_M)

# Legend ----------------------------------------
leg2 = ROOT.TLegend(0.6, 0.29, 0.87, 0.70)  # right positioning
leg2.SetHeader(" ")
leg2.SetNColumns(1)
leg2.SetFillColorAlpha(0, 0.)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)
leg2.AddEntry(cut_chi2_H,"#chi^{2} / ndof = " + cut_chi2_M + " / " + str(my_ndof),"brNDC")
leg2.Draw()

# Rho DCB parameters
box_M = draw_parameters(
    [dCB_pole_rho, dCB_width_rho, dCB_aL_rho, dCB_nL_rho,
     dCB_aR_rho, dCB_nR_rho],
    x1=0.52, y1=0.57, x2=0.93, y2=0.90
)

# --- Panel 3: 2D Heatmap (COLZ) ---
canvas.cd(3)
ROOT.gPad.SetRightMargin(0.15)

h2_heatmap = data.createHistogram("h2_heatmap", bosonMass, ROOT.RooFit.Binning(50), ROOT.RooFit.YVar(mesonMass, ROOT.RooFit.Binning(50)))
h2_heatmap.SetTitle("Scatter plot")
h2_heatmap.GetXaxis().SetTitle("m_{Trk^{+}Trk^{-}#gamma} [GeV]")
h2_heatmap.GetYaxis().SetTitle("m_{Trk^{+}Trk^{-}} [GeV]")

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
output_path_pdf = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/fit_projections_signal_2D.pdf"
output_path_png = "/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/fit_projections_signal_2D.png"
canvas.SaveAs(output_path_pdf)
canvas.SaveAs(output_path_png)
print("Plots saved successfully with 2D Heatmap.")


#create Workspace
#norm_H     = fileInput.Get("h_bosonMass").Integral() #get the normalization of ggH signal (area under ggH signal)
#sig_norm_H = ROOT.RooRealVar(signalPDF_H.GetName()+ "_norm", signalPDF_H.GetName()+ "_norm", norm_H) 
#norm_M     = fileInput.Get("h_mesonMass").Integral() #get the normalization of ggH signal (area under ggH signal)
#sig_norm_M = ROOT.RooRealVar(signalPDF_M.GetName()+ "_norm", signalPDF_M.GetName()+ "_norm", norm_M)

norm = fileInput.Get("h_bosonMass").Integral()
sig_norm = ROOT.RooRealVar(signalPDF_2D.GetName() + "_norm", signalPDF_2D.GetName() + "_norm",norm)

   
dCB_pole_H.setConstant(1)  
dCB_width_H.setConstant(1)  
dCB_aL_H.setConstant(1)     
dCB_aR_H.setConstant(1)    
dCB_nL_H.setConstant(1)    
dCB_nR_H.setConstant(1)     
#sig_norm_H.setConstant(1)
dCB_pole_rho.setConstant(1) 
dCB_width_rho.setConstant(1) 
dCB_aL_rho.setConstant(1)
dCB_nL_rho.setConstant(1) 
dCB_aR_rho.setConstant(1)
dCB_nR_rho.setConstant(1)
#sig_norm_M.setConstant(1)

sig_norm.setConstant(1)


fOut = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024.root","UPDATE")
fOut.cd()

workspace = fOut.Get("workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024")
getattr(workspace,'import')(signalPDF_2D)
#getattr(workspace,'import')(sig_norm_H)
#getattr(workspace,'import')(sig_norm_M)
getattr(workspace,'import')(sig_norm)


print('ggH Signal integral = ',sig_norm.getVal())

workspace.Print()

print('----------------------')
print('nParam_H   = ', nParam_H)
print('chi2_H     = ', chi2_H)
print('----------------------')
print('')

workspace.Write()
#xframe_ggH.Write("xframe_H")
#h_mH_ggH.Write("h_mH_H")


fOut.Close()	

f = ROOT.TFile.Open("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024.root")
w = f.Get("workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024")

# Numero di eventi nel dataset
observed_data = w.data("observed_data")
print("observed_data.numEntries():", observed_data.numEntries())

# Normalizzazione del background
bkg_norm = w.var("multipdf_Rho_GFcat_bdt1_bkg_norm")
print("Background norm:", bkg_norm.getVal())

# Normalizzazione del segnale
sig_norm = w.var("signalPDF_2D_Rho_GFcat_bdt1_ggH_norm")
print("Signal norm:", sig_norm.getVal())

print("h_bosonMass integral:", fileInput.Get("h_bosonMass").Integral())
print("Workspace signal norm:", sig_norm.getVal())