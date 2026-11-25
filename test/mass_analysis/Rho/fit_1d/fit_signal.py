import ROOT
import cmsstyle
import argparse

#bool initialization
includeResiduals   = False  # Flag to include the residuals subplot

parser = argparse.ArgumentParser()
parser.add_argument("--cat", type=int, choices=[0, 1], required=True)
args = parser.parse_args()

CAT = args.cat
CATEGORY = f"GFcat_bdt{CAT}"
CHANNEL = "Rho"
print ('H -> RhoGamma analysis')

#Signal input rootfiles ---------------------------------------------------------------
fileInput_ggH = ROOT.TFile(f"../../../histos/H/SR_H{CHANNEL}_ORTrigger_BDTCat{CAT}_Signal.root")
tree_ggH = fileInput_ggH.Get("tree_output")


#Supress the opening of many Canvas's ---------------------------------------------------
ROOT.gROOT.SetBatch(True)   

#CMS-style plotting ---------------------------------------------------------------
cmsstyle.setCMSStyle()
iPeriod = 4
iPos = 11
#CMS_lumi.lumiTextSize = 0.7
#CMS_lumi.lumiTextOffset = 0.25
#CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
#CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#prepare a rebinned TH1 for Higgs mass ---------------------------------------------------------------
xLowRange  = 110.
xHighRange = 140.

#h_mH_ggH = ROOT.TH1F("h_mH_ggH","h_mH_ggH", int(xHighRange - xLowRange)*4, xLowRange, xHighRange)
h_mH_ggH = ROOT.TH1D("h_mH_ggH", "h_mH_ggH",int(xHighRange - xLowRange)*4,xLowRange, xHighRange)
h_mH_ggH.Sumw2()

nentries_ggH = tree_ggH.GetEntriesFast()

for jentry in range(nentries_ggH):
    ientry = tree_ggH.LoadTree( jentry )
    if ientry < 0:
        break
    nb = tree_ggH.GetEntry(jentry)
    if nb <= 0:
        print ('nb < 0')
        continue

    h_mH_ggH.Fill(tree_ggH.bosonMass, tree_ggH._eventWeight)

#Define the observable ---------------------------------------------------------------
bosonMass = ROOT.RooRealVar("bosonMass", "bosonMass", 125., xLowRange, xHighRange, "GeV/c^2")

#Double Crystal Ball definition ---------------------------------------------------------------
dCB_pole_ggH  = ROOT.RooRealVar(f"dCB_pole_{CHANNEL}_{CATEGORY}_ggH", "Double CB pole", 125.,122.,127.)
dCB_width_ggH = ROOT.RooRealVar(f"dCB_width_{CHANNEL}_{CATEGORY}_ggH", "Double CB width",0.8,0.5,1.5)
dCB_aL_ggH    = ROOT.RooRealVar(f"dCB_aL_{CHANNEL}_{CATEGORY}_ggH", "Double CB alpha left", 1.2, 1., 2.)
dCB_aR_ggH    = ROOT.RooRealVar(f"dCB_aR_{CHANNEL}_{CATEGORY}_ggH", "Double CB alpha right", 1.5, 1., 2.)
dCB_nL_ggH    = ROOT.RooRealVar(f"dCB_nL_{CHANNEL}_{CATEGORY}_ggH", "Double CB n left", 3.1, 0.1, 10.)
dCB_nR_ggH    = ROOT.RooRealVar(f"dCB_nR_{CHANNEL}_{CATEGORY}_ggH", "Double CB n right", 3.1, 0.1, 10.)


signalPDF_ggH = ROOT.RooDoubleCBFast(f"signalPDF_{CHANNEL}_{CATEGORY}_ggH", "Double Crystal Ball", bosonMass, dCB_pole_ggH, dCB_width_ggH, dCB_aL_ggH, dCB_nL_ggH, dCB_aR_ggH, dCB_nR_ggH)


#Input tree ------------------------------------------------------------------------------------------------------------------------------
#fileInput_ggH.cd()
#tree_ggH = fileInput_ggH.Get("tree_output")
#Retrieve dataset from the tree, insert the variable also ---------------------------------------------------------------
#dataset_ggH = ROOT.RooDataSet("dataset_ggH","dataset_ggH", ROOT.RooArgSet(bosonMass), ROOT.RooFit.Import(tree_ggH))
dataset_ggH = ROOT.RooDataHist("dataset_ggH", "dataset_ggH", ROOT.RooArgList(bosonMass), h_mH_ggH)

#Do the fit ------------------------------------------------------------------------------------------------------------------------------
fitResult_ggH = signalPDF_ggH.fitTo(dataset_ggH,ROOT.RooFit.Save())

print("\n========== SIGNAL FIT QUALITY ==========")
print("Status  =", fitResult_ggH.status())
print("CovQual =", fitResult_ggH.covQual())
print("EDM     =", fitResult_ggH.edm())
print("minNll  =", fitResult_ggH.minNll())

print("\n========== SIGNAL PARAMETERS ==========")
print("Mean    =", dCB_pole_ggH.getVal())
print("Sigma   =", dCB_width_ggH.getVal())
print("Alpha L =", dCB_aL_ggH.getVal())
print("Alpha R =", dCB_aR_ggH.getVal())
print("n L     =", dCB_nL_ggH.getVal())
print("n R     =", dCB_nR_ggH.getVal())

print("\n========== SIGNAL NORMALIZATION ==========")
print("Histogram integral =", h_mH_ggH.Integral())
print("Dataset sumEntries =", dataset_ggH.sumEntries())

#Plot ------------------------------------------------------------------------------------------------------------------------------
# Custom canvas dimensions
canvas_width      = 800
canvas_height     = 800
canvas_margin_top = 100  # Empty space added at the top of the main canvas

xframe_ggH = bosonMass.frame(int(xHighRange - xLowRange)*4)
dataset_ggH.plotOn(xframe_ggH,ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
signalPDF_ggH.plotOn(xframe_ggH)
xframe_ggH.SetTitle("")
ROOT.gStyle.SetOptFit(1111)
xframe_ggH.SetMaximum(1.2 * xframe_ggH.GetMaximum())
signalPDF_ggH.paramOn(xframe_ggH, ROOT.RooFit.Layout(0.53, 0.94, 0.91), ROOT.RooFit.Format("NEU", ROOT.RooFit.AutoPrecision(1)))  # ,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_ggH.getAttText().SetLineWidth(0)
xframe_ggH.getAttText().SetTextSize(0.019)
xframe_ggH.GetXaxis().SetTitle("m_{M,#gamma} [GeV]")
xframe_ggH.GetXaxis().SetRangeUser(xLowRange, xHighRange)
#xframe_ggH.GetYaxis().SetMaxDigits(2)

# ChiSquare test ---------------------------------------------------------------------------------------------------------------------------
nParam_ggH = fitResult_ggH.floatParsFinal().getSize()
chi2_ggH = xframe_ggH.chiSquare()  # Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xHighRange-xLowRange)*4 - nParam_ggH
cut_chi2_ggH = "{:.2f}".format(chi2_ggH*my_ndof) #Crop the chi2 to 2 decimal digits

print('nParam_ggH =', nParam_ggH)
print('cut_chi2_ggH =', cut_chi2_ggH)
'''
# Count non-empty plotted bins
hist = xframe_ggH.getHist()
nBins_used = sum(
    1 for i in range(hist.GetN())
    if hist.GetPointY(i) != 0
)

# Degrees of freedom (RooPlot convention)
my_ndof = nBins_used - nParam_ggH

# Reduced chi-square
chi2_ggH = xframe_ggH.chiSquare(nParam_ggH)

# Total chi-square
cut_chi2_ggH = "{:.2f}".format(chi2_ggH * my_ndof)
'''
print("nParam_ggH =", nParam_ggH)
#print("nBins_used =", nBins_used)
print("ndof       =", my_ndof)
print("chi2       =", cut_chi2_ggH)
print("chi2/ndof  =", chi2_ggH)
# Create the first canvas (c1) for the main plot
c1 = ROOT.TCanvas("c1", "c1", canvas_width, canvas_height + canvas_margin_top)
c1.cd()
c1.SetTitle("")
c1.SetBottomMargin(0.15)  # Increase the empty space at the bottom of the main canvas

if includeResiduals:
    c1.Divide(1, 2, 0, 0.25)  # Divide canvas into two pads: 1x2 grid, with the first pad being 4 times larger
    c1.cd(1)  # Activate the first pad for the main plot
else:
    c1.Divide(1, 1)

xframe_ggH.Draw()

if includeResiduals:
    c1.cd(2)  # Activate the second pad for the residuals plot
    c1.GetPad(2).SetBottomMargin(0.15)  # Increase the empty space at the bottom of the residuals plot
    residuals = xframe_ggH.pullHist()
    residuals.SetFillColor(ROOT.kBlue)
    residuals.SetMarkerColor(ROOT.kBlack)
    residuals.GetYaxis().SetTitle("Resid")
    residuals.GetYaxis().SetTitleOffset(1.)
    residuals.GetYaxis().SetTitleSize(0.06)
    residuals.GetYaxis().SetLabelSize(0.05)
    residuals.GetXaxis().SetTitleSize(0.06)
    residuals.GetXaxis().SetLabelSize(0.05)
    residuals.GetYaxis().SetTitleOffset(1.)
    residuals.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
    residuals.GetXaxis().SetRangeUser(xLowRange, xHighRange)
    residuals.Draw("AP")

c1.Update()

# Legend ----------------------------------------
leg1 = ROOT.TLegend(0.54, 0.29, 0.81, 0.70)
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0, 0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
# leg1.AddEntry(histo_map["h_MesonMass"],"Data","elp")
# leg1.AddEntry("totPDFOffline","Fit","l")
# leg1.AddEntry("backgroundPDF","Bkg only fit","l")
leg1.AddEntry(cut_chi2_ggH,"#chi^{2} / ndof = " + cut_chi2_ggH + " / " + str(my_ndof),"brNDC")
leg1.Draw()                                                                                                                                         

#if includeResiduals:
#	c1.SaveAs(f"/eos/user/c/covarell/www/AN-26-106/fitsignal_{CHANNEL}_{CATEGORY}_ggH_resid.pdf")
#	c1.SaveAs(f"/eos/user/c/covarell/www/AN-26-106/fitsignal_{CHANNEL}_{CATEGORY}_ggH_resid.png")
#else:
c1.SaveAs(f"/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat{CAT}/Signal_{CHANNEL}_{CATEGORY}_DSCB.pdf")
c1.SaveAs(f"/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat{CAT}/Signal_{CHANNEL}_{CATEGORY}_DSCB.png")  
#######################################################################################################

                                      
#create Workspace
if CHANNEL == "Rho": #norm factor x2 to include the Dbar channel
    norm_ggH     = h_mH_ggH.Integral() #get the normalization of ggH signal (area under ggH signal)
    sig_norm_ggH = ROOT.RooRealVar(signalPDF_ggH.GetName()+ "_norm", signalPDF_ggH.GetName()+ "_norm", norm_ggH) 
    
else:
    norm_ggH     = fileInput_ggH.Get("h_bosonMass").Integral() #get the normalization of ggH signal (area under ggH signal)
    sig_norm_ggH = ROOT.RooRealVar(signalPDF_ggH.GetName()+ "_norm", signalPDF_ggH.GetName()+ "_norm", norm_ggH)
     

dCB_pole_ggH.setConstant(1)  
dCB_width_ggH.setConstant(1)  
dCB_aL_ggH.setConstant(1)     
dCB_aR_ggH.setConstant(1)    
dCB_nL_ggH.setConstant(1)    
dCB_nR_ggH.setConstant(1)     
sig_norm_ggH.setConstant(1)

fOut = ROOT.TFile(f"../workspaces/workspace_STAT_{CHANNEL}_{CATEGORY}_2024_1d.root","UPDATE")
fOut.cd()

workspace = fOut.Get(f"workspace_STAT_{CHANNEL}_{CATEGORY}_2024_1d")
getattr(workspace,'import')(signalPDF_ggH)
getattr(workspace,'import')(sig_norm_ggH)

print('ggH Signal integral = ',sig_norm_ggH.getVal())

workspace.Print()

print('----------------------')
print('nParam_ggH   = ', nParam_ggH)
print('chi2_ggH     = ', chi2_ggH)

print('----------------------')
print('')

workspace.Write()
xframe_ggH.Write("xframe_ggH")
h_mH_ggH.Write("h_mH_ggH")


fOut.Close()	