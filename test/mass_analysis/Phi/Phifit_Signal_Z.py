import ROOT
import cmsstyle
import argparse

#bool initialization
includeResiduals   = True  # Flag to include the residuals subplot

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Decay_channel_option', help='Type <<Phi>> for Phi') #flag for bkg estimation
args = p.parse_args()

CHANNEL = "Phi"
print ('Z -> PhiGamma analysis')

#Signal input rootfiles ---------------------------------------------------------------
fileInput_Z = ROOT.TFile("./SR_ZPhi_ORTrigger_BDT_Signal.root")
tree_Z = fileInput_Z.Get("tree_output")

#Import the Doube Crystal Ball PDF -----------------------------------------------------
#ROOT.gROOT.ProcessLineSync(".L ./dCB/RooDoubleCBFast.cc+")

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
xLowRange  = 80.
xHighRange = 105.

h_mH_Z = ROOT.TH1F("h_mH_Z","h_mH_Z", int(xHighRange - xLowRange)*4, xLowRange, xHighRange)

nentries_Z = tree_Z.GetEntriesFast()

for jentry in range(nentries_Z):
    ientry = tree_Z.LoadTree( jentry )
    if ientry < 0:
        break
    nb = tree_Z.GetEntry(jentry)
    if nb <= 0:
        print ('nb < 0')
        continue

    h_mH_Z.Fill(tree_Z.bosonMass, tree_Z._eventWeight)

#Define the observable ---------------------------------------------------------------
bosonMass = ROOT.RooRealVar("bosonMass", "bosonMass", 95., xLowRange, xHighRange, "GeV/c^2")

#Double Crystal Ball definition ---------------------------------------------------------------
dCB_pole_Z  = ROOT.RooRealVar("dCB_pole_"+CHANNEL+"_Z_bdt0_Z", "Double CB pole", 91.,89.,93.)
dCB_width_Z = ROOT.RooRealVar("dCB_width_"+CHANNEL+"_Z_bdt0_Z", "Double CB width",2.8,0.5,1.5)
dCB_aL_Z    = ROOT.RooRealVar("dCB_aL_"+CHANNEL+"_Z_bdt0_Z", "Double CB alpha left", 1.2, 1., 2.)
dCB_aR_Z    = ROOT.RooRealVar("dCB_aR_"+CHANNEL+"_Z_bdt0_Z", "Double CB alpha right", 1.5, 1., 2.)
dCB_nL_Z    = ROOT.RooRealVar("dCB_nL_"+CHANNEL+"_Z_bdt0_Z", "Double CB n left", 3.1, 0.1, 10.)
dCB_nR_Z    = ROOT.RooRealVar("dCB_nR_"+CHANNEL+"_Z_bdt0_Z", "Double CB n right", 3.1, 0.1, 10.)
#voi_pole_Z  = ROOT.RooRealVar("voi_pole_"+CHANNEL+"_Z_bdt0_Z", "Voightian pole", 91.,89.,93.)
#voi_width_Z = ROOT.RooRealVar("voi_width_"+CHANNEL+"_Z_bdt0_Z", "Voightian width",2.4,1.5,3.5)
#voi_sigma_Z = ROOT.RooRealVar("voi_sigma_"+CHANNEL+"_Z_bdt0_Z", "Voightian sigma",1.5,0.5,3.)
#gauss_mean_Z     = ROOT.RooRealVar("gauss_mean_"+CHANNEL+"_Z_bdt0_Z", "Gaussian width", 121., 118., 122.)
#gauss_width_Z    = ROOT.RooRealVar("gauss_width_"+CHANNEL+"_Z_bdt0_Z", "Gaussian width", 2., 1.5, 5.)
#f_gauss_Z        = ROOT.RooRealVar("f_gauss_"+CHANNEL+"_Z_bdt0_Z", "Gaussian fraction", 0.3, 0.2, 0.4)
#ellpis_axis_Z     = ROOT.RooRealVar("ellpis_axis_"+CHANNEL+"_Z_bdt0_Z", "Axis ellipse", 5.5, 3., 8.)
#ellpis_cent_Z     = ROOT.RooRealVar("ellpis_cent_"+CHANNEL+"_Z_bdt0_Z", "Center ellipse", 120.5, 118., 123.)
#f_ellpis_Z        = ROOT.RooRealVar("f_ellpis_"+CHANNEL+"_Z_bdt0_Z", "Ellipse fraction", 0.3, 0.1, 0.6)
#r_ellpis            = ROOT.RooRealVar("r_ellpis_"+CHANNEL+"_Z_bdt0_Z", "Ellipse normalization", 0.83, 0.0, 1.0)
#r_ellpis.setConstant(True)

signalPDF_Z = ROOT.RooDoubleCBFast("signalPDF_"+CHANNEL+"_Z_bdt0_Z", "Double Crystal Ball", bosonMass, dCB_pole_Z, dCB_width_Z, dCB_aL_Z, dCB_nL_Z, dCB_aR_Z, dCB_nR_Z)
#signalPDF_Z = ROOT.RooVoigtian("signalPDF_"+CHANNEL+"_Z_bdt0_Z", "Voigtian", bosonMass, voi_pole_Z, voi_width_Z, voi_sigma_Z)
#gauss_Z = ROOT.RooGaussian("gauss_"+CHANNEL+"_Z_bdt0_Z", "Gaussian", bosonMass, gauss_mean_Z, gauss_width_Z)
#halfEllipse_Z = ROOT.RooGenericPdf(
#    "halfEllipse_Z",
#    "((abs(bosonMass-ellpis_cent_Phi_Z_bdt0_Z) < r_ellpis_Phi_Z_bdt0_Z*ellpis_axis_Phi_Z_bdt0_Z) ? sqrt(1 - ((bosonMass-ellpis_cent_Phi_Z_bdt0_Z)/ellpis_axis_Phi_Z_bdt0_Z)^2) : sqrt(1-r_ellpis_Phi_Z_bdt0_Z^2)*exp(-2*(bosonMass-ellpis_cent_Phi_Z_bdt0_Z)^2/ellpis_axis_Phi_Z_bdt0_Z^2)* exp(2*r_ellpis_Phi_Z_bdt0_Z^2))",
#    ROOT.RooArgList(bosonMass, ellpis_cent_Z, ellpis_axis_Z, r_ellpis)
#)
#signalPDF_Z = ROOT.RooAddPdf("signalPDF_"+CHANNEL+"_Z_bdt0_Z", "total",ROOT.RooArgList(halfEllipse_Z,dcb_Z),ROOT.RooArgList(f_ellpis_Z))


#Input tree ------------------------------------------------------------------------------------------------------------------------------
fileInput_Z.cd()
tree_Z = fileInput_Z.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also ---------------------------------------------------------------
#dataset_Z = ROOT.RooDataSet("dataset_Z","dataset_Z", ROOT.RooArgSet(bosonMass), ROOT.RooFit.Import(tree_Z))
dataset_Z = ROOT.RooDataHist("dataset_Z", "dataset_Z", ROOT.RooArgList(bosonMass), h_mH_Z)

#Do the fit ------------------------------------------------------------------------------------------------------------------------------
fitResult_Z = signalPDF_Z.fitTo(dataset_Z,ROOT.RooFit.Save())

#Plot ------------------------------------------------------------------------------------------------------------------------------
# Custom canvas dimensions
canvas_width      = 800
canvas_height     = 800
canvas_margin_top = 100  # Empty space added at the top of the main canvas

xframe_Z = bosonMass.frame(int(xHighRange - xLowRange)*4)
dataset_Z.plotOn(xframe_Z)
signalPDF_Z.plotOn(xframe_Z)
xframe_Z.SetTitle("")
ROOT.gStyle.SetOptFit(1111)
xframe_Z.SetMaximum(1.2 * xframe_Z.GetMaximum())
signalPDF_Z.paramOn(xframe_Z, ROOT.RooFit.Layout(0.53, 0.94, 0.91), ROOT.RooFit.Format("NEU", ROOT.RooFit.AutoPrecision(1)))  # ,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_Z.getAttText().SetLineWidth(0)
xframe_Z.getAttText().SetTextSize(0.019)
xframe_Z.GetXaxis().SetTitle("m_{M,#gamma} [GeV]")
xframe_Z.GetXaxis().SetRangeUser(xLowRange, xHighRange)
#xframe_Z.GetYaxis().SetMaxDigits(2)

# ChiSquare test ---------------------------------------------------------------------------------------------------------------------------
nParam_Z = fitResult_Z.floatParsFinal().getSize()
chi2_Z = xframe_Z.chiSquare()  # Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xHighRange-xLowRange)*4 - nParam_Z
cut_chi2_Z = "{:.2f}".format(chi2_Z*my_ndof) #Crop the chi2 to 2 decimal digits

print('nParam_Z =', nParam_Z)
print('cut_chi2_Z =', cut_chi2_Z)

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

xframe_Z.Draw()

if includeResiduals:
    c1.cd(2)  # Activate the second pad for the residuals plot
    c1.GetPad(2).SetBottomMargin(0.15)  # Increase the empty space at the bottom of the residuals plot
    residuals = xframe_Z.pullHist()
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
leg1 = ROOT.TLegend(0.6, 0.29, 0.87, 0.70)  # right positioning
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
leg1.AddEntry(cut_chi2_Z,"#chi^{2} / ndof = " + cut_chi2_Z + " / " + str(my_ndof),"brNDC")
leg1.Draw()                                                                                                                                         

if includeResiduals:
	c1.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_Z_resid.pdf")
	c1.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_Z_resid.png")
else:
	c1.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_Z.pdf")
	c1.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_Z.png")  
#######################################################################################################

#the same for VBF
#fileInput_VBF.cd()
#tree_VBF = fileInput_VBF.Get("tree_output")

#Retrieve dataset from the tree, insert the variable also
#dataset_VBF = ROOT.RooDataSet("dataset_VBF","dataset_VBF",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree_VBF))
#dataset_VBF = ROOT.RooDataHist("dataset_VBF", "dataset_VBF", ROOT.RooArgList(mass), h_mH_VBF)

#Do the fit
#fitResult_VBF = signalPDF_VBF.fitTo(dataset_VBF,ROOT.RooFit.Save())

#Plot
#xframe_VBF = mass.frame(int(xHighRange - xLowRange)*10)
#dataset_VBF.plotOn(xframe_VBF)
#signalPDF_VBF.plotOn(xframe_VBF)
#xframe_VBF.SetTitle("")
#ROOT.gStyle.SetOptFit(1111)
#xframe_VBF.SetMaximum(1.2*xframe_VBF.GetMaximum())
#signalPDF_VBF.paramOn(xframe_VBF,ROOT.RooFit.Layout(0.53,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
#xframe_VBF.getAttText().SetLineWidth(0)
#xframe_VBF.getAttText().SetTextSize(0.019)
#xframe_VBF.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
#xframe_VBF.GetXaxis().SetRangeUser(xLowRange,xHighRange)
#xframe_VBF.GetYaxis().SetMaxDigits(2)

# ChiSquare test ---------------------------------------------------------------------------------------------------------------------------
#nParam_VBF = fitResult_VBF.floatParsFinal().getSize()
#chi2_VBF = xframe_VBF.chiSquare()  # Returns chi2/ndof. Remember to remove the option XErrorSize(0) from data.PlotOn
#cut_chi2_VBF = "{:.2f}".format(chi2_VBF)  # Crop the chi2 to 2 decimal digits
#
#print('nParam_VBF =", nParam_VBF
#print('cut_chi2_VBF =", cut_chi2_VBF
#
#c2 = ROOT.TCanvas("c2", "c2", canvas_width, canvas_height + canvas_margin_top)
#c2.cd()
#c2.SetTitle("")
#c2.SetBottomMargin(0.15)  # Increase the empty space at the bottom of the main canvas
#
#if includeResiduals:
#    c2.Divide(1, 2, 0, 0.25)  # Divide canvas into two pads: 1x2 grid, with the first pad being 4 times larger
#    c2.cd(1)  # Activate the first pad for the main plot
#else:
#    c2.Divide(1, 1)
#
#xframe_VBF.Draw()

#if includeResiduals:
#    c2.cd(2)  # Activate the second pad for the residuals plot
#    c2.GetPad(2).SetBottomMargin(0.35)  # Increase the empty space at the bottom of the residuals plot
#    residuals = xframe_VBF.residHist()
#    residuals.SetFillColor(ROOT.kBlue)
#    residuals.SetMarkerColor(ROOT.kBlack)
#    residuals.GetYaxis().SetTitle("Resid")
#    residuals.GetYaxis().SetTitleOffset(1.)
#    residuals.GetYaxis().SetTitleSize(0.06)
#    residuals.GetYaxis().SetLabelSize(0.05)
#    residuals.GetXaxis().SetTitleSize(0.06)
#    residuals.GetXaxis().SetLabelSize(0.05)
#    residuals.GetYaxis().SetTitleOffset(1.)
#    residuals.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
#    residuals.GetXaxis().SetRangeUser(xLowRange, xHighRange)
#    residuals.Draw("AP")
#
#c2.Update()

#if includeResiduals:
#	c2.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_VBF_resid.pdf")
#	c2.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_VBF_resid.png")
#else:
#	c2.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_VBF.pdf")
#	c2.SaveAs("/eos/user/c/covarell/www/AN-26-106/fitsignal_" + CHANNEL + "_BDT_VBF.png")
                                      
#create Workspace
if CHANNEL == "Phi": #norm factor x2 to include the Dbar channel
    norm_Z     = fileInput_Z.Get("h_bosonMass").Integral() #get the normalization of Z signal (area under Z signal)
    sig_norm_Z = ROOT.RooRealVar(signalPDF_Z.GetName()+ "_norm", signalPDF_Z.GetName()+ "_norm", norm_Z) 
    #norm_VBF     = fileInput_VBF.Get("h_bosonMass").Integral() #get the normalization of VBF signal (area under VBF signal)
    #sig_norm_VBF = ROOT.RooRealVar(signalPDF_VBF.GetName()+ "_norm", signalPDF_VBF.GetName()+ "_norm", 0.0394*2*norm_VBF)
else:
    norm_Z     = fileInput_Z.Get("h_bosonMass").Integral() #get the normalization of Z signal (area under Z signal)
    sig_norm_Z = ROOT.RooRealVar(signalPDF_Z.GetName()+ "_norm", signalPDF_Z.GetName()+ "_norm", norm_Z)
    #norm_VBF     = fileInput_VBF.Get("h_bosonMass").Integral() #get the normalization of VBF signal (area under VBF signal)
    #sig_norm_VBF = ROOT.RooRealVar(signalPDF_VBF.GetName()+ "_norm", signalPDF_VBF.GetName()+ "_norm", norm_VBF)    

dCB_pole_Z.setConstant(1)  
dCB_width_Z.setConstant(1)  
dCB_aL_Z.setConstant(1)     
dCB_aR_Z.setConstant(1)    
dCB_nL_Z.setConstant(1)    
dCB_nR_Z.setConstant(1)     
sig_norm_Z.setConstant(1)
#ellpis_axis_Z.setConstant(1)
#ellpis_cent_Z.setConstant(1)
#f_ellpis_Z.setConstant(1)

#dCB_pole_VBF.setConstant(1)  
#dCB_width_VBF.setConstant(1)  
#dCB_aL_VBF.setConstant(1)     
#dCB_aR_VBF.setConstant(1)    
#dCB_nL_VBF.setConstant(1)    
#dCB_nR_VBF.setConstant(1)     
#sig_norm_VBF.setConstant(1)

'''
#Systematics --------------------------------------------------------------------
fileInput_Z_TwoProngsSFUP  = ROOT.TFile("histos/latest_production/histos_SR_SignalZTwoProngsSFUp.root")
fileInput_Z_TwoProngsSFDW  = ROOT.TFile("histos/latest_production/histos_SR_SignalZTwoProngsSFDown.root")

signalZ_nominal = fileInput_Z.Get("h_bosonMass")
signalZ_nominal.SetName("signal")
signalZ_trigUP  = fileInput_Z_TwoProngsSFUP.Get("h_bosonMass")
signalZ_trigUP.SetName("signal_trigEffTwoProngsUp")
signalZ_trigDW  = fileInput_Z_TwoProngsSFDW.Get("h_bosonMass")
signalZ_trigDW.SetName("signal_trigEffTwoProngsDown")
#---------------------------------------------------------------------------------
'''
fOut = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_Z_bdt0_2024.root","UPDATE")
fOut.cd()

workspace = fOut.Get("workspace_STAT_"+CHANNEL+"_Z_bdt0_2024")
getattr(workspace,'import')(signalPDF_Z)
#getattr(workspace,'import')(signalPDF_VBF)
getattr(workspace,'import')(sig_norm_Z)
#getattr(workspace,'import')(sig_norm_VBF)
#getattr(workspace,'import')(dataset_Z)
#getattr(workspace,'import')(dataset_VBF)
print('Z Signal integral = ',sig_norm_Z.getVal())
#print('VBF Signal integral = ",sig_norm_VBF.Print()

workspace.Print()

print('----------------------')
print('nParam_Z   = ', nParam_Z)
print('chi2_Z     = ', chi2_Z)
#print('chi2 / ndof  = ', chi2_Z/float(nParam_Z))
#print('"
#print('nParam_VBF   = ", nParam_VBF
#print('chi2_VBF     = ", chi2_VBF
#print('chi2 / ndof  = ", chi2_VBF/float(nParam_VBF)
print('----------------------')
print('')

workspace.Write()
xframe_Z.Write("xframe_Z")
#xframe_VBF.Write("xframe_VBF")
h_mH_Z.Write("h_mH_Z")
#h_mH_VBF.Write("h_mH_VBF")
#signalZ_nominal.Write()
#signalZ_trigUP.Write()
#signalZ_trigDW.Write()

fOut.Close()	
