import ROOT
import cmsstyle
import argparse

# Suppress the opening of many Canvas's 
ROOT.gROOT.SetBatch(True) 

#input
#file_path = "../histos/H/SR_HRho2022Tau_BDT_Data.root"
#file_path = "../histos/H/SR_HRho2024_ORTrigger_BDTCat1_Data.root"
file_path = "../../histos/H/FullDataset_HRho2024_ORTrigger_BDTCat1_Data.root"
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

#Parameters of the PDF ---------------------------------------------------------------
xLowRange  = 110.
xHighRange = 140.
xSidLeft = 120.
xSidRight = 130.

mass = ROOT.RooRealVar("bosonMass","bosonMass",xLowRange,xHighRange,"GeV")
mass.setRange("LowSideband",xLowRange,120.)
mass.setRange("HighSideband",130.,xHighRange)
mass.setRange("full",xLowRange,xHighRange)

'''
#Initialize a Chebychev pdf
a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",-1.,-3.,0.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.3,-3.,1.)

#c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",-0.01,-1.,1.)
#d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",-0.05,-0.2,0.1)
#f_bkg = ROOT.RooRealVar("e_bkg","e_bkg",-0.05,-0.1,0.)

bkgPDF_chebychev  = ROOT.RooChebychev("chebychev","bkgPDF",mass,ROOT.RooArgList(a_bkg, b_bkg))
'''
#Initialize Chebychev pdfs
a_bkg = ROOT.RooRealVar("a_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","a_bkg",-1.,-3.,0.)
b_bkg = ROOT.RooRealVar("b_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","b_bkg",0.3,-3.,1.)
c_bkg = ROOT.RooRealVar("c_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","c_bkg",-0.01,-1.,1.)

a1_bkg = ROOT.RooRealVar("a1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","a1_bkg",-1.,-3.,0.)
b1_bkg = ROOT.RooRealVar("b1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","b1_bkg",0.3,-3.,1.)
c1_bkg = ROOT.RooRealVar("c1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","c1_bkg",-0.01,-1.,1.)
d1_bkg = ROOT.RooRealVar("d1_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","d1_bkg",-0.05,-0.2,0.1)

a2_bkg = ROOT.RooRealVar("a2_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","a2_bkg",-0.8,-15.,15.)
b2_bkg = ROOT.RooRealVar("b2_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","b2_bkg",-0.03,-1.,0.)
c2_bkg = ROOT.RooRealVar("c2_bkg_chebychev_"+CHANNEL+"_GFcat_bdt1","c2_bkg",-0.01,-1.,1.)

bkgPDF_chebychev  = ROOT.RooChebychev("chebychev_GFcat_bdt1_bkg","bkgPDF",mass,ROOT.RooArgList(a_bkg, b_bkg,c_bkg))
bkgPDF_chebychev_Nplus1 = ROOT.RooChebychev("chebychev_Nplus1_GFcat_bdt1_bkg","bkgPDF",mass,ROOT.RooArgList(a1_bkg,b1_bkg,c1_bkg,d1_bkg))
bkgPDF_chebychev_Nminus1 = ROOT.RooChebychev("chebychev_Nminus1_GFcat_bdt1_bkg","bkgPDF",mass,ROOT.RooArgList(a2_bkg, b2_bkg))

#Initialize an exponential pdf
e_bkg = ROOT.RooRealVar("e_bkg_exponential_"+CHANNEL+"_GFcat_bdt1","e_bkg",-0.031,-1.,0.)
f_bkg = ROOT.RooRealVar("f_bkg_exponential_"+CHANNEL+"_GFcat_bdt1","f_bkg",0.2,0.,0.7)
bkgPDF_exponential = ROOT.RooExponential("exponential_GFcat_bdt1_bkg","bkgPDF",mass,e_bkg)

bkgPDF_chebyplusexp = ROOT.RooAddPdf("chebyplusexp_GFcat_bdt1_bkg","bkgPDF",ROOT.RooArgSet(bkgPDF_exponential,bkgPDF_chebychev_Nminus1),ROOT.RooArgList(f_bkg))


# Create RooDataSet (Unbinned)
data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(mass), ROOT.RooFit.Import(tree))
print(f"Number of entries loaded in dataset: {data.numEntries()}")
nEntries = data.numEntries() 

#Give the blind range - just for plotting ------------------------------------------------------------------------------------------------------------------------
data_blinded = data.reduce("bosonMass < 120. || bosonMass > 130.")

#Do the fit ------------------------------------------------------------------------------------------------------------------------------
fitResult_chebychev   = bkgPDF_chebychev.fitTo(data, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Save())
fitResult_chebychev_Nplus1  = bkgPDF_chebychev_Nplus1.fitTo(data, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Save())
fitResult_exponential = bkgPDF_chebyplusexp.fitTo(data,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Save(),ROOT.RooFit.Strategy(2),ROOT.RooFit.PrintLevel(2))

#Do the F-test ------------------------------------------------------------------------------------------------------------------------------
print('################## F-TEST')
print('minNll = ', fitResult_chebychev.minNll())
print('2Delta_minNll = ', 2*(19265.2886079-fitResult_chebychev.minNll())) # If 2*(NLL(N)-NLL(N+1)) > 3.85 -> N+1 is significant improvement
print('##################')


#Plot ------------------------------------------------------------------------------------------------------------------------
canvas_chebychev = ROOT.TCanvas()
canvas_chebychev.cd()

#Chebychev frame
xframe_chebychev = mass.frame(30)

data_blinded.plotOn(xframe_chebychev)
bkgPDF_chebychev.plotOn(xframe_chebychev,ROOT.RooFit.NormRange("LowSideband,HighSideband"),ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_chebychev"),ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_chebychev.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
xframe_chebychev.GetXaxis().SetTitle("m_{M,#gamma} [GeV]")
xframe_chebychev.SetMaximum(1.3*xframe_chebychev.GetMaximum())
bkgPDF_chebychev.paramOn(xframe_chebychev,ROOT.RooFit.Layout(0.45,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_chebychev.getAttText().SetTextSize(0.02)
xframe_chebychev.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

#Calculate Chi square and parameters 
nParam_cheby = fitResult_chebychev.floatParsFinal().getSize()
chi2_cheby = xframe_chebychev.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xSidLeft-xLowRange+xHighRange-xSidRight) - nParam_cheby
cut_chi2_cheby = "{:.2f}".format(chi2_cheby*my_ndof) #Crop the chi2 to 2 decimal digits
print('Chi square cheby = ',chi2_cheby)
print('n param cheby = ',nParam_cheby)
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

#CMS_lumi.CMS_lumi(canvas_chebychev, iPeriod, iPos) #Print integrated lumi and energy information

canvas_chebychev.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/Rho_sidebands_GFcat_bdt1_chebychev.pdf")
canvas_chebychev.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/Rho_sidebands_GFcat_bdt1_chebychev.png")


#Exponential frame
canvas_exponential = ROOT.TCanvas()
canvas_exponential.cd()

xframe_exponential = mass.frame(30)

data_blinded.plotOn(xframe_exponential)
bkgPDF_chebyplusexp.plotOn(xframe_exponential,ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_chebyplusexp"), ROOT.RooFit.LineColor(ROOT.kRed))
xframe_exponential.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
xframe_exponential.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_exponential.SetMaximum(1.3*xframe_exponential.GetMaximum())
bkgPDF_chebyplusexp.paramOn(xframe_exponential,ROOT.RooFit.Layout(0.45,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_exponential.getAttText().SetTextSize(0.02)
xframe_exponential.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

nParam_exp = fitResult_exponential.floatParsFinal().getSize()
chi2_exp = xframe_exponential.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xSidLeft-xLowRange+xHighRange-xSidRight) - nParam_exp
cut_chi2_exp = "{:.2f}".format(chi2_exp*my_ndof) #Crop the chi2 to 2 decimal digits
print('Chi square exp = ',chi2_exp)
print('n param exp = ',nParam_exp)

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

#CMS_lumi.CMS_lumi(canvas_exponential, iPeriod, iPos) #Print integrated lumi and energy information

canvas_exponential.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/Rho_sidebands_GFcat_bdt1_chebyplusexp.pdf")
canvas_exponential.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/Rho_sidebands_GFcat_bdt1_chebyplusexp.png")


#Plot ------------------------------------------------------------------------------------------------------------------------
canvas_chebychev_Nplus1 = ROOT.TCanvas()
canvas_chebychev_Nplus1.cd()

#chebychev_Nplus1 frame
xframe_chebychev_Nplus1 = mass.frame(30)

data_blinded.plotOn(xframe_chebychev_Nplus1)
bkgPDF_chebychev_Nplus1.plotOn(xframe_chebychev_Nplus1,ROOT.RooFit.NormRange("LowSideband,HighSideband"),ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_chebychev_Nplus1"),ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_chebychev_Nplus1.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
xframe_chebychev_Nplus1.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_chebychev_Nplus1.SetMaximum(1.3*xframe_chebychev_Nplus1.GetMaximum())
bkgPDF_chebychev_Nplus1.paramOn(xframe_chebychev_Nplus1,ROOT.RooFit.Layout(0.45,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_chebychev_Nplus1.getAttText().SetTextSize(0.02)
xframe_chebychev_Nplus1.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

#Calculate Chi square and parameters 
nParam_cheby = fitResult_chebychev_Nplus1.floatParsFinal().getSize()
chi2_cheby = xframe_chebychev_Nplus1.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = int(xSidLeft-xLowRange+xHighRange-xSidRight) - nParam_cheby
cut_chi2_cheby = "{:.2f}".format(chi2_cheby*my_ndof) #Crop the chi2 to 2 decimal digits
print('Chi square cheby = ',chi2_cheby)
print('n param cheby = ',nParam_cheby)
print('')

leg1 = ROOT.TLegend(0.5,0.22,0.72,0.80) #right positioning
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

#CMS_lumi.CMS_lumi(canvas_chebychev_Nplus1, iPeriod, iPos) #Print integrated lumi and energy information

canvas_chebychev_Nplus1.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/Rho_sidebands_GFcat_bdt1_chebychev_Nplus1.pdf")
canvas_chebychev_Nplus1.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat1/Rho_sidebands_GFcat_bdt1_chebychev_Nplus1.png")


# Multipdf ------------------------------------------------------------------------------------------------------------------------------
a_bkg.setConstant(True);
b_bkg.setConstant(True);
c_bkg.setConstant(True);
e_bkg.setConstant(True);
a1_bkg.setConstant(True);
b1_bkg.setConstant(True);
c1_bkg.setConstant(True);
d1_bkg.setConstant(True);
a2_bkg.setConstant(True);
b2_bkg.setConstant(True);
c2_bkg.setConstant(True);
f_bkg.setConstant(True);


cat = ROOT.RooCategory("pdf_index_GFcat_bdt1","Index of Pdf which is active")
mypdfs = ROOT.RooArgList()
mypdfs.add(bkgPDF_chebychev)
#mypdfs.add(bkgPDF_bernstein)
mypdfs.add(bkgPDF_chebyplusexp)
mypdfs.add(bkgPDF_chebychev_Nplus1)

multipdf = ROOT.RooMultiPdf("multipdf_"+CHANNEL+"_GFcat_bdt1_bkg","All Pdfs",cat,mypdfs)

'''
#create Workspace ------------------------------------------------------------------------------------------------------------------------------
norm     = nEntries #fileInput.Get("h_InvMass_TwoTrk_Photon").Integral() #get the normalization of ggH signal (area under ggH signal)
bkg_norm = ROOT.RooRealVar(multipdf.GetName()+ "_norm", multipdf.GetName()+ "_norm", norm,0.5*norm, 2*norm)

print('************************************** n. events = ',nEntries)
workspace = ROOT.RooWorkspace("workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024")
#getattr(workspace,'import')(bkgPDF_chebychev)
#getattr(workspace,'import')(bkgPDF_bernstein)

#workspace = inputWS.Get("workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024")
#getattr(workspace,'import')(bkgPDF_chebychev)
#getattr(workspace,'import')(bkgPDF_bernstein)
#getattr(workspace,'import')(bkgPDF_chebyplusexp)
getattr(workspace,'import')(cat)
getattr(workspace,'import')(multipdf)
getattr(workspace,'import')(data)
getattr(workspace,'import')(bkg_norm)
print('integral BKG',bkg_norm.getVal())
#workspace.Print()

inputWS = ROOT.TFile("Rho/workspaces/workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024.root","RECREATE") #there's only one ws for both ggH and VBF 
inputWS.cd()
#fOut = ROOT.TFile("workspaces/ws_sidebands.root","RECREATE")
#fOut = ROOT.TFile("workspaces/workspace_STAT_"+CHANNEL+"_GFcat_bdt1_2024.root","UPDATE")
#fOut.cd()
workspace.Write()
#print('-------------------------------------------"
#print('Final print to check the workspace update:"
#workspace.Print()
inputWS.Close()
'''