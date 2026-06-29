import ROOT

# Suppress the opening of many Canvas's 
ROOT.gROOT.SetBatch(True) 

#input
file_path = "./histos/H/SR_HRho2022Tau_BDT_Data.root"
tree_name = "tree_output"
infile = ROOT.TFile.Open(file_path)
tree = infile.Get(tree_name)


#Parameters of the PDF ---------------------------------------------------------------
xLowRange  = 110.
xHighRange = 160.

mass = ROOT.RooRealVar("bosonMass","bosonMass",xLowRange,xHighRange,"GeV")
mass.setRange("LowSideband",xLowRange,120.)
mass.setRange("HighSideband",130.,xHighRange)
mass.setRange("full",xLowRange,xHighRange)


#Initialize a Chebychev pdf
a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",-1.,-3.,0.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.3,-3.,1.)
c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",-0.01,-1.,1.)
d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",-0.05,-0.2,0.1)
f_bkg = ROOT.RooRealVar("e_bkg","e_bkg",-0.05,-0.1,0.)

bkgPDF_chebychev  = ROOT.RooChebychev("chebychev","bkgPDF",mass,ROOT.RooArgList(a_bkg, b_bkg))


# Create RooDataSet (Unbinned)
data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(mass), ROOT.RooFit.Import(tree))
print(f"Number of entries loaded in dataset: {data.numEntries()}")

#Give the blind range - just for plotting ------------------------------------------------------------------------------------------------------------------------
data_blinded = data.reduce("bosonMass < 120. || bosonMass > 130.")

#Do the fit ------------------------------------------------------------------------------------------------------------------------------
fitResult_chebychev = bkgPDF_chebychev.fitTo(data, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Save())


#Plot ------------------------------------------------------------------------------------------------------------------------
canvas_chebychev = ROOT.TCanvas()
canvas_chebychev.cd()

#Chebychev frame
xframe_chebychev = mass.frame(60)

data_blinded.plotOn(xframe_chebychev)#, ROOT.RooFit.CutRange("LowSideband,HighSideband"))
bkgPDF_chebychev.plotOn(xframe_chebychev,ROOT.RooFit.NormRange("LowSideband,HighSideband"),ROOT.RooFit.Range("LowSideband,HighSideband"),ROOT.RooFit.Name("bkgPDF_chebychev"),ROOT.RooFit.LineColor(ROOT.kBlue))

#Calculate Chi square and parameters 
nParam_cheby = fitResult_chebychev.floatParsFinal().getSize()
chi2_cheby = xframe_chebychev.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2_cheby = "{:.2f}".format(chi2_cheby) #Crop the chi2 to 2 decimal digits
print("Chi square cheby = ",chi2_cheby)
print("n param cheby = ",nParam_cheby)
print()

xframe_chebychev.Draw()
canvas_chebychev.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/mass_analysis/fit1d_higgs.pdf")
canvas_chebychev.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/mass_analysis/fit1d_higgs.png")