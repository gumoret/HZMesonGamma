import ROOT

# Suppress the opening of many Canvas's 
ROOT.gROOT.SetBatch(True) 

#input
#file_path = "../histos/H/SR_HRho2022Tau_BDT_Data.root"
#file_path = "../histos/H/SR_HRho2024_ORTrigger_BDTCat1_Data.root" 
file_path = "../../histos/H/FullDataset_HRho2024_ORTrigger_BDTCat0_Data.root"
tree_name = "tree_output"
infile = ROOT.TFile.Open(file_path)
tree = infile.Get(tree_name)
nEntriesData = tree.GetEntriesFast()


mass = ROOT.RooRealVar("mesonMass","mesonMass",0.55, 0.95,"GeV")


# Meson pdf - voigtian
mass_min = 0.55
mass_max = 0.95
mean_central, mean_min, mean_max = 0.77, 0.74, 0.79
sigma_central, sigma_min, sigma_max = 0.055, 0.01, 0.2
width_central, width_min, width_max = 0.047, 0.0, 0.2

mean  = ROOT.RooRealVar("mean","The mean of the gaussian pdf", mean_central, mean_min, mean_max)
sigma = ROOT.RooRealVar("sigma","The sigma of the gaussian pdf", sigma_central, sigma_min, sigma_max)
width = ROOT.RooRealVar("width","The width of the gaussian pdf", width_central, width_min, width_max)

signalPDF = ROOT.RooGaussian("signalPDF","Gauss pdf",mass,mean,sigma)
#signalPDF = ROOT.RooVoigtian("signalPDF","Voigtian pdf",mass,mean,width,sigma)

# Construct background pdf
a1_central, a1_min, a1_max = 0.3,-0.5,0.5
a2_central, a2_min, a2_max = -0.1,-2.,2.
a3_central, a3_min, a3_max = -0.02,-2.,2.

a1_off = ROOT.RooRealVar("a1_off","The a1 of background", a1_central, a1_min, a1_max)
a2_off = ROOT.RooRealVar("a2_off","The a2 of background", a2_central, a2_min, a2_max)
a3_off = ROOT.RooRealVar("a3_off","The a3 of background", a3_central, a3_min, a3_max)

backgroundPDF = ROOT.RooChebychev("backgroundPDF","The background PDF",mass,ROOT.RooArgList(a1_off,a2_off,a3_off))

# N events
# ---------------------------------------------------------------
data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(mass), ROOT.RooFit.Import(tree))

nData = data.numEntries()
print(f"Number of entries loaded in dataset: {nData}")

nsig = ROOT.RooRealVar("nsig", "signal yield",0.2*nData, 0.,nData)
nbkg = ROOT.RooRealVar("nbkg", "background yield",0.8*nData,0.,nData)


# Create the model
# ---------------------------------------------------------------
model = ROOT.RooAddPdf("model", "Signal and Background PDF", ROOT.RooArgList(signalPDF, backgroundPDF), ROOT.RooArgList(nsig, nbkg))


# Create RooDataSet (Unbinned)
data = ROOT.RooDataSet("data", "Data", ROOT.RooArgSet(mass), ROOT.RooFit.Import(tree))
print(f"Number of entries loaded in dataset: {data.numEntries()}")


#Do the fit ------------------------------------------------------------------------------------------------------------------------------
#fitResult = model.fitTo(data, ROOT.RooFit.Save())
fitResult = model.fitTo(data,ROOT.RooFit.Save(),ROOT.RooFit.Strategy(2))


#Plot ------------------------------------------------------------------------------------------------------------------------
canvas = ROOT.TCanvas()
canvas.cd()

#frame
xframe = mass.frame(30)

data.plotOn(xframe)
model.plotOn(xframe)
xframe.SetTitle("#sqrt{s} = 13.6 TeV       lumi = 108.95/fb")
xframe.GetXaxis().SetTitle("m_{M} [GeV]")
xframe.SetMaximum(1.3*xframe.GetMaximum())
model.paramOn(xframe,ROOT.RooFit.Layout(0.45,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe.getAttText().SetTextSize(0.02)
xframe.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

#Calculate Chi square and parameters 
nParam = fitResult.floatParsFinal().getSize()
chi2 = xframe.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
my_ndof = 40 - nParam
cut_chi2 = "{:.2f}".format(chi2*my_ndof) #Crop the chi2 to 2 decimal digits
print('Chi square = ',chi2)
print('n param = ',nParam)
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
leg1.AddEntry(cut_chi2,"#chi^{2}/ndof = " + cut_chi2 + " / " + str(my_ndof),"brNDC")

leg1.Draw()

#CMS_lumi.CMS_lumi(canvas_chebychev, iPeriod, iPos) #Print integrated lumi and energy information

canvas.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit1d_meson_bdt0.pdf")
canvas.SaveAs("/eos/user/e/eferrand/www/eferrand/HZMesonGamma/HRho/2024/mass_analysis/cat0/fit1d_meson_bdt0.png")