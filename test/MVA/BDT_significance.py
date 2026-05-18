import ROOT
import math
import numpy as np
import argparse

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

#Choose an injected BR (since you moved to BR = 1 here you can set it at 10^-5 to make comparison with previous productions)
BR_inj = 0.00001

p = argparse.ArgumentParser(description='Select rootfile')
p.add_argument('meson_option', help='Type <<rho>> for rho, <<phi>> for phi') #flag for type of meson
args = p.parse_args()

#INPUT FILES
if args.meson_option == "HRho" :
    fInputSignal = ROOT.TFile("../histos/H/SR_HRho_OldTrigger_preselection_Signal.root")####modified for 2022Tau
    fInputData  = ROOT.TFile("../histos/H/SR_HRho2022Tau_preselection_Data.root")    

elif args.meson_option == "HPhi" :
    fInputSignal = ROOT.TFile("../histos/H/SR_HPhi_OldTrigger_preselection_Signal.root") #SR_HPhi_ORTrigger_preselection_Signal.root
    fInputData  = ROOT.TFile("../histos/H/SR_HPhi2022Tau_preselection_Data.root") #SR_HPhi2024_ORTrigger_preselection_Data.root

elif args.meson_option == "HKst" :
    fInputSignal = ROOT.TFile("../histos/H/SR_HKst_OldTrigger_preselection_Signal.root")
    fInputData  = ROOT.TFile("../histos/H/SR_HKst2022Tau_preselection_Data.root")


elif args.meson_option == "HDst" :
    fInputSignal = ROOT.TFile("../histos/H/SR_HDst_OldTrigger_preselection_Signal.root")
    fInputData  = ROOT.TFile("../histos/H/SR_HDst2022Tau_preselection_Data.root")

sig_tree    = fInputSignal.Get("tree_output")
data_tree   = fInputData.Get("tree_output")
h_mesonMass = fInputData.Get("h_mesonMass")


Nbkg_passed = 0
nEntriesData = int(h_mesonMass.GetEntries())

for dataEntry in range(nEntriesData):
    ientry = data_tree.LoadTree( dataEntry )  
    if ientry < 0: break
    nb = data_tree.GetEntry(dataEntry )
    if nb <= 0:
        print("nb < 0")
        continue  

    bosonMass = data_tree.bosonMass
    
    if not args.meson_option == "HDst":
        if (bosonMass < 120. or bosonMass > 130.): continue
    elif args.meson_option == "HDst":
        if (bosonMass < 114. or bosonMass > 126.): continue
    Nbkg_passed += 1

Nsig_passed = 0 #Initialize the number of signal events from the sum of the weights (before applying BDT cuts), take this number running the generate_histos for signal
nEntriesSig = sig_tree.GetEntriesFast()

for sigEntry in range(nEntriesSig):
    ientry = sig_tree.LoadTree( sigEntry )  
    if ientry < 0: break
    nb = sig_tree.GetEntry(sigEntry )
    if nb <= 0:
        print("nb < 0")
        continue
        
    bosonMass = sig_tree.bosonMass
    if not args.meson_option == "HDst":
        if (bosonMass < 120. or bosonMass > 130.): continue
    elif args.meson_option == "HDst":
        if (bosonMass < 114. or bosonMass > 126.): continue
    #print Zmass
    Nsig_passed += sig_tree._eventWeight * BR_inj
    #Nsig_passed += 1. * BR_inj


print("nEntriesSig  = ",nEntriesSig)
print("Nsig_passed  = ",Nsig_passed)
print("nEntriesData = ",nEntriesData)
print("Nbkg_passed  = ",Nbkg_passed)


BDT_file = ROOT.TFile("outputs/Nominal_training.root")

h_BDT_effB_effS = BDT_file.Get("default/Method_BDT/BDT/MVA_BDT_effBvsS")

canvas1 = ROOT.TCanvas()
sig_eff = []
bkg_eff = []
signif  = []
_effS   = 0

for jbin in range(1,h_BDT_effB_effS.GetNbinsX()+1):
    if h_BDT_effB_effS.GetBinCenter(jbin) > 0.2: #insert the number before whom the function fluctuates too much to estimate the maximum
        sig_eff.append(h_BDT_effB_effS.GetBinCenter(jbin))
        if h_BDT_effB_effS.GetBinContent(jbin) <= 0.:
            bkg_eff.append(0.)
            signif.append(0)
        else:
            bkg_eff.append(h_BDT_effB_effS.GetBinContent(jbin))
            signif.append(Nsig_passed*h_BDT_effB_effS.GetBinCenter(jbin)/math.sqrt(Nbkg_passed*h_BDT_effB_effS.GetBinContent(jbin)))


sig_eff_array = np.array(sig_eff)
bkg_eff_array = np.array(bkg_eff)
signif_array = np.array(signif)
#print "signif_len: ", len(signif_array)
sign = ROOT.TGraph(70,sig_eff_array,signif_array)
sign.SetTitle("")
sign.GetXaxis().SetTitle("#varepsilon_{S}^{BDT}")
sign.GetYaxis().SetTitle("Significance")
sign.SetMinimum(0.)
sign.SetMaximum(2.*ROOT.TMath.MaxElement(sign.GetN(),sign.GetY()))
sign.SetMarkerStyle(8)
sign.SetMarkerColor(4)
sign.GetXaxis().SetRangeUser(0.1,1.)

sign.Draw("AP")

canvas2 = ROOT.TCanvas()
sign_vs_bkg = ROOT.TGraph(70,bkg_eff_array,signif_array)
sign_vs_bkg.SetTitle("")
sign_vs_bkg.GetXaxis().SetTitle("#varepsilon_{B}^{BDT}")
sign_vs_bkg.GetYaxis().SetTitle("Significance")
sign_vs_bkg.GetXaxis().SetRangeUser(0.,0.8)
sign_vs_bkg.SetMinimum(0.)
sign_vs_bkg.SetMaximum(2.*ROOT.TMath.MaxElement(sign_vs_bkg.GetN(),sign_vs_bkg.GetY()))
sign_vs_bkg.SetMarkerStyle(8)
sign_vs_bkg.SetMarkerColor(4)
sign_vs_bkg.Draw("AP")

if args.meson_option == "HRho" :
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Rho/MVA/signif_vs_effS.pdf")
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Rho/MVA/signif_vs_effS.png")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Rho/MVA/signif_vs_effB.pdf")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Rho/MVA/signif_vs_effB.png")
elif args.meson_option == "HPhi" :
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Phi/MVA/signif_vs_effS.pdf")
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Phi/MVA/signif_vs_effS.png")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Phi/MVA/signif_vs_effB.pdf")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Phi/MVA/signif_vs_effB.png")
elif args.meson_option == "HKst" :
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Kst/MVA/signif_vs_effS.pdf")
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Kst/MVA/signif_vs_effS.png")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Kst/MVA/signif_vs_effB.pdf")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Kst/MVA/signif_vs_effB.png")
elif args.meson_option == "HDst" :
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Dst/MVA/signif_vs_effS.pdf")
    canvas1.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Dst/MVA/signif_vs_effS.png")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Dst/MVA/signif_vs_effB.pdf")
    canvas2.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/plots/H/Dst/MVA/signif_vs_effB.png")

#---- Now find the BDT output corresponding to the highest significance

h_BDT_effS = BDT_file.Get("default/Method_BDT/BDT/MVA_BDT_effS")
signif_maximizing_eff = sig_eff_array[np.argmax(signif_array)]
print("signif_maximizing_eff = ",signif_maximizing_eff)

BDT_output = 0. 
t = 1.
for entry in range(h_BDT_effS.GetNbinsX()):
    effS = h_BDT_effS.GetBinContent(entry)
    effS = float(format(effS, '.4f'))
    signif_maximizing_eff = float(format(signif_maximizing_eff, '.4f'))
    #print "effS: ", effS, "signif_max_eff: ", signif_maximizing_eff
    '''
    if effS == signif_maximizing_eff:
        BDT_output =  h_BDT_effS.GetBinCenter(entry)
        _effS = effS
    '''
    if abs(effS - signif_maximizing_eff) < t :#and not effS == 0:
        t = abs(effS - signif_maximizing_eff)
        BDT_output =  h_BDT_effS.GetBinCenter(entry)
        _effS = effS

print("t = " ,t)
print("RESULT")
print("-------------------------------------------------------------")
print("Nsig_passed       = " ,Nsig_passed)
print("Nbkg_passed       = " ,Nbkg_passed)
print("Signal efficiency = ", _effS)
print("Significance:   Z = ",signif_array[np.argmax(signif_array)])
print("BDT output        = ", BDT_output)
print("-------------------------------------------------------------")

# Create a ROOT file to store the tree
output_file = ROOT.TFile("BDToutput.root", "RECREATE")

# Create the tree and add a branch for the BDT output variable
tree = ROOT.TTree("BDTtree", "Tree with BDT output")

output_file.cd()

_BDT_output = np.zeros(1, dtype=float)
tree.Branch("_BDT_output",_BDT_output, "_BDT_output/D")
_BDT_output[0] = BDT_output
tree.Fill()
# Save the tree to the output file and close the file
output_file.Write()
output_file.Close()

#raw_input()