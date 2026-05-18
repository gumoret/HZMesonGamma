import ROOT
import os
import argparse

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile')
p.add_argument('meson_option', help='Type <<rho>> for rho, <<phi>> for phi') #flag for type of meson
args = p.parse_args()


#INPUT FILES
if args.meson_option == "HRho" :
    fIn_bkg  = ROOT.TFile("../histos/H/CR_HRho2022Tau_preselection_Sidebands.root")####modified for 2022Tau
    fIn_sig  = ROOT.TFile("../histos/H/SR_HRho_OldTrigger_preselection_Signal.root")

elif args.meson_option == "HPhi" :
    fIn_bkg  = ROOT.TFile("../histos/H/CR_HPhi2022Tau_preselection_Sidebands.root") #CR_HPhi2024_ORTrigger_preselection_Sidebands.root
    fIn_sig  = ROOT.TFile("../histos/H/SR_HPhi_OldTrigger_preselection_Signal.root") #SR_HPhi_ORTrigger_preselection_Signal.root

elif args.meson_option == "HKst" :
    fIn_bkg  = ROOT.TFile("../histos/H/CR_HKst2022Tau_preselection_Sidebands.root") #CR_HKst2024_ORTrigger_preselection_Sidebands.root
    fIn_sig  = ROOT.TFile("../histos/H/SR_HKst_OldTrigger_preselection_Signal.root") #SR_HKst_ORTrigger_preselection_Signal.root

elif args.meson_option == "HDst" :
    fIn_bkg  = ROOT.TFile("../histos/H/CR_HDst2022Tau_preselection_Sidebands.root")
    fIn_sig  = ROOT.TFile("../histos/H/SR_HDst_OldTrigger_preselection_Signal.root")

tree_bkg = fIn_bkg.Get("tree_output")
tree_sig = fIn_sig.Get("tree_output")   

#OUTPUT FILE
fOut = ROOT.TFile("outputs/Nominal_training.root","RECREATE")

#START MVA
ROOT.TMVA.Tools.Instance()

factory = ROOT.TMVA.Factory("TMVAClassification", fOut,":".join(["!V","Transformations=I;D;P;G,D","AnalysisType=Classification"]))

dataloader = ROOT.TMVA.DataLoader()

#VARIABLES FROM THE TREE -------------------------------------------------------------------------------

#first track:  firstTrkPt   firstTrkEta    firstTrkPhi
#second track: secondTrkPt  secondTrkEta   secondTrkPhi
#meson:   mesonPt   mesonIso (charged)   _bestCoupleEta   trksDeltaR  mesonMass (using mesonMass variable is not correct since you use it to define the bkg estimation, use it only to compute scatter plots to see correlation with other variables)
#photon:   photonEt  photonEta
#invariant mass: bosonMass

#Pay attention to the order, it must be the same in the function_smuggler.py

#dataloader.AddVariable("firstTrkIsoCh","F")###old variable
dataloader.AddVariable("mesonIsoCh","F")
dataloader.AddVariable("mesonIsoNeu","F")###old variable
dataloader.AddVariable("mesonPt/bosonMass","F")
dataloader.AddVariable("mesonEta","F")###old variable
#dataloader.AddVariable("mesonGammaMass","F") #used just to check the correlation of Hmass and other input variables
#dataloader.AddVariable("_HpT","F") #used just to check the correlation of Hmass and other input variables
dataloader.AddVariable("photonEt/bosonMass","F")
#dataloader.AddVariable("secondTrkPt","F")
dataloader.AddVariable("mesonGammaDeltaPhi","F")
if args.meson_option == "HDst": dataloader.AddVariable("lxy","F")

##-------------------------------------------------------------------------------------------------------

sig_weight = 1.
bkg_weight = 1.

print("before AddSignalTree")
dataloader.AddSignalTree(tree_sig, sig_weight)
print("before AddBackgroundTree")
dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

###############dataloader.SetWeightExpression("_BDTweight") #_BDTweight is the weight variable of the tree

mycutSig = ROOT.TCut("")
mycutBkg = ROOT.TCut("")


dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"]))
#dataloader.PrepareTrainingAndTestTree(mycutBkg, ":".join(["!V","nTrain_Background=0:nTest_Background=0"]))

method_btd    = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=600:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:UseRandomisedTrees:UseNvars=4:UseBaggedBoost:BaggedSampleFraction=0.5"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()

#weightfile_dir = "default/weights/TMVAClassification_BDT.weights.xml"

#if evaluate_BDT_systematic:
    # weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_shifted.xml"
    # weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_shifted.xml"
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_secondPart.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_secondPart.xml"
#elif test_on_signal_shifted_up:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_up.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_up.xml"
#elif test_on_signal_shifted_down:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_down.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_down.xml"
#elif test_on_signal_sin2:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_sin2.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_sin2.xml"
#elif test_on_signal_cos:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_cos_with_flat_theta.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_cos_with_flat_theta.xml"
#else:
    #weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu.xml"
    #weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele.xml"
    #weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Wmass.xml"
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_firstPart.xml"
    #weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Wmass.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_firstPart.xml"

#if isMuon:
 #   rename_weightfile = "mv " + weightfile_dir + " " + weightfile_mu
  #  os.system(rename_weightfile)
#else:
 #   rename_weightfile = "mv " + weightfile_dir + " " + weightfile_ele
  #  os.system(rename_weightfile)