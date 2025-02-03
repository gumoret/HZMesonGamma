import ROOT
import argparse
import math
import numpy as np
import sys
from array import array

#Following bools are given as input
#isDataBlind = False #Bool for blind analysis
isPhiAnalysis  = False # for Z -> Phi Gamma
isRhoAnalysis  = True # for Z -> Rho Gamma

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("tree_output")


#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = [ "h_MesonMassTrkTrk", "h_MesonTreeMass", "h_MesonTreeMassKin", "h_MesonPt", "h_MesonTreePt", "h_MesonEta", "h_MesonTreeEta","h_MesonPhi", "h_MesonTreePhi" ] 

if isPhiAnalysis:
    histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{TrkTrk}", 100, 1., 1.05)
    histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{mesonTree}", 100, 1., 1.05)
    histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"M_{mesonTreeKin}", 100, 1., 1.05)
    
    
elif isRhoAnalysis: 
    histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{mesonTrkTrk}", 100, 0.5, 1.) 
    histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{mesonTreeMass}", 100, 0.5, 1.) 
    histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"M_{mesonTreeMassKin}", 100, 0.5, 1.) 

histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"pT_{trktrk}", 100, 35., 110.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"pT_{mesonTree}", 100, 35., 110.)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta_{trktrk}", 100, -2.5, 2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#eta_{mesonTree}", 100, -2.5, 2.5) 
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi_{trktrk}", 100, -math.pi, math.pi)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"#phi_{mesonTree}", 100, -math.pi, math.pi)


#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()



#EVENTS LOOP ########################################################################################################
nentries = mytree.GetEntriesFast()
for jentry in range(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry)
    if nb <= 0:
        print("nb < 0")
        continue

    #Retrieve variables from the tree     
    mesonMass          = mytree.mesonMassTrkTrk
    mesonTreeMass      = mytree.mesonTreeMass
    mesonTreeMassKin   = mytree.mesonTreeMassKin
    mesonPt            = mytree.mesonPt
    mesonTreePt        = mytree.mesonTreePt
    mesonEta           = mytree.mesonEta 
    mesonTreeEta       = mytree.mesonTreeEta
    mesonPhi           = mytree.mesonPhi
    mesonTreePhi       = mytree.mesonTreePhi

    

    eventWeight = 1.
    #FILL HISTOS #####################################################################################################
    #if DATA -> Blind Analysis on Z inv mass plot    
    histo_map["h_MesonMassTrkTrk"].Fill(mesonMass, eventWeight)
    histo_map["h_MesonTreeMass"].Fill(mesonTreeMass, eventWeight)
    histo_map["h_MesonTreeMassKin"].Fill(mesonTreeMassKin, eventWeight)
    histo_map["h_MesonPt"].Fill(mesonPt, eventWeight)
    histo_map["h_MesonTreePt"].Fill(mesonTreePt, eventWeight)
    histo_map["h_MesonEta"].Fill(mesonEta, eventWeight)
    histo_map["h_MesonTreeEta"].Fill(mesonTreeEta, eventWeight)
    histo_map["h_MesonPhi"].Fill(mesonPhi, eventWeight)
    histo_map["h_MesonTreePhi"].Fill(mesonTreePhi, eventWeight)



#HISTO LABELS #########################################################################################################
histo_map["h_MesonMassTrkTrk"].GetXaxis().SetTitle("m_{TrkTrk } [GeV/c^2]")
histo_map["h_MesonMassTrkTrk"].SetTitle("Meson mass from tracks")

histo_map["h_MesonTreeMass"].GetXaxis().SetTitle("m_{tree } [GeV/c^2]")
histo_map["h_MesonTreeMass"].SetTitle("Meson mass from tree")

histo_map["h_MesonTreeMassKin"].GetXaxis().SetTitle("m_{treeKin } [GeV/c^2]")
histo_map["h_MesonTreeMassKin"].SetTitle("Meson kin mass from tree")

histo_map["h_MesonPt"].GetXaxis().SetTitle("pT_{trktrk } [GeV/c]")
histo_map["h_MesonPt"].SetTitle("Meson pT")

histo_map["h_MesonTreePt"].GetXaxis().SetTitle("pT_{treeKin } [GeV/c]")
histo_map["h_MesonTreePt"].SetTitle("Meson pT from tree")

histo_map["h_MesonEta"].GetXaxis().SetTitle("#eta_{trktrk } ")
histo_map["h_MesonEta"].SetTitle("Meson pseudorapidity") 

histo_map["h_MesonTreeEta"].GetXaxis().SetTitle("#eta_{treeKin } ")
histo_map["h_MesonTreeEta"].SetTitle("Meson pseudorapidity from tree")

histo_map["h_MesonPhi"].GetXaxis().SetTitle("#phi_{trktrk } [rad]")
histo_map["h_MesonPhi"].SetTitle("Meson Azimuthal angle")

histo_map["h_MesonTreePhi"].GetXaxis().SetTitle("#phi_{treeKin } [rad]")
histo_map["h_MesonTreePhi"].SetTitle("Meson Azimuthal angle from tree")



#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()