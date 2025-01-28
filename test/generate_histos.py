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
list_histos = [ "h_MesonMassTrkTrk", "h_MesonTreeMass" ] 

if isPhiAnalysis:
    histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{mesonTrkTrk}", 100, 1., 1.05)
    histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{mesonTreeMass}", 100, 1., 1.05)

elif isRhoAnalysis: 
    histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{mesonTrkTrk}", 100, 0.5, 1.) 
    histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{mesonTreeTrk}", 100, 0.5, 1.) 


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
    
    mesonMassTrkTrk    = mytree.mesonMassTrkTrk
    mesonTreeMass      = mytree.mesonTreeMass

    

    eventWeight = 1.
    #FILL HISTOS #####################################################################################################
    #if DATA -> Blind Analysis on Z inv mass plot
    
    histo_map["h_MesonMassTrkTrk"].Fill(mesonMassTrkTrk, eventWeight)
    histo_map["h_MesonTreeMass"].Fill(mesonTreeMass, eventWeight)



#HISTO LABELS #########################################################################################################
histo_map["h_MesonMassTrkTrk"].GetXaxis().SetTitle("m_{mesonTrkTrk } [GeV/c^2]")
histo_map["h_MesonMassTrkTrk"].SetTitle("Meson mass from tracks")

histo_map["h_MesonTreeMass"].GetXaxis().SetTitle("m_{mesonTree } [GeV/c^2]")
histo_map["h_MesonTreeMass"].SetTitle("Meson mass from tree")




#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()

    