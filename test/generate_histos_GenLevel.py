import ROOT
import argparse
import math
import numpy as np
import sys
from array import array

#Following bools are given as input
isPhiAnalysis = False # for Z -> Phi Gamma
isRhoAnalysis = False # for Z -> Rho Gamma
isKAnalysis   = False # for Z -> K* Gamma
isD0Analysis  = False # for Z -> D0* Gamma

isHAnalysis   = False
isZAnalysis   = False

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('boson_channel', help='type <<H>> or <<Z>>')
p.add_argument('meson_channel', help='type <<rho>> or <<phi>> or <<K>> or <<D0>>')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("tree_output")

if args.boson_channel == "H": isHAnalysis = True 
elif args.boson_channel == "Z": isZAnalysis = True 
if args.meson_channel == "phi": isPhiAnalysis = True 
if args.meson_channel == "rho": isRhoAnalysis = True 
if args.meson_channel == "K": isKAnalysis = True 
if args.meson_channel == "D0": isDAnalysis = True 


#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_bosonMass", "h_mesonMass", "h_firstTrkPt", "h_secondTrkPt", "h_firstTrkEta", "h_secondTrkEta", "h_firstTrkPhi", "h_secondTrkPhi", "h_mesonPt", "h_mesonEta", "h_mesonPhi", "h_photonPt", "h_photonEta", "h_photonPhi", "h_bosonPt", "h_bosonEta", "h_bosonPhi"] 

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{boson}", 300, 50., 200.) 
if   isPhiAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 1., 1.05) 
elif isRhoAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.5, 1.)
elif isKAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.6, 1.3) 
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T} of the 1st track", 100, 0.,70.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 0., 50.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"#eta of the 1st track", 100, -2.5, 2.5)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 2nd track", 100, -2.5, 2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#phi of the 1st track", 100, -math.pi, math.pi)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 2nd track", 100, -math.pi, math.pi)
histo_map[list_histos[8]] = ROOT.TH1F(list_histos[8],"p_{T} of the meson", 100, 0., 140.)
histo_map[list_histos[9]] = ROOT.TH1F(list_histos[9],"#eta_{meson}", 100, -2.5,2.5)
histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#phi_{meson}", 100, -math.pi, math.pi)
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"p_{T} of the #gamma", 100, 0., 140.)
histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"#eta_{#gamma}", 100, -2.5,2.5)
histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"#phi_{gamma}", 100, -math.pi, math.pi)
histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"p_{T} of boson", 100, 0., 250.)
histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"#eta_{boson}", 100, -2.5,2.5)
histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"#phi_{boson}", 100, -math.pi, math.pi)

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
        print( "nb < 0")
        continue

    #Retrieve variables from the tree 
    bosonMass      = mytree.genBoson_mass
    mesonMass      = mytree.genMeson_mass
    firstTrkPt     = mytree.genTrack1_pT      
    secondTrkPt    = mytree.genTrack2_pT
    firstTrkEta    = mytree.genTrack1_eta
    secondTrkEta   = mytree.genTrack2_eta
    firstTrkPhi    = mytree.genTrack1_phi
    secondTrkPhi   = mytree.genTrack2_phi
    mesonPt        = mytree.genMeson_pT      
    mesonEta       = mytree.genMeson_eta
    mesonPhi       = mytree.genMeson_phi
    photonPt       = mytree.genGamma_pT
    photonEta      = mytree.genGamma_eta
    photonPhi      = mytree.genGamma_phi
    bosonPt        = mytree.genBoson_pT
    bosonEta       = mytree.genBoson_eta
    bosonPhi       = mytree.genBoson_phi
    

    eventWeight = 1.
    #FILL HISTOS #####################################################################################################
    
    values_to_check = [bosonMass, mesonMass, mesonPt, mesonEta, mesonPhi, photonPt, photonEta, photonPhi, bosonPt, bosonEta, bosonPhi]
    isValid_values = False

    if all(value > -999 for value in values_to_check): isValid_values = True

    if isValid_values and not isKAnalysis:
        histo_map["h_bosonMass"].Fill(bosonMass, eventWeight)          
        histo_map["h_mesonMass"].Fill(mesonMass, eventWeight)
        histo_map["h_firstTrkPt"].Fill(firstTrkPt, eventWeight)
        histo_map["h_secondTrkPt"].Fill(secondTrkPt, eventWeight)
        histo_map["h_firstTrkEta"].Fill(firstTrkEta, eventWeight)    
        histo_map["h_secondTrkEta"].Fill(secondTrkEta, eventWeight)   
        histo_map["h_firstTrkPhi"].Fill(firstTrkPhi, eventWeight)    
        histo_map["h_secondTrkPhi"].Fill(secondTrkPhi, eventWeight)
        histo_map["h_mesonPt"].Fill(mesonPt, eventWeight)
        histo_map["h_mesonEta"].Fill(mesonEta, eventWeight)
        histo_map["h_mesonPhi"].Fill(mesonPhi, eventWeight)
        histo_map["h_photonPt"].Fill(photonPt, eventWeight)
        histo_map["h_photonEta"].Fill(photonEta, eventWeight)
        histo_map["h_photonPhi"].Fill(photonPhi, eventWeight)
        histo_map["h_bosonPt"].Fill(bosonPt, eventWeight)
        histo_map["h_bosonEta"].Fill(bosonEta, eventWeight)
        histo_map["h_bosonPhi"].Fill(bosonPhi, eventWeight)

    elif isKAnalysis:
        histo_map["h_bosonMass"].Fill(bosonMass, eventWeight)          
        histo_map["h_mesonMass"].Fill(mesonMass, eventWeight)
        histo_map["h_firstTrkPt"].Fill(firstTrkPt, eventWeight)
        histo_map["h_secondTrkPt"].Fill(secondTrkPt, eventWeight)
        histo_map["h_firstTrkEta"].Fill(firstTrkEta, eventWeight)    
        histo_map["h_secondTrkEta"].Fill(secondTrkEta, eventWeight)   
        histo_map["h_firstTrkPhi"].Fill(firstTrkPhi, eventWeight)    
        histo_map["h_secondTrkPhi"].Fill(secondTrkPhi, eventWeight)
        histo_map["h_mesonPt"].Fill(mesonPt, eventWeight)
        histo_map["h_mesonEta"].Fill(mesonEta, eventWeight)
        histo_map["h_mesonPhi"].Fill(mesonPhi, eventWeight)
        histo_map["h_photonPt"].Fill(photonPt, eventWeight)
        histo_map["h_photonEta"].Fill(photonEta, eventWeight)
        histo_map["h_photonPhi"].Fill(photonPhi, eventWeight)
        histo_map["h_bosonPt"].Fill(bosonPt, eventWeight)
        histo_map["h_bosonEta"].Fill(bosonEta, eventWeight)
        histo_map["h_bosonPhi"].Fill(bosonPhi, eventWeight)


        

#HISTO LABELS #########################################################################################################
histo_map["h_bosonMass"].GetXaxis().SetTitle("m_{boson} [GeV/c^2]")
histo_map["h_bosonMass"].SetTitle("boson invariant mass")
histo_map["h_mesonMass"].GetXaxis().SetTitle("m_{meson } [GeV/c^2]")
histo_map["h_mesonMass"].SetTitle("Meson mass")
histo_map["h_firstTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_firstTrkPt"].SetTitle("Transverse momentum of the first charged particle ")
histo_map["h_secondTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_secondTrkPt"].SetTitle("Transverse momentum of the second charged particle")
histo_map["h_firstTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_firstTrkEta"].SetTitle("Pseudorapidity of the first charged particle")
histo_map["h_secondTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_secondTrkEta"].SetTitle("Pseudorapidity of the second charged particle")
histo_map["h_firstTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_firstTrkPhi"].SetTitle("Azimuthal angle of the first charged particle")
histo_map["h_secondTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_secondTrkPhi"].SetTitle("Azimuthal angle of the second charged particle")
histo_map["h_mesonPt"].GetXaxis().SetTitle("pT_{meson} [GeV/c]") 
histo_map["h_mesonPt"].SetTitle("Transverse momentum of the meson")
histo_map["h_mesonEta"].GetXaxis().SetTitle("#eta")
histo_map["h_mesonEta"].SetTitle("Pseudorapidity of the meson")
histo_map["h_mesonPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_mesonPhi"].SetTitle("Azimuthal angle of the meson")
histo_map["h_photonPt"].GetXaxis().SetTitle("p_{T}^{#gamma} [GeV/c]")
histo_map["h_photonPt"].SetTitle("Transverse momentum of the photon")
histo_map["h_photonEta"].GetXaxis().SetTitle("#eta_{#gamma}")
histo_map["h_photonEta"].SetTitle("Eta of the photon")
histo_map["h_photonPhi"].GetXaxis().SetTitle("#phi_{#gamma} [rad]")
histo_map["h_photonPhi"].SetTitle("Phi of the photon")
histo_map["h_bosonPt"].GetXaxis().SetTitle("pT_{boson} [GeV/c]")
histo_map["h_bosonPt"].SetTitle("Transverse momentum of the boson")
histo_map["h_bosonEta"].GetXaxis().SetTitle("#eta_{boson}")
histo_map["h_bosonEta"].SetTitle("Eta of the Z")
histo_map["h_bosonPhi"].GetXaxis().SetTitle("#phi_{boson} [rad]")
histo_map["h_bosonPhi"].SetTitle("Phi of the boson")


#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()    