import ROOT
import argparse
import math
import numpy as np
import sys
from array import array

#Following bools are given as input
#isDataBlind = False #Bool for blind analysis
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
p.add_argument('meson_channel', help='type <<rho>> or <<phi>> or <<K*>> or <<D0*>>')
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
if args.meson_channel == "K*": isKAnalysis = True 
if args.meson_channel == "D0*": isDAnalysis = True 




#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_bosonMass", "h_mesonMass", "h_firstTrkPt", "h_secondTrkPt", "h_firstTrkEta", "h_secondTrkEta", "h_firstTrkPhi", "h_secondTrkPhi", "h_mesonPt", "h_mesonEta", "h_trksDeltaR","h_mesonIsoCh", "h_photonEnergy", "h_photonEta","h_nMuons","h_nElectrons"]  

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{boson}", 300, 50., 200.)   #120,60,120 
if isPhiAnalysis : histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 1., 1.05) 
elif isRhoAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.5, 1.)
elif isKAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.6, 1.3)
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T} of the 1st track", 100, 0.,70.)
histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 0.,55.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"#eta of the 1st track", 100, -2.5,2.5)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 2nd track", 100, -2.5,2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#phi of the 1st track", 100, -math.pi, math.pi)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 2nd track", 100, -math.pi, math.pi)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"p_{T} of the meson", 100, 0.,140.)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"#eta_{meson}", 100, -2.5,2.5)
histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#Delta R_{meson}", 100, 0.,0.07)
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"Iso_ch of the meson", 100, 0., 1.)
histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"E_{T} of the #gamma", 100, 0., 250.)
histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"#eta_{#gamma}", 100, -2.5,2.5)
histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"n. of muons", 6, -0.5, 5.5)
histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"n. of electrons", 5, -0.5, 5.5)


#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()


#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree ################################################################################
_bosonMass    = np.zeros(1, dtype=float)  
_mesonMass    = np.zeros(1, dtype=float)
_firstTrkPt   = np.zeros(1, dtype=float)
_secondTrkPt  = np.zeros(1, dtype=float)
_firstTrkEta  = np.zeros(1, dtype=float)  
_secondTrkEta = np.zeros(1, dtype=float)
_mesonPt      = np.zeros(1, dtype=float)
_mesonEta     = np.zeros(1, dtype=float)  
_mesonIsoCh   = np.zeros(1, dtype=float)
#_mesonIso     = np.zeros(1, dtype=float)
_photonEt     = np.zeros(1, dtype=float)
_photonEta    = np.zeros(1, dtype=float)  
_photonEt     = np.zeros(1, dtype=float)
_photonEta    = np.zeros(1, dtype=float) 
_trksDeltaR   = np.zeros(1, dtype=float)
_eventWeight  = np.zeros(1, dtype=float)


tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('bosonMass',_bosonMass,'_bosonMass/D')
tree_output.Branch('mesonMass',_mesonMass,'_MesonMass/D')
tree_output.Branch('firstTrkPt',_firstTrkPt,'_firstTrkPt/D')
tree_output.Branch('secondTrkPt',_secondTrkPt,'_secondTrkPt/D')
tree_output.Branch('firstTrkEta',_firstTrkEta,'_firstTrkEta/D')
tree_output.Branch('secondTrkEta',_secondTrkEta,'_secondTrkEta/D')
tree_output.Branch('mesonPt',_mesonPt,'_mesonPt/D')
tree_output.Branch('mesonEta',_mesonEta,'_mesonEta/D')
tree_output.Branch('mesonIsoCh',_mesonIsoCh,'_mesonIsoCh/D')
#tree_output.Branch('mesonIso',_mesonIso,'_mesonIso/D')
tree_output.Branch('photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('photonEta',_photonEta,'_photonEta/D')
tree_output.Branch('trksDeltaR',_trksDeltaR,'_trksDeltaR/D')
tree_output.Branch('eventWeight',_eventWeight,'_eventWeight/D')


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
    bosonMass    = mytree.bosonMassFrom2TrksPhoton     
    mesonMass    = mytree.mesonTreeMassKin 
    firstTrkPt   = mytree.firstTrkPt         
    secondTrkPt  = mytree.secondTrkPt 
    firstTrkEta  = mytree.firstTrkEta
    secondTrkEta = mytree.secondTrkEta
    firstTrkPhi  = mytree.firstTrkPhi
    secondTrkPhi = mytree.secondTrkPhi
    mesonPt      = mytree.mesonTreePt      
    mesonEta     = mytree.mesonTreeEta
    photonEt     = mytree.photon_pT
    photonEta    = mytree.photon_eta
    #MesonIso    = mytree.iso_couple
    mesonIsoCh   = mytree.mesonIsoCh
    nElectrons   = mytree.nElectrons20
    nMuons       = mytree.nMuons20
    isMeson      = mytree.isMeson 
    #phi angle folding
    trksDeltaPhi = math.fabs(mytree.firstTrkPhi - mytree.secondTrkPhi)
    if trksDeltaPhi > 3.14: trksDeltaPhi = 6.28 - trksDeltaPhi  
    deltaR = math.sqrt((mytree.firstTrkEta - mytree.secondTrkEta)**2 + (trksDeltaPhi)**2)

    eventWeight = 1.

    #FILL HISTOS #####################################################################################################
    #if DATA -> Blind Analysis on Z inv mass plot  

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
    histo_map["h_mesonIsoCh"].Fill(mesonIsoCh, eventWeight)
    histo_map["h_photonEnergy"].Fill(photonEt, eventWeight)
    histo_map["h_photonEta"].Fill(photonEta, eventWeight)
    histo_map["h_trksDeltaR"].Fill(deltaR, eventWeight)#/mesonMass
    histo_map["h_nMuons"].Fill(nMuons, eventWeight)
    histo_map["h_nElectrons"].Fill(nElectrons, eventWeight)

    #FILL TREE ########################################################################################################
    _bosonMass[0]    = bosonMass
    _mesonMass[0]    = mesonMass
    _firstTrkPt[0]   = firstTrkPt
    _secondTrkPt[0]  = secondTrkPt
    _firstTrkEta[0]  = firstTrkEta
    _secondTrkEta[0] = secondTrkEta
    _mesonPt[0]      = mesonPt
    _mesonEta[0]     = mesonEta
    _mesonIsoCh[0]   = mesonIsoCh
    _photonEt[0]     = photonEt
    _photonEta[0]    = photonEta 
    _trksDeltaR[0]   = deltaR
    _eventWeight[0]  = eventWeight

    tree_output.Fill()


#HISTO LABELS #########################################################################################################
histo_map["h_bosonMass"].GetXaxis().SetTitle("m_{meson#gamma} [GeV/c^2]")
histo_map["h_bosonMass"].SetTitle("Meson+Photon invariant mass (Cut on phi inv. mass)")

histo_map["h_mesonMass"].GetXaxis().SetTitle("m_{meson} [GeV/c^2]")
histo_map["h_mesonMass"].SetTitle("Meson kin mass")

histo_map["h_firstTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_firstTrkPt"].SetTitle("Transverse momentum of the first charged particle (pT_{trk1})")

histo_map["h_secondTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_secondTrkPt"].SetTitle("Transverse momentum of the second charged particle (pT_{trk2})")

histo_map["h_firstTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_firstTrkEta"].SetTitle("Pseudorapidity of the first charged particle (#eta_{trk1} )")

histo_map["h_secondTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_secondTrkEta"].SetTitle("Pseudorapidity of the second charged particle (#eta_{trk2})")

histo_map["h_firstTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_firstTrkPhi"].SetTitle("Azimuthal angle of the first charged particle (#phi_{trk1})")

histo_map["h_secondTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_secondTrkPhi"].SetTitle("Azimuthal angle of the second charged particle (#phi_{trk2})")

histo_map["h_mesonPt"].GetXaxis().SetTitle("pT_{meson} [GeV]") 
histo_map["h_mesonPt"].SetTitle("Transverse kin momentum of the meson")

histo_map["h_mesonEta"].GetXaxis().SetTitle("#eta")
histo_map["h_mesonEta"].SetTitle("Kin pseudorapidity of the meson")

histo_map["h_mesonIsoCh"].GetXaxis().SetTitle("sum pT_{trks}/pT_{meson}")
histo_map["h_mesonIsoCh"].SetTitle("Charged isolation of the meson")

histo_map["h_trksDeltaR"].GetXaxis().SetTitle("#DeltaR_{trk^{+}trk^{-} }")
histo_map["h_trksDeltaR"].SetTitle("Delta R of the couple")

histo_map["h_photonEnergy"].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
histo_map["h_photonEnergy"].SetTitle("Energy of the photon")

histo_map["h_photonEta"].GetXaxis().SetTitle("#eta_{#gamma}")
histo_map["h_photonEta"].SetTitle("Eta of the photon")

histo_map["h_nMuons"].GetXaxis().SetTitle("nMuons over selection")
histo_map["h_nMuons"].SetTitle("# of muons")

histo_map["h_nElectrons"].GetXaxis().SetTitle("nElectrons over selection")
histo_map["h_nElectrons"].SetTitle("# of electrons")


#Tree writing ##########################################################################################################
tree_output.Write()


#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()