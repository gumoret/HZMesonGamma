import ROOT
import argparse
import math
import numpy as np
import sys
from array import array
import os

#Following bools are given as input
debug         = False
verbose       = False
isBDT         = False #BDT bool
isDataBlind   = False #Bool for blind analysis

isPhiAnalysis = False # for Z -> Phi Gamma
isRhoAnalysis = False # for Z -> Rho Gamma
isKAnalysis   = False # for Z -> K* Gamma
isDAnalysis   = False # for Z -> D0* Gamma

isHAnalysis   = False
isZAnalysis   = False

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument("runningOnData_option", help="Type <<signal>> for signal, <<data>> for data")
p.add_argument('boson_channel', help='type <<H>> or <<Z>>')
p.add_argument('meson_channel', help='type <<rho>> or <<phi>> or <<K>> or <<D>>')
p.add_argument('region_option', help='Type <<SR>> for Signal Region, <<CR>> for Control Region') #flag for bkg estimation
p.add_argument('blind_option', help='Type <<blind>> for blind analysis, <<unblind>> for unblinded analysis') #flag for blindness
p.add_argument('isBDT_option', help='Type <<preselection>> or <<BDT>>') #flag for loose selection or tight selection (from BDT output)
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_name', help='Provide output file name')
args = p.parse_args()


fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_name
mytree = fInput.Get("tree_output")
h_Events = fInput.Get("nEvents")

if args.runningOnData_option == "signal": runningOnData = False
elif args.runningOnData_option == "data": runningOnData = True

if args.boson_channel == "H": 
    isHAnalysis = True
    print("Boson: Higgs")  
elif args.boson_channel == "Z":
    isZAnalysis = True
    print("Boson: Z")
if args.meson_channel == "phi":
    isPhiAnalysis = True
    print("Meson: phi") 
if args.meson_channel == "rho": 
    isRhoAnalysis = True
    print("Meson: rho") 
if args.meson_channel == "K": 
    isKAnalysis = True
    print("Meson: K*0") 
if args.meson_channel == "D": 
    isDAnalysis = True
    print("Meson: D*0") 

CRFlag = args.region_option
if CRFlag == "SR": print("Processing the signal region")
if CRFlag == "CR": print("Processing the control region")

if args.blind_option == "blind": isDataBlind = True
if args.blind_option == "unblind": isDataBlind = False

if (args.isBDT_option == "BDT"):
    isBDT = True

isWideRange = False  #bool for wide range or zoomed range

'''
#Combine luminosity
luminosity = 39.54 #total lumi delivered during the trigger activity: 39.54 #fb^-1
if isPhiAnalysis:
    normalization_weight = (1./1.) * (1928000./0.0336) * 0.49 
elif isRhoAnalysis:
    normalization_weight = (1./1.) * (1928000./0.0336) 
'''

#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_bosonMass", "h_mesonMass", "h_firstTrkPt", "h_secondTrkPt", "h_firstTrkEta", "h_secondTrkEta", 
               "h_firstTrkPhi", "h_secondTrkPhi", "h_mesonPt", "h_mesonEta", "h_trksDeltaR","h_mesonIso", 
               "h_photonEnergy", "h_photonEta","h_nMuons","h_nElectrons", "h_efficiency"]  

if isWideRange:
    if isZAnalysis: histo_map[list_histos[0]]   = ROOT.TH1F(list_histos[0],"M_{Z}", 300, 50., 150.) 
    elif isHAnalysis: histo_map[list_histos[0]] = ROOT.TH1F(list_histos[0],"M_{H}", 300, 90, 155) 
    if   isPhiAnalysis: histo_map[list_histos[1]] = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.95, 1.1) 
    elif isRhoAnalysis: histo_map[list_histos[1]] = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.5, 1.)
    elif isKAnalysis: histo_map[list_histos[1]]   = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.6, 1.3) 
    histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T} of the 1st track", 200, 0.,200.)
    histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 250, 0.,250.)
    histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"#eta of the 1st track", 100, -2.5,2.5)
    histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 2nd track", 100, -2.5,2.5)
    histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#phi of the 1st track", 100, -math.pi, math.pi)
    histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 2nd track", 100, -math.pi, math.pi)
    histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"p_{T} of the meson", 500, 0., 500.)
    histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"#eta_{meson}", 100, -2.5,2.5)
    histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#Delta R_{meson}", 100, 0.,0.07)
    histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"Iso_ch of the meson", 100, 0., 1.5)
    histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"E_{T} of the #gamma", 500, 0., 500.)
    histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"#eta_{#gamma}", 100, -2.5,2.5)
    histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"n. of muons", 6, -0.5, 5.5)
    histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"n. of electrons", 5, -0.5, 5.5)

else:
    if isZAnalysis: histo_map[list_histos[0]]   = ROOT.TH1F(list_histos[0],"M_{Z}", 300, 60., 200.) 
    elif isHAnalysis and not runningOnData: histo_map[list_histos[0]] = ROOT.TH1F(list_histos[0],"M_{H}", 150, 120, 130)
    elif isHAnalysis and runningOnData: histo_map[list_histos[0]] = ROOT.TH1F(list_histos[0],"M_{H}", 150, 100, 170) 
    if   isPhiAnalysis: histo_map[list_histos[1]] = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 1., 1.04) 
    elif isRhoAnalysis: histo_map[list_histos[1]] = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.55, 1.)
    elif isKAnalysis: histo_map[list_histos[1]]   = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.8, 0.99) 
    histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T} of the 1st track", 100, 0.,70.)
    histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 0.,70.)
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

if not isBDT:
    histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"Efficiency steps", 5, 0., 5.)
else :
    histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"Efficiency steps", 6, 0., 6.)

#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree 
_bosonMass    = np.zeros(1, dtype=float)  
_mesonMass    = np.zeros(1, dtype=float)
_firstTrkPt   = np.zeros(1, dtype=float)
_secondTrkPt  = np.zeros(1, dtype=float)
_firstTrkEta  = np.zeros(1, dtype=float)  
_secondTrkEta = np.zeros(1, dtype=float)
_mesonPt      = np.zeros(1, dtype=float)
_mesonEta     = np.zeros(1, dtype=float)  
#_mesonIsoCh   = np.zeros(1, dtype=float)
_mesonIso     = np.zeros(1, dtype=float)
_photonEt     = np.zeros(1, dtype=float)
_photonEta    = np.zeros(1, dtype=float)  
_photonEt     = np.zeros(1, dtype=float)
_photonEta    = np.zeros(1, dtype=float) 
_trksDeltaR   = np.zeros(1, dtype=float)
_eventWeight  = np.zeros(1, dtype=float)

#output tree
tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.SetDirectory(fOut)
tree_output.SetAutoSave(0)
tree_output.SetAutoFlush(0)

tree_output.Branch('bosonMass',_bosonMass,'_bosonMass/D')
tree_output.Branch('mesonMass',_mesonMass,'_MesonMass/D')
tree_output.Branch('firstTrkPt',_firstTrkPt,'_firstTrkPt/D')
tree_output.Branch('secondTrkPt',_secondTrkPt,'_secondTrkPt/D')
tree_output.Branch('firstTrkEta',_firstTrkEta,'_firstTrkEta/D')
tree_output.Branch('secondTrkEta',_secondTrkEta,'_secondTrkEta/D')
tree_output.Branch('mesonPt',_mesonPt,'_mesonPt/D')
tree_output.Branch('mesonEta',_mesonEta,'_mesonEta/D')
#tree_output.Branch('mesonIsoCh',_mesonIsoCh,'_mesonIsoCh/D')
tree_output.Branch('mesonIso',_mesonIso,'_mesonIso/D')
tree_output.Branch('photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('photonEta',_photonEta,'_photonEta/D')
tree_output.Branch('trksDeltaR',_trksDeltaR,'_trksDeltaR/D')
tree_output.Branch('eventWeight',_eventWeight,'_eventWeight/D')

#------------- counters -----------------
#nEventsMesonAnalysis  = 0
nEventsOverLeptonVeto = 0
nEventsOverDiPhotonVeto   = 0
nEventsAfterRegionDefiner = 0
nEventsInBosonMassRange       = 0
nEventsOverCuts           = 0
nEventsLeftSB             = 0
nEventsRightSB            = 0
nHmatched                 = 0
nMesonMatched             = 0
nTightSelection           = 0


#EVENTS LOOP ########################################################################################################
nentries = mytree.GetEntriesFast()
for jentry in range(nentries):
    ientry = mytree.LoadTree(jentry)
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry)
    if nb <= 0:
        print("nb < 0")
        continue

    if debug: print("Processing EVENT n.",jentry+1, " ...")

    #Retrieve variables from the tree
    bosonMass    = mytree.H_mass     
    mesonMass    = mytree.bestMeson_mass 
    firstTrkPt   = mytree.firstTrk_pt         
    secondTrkPt  = mytree.secondTrk_pt 
    firstTrkEta  = mytree.firstTrk_eta
    secondTrkEta = mytree.secondTrk_eta
    firstTrkPhi  = mytree.firstTrk_phi
    secondTrkPhi = mytree.secondTrk_phi
    mesonPt      = mytree.bestMeson_pt      
    mesonEta     = mytree.bestMeson_eta
    photonEt     = mytree.bestPhoton_pt
    photonEta    = mytree.bestPhoton_eta
    mesonIso     = mytree.isoMeson
    #mesonIsoCh   = mytree.mesonIsoCh
    nElectrons   = mytree.nElectrons20
    nMuons       = mytree.nMuons20
    #isMeson      = mytree.isMeson 
    #phi angle folding
    trksDeltaPhi = abs(firstTrkPhi -secondTrkPhi)
    if trksDeltaPhi > math.pi: trksDeltaPhi = 6.28 - trksDeltaPhi  
    deltaR       = math.sqrt((firstTrkEta - secondTrkEta)**2 + trksDeltaPhi*trksDeltaPhi)
    #eventWeight  = 1 ###placeholder 


    #Define Control and Signal regions: ------------------------------------------
    if isPhiAnalysis: 
        if CRFlag == "SR" and not (mesonMass > 1.008 and mesonMass < 1.032): continue
        if CRFlag == "CR" and (mesonMass > 1.008 and mesonMass < 1.032): continue

    if isRhoAnalysis:
        if CRFlag == "SR" and not (mesonMass > 0.62 and mesonMass < 0.92): continue
        if CRFlag == "CR" and (mesonMass > 0.62 and mesonMass < 0.92): continue

    if isKAnalysis:
        if CRFlag == "SR" and not (mesonMass > 0.842 and mesonMass < 0.942): continue
        if CRFlag == "CR" and (mesonMass > 0.842 and mesonMass < 0.942): continue

    nEventsAfterRegionDefiner+=1

    '''
    #NORMALIZATION for MC-------------------------------------------------------------------

    if not runningOnData:
        #PUWeight    = mytree.PU_Weight
        weight_sign = mytree.MC_Weight/abs(mytree.MC_Weight) #just take the sign of the MC gen weight

        eventWeight =  luminosity * normalization_weight * weight_sign # * PUWeight  
    '''
    #else:
    eventWeight = 1.


    #Lepton veto
    if nElectrons > 0: continue
    if nMuons     > 0: continue

    nEventsOverLeptonVeto += 1

    #-------------- n events in the sidebands -----------------------------
    if (bosonMass < 50. or bosonMass > 200.): continue
    nEventsInBosonMassRange+=1

    '''
    #TIGHT SELECTION from BDT output -------------------------------------------------  
    if isBDT: 
        BDT_out = myWF.get_BDT_output(firstTrkisoCh,MesonIso0,ZMass,mesonEta,MesonGammaDeltaPhi)#,mesonPt,photonEt,photonEta,nJets)#,JetNeutralEmEn,JetChargedHadEn,JetNeutralHadEn) 
        #histo_map["h_BDT_out"].Fill(BDT_out)

        if debug: print "BDT value before selection = ", BDT_out
        if args.isBDT_option == "BDT":
            if BDT_out < BDT_OUT: #Cut on BDT output
                if debug: print "BDT cut NOT passed"
                continue

        nTightSelection+=1
    '''
    if runningOnData == True:
        if isHAnalysis: #change limits 
            if (CRFlag == 'SR' and bosonMass > 50. and bosonMass < 120.) : nEventsLeftSB  += 1
            if (CRFlag == 'SR' and bosonMass > 130. and bosonMass < 200.) : nEventsRightSB += 1
        if isZAnalysis:
            if (CRFlag == 'SR' and bosonMass > 50. and bosonMass < 80.) : nEventsLeftSB  += 1
            if (CRFlag == 'SR' and bosonMass > 101. and bosonMass < 200.) : nEventsRightSB += 1

    '''
    if not runningOnData == True:
        weightSum += eventWeight
        #polWeightSum += polarizationWeight
    
        if isPhiAnalysis: eventWeight = eventWeight/(weightSum_eff)*(1928000/0.0336)*0.49*luminosity
        else: eventWeight = eventWeight/(weightSum_eff)*(1928000/0.0336)*luminosity
    '''

    #FILL HISTOS #####################################################################################################
    #if DATA -> Blind Analysis on boson inv mass plot  

    if runningOnData == True:
        if isDataBlind:
            if isZAnalysis:
                if bosonMass < 80. or bosonMass > 100: histo_map["h_bosonMass"].Fill(bosonMass, eventWeight)
            if isHAnalysis:
                if bosonMass < 120. or bosonMass > 130: histo_map["h_bosonMass"].Fill(bosonMass, eventWeight)
        else:
            histo_map["h_bosonMass"].Fill(bosonMass, eventWeight)
    else:
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
    histo_map["h_mesonIso"].Fill(mesonIso, eventWeight)
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
    _mesonIso[0]     = mesonIso
    _photonEt[0]     = photonEt
    _photonEta[0]    = photonEta 
    _trksDeltaR[0]   = deltaR
    _eventWeight[0]  = eventWeight

    tree_output.Fill()

    if debug:
        print("***********************************")
        print("*** EVENT RECORDED: tree filled ***")
        print("***********************************")
    
    #counters
    nEventsOverCuts += 1

#HISTO LABELS #########################################################################################################
histo_map["h_bosonMass"].GetXaxis().SetTitle("m_{meson#gamma} [GeV/c^2]")
histo_map["h_bosonMass"].SetTitle("Meson+Photon invariant mass")

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

histo_map["h_mesonIso"].GetXaxis().SetTitle("sum pT_{trks}/pT_{meson}")
histo_map["h_mesonIso"].SetTitle("Isolation of the meson")

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

histo_map["h_efficiency"].GetXaxis().SetTitle("")
histo_map["h_efficiency"].GetYaxis().SetTitle("#epsilon (%)")
#Tree writing ##########################################################################################################
#tree_output.Write()


#Variables for cut overflow
bin1content  = h_Events.GetBinContent(1)
bin2content  = h_Events.GetBinContent(2)
bin3content  = h_Events.GetBinContent(3)
bin4content  = h_Events.GetBinContent(4)
bin5content  = h_Events.GetBinContent(5)
if isBDT: bin6content = nTightSelection

#variables for final print
nEventsProcessed = bin1content
nEventsTriggered = bin2content
nEventsPhoton = bin3content
nEventsMesonMassSR = nEventsRightSB + nEventsLeftSB

nSignal      = bin1content
scale_factor = 100/nSignal

histo_map["h_efficiency"].Fill(0.5,bin1content*scale_factor)
histo_map["h_efficiency"].Fill(1.5,bin2content*scale_factor)
histo_map["h_efficiency"].Fill(2.5,bin3content*scale_factor)
histo_map["h_efficiency"].Fill(3.5,bin4content*scale_factor)
histo_map["h_efficiency"].Fill(4.5,bin5content*scale_factor)
if isBDT: histo_map["h_efficiency"].Fill(5.5,bin6content*scale_factor)    


histo_map["h_efficiency"].GetXaxis().SetBinLabel(1,"Events processed")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(2,"Events triggered")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(3,"Photon requested")
#histo_map["h_efficiency"].GetXaxis().SetBinLabel(4,"Iso selection")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(4,"Best couple found")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(5,"trk-cand pT selection")
if isBDT:  histo_map["h_efficiency"].GetXaxis().SetBinLabel(6,"Tight selection")

c11 = ROOT.TCanvas()
c11.cd()
histo_map["h_efficiency"].SetFillColor(1) 
histo_map["h_efficiency"].SetFillStyle(3003)
ROOT.gStyle.SetPaintTextFormat("4.2f %")
ROOT.gStyle.SetOptStat(0)
histo_map["h_efficiency"].SetMarkerSize(1.4)
histo_map["h_efficiency"].GetXaxis().SetRangeUser(0.,4.1)##modify
#histo_map["h_efficiency"].GetYaxis().SetRangeUser(0.,30.)
#histo_map["h_efficiency"].SetMaximum(max(histo_map["h_efficiency"].GetHistogram().GetMaximum(),30.))
histo_map["h_efficiency"].Draw("HIST TEXT0")

if not runningOnData:
    h_effinciencyName = os.path.basename(output_filename).replace(".root","") 
    if isHAnalysis:
        c11.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/histos/H/"+h_effinciencyName+"_hefficiency.pdf")       
        c11.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/histos/H/"+h_effinciencyName+"_hefficiency.png")
    if isZAnalysis:
        c11.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/histos/Z/"+h_effinciencyName+"_hefficiency.pdf")       
        c11.SaveAs("/eos/user/e/eferrand/Work/CMSSW_15_0_6/src/HZMesonGammaAnalysis/HZMesonGamma/test/histos/Z/"+h_effinciencyName+"_hefficiency.png")

#FINAL PRINTS ###########################################################
print("\n\nCUT OVERFLOW")
print("---------------------------------------")
print("CRflag                    = ",CRFlag)
print("nEventsProcessed          = ",nEventsProcessed," (n. events in dataset)")
print("nEventsTriggered          = ",nEventsTriggered," (n. events over HLT)")
print("nEventsPhoton             = ",nEventsPhoton," (n. events with a best photon found)")
print("nEventsAfterRegionDefiner = ",nEventsAfterRegionDefiner," (split in SR or CR of the meson inv mass)")
print("nEventsOverLeptonVeto     = ",nEventsOverLeptonVeto," (n. events without any electron or muon)")
print("nEventsInBosonMassRange   = ",nEventsInBosonMassRange," (n. events in 50 < bosonMass < 200 GeV)")
print("nEventsMesonMassSR        = ",nEventsMesonMassSR," (Events in SR of MesonMass and in SBs of bosonMass)")
print("\n----------- SUMMARY -----------------------")
if isBDT:
    print("BDT output used = ",BDT_OUT)
if not runningOnData:
    print("Signal MC sample")
if isBDT:
    print("n. events after preselection = ",jentry)
print("nEventsInBosonMassRange   = ",nEventsInBosonMassRange," (n. events in 50 < bosonMass < 200 GeV)")
print("n. events after cuts      = " , nEventsOverCuts)

if runningOnData and CRFlag == 'SR':
    print( "n. events in the left sideband counted = ",nEventsLeftSB)
    print("n. events in the right sideband counted = ",nEventsRightSB) 
    print( "Total events in the sidebands = ", nEventsRightSB + nEventsLeftSB)

if not runningOnData:
    #print("Signal weight sum   = ",float(weightSum))
    #if isPhiAnalysis: print "Signal integral     = ",histo_map["h_ZMass"].Integral()
    #else : print "Signal integral     = ",histo_map["h_ZMass"].Integral()#/(weightSum_eff)*(1928000/0.0336)*luminosity
    #print "Total signal weight = ",weightSum_eff
    print("Total signal events       = ",bin1content)
    print("Signal efficiency         = ",nEventsOverCuts/bin1content)
    #print "new eff             = ",weightSum/weightSum_eff
    
print("-------------------------------------------\n\n")

#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Write()
fOut.Close()