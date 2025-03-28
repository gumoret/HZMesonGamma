import ROOT
import argparse
import math
import numpy as np
import sys
from array import array
from ROOT import TLorentzVector


#Following bools are given as input
verbose       = False
isPhiAnalysis = False # for Z -> Phi Gamma
isRhoAnalysis = False # for Z -> Rho Gamma
isKAnalysis   = False # for Z -> K* Gamma
isD0Analysis  = False # for Z -> D0* Gamma



#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('meson_option', help='Type <<rho>> for rho, <<phi>> for phi, <<K>> for K*, <<D0>> for D0*') #flag for type of meson
p.add_argument('runningOnData_option', help='Type <<signal>> for signal, <<data>> for data') #flag for data or signal
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("Events")

if not mytree:
    print("Error: Tree 'Events' not found in file", args.rootfile_name)
    sys.exit(1)

if args.meson_option == "phi": isPhiAnalysis = True
elif args.meson_option == "rho": isRhoAnalysis = True
elif args.meson_option == "K": isKAnalysis = True
elif args.meson_option == "D0": isD0Analysis = True
else: print("meson_option must be <<phi>> or <<rho>> or <<K*>> or <<DO*>>")


if args.runningOnData_option == "signal": runningOnData = False
elif args.runningOnData_option == "data": runningOnData = True
else: print("runninOnData must be <<signal>> or <<data>>")


#print("isPhiAnalysis =", isPhiAnalysis, "isRhoAnalysis =", isRhoAnalysis, "runningOnData =", runningOnData)

#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["nEvents", "pileup"]

histo_map[list_histos[0]] = ROOT.TH1F(list_histos[0],"Event counting in different steps", 5, 0., 5.)
histo_map[list_histos[1]] = ROOT.TH1F(list_histos[1],"pileup", 130, 0, 130)


#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree ################################################################################
_MC_Weight                = np.zeros(1, dtype=float)
_isTwoProngTrigger        = np.zeros(1, dtype=bool)
_nMuons10                 = np.zeros(1, dtype=int)
_nMuons20                 = np.zeros(1, dtype=int)
_nElectrons10             = np.zeros(1, dtype=int)
_nElectrons20             = np.zeros(1, dtype=int)
_nPhotonsChosen           = np.zeros(1, dtype=int)
_nPhotons20WP90           = np.zeros(1, dtype=int)
_nPhotons35WP80           = np.zeros(1, dtype=int)
_photon_pT                = np.zeros(1, dtype=float)
_photon_eta               = np.zeros(1, dtype=float)
_photonEtaSC              = np.zeros(1, dtype=float)
_photonRegressionError    = np.zeros(1, dtype=float)
_firstTrkPt               = np.zeros(1, dtype=float)
_firstTrkEta              = np.zeros(1, dtype=float)
_firstTrkPhi              = np.zeros(1, dtype=float)
_secondTrkPt              = np.zeros(1, dtype=float)
_secondTrkEta             = np.zeros(1, dtype=float)
_secondTrkPhi             = np.zeros(1, dtype=float)
_isMeson                  = np.zeros(1, dtype=bool)
_mesonPt                  = np.zeros(1, dtype=float)
_mesonEta                 = np.zeros(1, dtype=float)
_mesonPhi                 = np.zeros(1, dtype=float)
_mesonIsoCh               = np.zeros(1, dtype=float)
_mesonMassTrkTrk          = np.zeros(1, dtype=float)
_bosonMassFrom2TrksPhoton = np.zeros(1, dtype=float)

_mesonTreeMass            = np.zeros(1, dtype=float)
_mesonTreeMassKin         = np.zeros(1, dtype=float)
_mesonTreePt              = np.zeros(1, dtype=float)
_mesonTreeEta             = np.zeros(1, dtype=float)
_mesonTreePhi             = np.zeros(1, dtype=float)



tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('MC_Weight', _MC_Weight, '_MC_Weight/D')
tree_output.Branch('isTwoProngTrigger', _isTwoProngTrigger, '_isTwoProngTrigger/O')
tree_output.Branch('nMuons10', _nMuons10, '_nMuons10/I')
tree_output.Branch('nMuons20', _nMuons20, '_nMuons20/I')
tree_output.Branch('nElectrons10', _nElectrons10, '_nElectrons10/I')
tree_output.Branch('nElectrons20', _nElectrons20, '_nElectrons20/I')
tree_output.Branch('nPhotonsChosen', _nPhotonsChosen, '_nPhotonsChosen/I')
tree_output.Branch('nPhotons20WP90', _nPhotons20WP90, '_nPhotons20WP90/I')
tree_output.Branch('nPhotons35WP80', _nPhotons35WP80, '_nPhotons35WP80/I')
tree_output.Branch('photon_pT', _photon_pT, '_photon_pT/D')
tree_output.Branch('photon_eta', _photon_eta, '_photon_eta/D')
tree_output.Branch('photonEtaSC', _photonEtaSC, '_photonEtaSC/D')
tree_output.Branch('photonRegressionError', _photonRegressionError, '_photonRegressionError/D')
tree_output.Branch('firstTrkPt', _firstTrkPt, '_firstTrkPt/D')
tree_output.Branch('firstTrkEta', _firstTrkEta, '_firstTrkEta/D')
tree_output.Branch('firstTrkPhi', _firstTrkPhi, '_firstTrkPhi/D')
tree_output.Branch('secondTrkPt', _secondTrkPt, '_secondTrkPt/D')
tree_output.Branch('secondTrkEta', _secondTrkEta, '_secondTrkEta/D')
tree_output.Branch('secondTrkPhi', _secondTrkPhi, '_secondTrkPhi/D')
tree_output.Branch('isMeson', _isMeson, '_isMeson/O')
tree_output.Branch('mesonPt', _mesonPt, '_mesonPt/D')
tree_output.Branch('mesonEta', _mesonEta, '_mesonEta/D')
tree_output.Branch('mesonPhi', _mesonPhi, '_mesonPhi/D')
tree_output.Branch('mesonIsoCh', _mesonIsoCh, '_mesonIsoCh/D')
tree_output.Branch('mesonMassTrkTrk', _mesonMassTrkTrk, '_mesonMassTrkTrk/D')
tree_output.Branch('bosonMassFrom2TrksPhoton', _bosonMassFrom2TrksPhoton, '_bosonMassFrom2TrksPhoton/D')

tree_output.Branch('mesonTreeMass', _mesonTreeMass, '_mesonTreeMass/D')
tree_output.Branch('mesonTreeMassKin', _mesonTreeMassKin, '_mesonTreeMassKin/D')
tree_output.Branch('mesonTreePt', _mesonTreePt, '_mesonTreePt/D')
tree_output.Branch('mesonTreeEta', _mesonTreeEta, '_mesonTreeEta/D')
tree_output.Branch('mesonTreePhi', _mesonTreePhi, '_mesonTreePhi/D')


########### MESON ANALYSIS FUNCTION ###############
def process_meson(mytree, meson_type, mass_cut, verbose=False):
    isBestMesonOfTheEventFound = False
    bestMesonPt = 0.
    # Initialize variables to return
    (deltaRTrks_chosen, isMeson_chosen, firstTrkP4_chosen, secondTrkP4_chosen, mesonP4_chosen, mesonIsoCh_chosen, mesonTreeMass_chosen, mesonTreeMassKin_chosen, mesonTreePt_chosen, mesonTreeEta_chosen, mesonTreePhi_chosen) = (None,) * 11


    if meson_type == "rho":
        n = len(mytree.rho_mass)
        trk1_pt, trk1_eta, trk1_phi = mytree.rho_trk1_pt, mytree.rho_trk1_eta, mytree.rho_trk1_phi
        trk2_pt, trk2_eta, trk2_phi = mytree.rho_trk2_pt, mytree.rho_trk2_eta, mytree.rho_trk2_phi
        iso_meson = mytree.rho_iso
        meson_mass, meson_mass_kin = mytree.rho_mass, mytree.rho_kin_mass
        meson_pt, meson_eta, meson_phi = mytree.rho_kin_pt, mytree.rho_kin_eta, mytree.rho_kin_phi
        mass_trk1 = 0.13957  #Pion mass
        mass_trk2 = 0.13957  

    elif meson_type == "phi":
        n = len(mytree.phi_mass)
        trk1_pt, trk1_eta, trk1_phi = mytree.phi_trk1_pt, mytree.phi_trk1_eta, mytree.phi_trk1_phi
        trk2_pt, trk2_eta, trk2_phi = mytree.phi_trk2_pt, mytree.phi_trk2_eta, mytree.phi_trk2_phi
        iso_meson = mytree.phi_iso
        meson_mass, meson_mass_kin = mytree.phi_mass, mytree.phi_kin_mass
        meson_pt, meson_eta, meson_phi = mytree.phi_kin_pt, mytree.phi_kin_eta, mytree.phi_kin_phi
        mass_trk1 = 0.4937  # Kaon mass
        mass_trk2 = 0.4937  

    elif meson_type == "K":
        n = len(mytree.K0Star_mass)
        trk1_pt, trk1_eta, trk1_phi = mytree.K0Star_pion_pt, mytree.K0Star_pion_eta, mytree.K0Star_pion_phi
        trk2_pt, trk2_eta, trk2_phi = mytree.K0Star_kaon_pt, mytree.K0Star_kaon_eta, mytree.K0Star_kaon_phi
        iso_meson = mytree.K0Star_iso
        meson_mass, meson_mass_kin = mytree.K0Star_mass, mytree.K0Star_kin_mass
        meson_pt, meson_eta, meson_phi = mytree.K0Star_kin_pt, mytree.K0Star_kin_eta, mytree.K0Star_kin_phi
        mass_trk1 = 0.13957  # Pion mass
        mass_trk2 = 0.4937  # Kaon mass

    else:
        return

    for i in range(n):
        firstTrkPt, firstTrkEta, firstTrkPhi = trk1_pt[i], trk1_eta[i], trk1_phi[i]
        secondTrkPt, secondTrkEta, secondTrkPhi = trk2_pt[i], trk2_eta[i], trk2_phi[i]
        isoMesonCh = iso_meson[i]

        mesonTreeMass, mesonTreeMassKin = meson_mass[i], meson_mass_kin[i]
        mesonTreePt, mesonTreeEta, mesonTreePhi = meson_pt[i], meson_eta[i], meson_phi[i]

        #Tracks pt cuts------------------------------------------------
        if firstTrkPt < 1. or secondTrkPt < 1.: 
            if verbose: print("if first trk pt or second trk pt<1 continue")
            continue
        if firstTrkPt < 10. and secondTrkPt < 10.: 
            if verbose: print("if first trk pt and second trk pt<10 continue")
            continue

        #DiTrks deltaR cut---------------------------------------------------
        deltaEta, deltaPhi = firstTrkEta - secondTrkEta, abs(firstTrkPhi - secondTrkPhi)
        if deltaPhi > math.pi: deltaPhi = 2 * math.pi - deltaPhi
        deltaRTrks = math.sqrt(deltaEta**2 + deltaPhi**2)
        if deltaRTrks > 0.07: 
            if verbose: print("if deltaRTrks>0.07 continue")
            continue

        #Quadrimomentum calculation -------------------------------------------
        firstTrkP4, secondTrkP4 = TLorentzVector(), TLorentzVector()
        firstTrkP4.SetPtEtaPhiM(firstTrkPt, firstTrkEta, firstTrkPhi, mass_trk1)
        secondTrkP4.SetPtEtaPhiM(secondTrkPt, secondTrkEta, secondTrkPhi, mass_trk2)
        pairP4 = firstTrkP4 + secondTrkP4
        if verbose: print("PiPi pT =", pairP4.Pt())

        #DiTrk pT cut----------------------------------------------------------
        if pairP4.Pt() < 38: continue

        #Meson inv mass---------------------------------------------------------
        mesonMassTrkTrk = pairP4.M()
        isMeson = mass_cut[0] < mesonMassTrkTrk < mass_cut[1]
        if not isMeson: continue

        #Isolation cut----------------------------------------------------------
        if isoMesonCh < 0.9: 
            print("No isolation cut passed, continue")
            continue

        #pT max of the event filter ----------------------------------------------
        if verbose: print("Current bestMeson_Pt =", bestMesonPt)
        if pairP4.Pt() <= bestMesonPt: 
            if verbose: print("Not passed: pT lower than the current best meson of the event. Continue")
            continue

        bestMesonPt = pairP4.Pt()
        if verbose:
                print(f"pairP4.Pt() = {bestMesonPt}")
                print("This is the best meson so far!")

        isBestMesonOfTheEventFound = True

        #Save variables if best meson has been found
        deltaRTrks_chosen    = deltaRTrks
        isMeson_chosen       = isMeson
        firstTrkP4_chosen    = firstTrkP4
        secondTrkP4_chosen   = secondTrkP4
        mesonP4_chosen       = pairP4
        mesonIsoCh_chosen    = isoMesonCh

        mesonTreeMass_chosen    = mesonTreeMass
        mesonTreeMassKin_chosen = mesonTreeMassKin
        mesonTreePt_chosen      = mesonTreePt
        mesonTreeEta_chosen     = mesonTreeEta
        mesonTreePhi_chosen     = mesonTreePhi

    return isBestMesonOfTheEventFound, deltaRTrks_chosen, isMeson_chosen, firstTrkP4_chosen, secondTrkP4_chosen, mesonP4_chosen, mesonIsoCh_chosen, mesonTreeMass_chosen, mesonTreeMassKin_chosen, mesonTreePt_chosen, mesonTreeEta_chosen, mesonTreePhi_chosen
##########################################

#EVENTS LOOP #######################################################################################################
nEventsProcessed            = 0
nEventsTriggered            = 0
nEventsIsPhoton             = 0
nEventsMesonIsolationFilter = 0
nEventsBestPairFound        = 0
nEventsTrkPtFilter          = 0
nEventsBosonMatched         = 0
nEventsBosonNotMatched      = 0

nentries = mytree.GetEntriesFast()

# Loop sugli eventi
for jentry in range(nentries):
    ientry = mytree.LoadTree(jentry)
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry)
    if nb <= 0:
        continue

    nEventsProcessed+=1
    histo_map["nEvents"].SetBinContent(1, nEventsProcessed)

    print("Processing event ", nEventsProcessed)



    #//*************************************************************//
    #//                                                             //
    #//---------------------------Vertices -------------------------//
    #//                                                             //
    #//*************************************************************//

    #----//----#

    
    #//*************************************************************//
    #//                                                             //
    #//--------------------------- Pile Up -------------------------//
    #//                                                             //
    #//*************************************************************//

    nPU = mytree.Pileup_nPU

    histo_map["pileup"].Fill(nPU)



    #//*************************************************************//
    #//                                                             //
    #//--------------------------MC weight -------------------------//
    #//                                                             //
    #//*************************************************************//

    MC_Weight = mytree.Generator_weight

    _MC_Weight[0] = MC_Weight



    #//*************************************************************//
    #//                                                             //
    #//--------------------------Trigger-- -------------------------//
    #//                                                             //
    #//*************************************************************//

    isTwoProngTrigger = mytree.HLT_Photon35_TwoProngs35

    if not isTwoProngTrigger:
        #if verbose: print("Event not triggered, RETURN.")
        continue 

    nEventsTriggered+=1
    histo_map["nEvents"].SetBinContent(2, nEventsTriggered)

    _isTwoProngTrigger[0] = isTwoProngTrigger

    

    #//*************************************************************//
    #//                                                             //
    #//--------------------------Muons------------------------------//
    #//                                                             //
    #//*************************************************************//

    nMuons10 = 0
    nMuons20 = 0    
    
    for i in range(len(mytree.Muon_pt)):        
        muon_pT = mytree.Muon_pt[i]
        isIdMedium = mytree.Muon_mediumId[i]
        muon_eta = mytree.Muon_eta[i]
        muon_dxy = mytree.Muon_dxy[i]
        muon_dz = mytree.Muon_dz[i]
        isPFIsoLoose = mytree.Muon_pfIsoId[i]

        #if verbose: print(f"Muon {i}: pT={muon_pT}, mediumId={isIdMedium}, eta={muon_eta}, dxy={muon_dxy}, dz={muon_dz}, pfIsoId={isPFIsoLoose}")

        if muon_pT < 10. or not isIdMedium or abs(muon_eta) > 2.4 or abs(muon_dxy) >= 0.2 or abs(muon_dz) >= 0.5: continue

        if isPFIsoLoose != 2: continue

        nMuons10+=1
                
        if muon_pT < 20.: continue

        nMuons20+=1


    _nMuons10[0] = nMuons10
    _nMuons20[0] = nMuons20



    #//*************************************************************//
    #//                                                             //
    #//----------------------Electrons------------------------------//
    #//                                                             //
    #//*************************************************************//

    nElectrons10 = 0
    nElectrons20 = 0 
    
    for i in range(len(mytree.Electron_pt)):        
        electron_pT = mytree.Electron_pt[i]
        electron_eta = mytree.Electron_eta[i]
        electron_dxy = mytree.Electron_dxy[i]
        electron_dz = mytree.Electron_dz[i]
        isMVAIso_WP80 = mytree.Electron_mvaIso_WP80[i]

        #if verbose: print(f"Electron {i}: pT={electron_pT}, eta={electron_eta}, dxy={electron_dxy}, dz={electron_dz}, MVAIso_WP80={isMVAIso_WP80}")

        if electron_pT < 10. or abs(electron_eta) > 2.5 or abs(electron_dxy) >= 0.2 or abs(electron_dz) >= 0.5: continue

        #-------------Conditions on loose/medium MVA electron ID-------------#
        if not isMVAIso_WP80: continue

        nElectrons10+=1

        if electron_pT < 20.: continue

        nElectrons20+=1


    _nElectrons10[0] = nElectrons10
    _nElectrons20[0] = nElectrons20



    #//*************************************************************//
    #//                                                             //
    #//------------------------Photons------------------------------//
    #//                                                             //
    #//*************************************************************//

    nPhotons20WP90    = 0
    nPhotons35WP80    = 0
    photonEtMax       = -1000.
    cand_photon_found = False


    for i in range(len(mytree.Photon_pt)):
        photon_pT              = mytree.Photon_pt[i]
        photon_eta             = mytree.Photon_eta[i]
        photon_phi             = mytree.Photon_phi[i]
        isMVAIso_WP80          = mytree.Photon_mvaID_WP80[i]
        isMVAIso_WP90          = mytree.Photon_mvaID_WP90[i]
        isElectronVeto         = mytree.Photon_electronVeto[i]
        photonIsoChargedHadron = mytree.Photon_pfChargedIso[i] 
        photonEtaSC            = mytree.Photon_superclusterEta[i]
        photonRegressionError  = mytree.Photon_energyErr[i]

        if verbose: print("Photon", i, ":", "pT=", photon_pT, "eta=", photon_eta, "WP80=", isMVAIso_WP80, "WP90=", isMVAIso_WP90, "veto=", isElectronVeto)

        if photon_pT < 20 or abs(photon_eta) > 2.5: continue

        if not isMVAIso_WP90 or not isElectronVeto: continue        
        nPhotons20WP90+=1

        if photon_pT < 35: continue
        if not isMVAIso_WP80: continue 
        nPhotons35WP80+=1
        
        if photon_pT < photonEtMax: continue

        photonEtMax = photon_pT
        
        photon_eta_chosen             = photon_eta
        photon_phi_chosen             = photon_phi
        photonIsoChargedHadron_chosen = photonIsoChargedHadron
        photonEtaSC_chosen            = photonEtaSC
        photonRegressionError_chosen  = photonRegressionError
        
        ph_p4 = TLorentzVector()
        ph_p4.SetPtEtaPhiM(photonEtMax, photon_eta_chosen, photon_phi_chosen, 0.)

        cand_photon_found = True


    if not cand_photon_found: continue 

    nEventsIsPhoton+=1
    histo_map["nEvents"].SetBinContent(3, nEventsIsPhoton) 

    if verbose: print("Chosen photon: pT=", photonEtMax, "eta=", photon_eta_chosen)
    #print(f"photons 20= {nPhotons20WP90}, photons 35 = {nPhotons35WP80}")


    _photon_pT[0]             = photonEtMax
    _photon_eta[0]            = photon_eta_chosen
    _photonEtaSC[0]           = photonEtaSC_chosen
    _photonRegressionError[0] = photonRegressionError_chosen
    _nPhotons20WP90[0]        = nPhotons20WP90
    _nPhotons35WP80[0]        = nPhotons35WP80


    
    #//*************************************************************//
    #//                                                             //
    #//-------------------------Mesons------------------------------//
    #//                                                             //
    #//*************************************************************//
    if isRhoAnalysis: isBestMesonOfTheEventFound, deltaRTrks_chosen, isMeson_chosen, firstTrkP4_chosen, secondTrkP4_chosen, mesonP4_chosen, mesonIsoCh_chosen, mesonTreeMass_chosen, mesonTreeMassKin_chosen, mesonTreePt_chosen, mesonTreeEta_chosen, mesonTreePhi_chosen = process_meson(mytree, "rho", (0.5, 1.0), verbose)
    if isPhiAnalysis: isBestMesonOfTheEventFound, deltaRTrks_chosen, isMeson_chosen, firstTrkP4_chosen, secondTrkP4_chosen, mesonP4_chosen, mesonIsoCh_chosen, mesonTreeMass_chosen, mesonTreeMassKin_chosen, mesonTreePt_chosen, mesonTreeEta_chosen, mesonTreePhi_chosen = process_meson(mytree, "phi", (1.0, 1.05), verbose)
    if isKAnalysis: isBestMesonOfTheEventFound, deltaRTrks_chosen, isMeson_chosen, firstTrkP4_chosen, secondTrkP4_chosen, mesonP4_chosen, mesonIsoCh_chosen, mesonTreeMass_chosen, mesonTreeMassKin_chosen, mesonTreePt_chosen, mesonTreeEta_chosen, mesonTreePhi_chosen   = process_meson(mytree, "K", (0.6, 1.0), verbose)
    
    if not isBestMesonOfTheEventFound:
        print("No best couple detected for current event, RETURN.")
        continue  

    nEventsBestPairFound+=1
    histo_map["nEvents"].SetBinContent(4, nEventsBestPairFound)


    #DATAMEMBER SAVING
    firstTrkPt   = firstTrkP4_chosen.Pt()
    firstTrkEta  = firstTrkP4_chosen.Eta()
    firstTrkPhi  = firstTrkP4_chosen.Phi()
    secondTrkPt  = secondTrkP4_chosen.Pt()       
    secondTrkEta = secondTrkP4_chosen.Eta()
    secondTrkPhi = secondTrkP4_chosen.Phi()
    mesonPt      = mesonP4_chosen.Pt()
    mesonEta     = mesonP4_chosen.Eta()
    mesonPhi     = mesonP4_chosen.Phi()


    #MESON MASS CALCULATION
    mesonMass = (firstTrkP4_chosen + secondTrkP4_chosen).M() 

    if verbose:
        print("m_TrkTrk =", mesonMass, "meson tree mass =", mesonTreeMass_chosen, "meson tree mass kin =", mesonTreeMassKin_chosen)
        print("pT trktrk =", mesonPt, "meson tree pT =", mesonTreePt_chosen, "eta trktrk =", mesonEta, "meson tree eta =", mesonTreeEta_chosen, "phi trktrk =", mesonPhi, "meson tree phi =", mesonTreePhi_chosen)


    #BOSON INV MASS CALCULATION
    bosonMassFrom2TrksPhoton = (firstTrkP4_chosen + secondTrkP4_chosen + ph_p4).M()

    #CANDIDATES SORTING
    if firstTrkPt < secondTrkPt:#swap-values loop, in order to fill the tree with the candidate with max pt of the couple in firstCand branches and one with min pt in secondCand branches
        a = firstTrkPt
        b = firstTrkEta
        c = firstTrkPhi
    
        firstTrkPt   = secondTrkPt
        firstTrkEta  = secondTrkEta
        firstTrkPhi  = secondTrkPhi
    
        secondTrkPt  = a
        secondTrkEta = b
        secondTrkPhi = c


    #CUTS ON CANDIDATES PT
    if firstTrkPt < 20. or secondTrkPt < 5.: 
        print("Final cut on candidates pT not passed, RETURN.")
        continue
  
    nEventsTrkPtFilter+=1
    histo_map["nEvents"].SetBinContent(5, nEventsTrkPtFilter)


    #ISOLATION DATAMEMBER FOR TREE FILLING
    mesonIsoCh = mesonIsoCh_chosen

    
  
    _firstTrkPt[0]               = firstTrkPt              
    _firstTrkEta[0]              = firstTrkEta              
    _firstTrkPhi[0]              = firstTrkPhi            
    _secondTrkPt[0]              = secondTrkPt            
    _secondTrkEta[0]             = secondTrkEta            
    _secondTrkPhi[0]             = secondTrkPhi           
    _isMeson[0]                  = isMeson_chosen              
    _mesonPt[0]                  = mesonPt             
    _mesonEta[0]                 = mesonEta              
    _mesonPhi[0]                 = mesonPhi            
    _mesonIsoCh[0]               = mesonIsoCh
    _mesonMassTrkTrk[0]          = mesonMass   
    _bosonMassFrom2TrksPhoton[0] = bosonMassFrom2TrksPhoton

    _mesonTreeMass[0]            = mesonTreeMass_chosen
    _mesonTreeMassKin[0]         = mesonTreeMassKin_chosen
    _mesonTreePt[0]              = mesonTreePt_chosen
    _mesonTreeEta[0]             = mesonTreeEta_chosen
    _mesonTreePhi[0]             = mesonTreePhi_chosen



    #//*************************************************************//
    #//                                                             //
    #//--------------------------MC Truth---------------------------//
    #//                                                             //
    #//*************************************************************//
    
      
    is_photon_matched = False
    is_meson_matched  = False
    is_boson_matched  = False    
    foundPhoton = False
    foundMeson = False
    genPhoton_eT = genPhoton_eta = genPhoton_phi = None  # Placeholder
    genMeson_pT = genMeson_eta = genMeson_phi = None

    if not runningOnData:
        for i in range(mytree.nGenPart):

            mother_idx = mytree.GenPart_genPartIdxMother[i]

            # Meson check (φ, ρ or K*)
            #if verbose: print(f"Indice: {i}, PDG ID: {mytree.GenPart_pdgId[i]}, madre: {mytree.GenPart_pdgId[mother_idx]}")

            if (mytree.GenPart_pdgId[i] in [333, 113, 313] and mother_idx >= 0 and mytree.GenPart_pdgId[mother_idx] in [23, 25]):             
                genMeson_pT = mytree.GenPart_pt[i]
                genMeson_eta = mytree.GenPart_eta[i]
                genMeson_phi = mytree.GenPart_phi[i]
                foundMeson = True
                if verbose: print(f"Trovato mesone! Indice: {i}, PDG ID: {mytree.GenPart_pdgId[i]}, madre: {mytree.GenPart_pdgId[mother_idx]}")

            # photon check with mother H or Z
            if (mytree.GenPart_pdgId[i] == 22 and mother_idx >= 0 and mytree.GenPart_pdgId[mother_idx] in [23, 25]):            
                genPhoton_eT = mytree.GenPart_pt[i]
                genPhoton_eta = mytree.GenPart_eta[i]
                genPhoton_phi = mytree.GenPart_phi[i]
                foundPhoton = True
                if verbose: print(f"Trovato fotone! Indice: {i}, madre: {mytree.GenPart_pdgId[mother_idx]}")


        if not foundPhoton:
            if verbose: print("Nessun fotone trovato in questo evento con madre Z o H")
            genPhoton_eT = genPhoton_eta = genPhoton_phi = -999  # Default value

        if not foundMeson:
            if verbose: print("Nessun mesone trovato in questo evento con madre Z o H")
            genMeson_pT = genMeson_eta = genMeson_phi = -999  # Default value


        # PHOTON MATCHING---------------------------------------------------
        if genPhoton_phi is not None:  
            deltaPhiPhoton = abs(photon_phi_chosen - genPhoton_phi)
            if deltaPhiPhoton > math.pi:
                deltaPhiPhoton = 2 * math.pi - deltaPhiPhoton

            deltaR_photonGenVsReco = math.sqrt((photon_eta_chosen - genPhoton_eta) ** 2 + deltaPhiPhoton ** 2)
            is_photon_matched = deltaR_photonGenVsReco < 0.2

        # MESON MATCHING----------------------------------------------------
        if genMeson_phi is not None: 
            deltaPhiMeson = abs(mesonTreePhi_chosen - genMeson_phi)
            if deltaPhiMeson > math.pi:
                deltaPhiMeson = 2 * math.pi - deltaPhiMeson

            deltaR_mesonGenVsReco = math.sqrt((mesonTreeEta_chosen - genMeson_eta) ** 2 + deltaPhiMeson ** 2)
            is_meson_matched = deltaR_mesonGenVsReco < 0.3

        # BOSON MATCHING---------------------------------------------------
        if is_photon_matched and is_meson_matched: 
            is_boson_matched = True
            print("************************* BOSON FOUND ****************************")
            nEventsBosonMatched += 1
        else:
            nEventsBosonNotMatched += 1
            print("THAT'S NOT A HIGGS or a Z!")

        # DEBUG PRINTS
        if verbose:
            print("Photon eT =", photonEtMax)
            print("Tracks pT =", firstTrkPt + secondTrkPt)
            print("Meson pT =", mesonTreePt_chosen)
            print("firstTrkPt =", firstTrkPt)
            print("secondTrkPt =", secondTrkPt)
            print("genMesonPt =", genMeson_pT)
            print("Tracks DeltaR =", deltaRTrks_chosen)
            print("Meson inv. mass from tracks =", mesonMass)
            print("Meson kin. mass =", mesonTreeMassKin_chosen)
            print("Boson inv. mass =", bosonMassFrom2TrksPhoton)
            print("---------------------------------------------------------------")
            print("MC H or Z FOUND =", nEventsBosonMatched, ", MC H or Z NOT FOUND =", nEventsBosonNotMatched)
            print("---------------------------------------------------------------")

            

    
    tree_output.Fill()
    
#Tree writing ##########################################################################################################
tree_output.Write()

#Set labels ############################################################################################################
histo_map["nEvents"].GetXaxis().SetBinLabel(1, "Processed")
histo_map["nEvents"].GetXaxis().SetBinLabel(2, "Triggered")
histo_map["nEvents"].GetXaxis().SetBinLabel(3, "Photon")
#histo_map["nEvents"].GetXaxis().SetBinLabel(4, "Meson iso")
histo_map["nEvents"].GetXaxis().SetBinLabel(4, "Best meson")
histo_map["nEvents"].GetXaxis().SetBinLabel(5, "Trks pT")



#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()