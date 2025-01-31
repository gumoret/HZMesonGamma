import ROOT
import argparse
import math
import numpy as np
import sys
from array import array
from ROOT import TLorentzVector


#Following bools are given as input
verbose       = True
isPhiAnalysis = False # for Z -> Phi Gamma
isRhoAnalysis = False # for Z -> Rho Gamma

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('meson_option', help='Type <<rho>> for rho, <<phi>> for phi') #flag for type of meson
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("Events")

if not mytree:
    print("Error: Tree 'Events' not found in file", args.rootfile_name)
    sys.exit(1)

if args.meson_option == "phi":
    isPhiAnalysis = True
elif args.meson_option == "rho":
    isRhoAnalysis = True
else:
    print("Error: meson_option must be 'phi' or 'rho'")
    sys.exit(1)


#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["nEvents", "pileup"]

histo_map[list_histos[0]] = ROOT.TH1F(list_histos[0],"Event counting in different steps", 6, 0., 6.)
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








#EVENTS LOOP #######################################################################################################
nEventsProcessed            = 0
nEventsTriggered            = 0
nEventsIsPhoton             = 0
nEventsMesonIsolationFilter = 0
nEventsBestPairFound        = 0
nEventsTrkPtFilter          = 0

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

        if electron_pT < 10. or abs(electron_eta) > 2.5 or abs(muon_dxy) >= 0.2 or abs(electron_dz) >= 0.5: continue

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

        if verbose: print(f"Photon {i}: pT={photon_pT}, eta={photon_eta}, WP80={isMVAIso_WP80}, WP90={isMVAIso_WP90}, veto={isElectronVeto}")

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

    if verbose: print(f"Chosen photon: pT={photonEtMax}, eta={photon_eta_chosen}")
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
    
    isBestMesonOfTheEventFound = False
    bestMesonPt = 0.

    for i in range(len(mytree.rho_mass)):
        if isRhoAnalysis:
            firstTrkPt   = mytree.rho_trk1_pt[i]
            firstTrkEta  = mytree.rho_trk1_eta[i]
            firstTrkPhi  = mytree.rho_trk1_phi[i]
            secondTrkPt  = mytree.rho_trk2_pt[i]
            secondTrkEta = mytree.rho_trk2_eta[i]
            secondTrkPhi = mytree.rho_trk2_phi[i]
            isoMesonCh   = mytree.rho_iso[i]

            mesonTreeMass    = mytree.rho_mass[i]
            mesonTreeMassKin = mytree.rho_kin_mass[i]
            mesonTreePt      = mytree.rho_kin_pt[i]
            mesonTreeEta     = mytree.rho_kin_eta[i] 
            mesonTreePhi     = mytree.rho_kin_phi[i]

        elif isPhiAnalysis:
            firstTrkPt   = mytree.phi_trk1_pt[i]
            firstTrkEta  = mytree.phi_trk1_eta[i]
            firstTrkPhi  = mytree.phi_trk1_phi[i]
            secondTrkPt  = mytree.phi_trk2_pt[i]
            secondTrkEta = mytree.phi_trk2_eta[i]
            secondTrkPhi = mytree.phi_trk2_phi[i]
            isoMesonCh   = mytree.phi_iso[i]

            mesonTreeMass    = mytree.phi_mass[i]
            mesonTreeMassKin = mytree.phi_kin_mass[i]
            mesonTreePt      = mytree.phi_kin_pt[i]
            mesonTreeEta     = mytree.phi_kin_eta[i] 
            mesonTreePhi     = mytree.phi_kin_phi[i] 
            

        #Tracks pt cuts------------------------------------------------
        if firstTrkPt < 1. or secondTrkPt < 1.:
            if verbose: print("if first trk pt or second trk pt<1 continue")
            continue

        if firstTrkPt < 10. and secondTrkPt < 10.:
            if verbose: print("if first trk pt and second trk pt<10 continue")
            continue

        #DiTrks deltaR cut---------------------------------------------------
        deltaEta = firstTrkEta - secondTrkEta
        deltaPhi = abs(firstTrkPhi - secondTrkPhi)
        if deltaPhi > math.pi:
            deltaPhi = 2*math.pi - deltaPhi

        deltaRTrks = math.sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi)        
        if deltaRTrks > 0.07:
            if verbose: print("if deltaRTrks>0.07 continue")
            continue

        #Quadrimomentum calculation -------------------------------------------
        firstTrkP4 = TLorentzVector() 
        secondTrkP4 = TLorentzVector()

        if isRhoAnalysis:
            firstTrkP4.SetPtEtaPhiM(firstTrkPt, firstTrkEta, firstTrkPhi, 0.13957) #Pi mass            
            secondTrkP4.SetPtEtaPhiM(secondTrkPt, secondTrkEta, secondTrkPhi, 0.13957)

        pairP4 = firstTrkP4  + secondTrkP4

        if verbose: print(f"PiPi pT = {pairP4.Pt()}")

        #DiTrk pT cut----------------------------------------------------------
        if pairP4.Pt() < 38:
            if verbose: print("couplePt cut NOT passed, continue")
            continue

        #Meson inv mass---------------------------------------------------------
        isMeson = False

        mesonMassTrkTrk = pairP4.M()        

        if isRhoAnalysis:
            if 0.5 < mesonMassTrkTrk < 1.: isMeson = True
            
        elif isPhiAnalysis:
            if 1.0 < mesonMassTrkTrk < 1.05: isMeson = True

        if not isMeson: continue

        #Isolation cut----------------------------------------------------------
        if isoMesonCh < 0.9:
            print("No isolation cut passed, continue")
            continue

        nEventsMesonIsolationFilter+=1
        histo_map["nEvents"].SetBinContent(4, nEventsMesonIsolationFilter)

        #pT max of the event filter ----------------------------------------------
        if verbose: print(f"Current bestMeson_Pt = {bestMesonPt}")

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


    #meson loop end------------------



    if not isBestMesonOfTheEventFound:
        print("No best couple detected for current event, RETURN.")
        continue  

    nEventsBestPairFound+=1
    histo_map["nEvents"].SetBinContent(5, nEventsBestPairFound)



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
        print(f"m_TrkTrk = {mesonMass}, meson tree mass = {mesonTreeMass_chosen}, meson tree mass kin = {mesonTreeMassKin_chosen} ")
        print(f"pT trktrk = {mesonPt}, meson tree pT = {mesonTreePt_chosen}, eta trktrk = {mesonEta}, meson tree eta = {mesonTreeEta_chosen}, phi trktrk = {mesonPhi}, meson tree phi = {mesonTreePhi_chosen}")


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
    histo_map["nEvents"].SetBinContent(6, nEventsTrkPtFilter)


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




    
    tree_output.Fill()
    
#Tree writing ##########################################################################################################
tree_output.Write()

#Set labels ############################################################################################################
histo_map["nEvents"].GetXaxis().SetBinLabel(1, "Processed")
histo_map["nEvents"].GetXaxis().SetBinLabel(2, "Triggered")
histo_map["nEvents"].GetXaxis().SetBinLabel(3, "Photon")
histo_map["nEvents"].GetXaxis().SetBinLabel(4, "Meson iso")
histo_map["nEvents"].GetXaxis().SetBinLabel(5, "Best meson")
histo_map["nEvents"].GetXaxis().SetBinLabel(6, "Trks pT")



#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()