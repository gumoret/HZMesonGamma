import ROOT
import argparse
import math
import numpy as np
import sys
from array import array


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
list_histos = ["h_MesonMass", "h_firstTrkPt", "h_secondTrkPt", "h_mesonPt", "h_photonPt", "h_photonEta", "h_photonPhi"] 

if   isPhiAnalysis: histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{meson}", 100, 1., 1.05) 
elif isRhoAnalysis: histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{meson}", 100, 0.5, 1.) 
histo_map[list_histos[1]] = ROOT.TH1F(list_histos[1],"p_{T} of the 1st track", 100, 0.1, 70.)
histo_map[list_histos[2]] = ROOT.TH1F(list_histos[2],"p_{T} of the 2nd track", 100, 0.1, 70.)
histo_map[list_histos[3]] = ROOT.TH1F(list_histos[3],"p_{T} of the meson", 100, 0.1, 140.)
histo_map[list_histos[4]] = ROOT.TH1F(list_histos[4],"p_{T} of the #gamma", 100, 0.1, 140.)
histo_map[list_histos[5]] = ROOT.TH1F(list_histos[5],"#eta_{#gamma}", 100, -2.5, 2.5)
histo_map[list_histos[6]] = ROOT.TH1F(list_histos[6],"#phi_{gamma}", 100, -math.pi, math.pi)


#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree ################################################################################
_mesonMass      = np.zeros(1, dtype=float)
_firstTrkPt     = np.zeros(1, dtype=float)
_secondTrkPt    = np.zeros(1, dtype=float)
_bestPairPt     = np.zeros(1, dtype=float)
_photonEt       = np.zeros(1, dtype=float)
_photonEta      = np.zeros(1, dtype=float) 


tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('mesonMass',_mesonMass,'_mesonMass/D')
tree_output.Branch('firstTrkPt',_firstTrkPt,'_firstTrkPt/D')
tree_output.Branch('secondTrkPt',_secondTrkPt,'_secondTrkPt/D')
tree_output.Branch('mesonPt',_bestPairPt,'_bestPairPt/D')
tree_output.Branch('photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('photonEta',_photonEta,'_photonEta/D')



#EVENTS LOOP #######################################################################################################
nentries = mytree.GetEntriesFast()

# Loop sugli eventi
for jentry in range(nentries):
    ientry = mytree.LoadTree(jentry)
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry)
    if nb <= 0:
        continue

    #Retrieve variables from the tree 
    if isPhiAnalysis:
        genMeson_ID   = mytree.phi_gen_pdgId
        genMeson_pT   = mytree.phi_gen_pt
        genMeson_mass = mytree.phi_gen_mass

        genTrk1_ID    = mytree.phi_gen_trk1_mpdgId
        genTrk1_pT    = mytree.phi_gen_trk1_pt

        genTrk2_ID    = mytree.phi_gen_trk2_pdgId
        genTrk2_pT    = mytree.phi_gen_trk2_pt

    if isRhoAnalysis:
        genMeson_ID   = mytree.rho_gen_pdgId
        genMeson_pT   = mytree.rho_gen_pt
        genMeson_mass = mytree.rho_gen_mass

        genTrk1_ID    = mytree.rho_gen_trk1_mpdgId
        genTrk1_pT    = mytree.rho_gen_trk1_pt

        genTrk2_ID    = mytree.rho_gen_trk2_pdgId
        genTrk2_pT    = mytree.rho_gen_trk2_pt

    photon_pT     = mytree.Photon_pt
    photon_eta    = mytree.Photon_eta
    photon_phi    = mytree.Photon_phi

    eventWeight = 1.

    #FILL HISTOS AND TREE #####################################################################################################
    if genMeson_mass and len(genMeson_mass) > 0:
        for i in range(len(genMeson_mass)):
            genMeson_mass_value = float(genMeson_mass[i])
            histo_map["h_MesonMass"].Fill(genMeson_mass_value, eventWeight)  # Riempie l'istogramma
            if verbose:
                print(f"Event {jentry}: genMeson_mass_value = {genMeson_mass_value} (type: {type(genMeson_mass_value)})")

            # Riempie il tree
            _mesonMass[0] = genMeson_mass_value
            tree_output.Fill()


    if genTrk1_pT and len(genTrk1_pT) > 0:
        for i in range(len(genTrk1_pT)):
            genTrk1_pT_value = float(genTrk1_pT[i])
            histo_map["h_firstTrkPt"].Fill(genTrk1_pT_value, eventWeight)  # Riempie l'istogramma
            if verbose:
                print(f"Event {jentry}: genTrk1_pT_value = {genTrk1_pT_value} (type: {type(genTrk1_pT_value)})")

            # Riempie il tree
            _firstTrkPt[0] = genTrk1_pT_value
            tree_output.Fill()

    
    if genTrk2_pT and len(genTrk2_pT) > 0:
        for i in range(len(genTrk2_pT)):
            genTrk2_pT_value = float(genTrk2_pT[i])
            histo_map["h_secondTrkPt"].Fill(genTrk2_pT_value, eventWeight)  # Riempie l'istogramma
            if verbose:
                print(f"Event {jentry}: genTrk2_pT_value = {genTrk2_pT_value} (type: {type(genTrk2_pT_value)})")            
            
            # Riempie il tree
            _secondTrkPt[0] = genTrk2_pT_value
            tree_output.Fill()
    

    if genMeson_pT and len(genMeson_pT) > 0:
        for i in range(len(genMeson_pT)):
            genMeson_pT_value = float(genMeson_pT[i])
            histo_map["h_mesonPt"].Fill(genMeson_pT_value, eventWeight)  # Riempie l'istogramma
            if verbose:
                print(f"Event {jentry}: genMeson_pT_value = {genMeson_pT_value} (type: {type(genMeson_pT_value)})")            
            
            # Riempie il tree
            _bestPairPt[0] = genMeson_pT_value
            tree_output.Fill()


    if photon_pT and len(photon_pT) > 0:
        for i in range(len(photon_pT)):
            photon_pT_value = float(photon_pT[i])
            histo_map["h_photonPt"].Fill(photon_pT_value, eventWeight)  # Riempie l'istogramma
            if verbose:
                print(f"Event {jentry}: photon_pT_value = {photon_pT_value} (type: {type(photon_pT_value)})")            
            
            # Riempie il tree
            _photonEt[0] = photon_pT_value
            tree_output.Fill()


    if photon_eta and len(photon_eta) > 0:
        for i in range(len(photon_eta)):
            photon_eta_value = float(photon_eta[i])
            histo_map["h_photonEta"].Fill(photon_eta_value, eventWeight)  # Riempie l'istogramma
            if verbose:
                print(f"Event {jentry}: photon_eta_value = {photon_eta_value} (type: {type(photon_eta_value)})")            
            
            # Riempie il tree
            _photonEta[0] = photon_eta_value
            tree_output.Fill()


    if photon_phi and len(photon_phi) > 0:
        for i in range(len(photon_phi)):
            photon_phi_value = float(photon_phi[i])
            histo_map["h_photonPhi"].Fill(photon_phi_value, eventWeight)  # Riempie l'istogramma
            if verbose:
                print(f"Event {jentry}: photon_phi_value = {photon_phi_value} (type: {type(photon_phi_value)})")            
            
            # Riempie il tree
            #_photon[0] = photon_eta_value
            #tree_output.Fill()




#HISTO LABELS #########################################################################################################
histo_map["h_MesonMass"].GetXaxis().SetTitle("m_{meson } [GeV/c^2]")
histo_map["h_MesonMass"].SetTitle("Meson mass")
histo_map["h_firstTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_firstTrkPt"].SetTitle("Transverse momentum of the first charged particle ")
histo_map["h_secondTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_secondTrkPt"].SetTitle("Transverse momentum of the second charged particle")
histo_map["h_mesonPt"].GetXaxis().SetTitle("pT_{meson} [GeV/c]") 
histo_map["h_mesonPt"].SetTitle("Transverse momentum of the meson")
histo_map["h_photonPt"].GetXaxis().SetTitle("p_{T}^{#gamma} [GeV/c]")
histo_map["h_photonPt"].SetTitle("Transverse momentum of the photon")
histo_map["h_photonEta"].GetXaxis().SetTitle("#eta_{#gamma}")
histo_map["h_photonEta"].SetTitle("Eta of the photon")
histo_map["h_photonPhi"].GetXaxis().SetTitle("#phi_{#gamma} [rad]")
histo_map["h_photonPhi"].SetTitle("phi of the photon")


            
#Tree writing ##########################################################################################################
tree_output.Write()



#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()