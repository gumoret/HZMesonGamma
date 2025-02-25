import ROOT
import argparse
import math
import numpy as np
import sys
from array import array
from ROOT import TLorentzVector


#Following bools are given as input
verbose = True

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("Events")

if not mytree:
    print("Error: Tree 'Events' not found in file", args.rootfile_name)
    sys.exit(1)


#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree ################################################################################
_genBoson_ID    = np.zeros(1, dtype=int)
_genBoson_pT    = np.zeros(1, dtype=float)
_genBoson_eta   = np.zeros(1, dtype=float)
_genBoson_phi   = np.zeros(1, dtype=float)
_genBoson_mass  = np.zeros(1, dtype=float)

_genMeson_ID    = np.zeros(1, dtype=int)
_genMeson_pT    = np.zeros(1, dtype=float)
_genMeson_eta   = np.zeros(1, dtype=float)
_genMeson_phi   = np.zeros(1, dtype=float)
_genMeson_mass  = np.zeros(1, dtype=float)

_genGamma_ID    = np.zeros(1, dtype=int)
_genGamma_pT    = np.zeros(1, dtype=float)
_genGamma_eta   = np.zeros(1, dtype=float)
_genGamma_phi   = np.zeros(1, dtype=float)

_genTrack1_ID   = np.zeros(1, dtype=int)
_genTrack1_pT   = np.zeros(1, dtype=float)
_genTrack1_eta  = np.zeros(1, dtype=float)
_genTrack1_phi  = np.zeros(1, dtype=float)
_genTrack1_mass = np.zeros(1, dtype=float)

_genTrack2_ID   = np.zeros(1, dtype=int)
_genTrack2_pT   = np.zeros(1, dtype=float)
_genTrack2_eta  = np.zeros(1, dtype=float)
_genTrack2_phi  = np.zeros(1, dtype=float)
_genTrack2_mass = np.zeros(1, dtype=float)


tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('genBoson_ID', _genBoson_ID, '_genBoson_ID/I')
tree_output.Branch('genBoson_pT', _genBoson_pT, '_genBoson_pT/D')
tree_output.Branch('genBoson_eta', _genBoson_eta, '_genBoson_eta/D')
tree_output.Branch('genBoson_phi', _genBoson_phi, '_genBoson_phi/D')
tree_output.Branch('genBoson_mass', _genBoson_mass, '_genBoson_mass/D')

tree_output.Branch('genMeson_ID', _genMeson_ID, '_genMeson_ID/I')
tree_output.Branch('genMeson_pT', _genMeson_pT, '_genMeson_pT/D')
tree_output.Branch('genMeson_eta', _genMeson_eta, '_genMeson_eta/D')
tree_output.Branch('genMeson_phi', _genMeson_phi, '_genMeson_phi/D')
tree_output.Branch('genMeson_mass', _genMeson_mass, '_genMeson_mass/D')

tree_output.Branch('genGamma_ID', _genGamma_ID, '_genGamma_ID/I')
tree_output.Branch('genGamma_pT', _genGamma_pT, '_genGamma_pT/D')
tree_output.Branch('genGamma_eta', _genGamma_eta, '_genGamma_eta/D')
tree_output.Branch('genGamma_phi', _genGamma_phi, '_genGamma_phi/D')

tree_output.Branch('genTrack1_ID', _genTrack1_ID, '_genTrack1_ID/I')
tree_output.Branch('genTrack1_pT', _genTrack1_pT, '_genTrack1_pT/D')
tree_output.Branch('genTrack1_eta', _genTrack1_eta, '_genTrack1_eta/D')
tree_output.Branch('genTrack1_phi', _genTrack1_phi, '_genTrack1_phi/D')
tree_output.Branch('genTrack1_mass', _genTrack1_mass, '_genTrack1_mass/D')

tree_output.Branch('genTrack2_ID', _genTrack2_ID, '_genTrack2_ID/I')
tree_output.Branch('genTrack2_pT', _genTrack2_pT, '_genTrack2_pT/D')
tree_output.Branch('genTrack2_eta', _genTrack2_eta, '_genTrack2_eta/D')
tree_output.Branch('genTrack2_phi', _genTrack2_phi, '_genTrack2_phi/D')
tree_output.Branch('genTrack2_mass', _genTrack2_mass, '_genTrack2_mass/D')

#EVENTS LOOP #######################################################################################################
nentries = mytree.GetEntriesFast()

genBoson_ID, genBoson_pT, genBoson_eta, genBoson_phi, genBoson_mass = [-999] * 5
genMeson_ID, genMeson_pT, genMeson_eta, genMeson_phi, genMeson_mass = [-999] * 5
genGamma_ID, genGamma_pT, genGamma_eta, genGamma_phi = [-999] * 4
genTrack1_ID, genTrack1_pT, genTrack1_eta, genTrack1_phi, genTrack1_mass = [-999] * 5
genTrack2_ID, genTrack2_pT, genTrack2_eta, genTrack2_phi, genTrack2_mass = [-999] * 5

boson_not_found_count = 0

for jentry in range(nentries):
    ientry = mytree.LoadTree(jentry)
    if ientry < 0: break
    nb = mytree.GetEntry(jentry)
    if nb <= 0: continue  
    print("EVENT", jentry, " NUMBER OF GEN PARTICLES = ", mytree.nGenPart)        

    # Create mother->daughters map for this event event
    mother_to_daughters = {}
    meson_to_daughters = {}  # meson map

    for d in range(mytree.nGenPart):
        mother_idx = mytree.GenPart_genPartIdxMother[d]
        if 0 <= mother_idx < mytree.nGenPart:
            if mother_idx not in mother_to_daughters: mother_to_daughters[mother_idx] = []
            mother_to_daughters[mother_idx].append(d)

            # Create mother->daughters map for mesons 
            if mytree.GenPart_pdgId[mother_idx] in [333, 113, 313]:
                if mother_idx not in meson_to_daughters: meson_to_daughters[mother_idx] = []
                meson_to_daughters[mother_idx].append(d)  


        print("   - particle number ", d, "  ID = ", mytree.GenPart_pdgId[d], "  mother ID = ", mytree.GenPart_pdgId[mother_idx]  )


    for i in range(mytree.nGenPart):
        daughters = mother_to_daughters.get(i, [])  # retrieve daughters

        #if it is not a Higgs or a Z with 2 daughters (sometimes Pythia sends a Higgs in itself, therefore only 1 daughter), continue
        if mytree.GenPart_pdgId[i] not in [23,25] or len(daughters) != 2: continue           
            
        meson_idx, gamma_idx = None, None

        #for each daughter
        for daughter in daughters:
            if mytree.GenPart_pdgId[daughter] == 22: gamma_idx = daughter
            elif mytree.GenPart_pdgId[daughter] in [333, 113, 313]: meson_idx = daughter

        #if daughters are not Phi or Rho or K0* and gamma, continue
        if meson_idx is None or gamma_idx is None: continue
            


        #save boson variables
        genBoson_ID = mytree.GenPart_pdgId[i]
        genBoson_pT = mytree.GenPart_pt[i]
        genBoson_eta = mytree.GenPart_eta[i]
        genBoson_phi = mytree.GenPart_phi[i]
        genBoson_mass = mytree.GenPart_mass[i]

        #if genBoson_ID == -999: boson_not_found_count += 1

        #save Gamma variables
        genGamma_ID = mytree.GenPart_pdgId[gamma_idx]
        genGamma_pT = mytree.GenPart_pt[gamma_idx]
        genGamma_eta = mytree.GenPart_eta[gamma_idx]
        genGamma_phi = mytree.GenPart_phi[gamma_idx]

        #save meson variables
        genMeson_ID = mytree.GenPart_pdgId[meson_idx]
        genMeson_pT = mytree.GenPart_pt[meson_idx]
        genMeson_eta = mytree.GenPart_eta[meson_idx]
        genMeson_phi = mytree.GenPart_phi[meson_idx]
        genMeson_mass = mytree.GenPart_mass[meson_idx] 


        meson_daughters = meson_to_daughters.get(meson_idx, [])#retrieve mesons daughters      

        #if meson has not two daughters, continue
        if len(meson_daughters) != 2: continue

        #for each Meson daughter
        for md in meson_daughters:
            #if meson daughter is a K- or a Pi-
            if mytree.GenPart_pdgId[md] in [-321, -211]:
                genTrack1_ID = mytree.GenPart_pdgId[md]
                genTrack1_pT = mytree.GenPart_pt[md]
                genTrack1_eta = mytree.GenPart_eta[md]
                genTrack1_phi = mytree.GenPart_phi[md]
                genTrack1_mass = mytree.GenPart_mass[md]

            #if meson daughter is a K+ or a Pi+
            elif mytree.GenPart_pdgId[md] in [321, 211]:
                genTrack2_ID = mytree.GenPart_pdgId[md]
                genTrack2_pT = mytree.GenPart_pt[md]
                genTrack2_eta = mytree.GenPart_eta[md]
                genTrack2_phi = mytree.GenPart_phi[md]
                genTrack2_mass = mytree.GenPart_mass[md]


    _genBoson_ID[0] = genBoson_ID
    _genBoson_pT[0] = genBoson_pT
    _genBoson_eta[0] = genBoson_eta
    _genBoson_phi[0] = genBoson_phi
    _genBoson_mass[0] = genBoson_mass

    _genMeson_ID[0] = genMeson_ID
    _genMeson_pT[0] = genMeson_pT
    _genMeson_eta[0] = genMeson_eta
    _genMeson_phi[0] = genMeson_phi
    _genMeson_mass[0] = genMeson_mass

    _genGamma_ID[0] = genGamma_ID
    _genGamma_pT[0] = genGamma_pT
    _genGamma_eta[0] = genGamma_eta
    _genGamma_phi[0] = genGamma_phi

    _genTrack1_ID[0] = genTrack1_ID
    _genTrack1_pT[0] = genTrack1_pT
    _genTrack1_eta[0] = genTrack1_eta
    _genTrack1_phi[0] = genTrack1_phi
    _genTrack1_mass[0] = genTrack1_mass

    _genTrack2_ID[0] = genTrack2_ID
    _genTrack2_pT[0] = genTrack2_pT
    _genTrack2_eta[0] = genTrack2_eta
    _genTrack2_phi[0] = genTrack2_phi
    _genTrack2_mass[0] = genTrack2_mass


    tree_output.Fill()

    if genBoson_ID == -999: boson_not_found_count+=1

    print("bosons with ID = -999:", boson_not_found_count)
    print("boson ID =", genBoson_ID, "boson mass =", genBoson_mass)
    print("meson ID =", genMeson_ID, "meson mass =", genMeson_mass)
    
#Tree writing ##########################################################################################################
tree_output.Write()

fOut.Close()