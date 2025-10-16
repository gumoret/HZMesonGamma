import ROOT
import argparse

# multithreading

# following bools are given as input
verbose       = True
isPhiAnalysis = False # for H -> Phi Gamma
isRhoAnalysis = False # for H -> Rho Gamma

# PARSER and INPUT 
p = argparse.ArgumentParser(description="RDataFrame analyzer for H→meson+γ")
p.add_argument("meson_option", help="Type <<rho>> for rho, <<phi>> for phi")
p.add_argument("runningOnData_option", help="Type <<signal>> for signal, <<data>> for data")
p.add_argument("rootfile_name", help="Input nanoAOD ROOT file")
p.add_argument("outputfile_option", help="Output ROOT file")
args = p.parse_args()

if args.meson_option == "phi": isPhiAnalysis = True
elif args.meson_option == "rho": isRhoAnalysis = True
else: print("meson_option must be <<phi>> or <<rho>> or <<K*>> or <<DO*>>")


if args.runningOnData_option == "signal": runningOnData = False
elif args.runningOnData_option == "data": runningOnData = True
else: print("runninOnData must be <<signal>> or <<data>>")


input_file = args.rootfile_name
output_file = args.outputfile_option
input_tree_name = "Events"


# ------------------------------------------------------------
# C++ functions 
# ------------------------------------------------------------
ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;

// muons
RVec<int> muon_counts(const RVec<float>& pt, const RVec<float>& eta, 
                      const RVec<bool>& mediumId, const RVec<float>& dxy, 
                      const RVec<float>& dz, const RVec<int>& pfIsoId) {

    auto mask10 = (pt > 10) && (abs(eta) < 2.4) && mediumId && (abs(dxy) < 0.2) && (abs(dz) < 0.5) && (pfIsoId == 2);
    auto mask20 = (pt > 20) && (abs(eta) < 2.4) && mediumId && (abs(dxy) < 0.2) && (abs(dz) < 0.5) && (pfIsoId == 2);

    return {Sum(mask10), Sum(mask20)};
}

// ------------------------------------------------------------------------------------------------

// electrons
RVec<int> electron_counts(const RVec<float>& pt, const RVec<float>& eta,
                          const RVec<float>& dxy, const RVec<float>& dz,
                          const RVec<bool>& mvaIso_WP80) {

    auto mask10 = (pt > 10) && (abs(eta) < 2.5) && (abs(dxy) < 0.2) && (abs(dz) < 0.5) && mvaIso_WP80;
    auto mask20 = (pt > 20) && (abs(eta) < 2.5) && (abs(dxy) < 0.2) && (abs(dz) < 0.5) && mvaIso_WP80;
    
    return {Sum(mask10), Sum(mask20)};
}

//---------------------------------------------------------------------------------------------------

// photons
struct PhotonSelection { //struct to define which values to give as output
    int nGood;
    int bestIdx;
};

PhotonSelection select_photons(const RVec<float>& pt, const RVec<float>& eta, 
                               const RVec<bool>& mvaID_WP80, const RVec<bool>& electronVeto) { //function of type "struct"
    
    auto mask = (pt > 35) && (abs(eta) < 2.5) && mvaID_WP80 && electronVeto; //photon cuts

    PhotonSelection output;
    output.nGood = Sum(mask); //first output to return

    // if no good photon, bestIdx = -1
    if (output.nGood == 0) {
        output.bestIdx = -1;
        return output;
    }

    // max pT photon index
    auto good_pts = pt[mask]; // subvector with good photons only
    auto best_pt_local = ArgMax(good_pts); //best photon index in the best photons collection
    auto good_idx = Nonzero(mask);  // original photon indexes
    output.bestIdx = good_idx[best_pt_local]; // best photon index in the original collection. second output to return
    return output;
}

#include "TLorentzVector.h"

TLorentzVector make_photon_p4(const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi, int idx) {

    TLorentzVector p4;
    if (idx < 0 || idx >= (int)pt.size()) return p4; // returns (0,0,0,0)
    p4.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], 0.);
    return p4;
}
""")

# Create the dataframe
df = ROOT.RDataFrame(input_tree_name, input_file)
n_total = df.Count().GetValue()

# ---------------------------------------------------------------------
# Trigger
# ---------------------------------------------------------------------
df = df.Filter("HLT_Photon35_TwoProngs35", "Two-prong photon trigger")
n_trigger = df.Filter("HLT_Photon35_TwoProngs35").Count().GetValue()

# ------------------------------------------------------------
# Muons
# ------------------------------------------------------------
df = df.Define("muon_counts", "muon_counts(Muon_pt, Muon_eta, Muon_mediumId, Muon_dxy, Muon_dz, Muon_pfIsoId)")
df = df.Define("nMuons10", "muon_counts[0]")
df = df.Define("nMuons20", "muon_counts[1]")

# ------------------------------------------------------------
# Electrons
# ------------------------------------------------------------
df = df.Define("electron_counts", "electron_counts(Electron_pt, Electron_eta, Electron_dxy, Electron_dz, Electron_mvaIso_WP80)")
df = df.Define("nElectrons10", "electron_counts[0]")
df = df.Define("nElectrons20", "electron_counts[1]")

# ------------------------------------------------------------
# Photons
# ------------------------------------------------------------
df = df.Define("photon_sel", "select_photons(Photon_pt, Photon_eta, Photon_mvaID_WP80, Photon_electronVeto)")
df = df.Define("nGoodPhotons", "photon_sel.nGood")
df = df.Define("bestPhotonIdx", "photon_sel.bestIdx")
df = df.Filter("nGoodPhotons > 0", "At least one good photon")
n_photon = df.Filter("nGoodPhotons").Count().GetValue()

# pT, eta, phi, p4 of the selected photon (if it exists)
df = df.Define("bestPhoton_pt",  "(float)bestPhotonIdx >= 0 ? Photon_pt[bestPhotonIdx]  : -1.f")
df = df.Define("bestPhoton_eta", "(float)bestPhotonIdx >= 0 ? Photon_eta[bestPhotonIdx] : -999.f")
df = df.Define("bestPhoton_phi", "(float)bestPhotonIdx >= 0 ? Photon_phi[bestPhotonIdx] : -999.f")
df = df.Define("photon_p4", "make_photon_p4(Photon_pt, Photon_eta, Photon_phi, bestPhotonIdx)")


# ---------------------------------------------------------------------
# TTree writing
# ---------------------------------------------------------------------
columns_to_save = ["HLT_Photon35_TwoProngs35", "nMuons10", "nMuons20", "nElectrons10", "nElectrons20", 
                   "nGoodPhotons", "bestPhoton_pt", "bestPhoton_eta", "bestPhoton_phi"]

df.Snapshot("tree_output", output_file, columns_to_save)

# ---------------------------------------------------------------------
# efficiency histogram
# --------------------------------------------------------------------
cut_names = ["All events", "Triggered", "Best photon"]

h_cutflow = ROOT.TH1F("cutflow", "Cutflow;Selection;Events", len(cut_names), 0, len(cut_names))
for i, (label, val) in enumerate(zip(cut_names, [n_total, n_trigger, n_photon]), start=1):
    h_cutflow.SetBinContent(i, val)
    h_cutflow.GetXaxis().SetBinLabel(i, label)

# histo writing
output = ROOT.TFile(output_file, "UPDATE")
output.cd()
h_cutflow.Write()
output.Close()

print("Output saved in:", output_file)

# ---------------------------------------------------------------------
# DEBUG PRINTS
# ---------------------------------------------------------------------
if verbose:
    all_photon_pts  = df.Take["ROOT::VecOps::RVec<float>"]("Photon_pt")
    best_photon_pts = df.Take["float"]("bestPhoton_pt")

    for i, (a, b) in enumerate(zip(all_photon_pts, best_photon_pts)):
        print(f"Event {i}: Photon pTs = {list(a)},  best = {b}")