import ROOT
import argparse

# multithreading
ROOT.ROOT.EnableImplicitMT()

# following bools are given as input
verbose       = True
debug         = False
isPhiAnalysis = False # for H -> Phi Gamma
isRhoAnalysis = False # for H -> Rho Gamma
ismesonFromTracks = False # debug for reconstructing meson from tracks

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
#include "TLorentzVector.h"

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
                               const RVec<bool>& mvaID_WP80, const RVec<bool>& electronVeto) { //function of type "PhotonSelection"
    
    auto mask = (pt > 35) && (abs(eta) < 2.5) && mvaID_WP80 && electronVeto; //photon cuts

    PhotonSelection output; //define output of type PhotonSelection to be returned

    output.nGood = Sum(mask); //first output to return

    // if no good photon, bestIdx = -1
    if (output.nGood == 0) {
        output.bestIdx = -1;
        return output;
    }

    // Choose best photon = highest pT among selected
    auto good_pts = pt[mask]; // subvector with good photons only
    auto best_pt_local = ArgMax(good_pts); //best photon index in the best photons collection
    auto good_idx = Nonzero(mask);  // original photons indexes
    output.bestIdx = good_idx[best_pt_local]; // best photon index in the original collection. second output to return
    return output;
}

TLorentzVector make_photon_p4(const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi, int idx) {
    TLorentzVector p4;
    if (idx < 0 || idx >= (int)pt.size()) return p4; // returns (0,0,0,0)
    p4.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], 0.);
    return p4;
}

//-------------------------------------------------------------------------------------------------------

// mesons
struct MesonSelection {
    int nGood;
    int bestIdx;
};

MesonSelection select_mesons_kin(const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi,
                             const RVec<float>& mass, const RVec<float>& iso, float mass_low, float mass_high,
                             const RVec<float>& trk1_pt, const RVec<float>& trk1_eta, const RVec<float>& trk1_phi,
                             const RVec<float>& trk2_pt, const RVec<float>& trk2_eta, const RVec<float>& trk2_phi) {

    // ΔR(trk1, trk2)
    auto dEta = trk1_eta - trk2_eta;
    auto dPhi = ROOT::VecOps::Map(trk1_phi - trk2_phi, [](float dphi) {
        if (fabs(dphi) > M_PI) dphi = 2 * M_PI - fabs(dphi);
        return dphi;
    });
    auto deltaR = sqrt(dEta * dEta + dPhi * dPhi);
    
    // Track-level cuts
    auto trk_mask = !((trk1_pt < 1.0f) || (trk2_pt < 1.0f)) && !((trk1_pt < 10.0f) && (trk2_pt < 10.0f)) && (deltaR < 0.07f);
        
    // Meson cuts: mass window, iso > 0.9, pt > 38, |eta| < 2.5
    auto meson_mask = (pt > 38) && (abs(eta) < 2.5) && (mass > mass_low) && (mass < mass_high) && (iso > 0.9);

    auto mask = trk_mask && meson_mask;

    MesonSelection out; //define output
    out.nGood = Sum(mask); //first output to return

    if (out.nGood == 0) {
        out.bestIdx = -1;
        return out;
    }

    // Choose best meson = highest pT among selected
    auto good_pts = pt[mask]; // subvector with good mesons only
    auto best_pt_local = ArgMax(good_pts); //best meson index in the best mesons collection
    auto good_idx = Nonzero(mask); // original mesons indexes
    out.bestIdx = good_idx[best_pt_local]; // best meson index in the original collection. second output to return
    return out;
}

TLorentzVector make_meson_p4(const RVec<float>& pt, const RVec<float>& eta,
                             const RVec<float>& phi, const RVec<float>& mass, int idx) {

    TLorentzVector p4;
    if (idx < 0 || idx >= (int)pt.size()) return p4; // returns (0,0,0,0)
    p4.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], mass[idx]);
    return p4;
}

// debug: meson reconstruction using tracks instead of nanoaod variables

MesonSelection select_mesons_tracks(const RVec<float>& iso, 
                                    const RVec<float>& trk1_pt, const RVec<float>& trk1_eta, const RVec<float>& trk1_phi,
                                    const RVec<float>& trk2_pt, const RVec<float>& trk2_eta, const RVec<float>& trk2_phi,
                                    float mass_trk1, float mass_trk2, float mass_low, float mass_high, float dR_max, float pt_pair_min){

    MesonSelection out{0, -1};
    const int N = trk1_pt.size();
    float bestPt = -1.f;

    for (int i = 0; i < N; ++i) {

        // Track pT cuts
        if (trk1_pt[i] < 1.f || trk2_pt[i] < 1.f) continue;
        if (trk1_pt[i] < 10.f && trk2_pt[i] < 10.f) continue;

        // ΔR
        float dphi = fabs(trk1_phi[i] - trk2_phi[i]);
        if (dphi > M_PI) dphi = 2.f * M_PI - dphi;
        float deta = trk1_eta[i] - trk2_eta[i];
        float dR = sqrt(deta*deta + dphi*dphi);
        if (dR > dR_max) continue;

        // Build p4 from track
        TLorentzVector t1, t2;
        t1.SetPtEtaPhiM(trk1_pt[i], trk1_eta[i], trk1_phi[i], mass_trk1);
        t2.SetPtEtaPhiM(trk2_pt[i], trk2_eta[i], trk2_phi[i], mass_trk2);
        TLorentzVector pair = t1 + t2;

        // pT cut
        float ptPair = pair.Pt();
        if (ptPair < pt_pair_min) continue;

        // Mass window
        float mPair = pair.M();
        if (mPair <= mass_low || mPair >= mass_high) continue;

        // Isolation
        if (iso[i] < 0.9f) continue;

        out.nGood++;
        if (ptPair > bestPt) {
            bestPt = ptPair;
            out.bestIdx = i;
        }
    }

    return out;
}

// build meson 4-vector from 2 tracks
TLorentzVector build_meson_from_tracks(float pt1, float eta1, float phi1, float mass1,
                                      float pt2, float eta2, float phi2, float mass2)
{
    TLorentzVector t1, t2;
    t1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
    t2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    return (t1 + t2);
}

//--------------------------------------------------------------------------------------

// MC truth matching
struct MCMatching {
    bool isPhotonMatched;
    bool isMesonMatched;
    bool isBosonMatched;
};

MCMatching match_mc(const RVec<int>& pdgId, const RVec<int>& motherIdx, const RVec<int>& motherPdgId,
                    const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi,
                    float recoPhoton_eta, float recoPhoton_phi, float recoMeson_eta,  float recoMeson_phi){

    MCMatching out{false,false,false};

   // boolean vector for photon and meson candidates
    auto isPhoton = (pdgId == 22) && (Take(motherPdgId, motherIdx) == 25 || Take(motherPdgId, motherIdx) == 23); // take(v, indices) returns a vector with the elements of v corresponding to the position of indices
    auto isMeson  = ((abs(pdgId) == 333) || (abs(pdgId) == 113) || (abs(pdgId) == 313)) && (Take(motherPdgId, motherIdx) == 25 || Take(motherPdgId, motherIdx) == 23);
 
    
    // gen eta/phi vectors of the candidates
    auto genPhoton_eta = eta[isPhoton]; // vector with the elements fulfilling the condition isPhoton==true
    auto genPhoton_phi = phi[isPhoton];
    auto genMeson_eta  = eta[isMeson];
    auto genMeson_phi  = phi[isMeson];

    // if they exist, match in ΔR
    if (genPhoton_eta.size() > 0) {
        auto dphi = abs(recoPhoton_phi - genPhoton_phi); // vector, dphi[i] = abs(recoPhoton_phi - genPhoton_phi[i])  
        dphi = Map(dphi, [](float x){ return x > M_PI ? 2*M_PI - x : x; }); // 1st argument: input vector (dphi), 2nd arg: lambda function (return 2pi-dphi if dphi > pi, else return dphi)
        auto dR = sqrt((recoPhoton_eta - genPhoton_eta)*(recoPhoton_eta - genPhoton_eta) + dphi*dphi); //vector
        out.isPhotonMatched = Any(dR < 0.2); // return True if any element of the dR vector is < 0.2 -> if at least 1 gen photon is at dR<0.2 from the reco one, match done 
    }
    if (genMeson_eta.size() > 0) {
        auto dphi = abs(recoMeson_phi - genMeson_phi);
        dphi = Map(dphi, [](float x){ return x > M_PI ? 2*M_PI - x : x; });
        auto dR = sqrt((recoMeson_eta - genMeson_eta)*(recoMeson_eta - genMeson_eta) + dphi*dphi);
        out.isMesonMatched = Any(dR < 0.3);
    }

    out.isBosonMatched = out.isPhotonMatched && out.isMesonMatched;
    return out;
}

""")

# Create the dataframe
df = ROOT.RDataFrame(input_tree_name, input_file)
n_total = df.Count().GetValue()

# ------------------------------------------------------------
# Trigger
# ------------------------------------------------------------
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
df = df.Define("bestPhoton_pt",  "bestPhotonIdx >= 0 ? Photon_pt[bestPhotonIdx]  : -1.f")
df = df.Define("bestPhoton_eta", "bestPhotonIdx >= 0 ? Photon_eta[bestPhotonIdx] : -999.f")
df = df.Define("bestPhoton_phi", "bestPhotonIdx >= 0 ? Photon_phi[bestPhotonIdx] : -999.f")
df = df.Define("photon_p4", "make_photon_p4(Photon_pt, Photon_eta, Photon_phi, bestPhotonIdx)")

# ------------------------------------------------------------
# Mesons
# ------------------------------------------------------------
if isPhiAnalysis:
    mass_low, mass_high = 1.00, 1.05
    mass_trk1, mass_trk2 = 0.4937, 0.4937  # Kaon mass
    meson_prefix = "phi"
elif isRhoAnalysis:
    mass_low, mass_high = 0.50, 1.00
    mass_trk1, mass_trk2 = 0.13957, 0.13957  #Pion mass
    meson_prefix = "rho"
else:
    print("Unknown meson type!")
    mass_low, mass_high = 0.0, 99.0
    meson_prefix = "phi"

if not ismesonFromTracks:
    if verbose:
        print("### Mode: meson from NanoAOD (kinematics branches) ###")

    df = df.Define("meson_sel", 
                f"select_mesons_kin({meson_prefix}_kin_pt, {meson_prefix}_kin_eta, {meson_prefix}_kin_phi, "
                f"{meson_prefix}_kin_mass, {meson_prefix}_iso, {mass_low}, {mass_high}, "
                f"{meson_prefix}_trk1_pt, {meson_prefix}_trk1_eta, {meson_prefix}_trk1_phi, "
                f"{meson_prefix}_trk2_pt, {meson_prefix}_trk2_eta, {meson_prefix}_trk2_phi)")
    df = df.Define("nGoodMesons", "meson_sel.nGood")
    df = df.Define("bestMesonIdx", "meson_sel.bestIdx")


    df = df.Filter("nGoodMesons > 0", "At least one good meson")
    n_meson = df.Filter("nGoodMesons").Count().GetValue()

    # pT, eta, phi, p4 of the selected meson (if it exists)
    df = df.Define("meson_p4", f"make_meson_p4({meson_prefix}_kin_pt, {meson_prefix}_kin_eta, {meson_prefix}_kin_phi, {meson_prefix}_kin_mass, bestMesonIdx)")
    df = df.Define("bestMeson_pt",  f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_pt[bestMesonIdx] : -1.f")
    df = df.Define("bestMeson_eta", f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_eta[bestMesonIdx] : -999.f")
    df = df.Define("bestMeson_phi", f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_phi[bestMesonIdx] : -999.f")
    df = df.Define("bestMeson_mass",f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_mass[bestMesonIdx] : -1.f")

    # final tracks pT selection
    df = df.Define("trk1_pt_best", f"bestMesonIdx  >= 0 ? {meson_prefix}_trk1_pt[bestMesonIdx] : -1.f")
    df = df.Define("trk2_pt_best", f"bestMesonIdx  >= 0 ? {meson_prefix}_trk2_pt[bestMesonIdx] : -1.f")
    df = df.Define("firstTrk_pt", "trk1_pt_best > trk2_pt_best ? trk1_pt_best : trk2_pt_best")
    df = df.Define("secondTrk_pt", "trk1_pt_best > trk2_pt_best ? trk2_pt_best : trk1_pt_best")
    df = df.Filter("firstTrk_pt >= 20 && secondTrk_pt >= 5", "Final track pT selection")
    n_meson_trks = df.Filter("firstTrk_pt").Count().GetValue()

    # eta, phi of the best meson tracks
    df = df.Define("trk1_eta_best", f"bestMesonIdx  >= 0 ? {meson_prefix}_trk1_eta[bestMesonIdx] : -999.f")
    df = df.Define("trk1_phi_best", f"bestMesonIdx  >= 0 ? {meson_prefix}_trk1_phi[bestMesonIdx] : -999.f")
    df = df.Define("trk2_eta_best", f"bestMesonIdx  >= 0 ? {meson_prefix}_trk2_eta[bestMesonIdx] : -999.f")
    df = df.Define("trk2_phi_best", f"bestMesonIdx  >= 0 ? {meson_prefix}_trk2_phi[bestMesonIdx] : -999.f")
    df = df.Define("firstTrk_eta", "trk1_pt_best  > trk2_pt_best ? trk1_eta_best : trk2_eta_best")
    df = df.Define("firstTrk_phi", "trk1_pt_best  > trk2_pt_best ? trk1_phi_best : trk2_phi_best")
    df = df.Define("secondTrk_eta", "trk1_pt_best  > trk2_pt_best ? trk2_eta_best : trk1_eta_best")
    df = df.Define("secondTrk_phi", "trk1_pt_best  > trk2_pt_best ? trk2_phi_best : trk1_phi_best")

else:
    if verbose:
        print("### Mode: meson reconstructed from tracks (vectorized) ###")

    # meson selection 
    df = (df         
            .Define("meson_sel",  f"select_mesons_tracks({meson_prefix}_iso, "
                                                       f"{meson_prefix}_trk1_pt, {meson_prefix}_trk1_eta, {meson_prefix}_trk1_phi, "
                                                       f"{meson_prefix}_trk2_pt, {meson_prefix}_trk2_eta, {meson_prefix}_trk2_phi, "
                                                       f"{mass_trk1}f, {mass_trk2}f, {mass_low}f, {mass_high}f, 0.07f, 38.f)")
            .Define("nGoodMesons", "meson_sel.nGood")
            .Define("bestMesonIdx", "meson_sel.bestIdx")
            .Filter("nGoodMesons > 0", "At least one good meson (tracks-based)")
         )               
    n_meson = df.Filter("nGoodMesons").Count().GetValue()
    
    # tracks of the selected meson
    df = (df
            .Define("trk1_pt_best",  f"{meson_prefix}_trk1_pt[bestMesonIdx]")
            .Define("trk1_eta_best", f"{meson_prefix}_trk1_eta[bestMesonIdx]")
            .Define("trk1_phi_best", f"{meson_prefix}_trk1_phi[bestMesonIdx]")
            .Define("trk2_pt_best",  f"{meson_prefix}_trk2_pt[bestMesonIdx]")
            .Define("trk2_eta_best", f"{meson_prefix}_trk2_eta[bestMesonIdx]")
            .Define("trk2_phi_best", f"{meson_prefix}_trk2_phi[bestMesonIdx]")
            # meson p4 from tracks 
            .Define("meson_p4", f"build_meson_from_tracks(trk1_pt_best, trk1_eta_best, trk1_phi_best, {mass_trk1}f, "
                                                          f"trk2_pt_best, trk2_eta_best, trk2_phi_best, {mass_trk2}f)")
            .Define("bestMeson_pt",   "float(meson_p4.Pt())")
            .Define("bestMeson_eta",  "float(meson_p4.Eta())")
            .Define("bestMeson_phi",  "float(meson_p4.Phi())")
            .Define("bestMeson_mass", "float(meson_p4.M())")
            # tracks sorting
            .Define("firstTrk_pt",  "trk1_pt_best  > trk2_pt_best ? trk1_pt_best  : trk2_pt_best")
            .Define("firstTrk_eta", "trk1_pt_best  > trk2_pt_best ? trk1_eta_best : trk2_eta_best")
            .Define("firstTrk_phi", "trk1_pt_best  > trk2_pt_best ? trk1_phi_best : trk2_phi_best")
            .Define("secondTrk_pt",  "trk1_pt_best  > trk2_pt_best ? trk2_pt_best  : trk1_pt_best")
            .Define("secondTrk_eta", "trk1_pt_best  > trk2_pt_best ? trk2_eta_best : trk1_eta_best")
            .Define("secondTrk_phi", "trk1_pt_best  > trk2_pt_best ? trk2_phi_best : trk1_phi_best")
            # final tracks cut
            .Filter("firstTrk_pt >= 20 && secondTrk_pt >= 5", "Final track pT selection") 
        )       
    n_meson_trks = df.Filter("firstTrk_pt").Count().GetValue()


# ------------------------------------------------------------
# Boson
# ------------------------------------------------------------
df = df.Define("H_p4", "photon_p4 + meson_p4")
df = df.Define("H_mass", "H_p4.M()")
df = df.Define("H_pt", "H_p4.Pt()")
df = df.Define("H_eta", "H_p4.Eta()")
df = df.Define("H_phi", "H_p4.Phi()")


# ------------------------------------------------------------
# MC truth
# ------------------------------------------------------------
if not runningOnData:
    df = df.Define("mc_match", "match_mc(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pdgId, "
                                         "GenPart_pt, GenPart_eta, GenPart_phi, "
                                         "bestPhoton_eta, bestPhoton_phi, bestMeson_eta, bestMeson_phi)")

    df = df.Define("isPhotonMatched", "mc_match.isPhotonMatched")
    df = df.Define("isMesonMatched",  "mc_match.isMesonMatched")
    df = df.Define("isBosonMatched",  "mc_match.isBosonMatched")
    n_matched_boson = df.Filter("isBosonMatched").Count().GetValue()
    #print(f"n_matched_boson = {n_matched_boson}")

else:
    df = df.Define("isPhotonMatched", "false") \
           .Define("isMesonMatched",  "false") \
           .Define("isBosonMatched",  "false")

# ---------------------------------------------------------------------
# TTree writing
# ---------------------------------------------------------------------
columns_to_save = ["HLT_Photon35_TwoProngs35", "nMuons10", "nMuons20", "nElectrons10", "nElectrons20", 
                   "nGoodPhotons", "bestPhoton_pt", "bestPhoton_eta", "bestPhoton_phi",
                   "bestMeson_pt", "bestMeson_eta", "bestMeson_phi", "bestMeson_mass",
                   "firstTrk_pt", "firstTrk_eta", "firstTrk_phi", "secondTrk_pt", "secondTrk_eta", "secondTrk_phi",
                   "H_mass", "H_pt", "H_eta", "H_phi",
                   "isPhotonMatched", "isMesonMatched", "isBosonMatched"]

df.Snapshot("tree_output", output_file, columns_to_save)

# ---------------------------------------------------------------------
# efficiency histogram
# --------------------------------------------------------------------
cut_names = ["All events", "Triggered", "Best photon", "Best meson", "Trks cut"]

h_cutflow = ROOT.TH1F("nEvents", "Event counting in different steps; ; Events", len(cut_names), 0, len(cut_names))
for i, (label, val) in enumerate(zip(cut_names, [n_total, n_trigger, n_photon, n_meson, n_meson_trks]), start=1):
    h_cutflow.SetBinContent(i, val)
    h_cutflow.GetXaxis().SetBinLabel(i, label)

# histo writing
output = ROOT.TFile(output_file, "UPDATE")
output.cd()
h_cutflow.Write()
output.Close()

if verbose:
    print("Output saved in:", output_file)

# ---------------------------------------------------------------------
# DEBUG PRINTS
# ---------------------------------------------------------------------
if debug:
    # --- Photons ---
    all_photon_pts   = df.Take["ROOT::VecOps::RVec<float>"]("Photon_pt")
    best_photon_pts  = df.Take["float"]("bestPhoton_pt")
    n_good_photons   = df.Take["int"]("nGoodPhotons")

    # --- Mesons ---
    all_meson_pts    = df.Take["ROOT::VecOps::RVec<float>"](f"{meson_prefix}_kin_pt")
    best_meson_pts   = df.Take["float"]("bestMeson_pt")
    n_good_mesons    = df.Take["int"]("nGoodMesons")

    print("\n==================== EVENT DEBUG ====================")
    for i, (ph_all, n_ph, ph_best, mes_all, n_mes, mes_best) in enumerate(zip(all_photon_pts, n_good_photons, best_photon_pts, all_meson_pts, n_good_mesons, best_meson_pts)):
        print(f"\nEvent {i}:")
        print(f"  Photons -> all pTs = {list(ph_all)},  nGood = {n_ph},  best = {ph_best:.2f}")
        print(f"  Mesons  -> all pTs = {list(mes_all)},  nGood = {n_mes},  best = {mes_best:.2f}")
    print("=====================================================\n")       

# ------------------------------------------------------------
# Summary prints
# ------------------------------------------------------------
print("\nEVENT SUMMARY REPORT\n")
print(f"Total events:            {n_total}")
print(f"After trigger:           {n_trigger}")
print(f"After photon selection:  {n_photon}")
print(f"After meson selection:   {n_meson}")
print(f"After track-level cuts:  {n_meson_trks}")

if not runningOnData:
    print(f"Selection efficiency:    {100.0 * n_meson_trks / n_total:.1f}%\n")
    print(f"Bosons matched to MC truth:  {n_matched_boson}")
    print(f"Matching efficiency:         {100.0 * n_matched_boson / n_meson_trks:.1f}%")