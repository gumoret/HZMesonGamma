import ROOT
import argparse

# multithreading
ROOT.ROOT.EnableImplicitMT(8)

# following bools are given as input
verbose       = True
debug         = False
isPhiAnalysis = False # for H -> Phi Gamma
isRhoAnalysis = False # for H -> Rho Gamma
isKAnalysis   = False # for H -> K*0 Gamma
isDAnalysis   = False # for H -> D*0 Gamma
ismesonFromTracks = False # debug for reconstructing meson from tracks

# PARSER and INPUT 
p = argparse.ArgumentParser(description="RDataFrame analyzer for H→meson+γ")
p.add_argument("meson_option", help="Type <<rho>> for rho, <<phi>> for phi, <<K>> for K*")
p.add_argument("runningOnData_option", help="Type <<signal>> for signal, <<data>> for data")
p.add_argument("rootfile_name", help="Input nanoAOD ROOT file")
p.add_argument("outputfile_option", help="Output ROOT file")
args = p.parse_args()

if args.meson_option == "phi": isPhiAnalysis = True
elif args.meson_option == "rho": isRhoAnalysis = True
elif args.meson_option == "K": isKAnalysis = True
elif args.meson_option == "D": isDAnalysis = True
else: print("meson_option must be <<phi>> or <<rho>> or <<K*>> or <<DO*>>")


if args.runningOnData_option == "signal": runningOnData = False
elif args.runningOnData_option == "data": runningOnData = True
else: print("runninOnData must be <<signal>> or <<data>>")

read_list = True
input_file = args.rootfile_name

import glob

#Gamma+Jets MC background
#base_dir = "/ceph/submit/data/group/cms/store/user/mariadlf/D07/GJ-4Jets_Bin-HT-1000-PTG-100to200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8+RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2+MINIAODSIM/"
#files = sorted(glob.glob(f"{base_dir}/*.root")) 


#Tau2022
'''
base_dir = "/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D05" 

files = []

for era in ["12022", "22022"]:
    files.extend(glob.glob(f"{base_dir}/{era}/Tau+Run*/*.root"))
'''
#EGamma 2024
'''
base_dir = "/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D07/2024"
files = []

files.extend(glob.glob(f"{base_dir}/EGamma*/*.root"))

files = sorted(files)
'''

#HRhogamma signal
#base_dir = "/ceph/submit/data/group/cms/store/user/mariadlf/D07/GluGluHtoRhoG_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8-evtgen+RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2+MINIAODSIM"

#HPhiGamma signal
base_dir = "/ceph/submit/data/group/cms/store/user/mariadlf/D07/GluGluHtoPhiG_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8-evtgen+RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2+MINIAODSIM"

#HKstGamma signal
#base_dir = "/ceph/submit/data/group/cms/store/user/mariadlf/D07/GluGluHtoKStar0G_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8-evtgen+RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2+MINIAODSIM"

#ZPhiGamma signal
#base_dir = "/ceph/submit/data/group/cms/store/user/mariadlf/D07/ZtoPhiG_TuneCP5_13p6TeV_madgraphMLM-pythia8+RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v3+MINIAODSIM"

#ZRhoGamma signal
#base_dir = "/ceph/submit/data/group/cms/store/user/mariadlf/D07/ZtoRhoG_TuneCP5_13p6TeV_madgraphMLM-pythia8+RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v3+MINIAODSIM"


files = sorted(glob.glob(f"{base_dir}/*.root"))




if read_list: input_file = files
output_file = args.outputfile_option
input_tree_name = "Events"

print("running...")

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
    float deltaMesonMass;
};

MCMatching match_mc(const RVec<int>& pdgId, const RVec<int>& motherIdx, const RVec<int>& motherPdgId,
                    const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi, const RVec<float>& mass,
                    float recoPhoton_eta, float recoPhoton_phi, float recoMeson_eta,  float recoMeson_phi, float recoMeson_mass){

    MCMatching out{false,false,false,-1.};

   // boolean vector for photon and meson candidates
    auto isPhoton = (pdgId == 22) && (Take(motherPdgId, motherIdx) == 25 || Take(motherPdgId, motherIdx) == 23); // take(v, indices) returns a vector with the elements of v corresponding to the position of indices
    auto isMeson  = ((abs(pdgId) == 333) || (abs(pdgId) == 113) || (abs(pdgId) == 313) || (abs(pdgId) == 423)) && (Take(motherPdgId, motherIdx) == 25 || Take(motherPdgId, motherIdx) == 23);
 
    
    // gen eta/phi vectors of the candidates
    auto genPhoton_eta = eta[isPhoton]; // vector with the elements fulfilling the condition isPhoton==true
    auto genPhoton_phi = phi[isPhoton];
    auto genMeson_eta  = eta[isMeson];  // vector with the elements fulfilling the condition isMeson==true
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
        
        auto idx_min = ArgMin(dR); //take the index of the minimum dR, to select the best matched meson
        if (dR[idx_min] < 0.3) {
            out.isMesonMatched = true;
            auto genMeson_mass = mass[isMeson];
            float matchedGenMass = genMeson_mass[idx_min];
            out.deltaMesonMass = recoMeson_mass - matchedGenMass;
        }
    }

    out.isBosonMatched = out.isPhotonMatched && out.isMesonMatched;
    return out;
}

""")

# Create the dataframe
df = ROOT.RDataFrame(input_tree_name, input_file)
n_total = df.Count().GetValue()

# ------------------------------------------------------------
# Pileup
# ------------------------------------------------------------
if not runningOnData: df = df.Define("nPU", "Pileup_nPU")
else: df = df.Define("nPU", "0")
h_pu = df.Histo1D(("pileup", "Pileup distribution", 130, 0, 130), "nPU")

# ------------------------------------------------------------
# MC weight
# ------------------------------------------------------------
if not runningOnData: df = df.Define("MC_Weight", "Generator_weight")
else: df = df.Define("MC_Weight", "1.0")


# ------------------------------------------------------------
# Trigger
# ------------------------------------------------------------
df = df.Filter("HLT_Photon35_TwoProngs35", "Two-prong photon trigger")
n_trigger = df.Filter("HLT_Photon35_TwoProngs35").Count().GetValue()

#df = df.Filter("HLT_Photon50EB_TightID_TightIso", "Single photon trigger")
#n_trigger = df.Filter("HLT_Photon50EB_TightID_TightIso").Count().GetValue() ##modify colums_to_save also

#trigger_expr = "HLT_Photon35_TwoProngs35 || HLT_Photon50EB_TightID_TightIso"

#df = df.Define("passTrigger", trigger_expr)
#df = df.Filter("passTrigger", "OR trigger") ##modify colums_to_save also
#n_trigger = df.Filter("passTrigger").Count().GetValue()

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
    #dictionary for tracks names
    trk1_pt   = f"{meson_prefix}_trk1_pt"
    trk1_eta  = f"{meson_prefix}_trk1_eta"
    trk1_phi  = f"{meson_prefix}_trk1_phi"
    trk2_pt   = f"{meson_prefix}_trk2_pt"
    trk2_eta  = f"{meson_prefix}_trk2_eta"
    trk2_phi  = f"{meson_prefix}_trk2_phi"
elif isRhoAnalysis:
    mass_low, mass_high = 0.50, 1.00
    mass_trk1, mass_trk2 = 0.13957, 0.13957  #Pion mass
    meson_prefix = "rho"
    #dictionary for tracks names
    trk1_pt   = f"{meson_prefix}_trk1_pt"
    trk1_eta  = f"{meson_prefix}_trk1_eta"
    trk1_phi  = f"{meson_prefix}_trk1_phi"
    trk2_pt   = f"{meson_prefix}_trk2_pt"
    trk2_eta  = f"{meson_prefix}_trk2_eta"
    trk2_phi  = f"{meson_prefix}_trk2_phi"
elif isKAnalysis:
    mass_low, mass_high = 0.60, 1.00
    mass_trk1, mass_trk2 = 0.4937, 0.13957  #Kaon mass, Pion mass
    meson_prefix = "K0Star" 
    #dictionary for tracks names
    trk1_pt  = f"{meson_prefix}_kaon_pt"
    trk1_eta = f"{meson_prefix}_kaon_eta"
    trk1_phi = f"{meson_prefix}_kaon_phi"
    trk2_pt  = f"{meson_prefix}_pion_pt"
    trk2_eta = f"{meson_prefix}_pion_eta"
    trk2_phi = f"{meson_prefix}_pion_phi"  
elif isDAnalysis:
    mass_low, mass_high = 2.00, 2.05
    mass_trk1, mass_trk2 = 0.4937, 0.13957  # Kaon, Pion
    meson_prefix = "d0pi0"
    #dictionary for tracks names
    trk1_pt  = f"{meson_prefix}_kaon_pt"
    trk1_eta = f"{meson_prefix}_kaon_eta"
    trk1_phi = f"{meson_prefix}_kaon_phi"
    trk2_pt  = f"{meson_prefix}_pion_pt"
    trk2_eta = f"{meson_prefix}_pion_eta"
    trk2_phi = f"{meson_prefix}_pion_phi" 
else:
    print("Unknown meson type!")
    mass_low, mass_high = 0.0, 99.0
    meson_prefix = "phi"

if not ismesonFromTracks:
    if verbose:
        print("### Mode: meson from NanoAOD ###")

    df = df.Define("meson_sel", 
                f"select_mesons_kin({meson_prefix}_kin_pt, {meson_prefix}_kin_eta, {meson_prefix}_kin_phi, "
                f"{meson_prefix}_kin_mass, {meson_prefix}_iso, {mass_low}, {mass_high}, "
                f"{trk1_pt}, {trk1_eta}, {trk1_phi}, "
                f"{trk2_pt}, {trk2_eta}, {trk2_phi})")
    df = df.Define("nGoodMesons", "meson_sel.nGood")
    df = df.Define("bestMesonIdx", "meson_sel.bestIdx")


    df = df.Filter("nGoodMesons > 0", "At least one good meson")
    n_meson = df.Filter("nGoodMesons").Count().GetValue()

    # pT, eta, phi, p4, iso of the selected meson (if it exists)
    df = df.Define("meson_p4",       f"make_meson_p4({meson_prefix}_kin_pt, {meson_prefix}_kin_eta, {meson_prefix}_kin_phi, {meson_prefix}_kin_mass, bestMesonIdx)")
    df = df.Define("bestMeson_pt",   f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_pt[bestMesonIdx] : -1.f")
    df = df.Define("bestMeson_eta",  f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_eta[bestMesonIdx] : -999.f")
    df = df.Define("bestMeson_phi",  f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_phi[bestMesonIdx] : -999.f")
    df = df.Define("bestMeson_mass", f"bestMesonIdx  >= 0 ? {meson_prefix}_kin_mass[bestMesonIdx] : -1.f")
    df = df.Define("isoMeson",       f"bestMesonIdx  >= 0 ? {meson_prefix}_iso[bestMesonIdx] : -1.f")

    # final tracks pT selection
    df = df.Define("trk1_pt_best", f"bestMesonIdx  >= 0 ? {trk1_pt}[bestMesonIdx] : -1.f")
    df = df.Define("trk2_pt_best", f"bestMesonIdx  >= 0 ? {trk2_pt}[bestMesonIdx] : -1.f")
    df = df.Define("firstTrk_pt",  "trk1_pt_best > trk2_pt_best ? trk1_pt_best : trk2_pt_best")
    df = df.Define("secondTrk_pt", "trk1_pt_best > trk2_pt_best ? trk2_pt_best : trk1_pt_best")
    df = df.Filter("firstTrk_pt >= 20 && secondTrk_pt >= 5", "Final track pT selection")
    n_meson_trks = df.Filter("firstTrk_pt").Count().GetValue()

    # eta, phi of the best meson tracks
    df = df.Define("trk1_eta_best", f"bestMesonIdx  >= 0 ? {trk1_eta}[bestMesonIdx] : -999.f")
    df = df.Define("trk1_phi_best", f"bestMesonIdx  >= 0 ? {trk1_phi}[bestMesonIdx] : -999.f")
    df = df.Define("trk2_eta_best", f"bestMesonIdx  >= 0 ? {trk2_eta}[bestMesonIdx] : -999.f")
    df = df.Define("trk2_phi_best", f"bestMesonIdx  >= 0 ? {trk2_phi}[bestMesonIdx] : -999.f")
    df = df.Define("firstTrk_eta",  "trk1_pt_best  > trk2_pt_best ? trk1_eta_best : trk2_eta_best")
    df = df.Define("firstTrk_phi",  "trk1_pt_best  > trk2_pt_best ? trk1_phi_best : trk2_phi_best")
    df = df.Define("secondTrk_eta", "trk1_pt_best  > trk2_pt_best ? trk2_eta_best : trk1_eta_best")
    df = df.Define("secondTrk_phi", "trk1_pt_best  > trk2_pt_best ? trk2_phi_best : trk1_phi_best")

else:
    if verbose:
        print("### Mode: meson reconstructed from tracks ###")

    # meson selection 
    df = (df         
            .Define("meson_sel",  f"select_mesons_tracks({meson_prefix}_iso, "
                                                       f"{trk1_pt}, {trk1_eta}, {trk1_phi}, "
                                                       f"{trk2_pt}, {trk2_eta}, {trk2_phi}, "
                                                       f"{mass_trk1}f, {mass_trk2}f, {mass_low}f, {mass_high}f, 0.07f, 38.f)")
            .Define("nGoodMesons", "meson_sel.nGood")
            .Define("bestMesonIdx", "meson_sel.bestIdx")
            .Filter("nGoodMesons > 0", "At least one good meson (tracks-based)")
         )               
    n_meson = df.Filter("nGoodMesons").Count().GetValue()
    
    # tracks of the selected meson
    df = (df
            .Define("trk1_pt_best",  f"{trk1_pt}[bestMesonIdx]")
            .Define("trk1_eta_best", f"{trk1_eta}[bestMesonIdx]")
            .Define("trk1_phi_best", f"{trk1_phi}[bestMesonIdx]")
            .Define("trk2_pt_best",  f"{trk2_pt}[bestMesonIdx]")
            .Define("trk2_eta_best", f"{trk2_eta}[bestMesonIdx]")
            .Define("trk2_phi_best", f"{trk2_phi}[bestMesonIdx]")
            
            # meson p4 from tracks 
            .Define("meson_p4", f"build_meson_from_tracks(trk1_pt_best, trk1_eta_best, trk1_phi_best, {mass_trk1}f, "
                                                          f"trk2_pt_best, trk2_eta_best, trk2_phi_best, {mass_trk2}f)")
            .Define("bestMeson_pt",   "float(meson_p4.Pt())")
            .Define("bestMeson_eta",  "float(meson_p4.Eta())")
            .Define("bestMeson_phi",  "float(meson_p4.Phi())")
            .Define("bestMeson_mass", "float(meson_p4.M())")

            # meson isolation
            .Define("isoMeson", f"{meson_prefix}_iso[bestMesonIdx]")

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
                                         "GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, "
                                         "bestPhoton_eta, bestPhoton_phi, bestMeson_eta, bestMeson_phi, bestMeson_mass)")

    df = df.Define("isPhotonMatched", "mc_match.isPhotonMatched")
    df = df.Define("isMesonMatched",  "mc_match.isMesonMatched")
    df = df.Define("isBosonMatched",  "mc_match.isBosonMatched")
    df = df.Define("delta_meson_mass","mc_match.deltaMesonMass")
    n_matched_boson = df.Filter("isBosonMatched").Count().GetValue()
    #print(f"n_matched_boson = {n_matched_boson}")

    #for m_meson_reco - m_meson_gen
    #df = df.Define("delta_meson_mass", f"bestMeson_mass-{meson_prefix}_gen_mass")

else:
    df = df.Define("isPhotonMatched", "false") \
           .Define("isMesonMatched",  "false") \
           .Define("isBosonMatched",  "false") \
           .Define("delta_meson_mass", "0.")


# ---------------------------------------------------------------------
# TTree writing
# ---------------------------------------------------------------------
columns_to_save = ["nPU", "MC_Weight", "HLT_Photon35_TwoProngs35", 
                   "nMuons10", "nMuons20", "nElectrons10", "nElectrons20", 
                   "nGoodPhotons", "bestPhoton_pt", "bestPhoton_eta", "bestPhoton_phi",
                   "bestMeson_pt", "bestMeson_eta", "bestMeson_phi", "bestMeson_mass", "isoMeson",
                   "firstTrk_pt", "firstTrk_eta", "firstTrk_phi", "secondTrk_pt", "secondTrk_eta", "secondTrk_phi",
                   "H_mass", "H_pt", "H_eta", "H_phi",
                   "isPhotonMatched", "isMesonMatched", "isBosonMatched", "delta_meson_mass"]

df.Snapshot("tree_output", output_file, columns_to_save)

# pileup histogram construction -----------------
h_pu.GetValue()

# ---------------------------------------------------------------------
# efficiency histogram
# --------------------------------------------------------------------
cut_names = ["All events", "Triggered", "Best photon", "Best meson", "Trks cut"]

h_cutflow = ROOT.TH1F("nEvents", "Event counting in different steps; ; Events", len(cut_names), 0, len(cut_names))

for i, (label, val) in enumerate(zip(cut_names, [n_total, n_trigger, n_photon, n_meson, n_meson_trks]), start=1):
    h_cutflow.SetBinContent(i, val)
    h_cutflow.GetXaxis().SetBinLabel(i, label)

# ---------------------------------------------------------------------
# histos writing
# ---------------------------------------------------------------------
output = ROOT.TFile(output_file, "UPDATE")
output.cd()
h_pu.Write()
h_cutflow.Write()
output.Close()

if verbose:
    print("Output saved in:", output_file)

# ---------------------------------------------------------------------
# DEBUG PRINTS
# ---------------------------------------------------------------------
if debug:
    # --- Photons ---
    all_photon_pts   = df.Take["ROOT::VecOps::RVec<float>"]("Photon_pt").GetValue()
    best_photon_pt   = df.Take["float"]("bestPhoton_pt").GetValue()
    n_good_photons   = df.Take["int"]("nGoodPhotons").GetValue()

    # --- Mesons ---
    all_meson_pts    = df.Take["ROOT::VecOps::RVec<float>"](f"{meson_prefix}_kin_pt").GetValue()
    best_meson_pt    = df.Take["float"]("bestMeson_pt").GetValue()
    n_good_mesons    = df.Take["int"]("nGoodMesons").GetValue()
    best_meson_mass  = df.Take["float"]("bestMeson_mass").GetValue()

    # --- Boson ---
    boson_pt   = df.Take["double"]("H_pt").GetValue()
    boson_mass = df.Take["double"]("H_mass").GetValue()

    # --- MC matching ---
    isPhotonMatched = df.Take["bool"]("isPhotonMatched").GetValue()
    isMesonMatched  = df.Take["bool"]("isMesonMatched").GetValue()
    is_boson_matched  = df.Take["bool"]("isBosonMatched").GetValue()

    '''
    print("\n==================== EVENT DEBUG ====================")
    for i, (ph_all, n_ph, ph_best, mes_all, n_mes, mes_best) in enumerate(zip(all_photon_pts, n_good_photons, best_photon_pt, all_meson_pts, n_good_mesons, best_meson_pt)):
        print(f"\nEvent {i}:")
        print(f"  Photons -> all pTs = {list(ph_all)},  nGood = {n_ph},  best = {ph_best:.2f}")
        print(f"  Mesons  -> all pTs = {list(mes_all)},  nGood = {n_mes},  best = {mes_best:.2f}")
    print("=====================================================\n")       
    '''

    print("\n==================== FULL EVENT DEBUG ====================\n")

    # counters
    nEventsBosonMatched = 0
    nEventsBosonNotMatched = 0

    for i in range(len(best_photon_pt)):

        print(f"Event {i}")
        #print("-----------------------------------------")

        # RECO INFO
        print(f"Reco Photon:")
        print(f"  pT  = {best_photon_pt[i]:.2f}")

        print(f"\nReco Meson (NanoAOD kin):")
        print(f"  pT  = {best_meson_pt[i]:.2f}")
        print(f"  mass     = {best_meson_mass[i]:.2f}")

        print(f"\nReco Boson:")
        print(f"  pT   = {boson_pt[i]:.2f}")
        print(f"  mass = {boson_mass[i]:.3f}")

        # MC TRUTH MATCHING
        #print("\nMC Truth Matching:")
        #print(f"  Photon matched = {isPhotonMatched[i]}")
        #print(f"  Meson  matched = {isMesonMatched[i]}")

        if not runningOnData:
            if is_boson_matched[i]:
                nEventsBosonMatched += 1
                print("\n**BOSON FOUND **")
            else:
                nEventsBosonNotMatched += 1
                print("THAT'S NOT A HIGGS or a Z!")

            print(f"MC H or Z FOUND = {nEventsBosonMatched}, MC H or Z NOT FOUND = {nEventsBosonNotMatched}")
        print("---------------------------------------------------------------")

print("===============================================================\n")

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