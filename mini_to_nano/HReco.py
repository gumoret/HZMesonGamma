import ROOT
import os
import sys
import json
from utilsHrare import getMClist, getDATAlist, getSkims, computeWeigths
from utilsAna import getMesonFromJson, pickTRG, getMVAFromJson, getTriggerFromJson
from utilsAna import loadCorrectionSet, loadPolarizationTree, loadSFhisto, loadUserCode, loadtmva_helper
import tmva_helper_xml

#ROOT.gSystem.Load("libGenVector")##
ROOT.gROOT.SetBatch()
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()

doSyst = True
doPol = True
doMVA = True
if sys.argv[1]=='isZtag': doMVA = False
if sys.argv[1]=='isWtag': doMVA = False

if sys.argv[2]=='isD0Pi0StarCat': doMVA = False
if sys.argv[2]=='isOmega3PiCat': doPol = False
if sys.argv[2]=='isPhi3PiCat': doPol = False
if sys.argv[2]=='isD0StarCat': doPol = False
if sys.argv[2]=='isD0Pi0StarCat': doPol = False


useD03=False
if sys.argv[2]=='isOmega3PiCat' or sys.argv[2]=='isPhi3PiCat' or sys.argv[2]=='isD0StarCat' or sys.argv[2]=='isD0Pi0StarCat': useD03=True



isGF = False
isZinv = False
isZ = False
isW = False
isVBF = False
isVBFlow = False
isBPH = False######

isPhiCat = "false"
isRhoCat = "false"
isK0StarCat = "false"
isOmega3PiCat = "false"
isPhi3PiCat = "false"
isD0StarCat = "false"
isD0Pi0StarCat = "false"



if sys.argv[1]=='isGFtag': isGF = True
if sys.argv[1]=='isZinvtag': isZinv = True
if sys.argv[1]=='isZtag': isZ = True
if sys.argv[1]=='isWtag': isW = True
if sys.argv[1]=='isVBFtag': isVBF = True
if sys.argv[1]=='isVBFtaglow': isVBFlow = True

if sys.argv[2]=='isPhiCat': isPhiCat = "true"


if sys.argv[4]=='2018': year = 2018
if sys.argv[4]=='2017': year = 2017
if sys.argv[4]=='22016': year = 22016 #F-H
if sys.argv[4]=='12016': year = 12016 #B-F

lumis={
    '12016': 19.52, #APV #(B-F for 2016 pre)
    '22016': 16.80, #postVFP
    '2016': 35.9,
    '2017': 41.5,
    '12017': 7.7, #(F for 2017)
    '2018': 59.70,
    '12018': 39.54,
    'all': 86.92,      #19.52 + 7.7 + 59.70
}





PRESELECTION = "(nPhoton>0 && (nphi>0 or nrho>0 or nK0Star>0 or nomega>0) && PV_npvsGood>0)"

with open("../analysis/config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

with open("../analysis/config/trigger.json") as trgJsonFile:
    trgObject = json.load(trgJsonFile)
    trgJsonFile.close()


GOODJETS = jsonObject['GOODJETS']
LOOSEmuons = jsonObject['LOOSEmuons']
LOOSEelectrons = jsonObject['LOOSEelectrons']
GOODMUON = jsonObject['GOODMUON']
GOODELE = jsonObject['GOODELE']
JSON = jsonObject['JSON']
BARRELphotons = jsonObject['BARRELphotons']
ENDCAPphotons = jsonObject['ENDCAPphotons']
ENDCAPphotonsLoose = jsonObject['ENDCAPphotonsLoose']
photonsLoose = jsonObject['photonsLOOSE']

METFLAG = jsonObject['METFLAG']

MVA = jsonObject['MVAweights']
TRIGGERS = trgObject['triggers']
mesons = jsonObject['mesons']


def selectionTAG(df, doSyst, isData):

    if isZ:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_mediumId and Muon_pfRelIso04_all < 0.25") # iso same as the loose
                 .Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("goodElectrons","{}".format(GOODELE)+" and Electron_mvaFall17V2Iso_WP90") # medium
                 .Define("looseMu","{}".format(LOOSEmuons))
                 .Define("looseEle","{}".format(LOOSEelectrons))
                 .Filter("Sum(looseEle)+Sum(looseMu)==2", "at least two muons or electrons, and no extra loose leptons")
                 .Define("isMuorEle","(Sum(looseMu)==2 and Sum(goodMuons)>=1)?1:(Sum(looseEle)==2 and Sum(goodElectrons)>=1)?2 :0")
                 .Filter("isMuorEle>0","at least one leading lepton")
                 .Define("Z_veto1", "Sum(looseEle)==2 ? Minv2(Electron_pt[looseEle][0], Electron_eta[looseEle][0], Electron_phi[looseEle][0], Electron_mass[looseEle][0], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]).first: -1")
                 .Define("Z_veto2", "Sum(looseEle)==2 ? Minv2(Electron_pt[looseEle][1], Electron_eta[looseEle][1], Electron_phi[looseEle][1], Electron_mass[looseEle][1], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]).first: -1")
                 .Filter("abs(Z_veto1-91) > 5 and abs(Z_veto2-91) > 5","kill the Z recontructed as gamma + electron")
                 .Define("V_mass", "(Sum(looseMu)==2 and Sum(Muon_charge[looseMu])==0)? Minv(Muon_pt[looseMu], Muon_eta[looseMu], Muon_phi[looseMu], Muon_mass[looseMu]) : (Sum(looseEle)==2 and Sum(Electron_charge[looseEle])==0) ? Minv(Electron_pt[looseEle], Electron_eta[looseEle], Electron_phi[looseEle], Electron_mass[looseEle]): 0.")
                 .Filter("(V_mass>(91-10) and V_mass<(91+15))","At least one good Z")
                 .Define("Visr_mass", "(Sum(looseMu)==2 and Sum(Muon_charge[looseMu])==0)? Minv3(Muon_pt[looseMu], Muon_eta[looseMu], Muon_phi[looseMu], Muon_mass[looseMu], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]) : (Sum(looseEle)==2 and Sum(Electron_charge[looseEle])==0) ? Minv3(Electron_pt[looseEle], Electron_eta[looseEle], Electron_phi[looseEle], Electron_mass[looseEle], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]): 0.")
                 .Filter("Visr_mass>91+5")
#                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
#                 .Define("Mu2_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][1], Muon_phi[goodMuons][1], TrigObj_eta, TrigObj_phi)")
#                 .Define("Ele1_hasTriggerMatch", "Sum(goodElectrons)>1 ? hasTriggerMatch(Electron_eta[goodElectrons][0], Electron_phi[goodElectrons][0], TrigObj_eta, TrigObj_phi) : 0")
#                 .Define("Ele2_hasTriggerMatch", "Sum(goodElectrons)>1 ? hasTriggerMatch(Electron_eta[goodElectrons][1], Electron_phi[goodElectrons][1], TrigObj_eta, TrigObj_phi) : 0")
#                 .Filter("trigger>0 and ((Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26) or (Mu2_hasTriggerMatch and Muon_pt[goodMuons][1]>26))","pass trigger")
                 .Define("LeadingLeptonPt", "isMuorEle==1 ? Muon_pt[looseMu][0]: isMuorEle==2 ? Electron_pt[looseEle][0] :0")
                 .Define("SubLeadingLeptonPt", "isMuorEle==1 ? Muon_pt[looseMu][1]: isMuorEle==2 ? Electron_pt[looseEle][1] :0")
                 .Define("LeadingLeptonEta", "isMuorEle==1 ? Muon_eta[looseMu][0]: isMuorEle==2 ? Electron_eta[looseEle][0] :0")
                 .Define("SubLeadingLeptonEta", "isMuorEle==1 ? Muon_eta[looseMu][1]: isMuorEle==2 ? Electron_eta[looseEle][1] :0")
                 )
        return dftag

    if isW:

        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_tightId and Muon_pfRelIso04_all < 0.15") ## tight
                 .Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("goodElectrons","{}".format(GOODELE)+" and Electron_mvaFall17V2Iso_WP80 and Electron_pt>30") ## tight
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(goodMuons)+Sum(goodElectrons))==1 and (Sum(vetoEle)+Sum(vetoMu))==1","one lepton")
                 .Define("isMuorEle","Sum(goodMuons)==1?1: Sum(goodElectrons)==1?2 :0")
                 .Define("V_mass","Sum(goodMuons)>0 ? mt(Muon_pt[goodMuons][0], Muon_phi[goodMuons][0], DeepMETResolutionTune_pt, DeepMETResolutionTune_phi) : mt(Electron_pt[goodElectrons][0], Electron_phi[goodElectrons][0], DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")
                 .Define("Z_veto", "Sum(goodElectrons)==1 ? Minv2(Electron_pt[goodElectrons][0], Electron_eta[goodElectrons][0], Electron_phi[goodElectrons][0], Electron_mass[goodElectrons][0], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]]).first: -1")
                 .Filter("abs(Z_veto-91) > 10","kill the Z recontructed as gamma + electron")
#                 .Filter("DeepMETResolutionTune_pt>15","MET>15")
#                 .Filter("V_mass>15","MT>15")
#                 .Define("trigger","{}".format(TRIGGER))
#                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
#                 .Filter("trigger>0 and Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26","pass trigger")
                 .Define("deltaLepMeson","{}".format(CLEAN_LepMes))
                 .Filter("deltaLepMeson>0.5","kill the muon reconstructed as meson")
                 .Define("dPhiGammaMET","abs(deltaPhi(goodPhotons_phi[index_pair[1]], DeepMETResolutionTune_phi))")
                 .Define("dPhiMesonMET","abs(deltaPhi(goodMeson_phi[index_pair[0]], DeepMETResolutionTune_phi))")
                 .Define("dPhiLeptonMET","isMuorEle==1 ? abs(deltaPhi(Muon_phi[goodMuons][0], DeepMETResolutionTune_phi)): isMuorEle==2 ? abs(deltaPhi(Electron_phi[goodElectrons][0], DeepMETResolutionTune_phi)): 999.")
                 .Define("LeadingLeptonPt", "isMuorEle==1 ? Muon_pt[goodMuons][0]: isMuorEle==2 ? Electron_pt[goodElectrons][0] :0")
                 .Define("LeadingLeptonEta", "isMuorEle==1 ? Muon_eta[goodMuons][0]: isMuorEle==2 ? Electron_eta[goodElectrons][0] :0")
                 )
        return dftag

    if isVBF or isVBFlow:

        VBFcut = "mJJ>300 and dEtaJJ>3 and Y1Y2<0 and Jet_pt[goodJets][0]>40"
        if isVBFlow: VBFcut = "mJJ>250 and dEtaJJ>3. and Y1Y2<0"

# tight means less PU
## default is the medium (Jet_puId == 2)
## use the tight PU id (Jet_puId == 1) for 2017-2018: note those are swapped for the 2016
## use the tight PU id (Jet_puId == 4) for 2016-12016
## https://github.com/cms-nanoAOD/cmssw/issues/583
## flag means passlooseID*4+passmediumID*2+passtightID*1.
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID#miniAOD_and_nanoAOD

        PUjetID = "true"
#        if year == 2018 or year == 2017: PUjetID = "(((Jet_puId & 1) and abs(Jet_eta)>2.75) or ((Jet_puId & 2) and abs(Jet_eta)<=2.75))"
#        if year == 2016 or year == 12016: PUjetID = "(((Jet_puId & 4) and abs(Jet_eta)>2.75 or ((Jet_puId & 2) and abs(Jet_eta)<=2.75))"

        if (doSyst and isData == "false"):
            df = df.Define("Jet_delta",'computeJECuncertainties(corr_sf, Jet_pt, Jet_eta)').Vary("Jet_pt", "ROOT::RVec<ROOT::RVecF>{Jet_pt*(1-Jet_delta),Jet_pt*(1+Jet_delta)}", variationTags=["dn","up"], variationName = "JetSYST")

        dftag = (df.Define("goodJets","{}".format(GOODJETS)+" and {}".format(PUjetID))
                 .Define("nGoodJets","Sum(goodJets)*1.0f").Filter("Sum(goodJets)>1","two jets for VBF")
                 .Define("mJJ","Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets])")
                 .Define("dEtaJJ","abs(Jet_eta[goodJets][0] - Jet_eta[goodJets][1])")
                 .Define("dPhiJJ","abs(deltaPhi(Jet_phi[goodJets][0],Jet_phi[goodJets][1]))")
                 .Define("Y1Y2","Jet_eta[goodJets][0]*Jet_eta[goodJets][1]")
                 .Filter("{}".format(VBFcut),"Filter on MJJ , Deta, Y1Y2")
                 .Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(vetoEle)+Sum(vetoMu))==0", "no leptons")
#                 .Define("trigger","{}".format(TRIGGER))
#                 .Filter("trigger>0", "pass triggers")
                 .Define("SoftActivityJetNjets5F","SoftActivityJetNjets5*1.0f")
                 .Define("jet1Pt","Jet_pt[goodJets][0]")
                 .Define("jet2Pt","Jet_pt[goodJets][1]")
                 .Define("jet1Eta","Jet_eta[goodJets][0]")
                 .Define("jet2Eta","Jet_eta[goodJets][1]")
                 .Define("jet1hfsigmaPhiPhi","Jet_hfsigmaPhiPhi[goodJets][0]")
                 .Define("jet1hfsigmaEtaEta","Jet_hfsigmaEtaEta[goodJets][0]")
                 .Define("jet2hfsigmaPhiPhi","Jet_hfsigmaPhiPhi[goodJets][1]")
                 .Define("jet2hfsigmaEtaEta","Jet_hfsigmaEtaEta[goodJets][1]")
                 .Define("deltaJetMeson","{}".format(CLEAN_JetMes))
                 .Define("deltaJetPhoton","{}".format(CLEAN_JetPH))
                 .Define("zepVar","compute_jet_HiggsVars_var(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]], goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]], goodMeson_mass[index_pair[0]], 0)")
                 .Define("detaHigJet1","compute_jet_HiggsVars_var(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]], goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]], goodMeson_mass[index_pair[0]], 1)")
                 .Define("detaHigJet2","compute_jet_HiggsVars_var(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]], goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]], goodMeson_mass[index_pair[0]], 2)")
#                 .Filter("DeepMETResolutionTune_pt<75","DeepMETResolutionTune_pt<75") # not doing Zinv as nominal
                 )
        return dftag

    if isZinv:
        dftag = (df.Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(vetoEle)+Sum(vetoMu))==0", "no leptons")
#                 .Define("trigger","{}".format(TRIGGER))
#                 .Filter("trigger>0", "pass triggers")
                 .Filter("DeepMETResolutionTune_pt>50","MET>50")
                 .Define("metFilter","{}".format(METFLAG))
                 .Filter("metFilter", "pass METfilter")
                 .Define("dPhiGammaMET","abs(deltaPhi(goodPhotons_phi[index_pair[1]], DeepMETResolutionTune_phi))")
                 .Define("dPhiMesonMET","abs(deltaPhi(goodMeson_phi[index_pair[0]], DeepMETResolutionTune_phi))")
                 .Define("ptRatioMEThiggs","abs(DeepMETResolutionTune_pt-HCandPT)/HCandPT")
#                 .Filter("ptRatioMEThiggs<0.8","ptRatioMEThiggs<0.8")
                 .Filter("dPhiGammaMET>1","dPhiGammaMET>1")
                 .Filter("dPhiMesonMET>1","dPhiMesonMET>1")
                 ##
                 .Define("goodJets","{}".format(GOODJETS))
                 .Define("nGoodJets","Sum(goodJets)*1.0f")
                 .Define("bjet", "Jet_btagDeepB[goodJets] > {}".format(DEEP_B_MEDIUM['2018']))
                 .Define("nbtag", "Sum(bjet)*1.0f")
                 .Define("PV_npvsGoodF","PV_npvsGood*1.0f")
                 .Define("SoftActivityJetNjets5F","SoftActivityJetNjets5*1.0f")
                 )
        return dftag

    if isGF:

        if (doSyst and isData == "false"):
            df = df.Define("Jet_delta",'computeJECuncertainties(corr_sf, Jet_pt, Jet_eta)').Vary("Jet_pt", "ROOT::RVec<ROOT::RVecF>{Jet_pt*(1-Jet_delta),Jet_pt*(1+Jet_delta)}", variationTags=["dn","up"], variationName = "JetSYST")
            #df = df.Define("Jet_delta", f'computeJECuncertainties({corr_sf}, Jet_pt, Jet_eta)').Vary("Jet_pt", "ROOT::RVec<ROOT::RVecF>{Jet_pt*(1-Jet_delta),Jet_pt*(1+Jet_delta)}", variationTags=["dn","up"], variationName = "JetSYST")

        dftag = (df.Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(vetoEle)+Sum(vetoMu))==0", "no leptons")
                 #                 .Define("trigger","{}".format(TRIGGER))
                 #                 .Filter("trigger>0", "pass triggers")
                 #.Filter("DeepMETResolutionTune_pt<75","DeepMETResolutionTune_pt<75") # not doing Zinv as nominal
                 .Define("goodJets","{}".format(GOODJETS))
                 .Define("nGoodJets","Sum(goodJets)*1.0f")
                 .Define("SoftActivityJetNjets5F","SoftActivityJetNjets5*1.0f")
#                 .Filter("Sum(goodJets)<2 or (Sum(goodJets)>=2 and Jet_pt[goodJets][0]<30) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]>0) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]<0 and abs(Jet_eta[goodJets][0] - Jet_eta[goodJets][1])<3 )","0 or 1 jets (pt20, |eta|<4.7) or >=2 with dEta<3")
                 .Define("mJJ","Sum(goodJets)>=2?Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets]):0.")
                 .Filter("Sum(goodJets)<2 or (Sum(goodJets)>=2 and Jet_pt[goodJets][0]<30) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]>0) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]<0 and abs(Jet_eta[goodJets][0] - Jet_eta[goodJets][1])<3 ) or (Sum(goodJets)>=2 and mJJ<300) ","0 or 1 jets (pt20, |eta|<4.7) or >=2 with dEta<3 or >=2 with mJJ<300")
                 )
        return dftag




def dfGammaMeson(df,PDType):

    TRIGGER=pickTRG(TRIGGERS,year,PDType,isVBF,isW,isZ,(isZinv or isVBFlow or isGF),isBPH)####isBPH argument added

    GOODphotons = ""
    GOODphotons = "{} and Photon_pt>75 and Photon_electronVeto".format(BARRELphotons) #90    
    print("PHOTONS = ", GOODphotons)

    dfOBJ = (df.Filter("nPhoton>0 and PV_npvsGood>0","photon from nano >0 and PV_npvsGood > 0")
             .Define("triggerAna","{}".format(TRIGGER))
             .Filter("triggerAna>0", "pass triggers")  ## comment while doing trigger studies
             .Define("loosePhotons", "{}".format(photonsLoose))
             .Define("nPhotonsVeto","Sum(loosePhotons)")
             .Define("goodPhotons", "{}".format(GOODphotons))
             .Define("nGoodPhotons","Sum(goodPhotons)*1.0f")
             .Filter("Sum(goodPhotons)>0", "At least one good Photon")
             .Define("goodPhotons_pt", "Photon_pt[goodPhotons]")
             .Define("goodPhotons_eta", "Photon_eta[goodPhotons]")
             .Define("goodPhotons_phi", "Photon_phi[goodPhotons]")
             .Define("goodPhotons_pfRelIso03_all", "Photon_pfRelIso03_all_quadratic[goodPhotons]")###quadratic
             .Define("goodPhotons_pfRelIso03_chg", "Photon_pfRelIso03_chg_quadratic[goodPhotons]")###quadratic
             .Define("goodPhotons_hoe", "Photon_hoe[goodPhotons]")
             .Define("goodPhotons_mvaID", "Photon_mvaID[goodPhotons]")
             .Define("goodPhotons_energyErr", "Photon_energyErr[goodPhotons]")
             .Define("goodPhotons_isScEtaEB", "Photon_isScEtaEB[goodPhotons]")
             .Define("jet_mask", "cleaningMask(Photon_jetIdx[goodPhotons],nJet)")
             )

    return dfOBJ




def dfHiggsCand(df, isData):

    GOODPHI = ""    
    GOODPHI = "{}".format(getMesonFromJson(mesons, "isVBF", "isPhiCat"))
    GOODRHO = ""
    GOODRHO = "{}".format(getMesonFromJson(mesons, "isVBF", "isRhoCat"))
    GOODK0STAR = ""
    GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isVBF", "isK0StarCat"))

    print("PHI = ", GOODPHI)
    dfbase = (df.Filter("nphi>0").Define("goodMeson","({}".format(GOODPHI)+" && {}".format(isPhiCat)+")")
                .Filter("Sum(goodMeson)>0", "one good Phi (ptPhi, validfit, ptTracks)")
                .Define("goodMeson_pt", "phi_kin_pt[goodMeson]")
                .Define("goodMeson_eta", "phi_kin_eta[goodMeson]")
                .Define("goodMeson_phi", "phi_kin_phi[goodMeson]")
                .Define("goodMeson_mass", "phi_kin_mass[goodMeson]")
                .Define("goodMeson_iso", "phi_iso[goodMeson]")
                .Define("goodMeson_isoNeuHad", "phi_isoNeuHad[goodMeson]")
                .Define("goodMeson_vtx_chi2dof", "phi_kin_vtx_chi2dof[goodMeson]")
                .Define("goodMeson_vtx_prob", "phi_kin_vtx_prob[goodMeson]")
                .Define("goodMeson_sipPV", "phi_kin_sipPV[goodMeson]")
                .Define("goodMeson_bestVtx_idx", "phi_bestVtx_idx[goodMeson]")
                .Define("goodMeson_bestVtx_X", "phi_bestVtx_X[goodMeson]")
                .Define("goodMeson_bestVtx_Y", "phi_bestVtx_Y[goodMeson]")
                .Define("goodMeson_bestVtx_Z", "phi_bestVtx_Z[goodMeson]")
                .Define("goodMeson_massErr", "phi_kin_massErr[goodMeson]")
                .Define("goodMeson_trk1_pt", "phi_trk1_pt[goodMeson]")
                .Define("goodMeson_trk2_pt", "phi_trk2_pt[goodMeson]")
                .Define("goodMeson_trk1_eta", "phi_trk1_eta[goodMeson]")
                .Define("goodMeson_trk2_eta", "phi_trk2_eta[goodMeson]")
                .Define("goodMeson_DR","DeltaR(phi_trk1_eta[goodMeson],phi_trk2_eta[goodMeson],phi_trk1_phi[goodMeson],phi_trk2_phi[goodMeson])")
                .Define("wrongMeson","({}".format(GOODRHO)+")")
                .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                .Define("wrongMeson2","({}".format(GOODK0STAR)+")")
                .Define("wrongMeson2_pt","Sum(wrongMeson2) > 0 ? K0Star_kin_pt[wrongMeson2]: ROOT::VecOps::RVec<float>(0.f)")
                ) 


    genMatchPDFNum='-1'
    if isPhiCat=="true": genMatchPDFNum='333'
    if isRhoCat=="true": genMatchPDFNum='113'
    if isK0StarCat=="true": genMatchPDFNum='313'
    if isOmega3PiCat=="true": genMatchPDFNum='223'
    if isPhi3PiCat=="true": genMatchPDFNum='333'
    
    
    dfHig = (dfbase.Define("index_pair","HiggsCandFromRECO(goodMeson_pt, goodMeson_eta, goodMeson_phi, goodMeson_mass, goodMeson_trk1_pt, goodMeson_trk2_pt, wrongMeson_pt, wrongMeson2_pt, goodPhotons_pt, goodPhotons_eta, goodPhotons_phi)").Filter("index_pair[0]!= -1", "at least a good meson candidate")
             .Define("jet_mask2", "cleaningJetFromOBJ(Jet_eta, Jet_phi, goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]])")
             .Define("photon_pt", "(index_pair[1]!= -1) ? goodPhotons_pt[index_pair[1]]: 0.f")
             .Define("meson_isoNeuHad", "(index_pair[0]!= -1) ? goodMeson_isoNeuHad[index_pair[0]]: 0.f")
             )

    if(isOmega3PiCat=="true" or isPhi3PiCat=="true" or isD0StarCat=="true" or isD0Pi0StarCat=="true"): dfHig_ = callMVAregress(dfHig,isVBF,isVBFlow,isGF,isZinv)
    if(isRhoCat=="true" or isPhiCat=="true" or isK0StarCat=="true"): dfHig_ = dfHig.Define("meson_pt", "(index_pair[0]!= -1) ? goodMeson_pt[index_pair[0]]: 0.f")


    if (isData == "false"):
        dfFinal = dfHig_.Define("HCandMass", "compute_HiggsVars_var(meson_pt,goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],0)")
               
    return dfFinal




def dfCommon(df,year,isData,mc,sumw,isVBF,isVBFlow,isGF,isZinv):

    lumi = 1.
    weight = "{0}".format(1.)
    if mc>=0: weight = "{0}*genWeight*{1}".format(lumi,sumw)

    lumiIntegrated = 1.
    print('isData = ',isData)
    if (isData == "false"):
        if((isVBF or isW or isZ) and year == 2018): lumiIntegrated = lumis['2018']
        if((isW or isZ) and year == 2017): lumiIntegrated = lumis['2017']
        if((isVBF) and year == 2017): lumiIntegrated = lumis['12017']
        if((isVBF or isW or isZ) and year == 12016): lumiIntegrated = lumis['12016']
        if((isW or isZ) and year == 22016): lumiIntegrated = lumis['22016']
        if((isVBFlow or isGF or isZinv) and year == 2018): lumiIntegrated = lumis['12018']
        print('lumiIntegrated=',lumiIntegrated, ' year=',year)

    dfComm = (df
              .Define("mc","{}".format(mc))
              .Define("isData","{}".format(isData))
              .Define("applyJson","{}".format(JSON)).Filter("applyJson","pass JSON")
              .Define("w","{}".format(weight))
              .Define("wraw","{}".format(weight))
              .Define("lumiIntegrated","{}".format(lumiIntegrated))
              .Filter("PV_npvsGood>0","one good PV")
              )

    return dfComm



    
def analysis(df,year,mc,sumw,isData,PDType):

    dfCom = dfCommon(df,year,isData,mc,sumw,isVBF,isVBFlow,isGF,isZinv)
    dfOBJ = dfGammaMeson(dfCom,PDType)
    dfbase = dfHiggsCand(dfOBJ,isData)
    #dfFINAL = selectionTAG(dfbase,doSyst,isData)

    output_file = "/eos/user/e/eferrand/Work/CMSSW_13_0_13/src/Hrare/Snapshot_HPhiGamma.root"
    
    snapshot_tdf = dfbase.Snapshot("Events", output_file)
    #snapshot_tdf = dfFINAL.Snapshot("Events", output_file)



    hists = {
             "HCandMass":  {"name":"HCandMass", "title":"H mass; m_{k^{+}k^{-}#gamma} (GeV); N_{Events}", "bin":70, "xmin":100, "xmax":170}
            }

    histos = []

    for h in hists:
        # 1D is for nom only
        model = (hists[h]["name"], hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
#       h1d = dfFINAL.Histo1D(model, hists[h]["name"], "w_allSF")
        h1d = dfbase.Histo1D(model, hists[h]["name"])
        histos.append(h1d)

    for h in histos:
        canv = ROOT.TCanvas("canvas","canvas",800,700)
        h.Draw("hist")

        canv.SaveAs("/eos/user/e/eferrand/Work/CMSSW_13_0_13/src/Hrare/"+h.GetName()+".png")




def runTest():
   
    fileName = "ntuple_HPhiGamma.root"
    treeName = "Events"
    df = ROOT.RDataFrame(treeName, fileName)
    
    w=1.
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)
    #########loadCorrectionSet(year)
    sampleNOW=-1
    analysis(df,2018,sampleNOW,w,"false","NULL")

    

def readMCSample(year,sampleNOW,useD03):

    fileName = "ntuple_HPhiGamma.root"
    treeName = "Events"
    runTreeName = "Runs"
    df = ROOT.RDataFrame(treeName, fileName)
    rdf = ROOT.RDataFrame(runTreeName, fileName)

    sumW = computeWeigths(df, rdf, sampleNOW, year, True, useD03)
    if (doPol and sampleNOW>1000 and sampleNOW<1039): loadPolarizationTree(sampleNOW,year)
    if doSyst: loadSFhisto(sampleNOW,year)
    loadCorrectionSet(year)
    analysis(df,year,sampleNOW,sumW,"false","NULL")



loadUserCode()

runTest()
#readMCSample(int(sys.argv[4]),int(sys.argv[3]),useD03)



print("Analysis done")

#to run: python3 -i HReco.py isGFtag isPhiCat 12 2018