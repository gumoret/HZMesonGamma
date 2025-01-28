#define _USE_MATH_DEFINES
#include <cmath> 
#include <iostream>

//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>


#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "FWCore/Framework/interface/stream/EDAnalyzer.h"////////
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

//Vertex inclusions
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/BeamSpot/interface/BeamSpot.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Electron ID stuff
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"//////
#include "CommonTools/Egamma/interface/ConversionTools.h"/////
//#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

//Photon ID stuff
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//Parton distribution and QCD scale variations stuff 
#include "FWCore/Framework/interface/Run.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h" //LHE reader
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h" //LHE reader

//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" //JEC uncertainties
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h" //JEC uncertainties
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

//JEC uncertainties
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
//#include "CondFormats/JetMETObjects/interface/JetResolution.h"
//#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//#include "CondFormats/DataRecord/interface/JetCorrectionsRecord.h"

//#include "CondFormats/JetMETObjects/interface/JetCorrectionsRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


typedef math::XYZTLorentzVector LorentzVector;

using namespace std;  

#include "ZMesonGamma.h"

// constructors and destructor
ZMesonGamma::ZMesonGamma(const edm::ParameterSet& iConfig) :
runningOnData_(iConfig.getParameter<bool>("runningOnData")),
verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose"))
{
  packedPFCandidatesToken_            = consumes<std::vector<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 
  slimmedMuonsToken_                  = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  prunedGenParticlesToken_            = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  photonsMiniAODToken_                = consumes<std::vector<pat::Photon> > (edm::InputTag("slimmedPhotons"));
  electronsMiniAODToken_              = consumes<std::vector<pat::Electron> > (edm::InputTag("slimmedElectrons"));
  slimmedJetsToken_                   = consumes<std::vector<pat::Jet> >(edm::InputTag("slimmedJets"));
  slimmedMETsToken_                   = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));
  slimmedMETsPuppiToken_              = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETsPuppi"));
  offlineSlimmedPrimaryVerticesToken_ = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));  
  offlineBeamSpotToken_               = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));
  pileupSummaryToken_                 = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  GenInfoToken_                       = consumes<GenEventInfoProduct> (edm::InputTag("generator"));
  triggerBitsToken_                   = consumes<edm::TriggerResults> (edm::InputTag("TriggerResults","","HLT"));
  rhoToken_                           = consumes<double> (iConfig.getParameter <edm::InputTag>("rho"));
  //packedGenParticlesToken_            = consumes<std::vector<pat::GenParticle>>(edm::InputTag("packedGenParticles", "", "PAT"));

  hEvents = fs->make<TH1F>("hEvents", "Event counting in different steps", 6, 0., 6.);

  nEventsProcessed           = 0;
  nEventsTriggered           = 0;
  nEventsIsTwoKaons          = 0;
  nEventsIsPhoton            = 0;
  nEventsZMatched            = 0;
  nEventsZNotMatched         = 0;
  nEventsMesonPtNotMatched   = 0;
  nEventsBestPairFound       = 0;
  nEventsTrkPtFilter         = 0;
  nEventsPairIsolationFilter = 0;


  debug=false;  //DEBUG datamember 
  verbose=false; 
  

  hPileup   = fs->make<TH1F>("pileup", "pileup", 130, 0, 130);

  create_trees();
}


ZMesonGamma::~ZMesonGamma()
{
}


// ------------ method called for each event  ------------
void ZMesonGamma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::PackedCandidate> > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon> > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle> > prunedGenParticles;
  if(!runningOnData_)iEvent.getByToken(prunedGenParticlesToken_, prunedGenParticles);

  edm::Handle<std::vector<pat::Photon> > slimmedPhotons;
  iEvent.getByToken(photonsMiniAODToken_,slimmedPhotons);

  edm::Handle<std::vector<pat::Electron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<pat::Jet> > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_, slimmedJets);

  edm::Handle<std::vector<pat::MET> > slimmedMETs;
  iEvent.getByToken(slimmedMETsToken_, slimmedMETs);

  edm::Handle<std::vector<pat::MET> > slimmedMETsPuppi;
  iEvent.getByToken(slimmedMETsPuppiToken_, slimmedMETsPuppi);

  edm::Handle<std::vector<reco::Vertex> > slimmedPV;
  iEvent.getByToken(offlineSlimmedPrimaryVerticesToken_, slimmedPV);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits); 


  nEventsProcessed++; //This will be saved in the output tree, giving the number of processed events

  //Retrieve the run number
  if(runningOnData_){
    runNumber = iEvent.id().run();
  }

  eventNumber = iEvent.id().event();


  //*************************************************************//
  //                                                             //
  //-------------------------- Vertices -------------------------//
  //                                                             //
  //*************************************************************//

  //Count the number of vertices and return if there's no vertex
  nPV = 0;
  if(slimmedPV->size()<=0){
    if(verbose) cout<<"No primary vertex found, RETURN."<<endl;
    return;
  }

  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      nPV++;
    }
  } 
  

  //*************************************************************//
  //                                                             //
  //--------------------------- Pile Up -------------------------//
  //                                                             //
  //*************************************************************//

  PU_Weight = -1.;
  float npT = -1.;

  if(!runningOnData_){
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI; 

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      const int BX = PVI->getBunchCrossing();
      if(BX == 0) {
        npT  = PVI->getTrueNumInteractions();
      }
    }

    if(npT == -1) {
      std::cout << "!!!! npT = -1 !!!!" << std::endl;
      abort();
    }

    // Calculate weight using above code
    //PU_Weight = Lumiweights_.weight(npT);######################

    // Fill histogram with PU distribution
    hPileup->Fill(npT);
  }



  //*************************************************************//
  //                                                             //
  //-------------------------- MC Weight ------------------------//
  //                                                             //
  //*************************************************************//

  MC_Weight = -10000000.;

  if(!runningOnData_){
    edm::Handle<GenEventInfoProduct> GenInfo;
    iEvent.getByToken(GenInfoToken_, GenInfo);
    
    float _aMCatNLOweight = GenInfo->weight();
    MC_Weight = _aMCatNLOweight;

    if(MC_Weight == -10000000.) {
      std::cout << "!!!! MC_Weight = -10000000 !!!!" << std::endl;
      abort();
    }
  }



  //*************************************************************//
  //                                                             //
  //--------------------------- Trigger -------------------------//
  //                                                             //
  //*************************************************************//

  //Examine the trigger information, return if the trigger doesn't switch on and count the number of events where the trigger has switched on
  isTwoProngTrigger = false;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for(unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);

    if( tmp_triggername.find("HLT_Photon35_TwoProngs35_v") != std::string::npos ){
      isTwoProngTrigger = true;
    }
  }

  if(!isTwoProngTrigger){
   if(verbose) cout<<"Event not triggered, RETURN."<<endl;
   return;
  }
  nEventsTriggered++;


  
  //*************************************************************//
  //                                                             //
  //------------------ Variable initialization ------------------//
  //                                                             //
  //*************************************************************//

  nMuons10              = 0;
  nElectrons10          = 0;
  nMuons20              = 0;
  nElectrons20          = 0;
  nPhotons38WP80        = 0;
  nPhotons20WP90        = 0;
  nPhotonsChosen        = 0;
  nJets30               = 0;
  nJets25               = 0;
  nJets20               = 0;
 
  //These variables will go in the tree
  photonEt               = 0.;
  photonEta              = 0.;
  photonEtaSC            = 0.;
  photonPhi              = 0.;
  LorentzVector ph_p4;
  photonIsoChargedHadron = 0.;
  photonIsoNeutralHadron = 0.;
  photonIsoPhoton        = 0.;
  photonIsoEArho         = 0.;
  photonRegressionError  = 0.;
  photonEtMax        = -1000.;

  jetPhotonInvMass  = -1.;
  mesonMass         = -1.;
  ZMassFrom2KPhoton = -1.;
 
  metPt      = 0.;
  metpuppiPt = 0.;

  bestJetPt                       =-1.;
  bestJetEta                      =-1.;
  bestJetPhi                      =-1.;
  bestJetnDaughters               = 0;
  bestJetPtMax                    =-1.;
  bestJetChargedEmEnergy          = 0.;
  bestJetNeutralEmEnergy          = 0.;
  bestJetChargedHadEnergy         = 0.;
  bestJetNeutralHadEnergy         = 0.;
  bestJetChargedEmEnergyFraction  = 0.;
  bestJetNeutralEmEnergyFraction  = 0.;
  bestJetChargedHadEnergyFraction = 0.;
  bestJetNeutralHadEnergyFraction = 0.;
  bestJetChargedHadMultiplicity   = 0.;
  bestJetInvMass                  = 0.;
  bestJetPhotonInvMass            = 0.;
  bestJetJECunc                   = 0.;
  firstTrkPt                      = 0.;
  firstTrkEta                     = 0.;
  firstTrkPhi                     = 0.;
  firstTrkCharge                  = 0.;
  secondTrkPt                     = 0.;
  secondTrkEta                    = 0.;
  secondTrkPhi                    = 0.;
  secondTrkCharge                 = 0.;
  bestPairPt                      = 0.;
  bestPairEta                     = 0.;
  bestPairPhi                     = 0.;  


  //K-candidates and PHI ISOLATION
  

  isoK1     = 0.;
  isoK1Ch   = 0.;
  isoK2     = 0.;
  isoK2Ch   = 0.;
  isoPair   = 0.;
  isoPairCh = 0.;
  


  //*************************************************************//
  //                                                             //
  //----------------------------- MET ---------------------------//
  //                                                             //
  //*************************************************************//
  for(auto met = slimmedMETs->begin(); met != slimmedMETs->end(); ++met){
  metPt = met->pt();
  }

  for(auto metpuppi = slimmedMETsPuppi->begin(); metpuppi != slimmedMETsPuppi->end(); ++metpuppi){
  metpuppiPt = metpuppi->pt();
  }



  //*************************************************************//
  //                                                             //
  //---------------------------- Muons --------------------------//
  //                                                             //
  //*************************************************************//
  //Count muons for each event
  for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
    if(mu->pt() < 10. || !mu->CutBasedIdMedium || fabs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
  if(!mu->PFIsoLoose) continue;
  nMuons10++;
  if(mu->pt() < 20.) continue;
  nMuons20++;
  }



  //*************************************************************//
  //                                                             //
  //-------------------------- Electrons ------------------------//
  //                                                             //
  //*************************************************************//
  //Count the number of electrons
  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  float corr_pt = 0.;


  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    //Calculate electron p4, correct it with the Scale&Smearing correction and extract the pT
    LorentzVector el_p4 = el->p4();// * el->userFloat("ecalTrkEnergyPostCorr")/el->energy();
    corr_pt = el_p4.pt();

    if(corr_pt < 10. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

    //float abseta = fabs(el->superCluster()->eta());
    //float eA     = effectiveAreas_el_.getEffectiveArea(abseta);
    //float el_iso   = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/corr_pt;
    //if(el_iso > 0.35) continue;

    //-------------Conditions on loose/medium MVA electron ID-------------//
    if(el->electronID("mvaEleID-Fall17-iso-V2-wp80") == 0) continue;
    nElectrons10++;
    if (corr_pt < 20.) continue;
    nElectrons20++;
  }

  //std::cout << "Nelectrons " << nElectrons << " Nmuons " << nMuons << std::endl;



  //*************************************************************//
  //                                                             //
  //--------------------------- Photons -------------------------//
  //                                                             //
  //*************************************************************//
  if(verbose) cout<< "PHOTONs"<<" --------------------------------"<<endl;
  
  bool cand_photon_found = false; //initialize this bool to false, return if it doesn't turn into true
  float corr_et = -1.;

  for(auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){ //PHOTON FORLOOP START --------------------------------
    
    // Apply energy scale corrections, from 18Apr2023 the correction are embedded in the config file with the postReco tool
    corr_et   = photon->et();// * photon->userFloat("ecalEnergyPostCorr") / photon->energy(); 

    if(corr_et < 20. || fabs(photon->eta()) > 2.5) continue; //loose selection to reject diphoton bkg 
    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp90") == 0) continue; //WP90
    if(!photon->passElectronVeto()) continue; 

    nPhotons20WP90++;

    if(corr_et < 35.) continue; /////////////////////////////
    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp80") == 0) continue; //WP80

    //float abseta = fabs(photon->superCluster()->eta());
    //float eA = effectiveAreas_ph_.getEffectiveArea(abseta);

    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();

    //if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6

    nPhotons38WP80++;

    // Apply energy scale corrections
    ph_p4 = photon->p4();// * photon->userFloat("ecalEnergyPostCorr") / photon->energy();

    if(corr_et < photonEtMax) continue; //choose as best photon the one with highest eT
    
    photonEtMax = corr_et;
    photonIsoChargedHadron = photon->chargedHadronIso();
    photonIsoNeutralHadron = photon->neutralHadronIso();
    photonIsoPhoton        = photon->photonIso();
    //photonIsoEArho         = eA * rho_;
    photonEt     = corr_et;
    //ph_energy = photon->energy();
    photonEta    = photon->eta();
    photonEtaSC  = photon->superCluster()->eta();
    photonPhi    = photon->phi();

    photonRegressionError = photon->getCorrectedEnergyError(reco::Photon::P4type::regression2);
    if(debug) cout << "Regression2 Energy Error: " << photonRegressionError << endl;
    
    cand_photon_found = true;
    nPhotonsChosen++;


    
  }//PHOTON FORLOOP END -------------------------------------------------------------------------------------------------------------------------

  //Return if there are no photons chosen
  if(!cand_photon_found) {
    cout<<"No best photon found, RETURN."<<endl;
    return;
  }

  nEventsIsPhoton++;



  //*************************************************************//
  //                                                             //
  //--------------------------- N-jets --------------------------//
  //                                                             //
  //*************************************************************//

  //int nJet=1;
  int jetIndex      = -1;
  int bestJet_Index = -1;
  //int MCtruthIndex = -1;
  float deltaR      = -1;   
  int nDaughters    = 0;
  //bool isBestJetFound = false; 

  //daughters forloop variable
  int firstTrkCharge;
  int secondTrkCharge;
  //float firstTrkPt;
  //float secondTrkPt;
  LorentzVector firstTrkP4;
  LorentzVector secondTrkP4;
  LorentzVector firstTrkP4K;
  LorentzVector secondTrkP4K;
  LorentzVector firstTrkP4Pi;
  LorentzVector secondTrkP4Pi;
  LorentzVector pairP4;
  LorentzVector pairP4K;
  LorentzVector PairP4Pi;
  LorentzVector bestFirstTrkP4;
  LorentzVector bestSecondTrkP4;
  LorentzVector bestPairP4;
  float bestCoupleOfTheJetPt = 0.;
  float deltaRKChosen        = 0.;
  float deltaRK              = 0.;
  firstTrkEnergyK            = 0.;
  secondTrkEnergyK           = 0.;
  firstTrkEnergyPi           = 0.;
  secondTrkEnergyPi          = 0.;
  firstTrkPx                 = 0.;
  firstTrkPy                 = 0.;
  firstTrkPz                 = 0.;
  float firstTrkDxy          = -999.;
  float firstTrkDxyErr       = -999.;
  float firstTrkDz           = -999.;
  float firstTrkDzErr        = -999.;
  bestFirstTrkDxy            = -999.;
  bestFirstTrkDz             = -999.;
  bestFirstTrkDxyErr         = -999.;
  bestFirstTrkDzErr          = -999.;
  secondTrkPx                = 0.;
  secondTrkPy                = 0.;
  secondTrkPz                = 0.;
  float secondTrkDxy         = -999.;
  float secondTrkDxyErr      = -999.;
  float secondTrkDz          = -999.;
  float secondTrkDzErr       = -999.;
  bestSecondTrkDxy           = -999.;
  bestSecondTrkDz            = -999.;
  bestSecondTrkDxyErr        = -999.;
  bestSecondTrkDzErr         = -999.;
  firstTrkEta                = 0.;
  firstTrkPhi                = 0.;
  secondTrkEta               = 0.;
  secondTrkPhi               = 0.;
  float PhiMass              = 0.;
  float RhoMass              = 0.;
  float kMass                = 0.4937;
  float PiMass               = 0.13957;
  bool isBestCoupleOfTheEventFound = false;
  bool isPhi                       = false;
  bool isRho                       = false;


  bool isJet = false;
  bool jet_verbose=false;

  
  if(isJet){
  //JET LOOP
  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet) { //JET LOOP START -------------------------------------------------------- 
    if(verbose||jet_verbose) cout << "jet loop starts" << endl;
    jetIndex++;

    jetPhotonInvMass=(jet->p4()+ph_p4).M(); //calculate inv mass
    nDaughters= jet->numberOfDaughters(); //calculate number of daughters

    //----------------------------- Pre-Filters --------------------------------------------------------
    float neutralHadEnergyFrac = jet->neutralHadronEnergyFraction();
    float neutralEmEnergyFrac  = jet->neutralEmEnergyFraction();
    float muonEnergyFrac       = jet->muonEnergyFraction();    
    float chargedHadEnergyFrac = jet->chargedHadronEnergyFraction();
    float chargedHadMult       = jet->chargedHadronMultiplicity();
    float chargedEmEnergyFrac  = jet->chargedEmEnergyFraction();
    float pt                   = jet->pt();
    float eta                  = abs(jet->eta());

    if(neutralHadEnergyFrac > 0.9 || neutralEmEnergyFrac > 0.9 || nDaughters < 2. || muonEnergyFrac > 0.8 || chargedHadEnergyFrac <= 0. || chargedHadMult == 0. || chargedEmEnergyFrac > 0.8 || pt < 20. || eta > 4.7) {
      if(verbose) cout << "if .... continue" << endl;
      continue;
      }
    
    if(jet->pt() < 38. || abs(jet->eta()) > 2.5) {
      if(verbose||jet_verbose) cout << "if jet pt<38 continue" << endl;
      continue;
    }
    if(jetPhotonInvMass < 30.) {
      if(verbose||jet_verbose) cout << "if jet mass<30 continue" << endl;
      continue; //reject jets with inv mass lower then 30 GeV
    }

                           
     //-------------------------------------------------------------------------------------------------      
    
    if (verbose) cout<<" Jet at index = "<<jetIndex<<" passed the cuts:"<<endl; 
       

    //-------------------------------------daughters forloop----------------------------

    for(int firstTrkIndex=0; firstTrkIndex < nDaughters; firstTrkIndex++){ //1ST LOOP STARTS
      if(verbose) cout << "1st particle loop starts" << endl;

      if (verbose) cout<<"Daughter n."<<firstTrkIndex+1<<" pT = "<<slimmedJets->at(jetIndex).daughter(firstTrkIndex)->pt()<<endl;

      firstTrkCharge = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->charge();  //take firstCand charge
      firstTrkPt    = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->pt();  //take firstCand pt
      firstTrkEta   = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->eta(); //take firstCand eta
      firstTrkPhi   = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->phi(); //take firstCand phi
      if(slimmedJets->at(jetIndex).daughter(firstTrkIndex)->bestTrack() == NULL) {
        if(verbose) cout << "if bestTrack null continue" << endl;
        continue;//loop only over charged daughters
        }
      firstTrkDxy    = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->bestTrack()->dxy((&slimmedPV->at(0))->position()); //take firstCand dxy
      firstTrkDxyErr = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->bestTrack()->dxyError();        
      firstTrkDz     = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->bestTrack()->dz((&slimmedPV->at(0))->position()); //take firstCand dz
      firstTrkDzErr  = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->bestTrack()->dzError();

        
      if(firstTrkCharge == 0 || abs(firstTrkDxy) >= 0.2 || abs(firstTrkDz) >= 0.5 || !(slimmedJets->at(jetIndex).daughter(firstTrkIndex)->bestTrack()->quality(reco::Track::highPurity))) {
          if(verbose) cout << "if 1st trk...... continue" << endl;
          continue;
          }


      if(firstTrkPt < 1.){
        if(verbose) cout << "if 1st trk pt <1 continue" << endl;
        continue; //firstCand filter if pT < 1 GeV
        }


      for(int secondTrkIndex=firstTrkIndex+1; secondTrkIndex < nDaughters; secondTrkIndex++){ //2ND LOOP STARTS
        if(verbose) cout << "2nd particle loop starts" << endl;

        secondTrkCharge = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->charge();
        secondTrkPt     = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->pt();
        secondTrkEta    = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->eta();
        secondTrkPhi    = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->phi();
        if(slimmedJets->at(jetIndex).daughter(secondTrkIndex)->bestTrack() == NULL) {
          if(verbose) cout << "if bestTrack null continue" << endl;
          continue;
          }
        secondTrkDxy    = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->bestTrack()->dxy((&slimmedPV->at(0))->position());
        secondTrkDxyErr = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->bestTrack()->dxyError();          
        secondTrkDz     = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->bestTrack()->dz((&slimmedPV->at(0))->position());
        secondTrkDzErr  = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->bestTrack()->dzError();

        //if (slimmedJets->at(jetIndex).daughter(secondTrkIndex)->charge() == 0) continue;
        if(secondTrkCharge == 0 || abs(secondTrkDxy) >= 0.2 || abs(secondTrkDz) >= 0.5 || !(slimmedJets->at(jetIndex).daughter(secondTrkIndex)->bestTrack()->quality(reco::Track::highPurity))) {
          if(verbose) cout << "if 2nd trk...... continue" << endl;
          continue;
          }


        //TRKs PT CUT --------------------------------------------------------------------------
        //secondTrkPt  = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->pt();
        if(secondTrkPt < 1.) {
          if(verbose) cout << "if second trk pt<1 continue" << endl;
          continue; //firstCand filter if pT < 1 GeV
          }

        if(firstTrkPt < 10. && secondTrkPt < 10.) {
          if(verbose) cout << "if firstTrkPt < 10. && secondTrkPt < 10 continue" << endl;
          continue;  //filter if both cand pT are < 10GeV
          }

        //DITRK DELTA R CUT --------------------------------------------------------------------------
        //secondTrkEta = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->eta();
        //secondTrkPhi = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->phi();
        
        float deltaEta= firstTrkEta - secondTrkEta;

        float deltaPhi = fabs(firstTrkPhi - secondTrkPhi);  //phi folding 
        if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

        deltaRK= sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
        if(deltaRK > 0.07) {
          if(verbose) cout << "deltaRK<0.07 continue" << endl;
          continue;
          }

        //OPPOSITE CHARGE - FILTER ------------------------------------------------------------
        //firstTrkCharge  = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->charge(); //take firstCand charge
        //secondTrkCharge = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->charge(); //take secondCand charge
        if(firstTrkCharge * secondTrkCharge >= 0) {
          if(verbose) cout << "firstTrkCharge * secondTrkCharge >= 0" << endl;
          continue; //choose only opposite charges
          }


        //QUADRIMOMENTUM CALCULATION ------------------------------------------------------------
        firstTrkP4  = slimmedJets->at(jetIndex).daughter(firstTrkIndex)->p4(); //take quadrimomentum
        secondTrkP4 = slimmedJets->at(jetIndex).daughter(secondTrkIndex)->p4();

        firstTrkPx  = firstTrkP4.px(); //take px, py, pz of the first candidate
        firstTrkPy  = firstTrkP4.py();
        firstTrkPz  = firstTrkP4.pz();
        secondTrkPx = secondTrkP4.px(); //take px, py, pz of the second candidate
        secondTrkPy = secondTrkP4.py();
        secondTrkPz = secondTrkP4.pz();

        //PIONS OR KAONS HYPOTHESIS -----------------------------------------------------------------------------------------------------------------------------------------------            
        firstTrkEnergyK   = sqrt(firstTrkPx  * firstTrkPx  + firstTrkPy  * firstTrkPy  + firstTrkPz  * firstTrkPz  + kMass  * kMass ); //Kaon hypothesis energy recalculation
        secondTrkEnergyK  = sqrt(secondTrkPx * secondTrkPx + secondTrkPy * secondTrkPy + secondTrkPz * secondTrkPz + kMass  * kMass ); //Kaon hypothesis energy recalculation
        firstTrkEnergyPi  = sqrt(firstTrkPx  * firstTrkPx  + firstTrkPy  * firstTrkPy  + firstTrkPz  * firstTrkPz  + PiMass * PiMass); //Pion hypothesis energy recalculation
        secondTrkEnergyPi = sqrt(secondTrkPx * secondTrkPx + secondTrkPy * secondTrkPy + secondTrkPz * secondTrkPz + PiMass * PiMass); //Pion hypothesis energy recalculation
        
        /*
        if (verbose) {
          cout<<"firstTrkEnergyK   = "<<firstTrkEnergyK<<endl;
          cout<<"firstTrkEnergyPi  = "<<firstTrkEnergyPi<<endl;
          cout<<"secondTrkEnergyK  = "<<secondTrkEnergyK<<endl;
          cout<<"secondTrkEnergyPi = "<<secondTrkEnergyPi<<endl;
        }
        */

        firstTrkP4K   = firstTrkP4.SetE(firstTrkEnergyK); //Kaon hypothesis quadrimomentum correction
        secondTrkP4K  = secondTrkP4.SetE(secondTrkEnergyK); //Kaon hypothesis quadrimomentum correction
        firstTrkP4Pi  = firstTrkP4.SetE(firstTrkEnergyPi); //Kaon hypothesis quadrimomentum correction
        secondTrkP4Pi = secondTrkP4.SetE(secondTrkEnergyPi); //Kaon hypothesis quadrimomentum correction

        pairP4K  = firstTrkP4K  + secondTrkP4K; //calculation of the couple-quadrimomentum after the correction
        PairP4Pi = firstTrkP4Pi + secondTrkP4Pi; //calculation of the couple-quadrimomentum after the correction
        
        if (verbose) {
          cout<<"KK pT = "<<pairP4K.pt()<<endl;
          cout<<"PiPi pT = "<<PairP4Pi.pt()<<endl;
        }

        //DITRK PT CUT -------------------------------------------------------------------------
        if(pairP4K.pt() < 38.) {
          if(verbose || jet_verbose) cout<<"couplePt cut NOT passed, continue"<<endl;
          continue;
        }  
        
        //MESON INV MASS CUT -------------------------------------------------------------------------
        isPhi = false;
        isRho = false;

        PhiMass = (pairP4K).M(); //calculate inv mass of the Phi candidate  
        if (verbose || jet_verbose) cout<<"mKK (before the meson mass selection) =  "<<PhiMass<<endl;
        if(PhiMass > 1. && PhiMass < 1.05) isPhi = true; //filter on Phi invariant mass  

        RhoMass = (PairP4Pi).M(); //calculate inv mass of the Rho candidate  
        if (verbose) cout<<"mPiPi (before the meson mass selection) =  "<<RhoMass<<endl;
        if(RhoMass > 0.5 && RhoMass < 1.) isRho = true; //filter on Rho invariant mass   

        if (!isPhi && !isRho) {
          if (verbose || jet_verbose) cout<<"the pair mass doesn't match any of the two mass hypothesis, continue "<<endl;
          continue; //continue if the pair mass doesn't match any of the two mass hypothesis
          }

        if (isPhi && isRho){ //if both hypothesis are true, mark it as a Phi candidate (this is done because the Phi mass window is tighter)
          isPhi = true;
          isRho = false;
        }

        //update values of quadrimomenta
        if(isPhi){
          firstTrkP4  = firstTrkP4K; 
          secondTrkP4 = secondTrkP4K;
          pairP4      = pairP4K;
        }  
        if(isRho){
          firstTrkP4  = firstTrkP4Pi;
          secondTrkP4 = secondTrkP4Pi;
          pairP4      = PairP4Pi;
        }


        K1SumPt05     = 0.;
        K1SumPt05Ch   = 0.;

        K2SumPt05     = 0.;
        K2SumPt05Ch   = 0.;

        pairSumPt05   = 0.;
        pairSumPt05Ch = 0.;

        // ISOLATION CUT -------------------------------------------------------------------------  
        for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){ //ISOLATION FORLOOP START
          if(verbose) cout << "isolation loop starts" << endl;
          if(debug){
            cout <<endl<<"ISO CALC DETAILS ---------------------"<<endl;
            cout << "pt cand_iso = "<<cand_iso->pt()<<endl;
          }
          if(cand_iso->pt() < 0.5) {
            if(verbose) cout << "cand_iso->pt() < 0.5, continue" << endl;
            continue; //do not consider tracks with pT < 500MeV
            }

          //calculate the deltaR between the track and the first candidate ---------------------------------------
          float deltaPhi_K1 = fabs(firstTrkP4.phi()-cand_iso->phi());  //phi folding 
          if (deltaPhi_K1 > M_PI) deltaPhi_K1 = 2*M_PI - deltaPhi_K1;

          float deltaRK1 = sqrt((firstTrkP4.eta()-cand_iso->eta())*(firstTrkP4.eta()-cand_iso->eta()) + deltaPhi_K1*deltaPhi_K1);
          if (debug) cout << "deltaRK1 = "<<deltaRK1<<endl;
          if(deltaRK1 < 0.0005) {
            if(verbose) cout<<"deltaRK1 < 0.0005, continue" << endl;
            continue; //remove first candidate from the sum
            }

          //calculate the deltaR between the track and the second candidate ---------------------------------------
          float deltaPhi_K2 = fabs(secondTrkP4.phi()-cand_iso->phi());  //phi folding  
          if (deltaPhi_K2 > M_PI) deltaPhi_K2 = 2*M_PI - deltaPhi_K2;

          float deltaRK2 = sqrt((secondTrkP4.eta()-cand_iso->eta())*(secondTrkP4.eta()-cand_iso->eta()) + deltaPhi_K2*deltaPhi_K2);
          if (debug) cout << "deltaRK2 = "<<deltaRK2<<endl;
          if(deltaRK2 < 0.0005) {
            if(verbose) cout << "deltaRK2 < 0.0005, continue" << endl;
            continue; //remove second candidate from the sum
            }

          //calculate the deltaR between the track and the best pair ---------------------------------------
          float deltaPhi_Couple = fabs(pairP4.phi()-cand_iso->phi());  //phi folding  
          if (deltaPhi_Couple > M_PI) deltaPhi_Couple = 2*M_PI - deltaPhi_Couple;

          float deltaR_Couple = sqrt((pairP4.eta()-cand_iso->eta())*(pairP4.eta()-cand_iso->eta()) + deltaPhi_Couple*deltaPhi_Couple);

          //sum pT of the tracks inside a cone of deltaR = 0.3 ---------------------------------------
          if(deltaRK1 <= 0.3) K1SumPt05 += cand_iso->pt();
          if(deltaRK2 <= 0.3) K2SumPt05 += cand_iso->pt();
          if(deltaR_Couple <= 0.3) pairSumPt05 += cand_iso->pt();
          //cout<< "charge before = "<<cand_iso->charge()<<endl;

          //sum pT of the charged tracks inside a cone of deltaR = 0.3 ---------------------------------------
          if (debug) cout << "Charge = "<< cand_iso->charge()<<endl;
          if(cand_iso->charge() == 0) {
            if(verbose) cout << "cand_iso->charge() == 0, continue" << endl;
            continue;
            }
          // cout << "particle charge = "<<cand_iso->charge()<<endl;
          if (debug) cout << "dxy = "<<fabs(cand_iso->dxy())<<" and dz = "<< fabs(cand_iso->dz())<<endl;
          if(fabs(cand_iso->dxy()) >= 0.2 || fabs(cand_iso->dz()) >= 0.5) {
            if(verbose) cout << "cand_iso->dxy()) >= 0.2 || cand_iso->dz() >= 0.5, continue" << endl;
            continue; // Requesting charged particles to come from PV
            }
          //cout<< "charge after = "<<cand_iso->charge()<<endl;
          if(deltaRK1 <= 0.3) K1SumPt05Ch += cand_iso->pt();
          if(deltaRK2 <= 0.3) K2SumPt05Ch += cand_iso->pt();
          if (debug) cout <<"deltaR_Couple = "<<deltaR_Couple<<endl;
          if(deltaR_Couple <= 0.3){
            pairSumPt05Ch += cand_iso->pt();
            if (verbose) cout<<"Particle in the cone: SumPt = "<<pairSumPt05Ch<<endl;
          }
          if(verbose) cout << "isolations loop ends" << endl;
        } //ISOLATION FORLOOP END

        float isoCoupleCh = pairP4.pt()/(pairSumPt05Ch + pairP4.pt());
        if(verbose || jet_verbose) cout << "pairP4.pt() = " << pairP4.pt() << " pairSumPt05Ch = " << pairSumPt05Ch << " isoCoupleCh = " << isoCoupleCh << endl;
        if(isoCoupleCh < 0.9) {
          cout<<"No isolation cut passed, continue."<<endl;
          continue; 
        }

        nEventsPairIsolationFilter++;


        //PT MAX OF THE JET - FILTER -------------------------------------------------
        if (verbose || jet_verbose) cout<<"Current bestCoupleOfTheEvent_Pt = "<<bestCoupleOfTheJetPt<<endl;
        
        if(pairP4.pt() <= bestCoupleOfTheJetPt) {
          if(verbose || jet_verbose) cout<<"Not passed: pT lower than the current best pair of the event. Continue"<<endl;
          continue; //choose the couple with greatest pt
        }

        //If passed, this is the pair with the largest pT of the event so far
        bestCoupleOfTheJetPt = pairP4.pt();     
        if (verbose || jet_verbose) cout<<"pairP4.pt() = "<<bestCoupleOfTheJetPt << endl; //", isoCoupleCh = " << isoCoupleCh << endl;

        if(verbose || jet_verbose) cout<<"This is the best pair so far!"<<endl<<"-------------------------"<<endl;
        isBestCoupleOfTheEventFound = true;

        //Save if best pair has been found
        bestJet_Index        = jetIndex; //note the position of the chosen jet inside the vector   
        deltaRKChosen        = deltaRK;
        bestJetPhotonInvMass = jetPhotonInvMass;
        _isPhi               = isPhi;
        _isRho               = isRho;
        //bestJetJECunc       = unc;
        firstTrkCharge       = firstTrkCharge;
        secondTrkCharge      = secondTrkCharge;
        bestFirstTrkDxy      = firstTrkDxy;
        bestFirstTrkDz       = firstTrkDz;
        bestSecondTrkDxy     = secondTrkDxy;
        bestSecondTrkDz      = secondTrkDz;
        bestFirstTrkDxyErr   = firstTrkDxyErr;
        bestFirstTrkDzErr    = firstTrkDzErr;
        bestSecondTrkDxyErr  = secondTrkDxyErr;
        bestSecondTrkDzErr   = secondTrkDzErr;
        

        bestFirstTrkP4       = firstTrkP4; 
        bestSecondTrkP4      = secondTrkP4;
        bestPairP4           = pairP4;
        bestPairIsoCh        = isoCoupleCh;
        bestPairSumPt05Ch    = pairSumPt05Ch;
          
        if(verbose) cout << "2nd loop ends" << endl;
      } //2ND LOOP ENDS
      if(verbose) cout << "1st loop ends" << endl;
    } //1ST LOOP ENDS

    if(jet->pt() < 25.) {
      if(verbose) cout << "jet->pt() < 25., continue" << endl;
      continue;
      }
    nJets25++;
    if(jet->pt() < 30.) {
      if(verbose) cout << "jet->pt() < 30., continue" << endl;
      continue;
      }  
    nJets30++;
    if(verbose) cout << "jet loop ends" << endl;
  } //JET LOOP END
  } //isJet BOOL END


  //if(jet_verbose) cout << "n_jet = " << nJets25 << endl;


  else{
    //FIRST TRACK LOOP STARTS
for(auto cand1 = PFCandidates->begin(); cand1 != PFCandidates->end(); ++cand1){
  if(verbose) cout << "1st particle loop starts" << endl;

  firstTrkCharge = cand1->charge();  //take firstCand charge
    firstTrkPt    = cand1->pt();  //take firstCand pt
    firstTrkEta   = cand1->eta(); //take firstCand eta
    firstTrkPhi   = cand1->phi(); //take firstCand phi

    if(cand1->bestTrack() == NULL) {
       if(verbose) cout << "if bestTrack null continue" << endl;
       continue;//loop only over charged daughters
    }

    firstTrkDxy    = cand1->bestTrack()->dxy((&slimmedPV->at(0))->position()); //take firstCand dxy
    firstTrkDxyErr = cand1->bestTrack()->dxyError();        
    firstTrkDz     = cand1->bestTrack()->dz((&slimmedPV->at(0))->position()); //take firstCand dz
    firstTrkDzErr  = cand1->bestTrack()->dzError();

    if(firstTrkCharge == 0 || abs(firstTrkDxy) >= 0.2 || abs(firstTrkDz) >= 0.5 || !(cand1->bestTrack()->quality(reco::Track::highPurity))) {
          if(verbose) cout << "if 1st trk...... continue" << endl;
          continue;
    }

    if(firstTrkPt < 1.){
       if(verbose) cout << "if 1st trk pt <1 continue" << endl;
       continue; //firstCand filter if pT < 1 GeV
    }


    //SECOND TRACK LOOP STARTS----------------------------------------------------------------------
    for(auto cand2 = PFCandidates->begin(); cand2 != PFCandidates->end(); ++cand2){
      if(verbose) cout << "2nd particle loop starts" << endl;

      secondTrkCharge = cand2->charge();
        secondTrkPt     = cand2->pt();
        secondTrkEta    = cand2->eta();
        secondTrkPhi    = cand2->phi();

        if(cand2->bestTrack() == NULL) {
          if(verbose) cout << "if bestTrack null continue" << endl;
          continue;
        }

        secondTrkDxy    = cand2->bestTrack()->dxy((&slimmedPV->at(0))->position());
        secondTrkDxyErr = cand2->bestTrack()->dxyError();          
        secondTrkDz     = cand2->bestTrack()->dz((&slimmedPV->at(0))->position());
        secondTrkDzErr  = cand2->bestTrack()->dzError();

        if(secondTrkCharge == 0 || abs(secondTrkDxy) >= 0.2 || abs(secondTrkDz) >= 0.5 || !(cand2->bestTrack()->quality(reco::Track::highPurity))) {
          if(verbose) cout << "if 2nd trk...... continue" << endl;
          continue;
        }

        //TRKs PT CUT --------------------------------------------------------------------------
        if(secondTrkPt < 1.) {
          if(verbose) cout << "if second trk pt<1 continue" << endl;
          continue; //firstCand filter if pT < 1 GeV
        }

        if(firstTrkPt < 10. && secondTrkPt < 10.) {
          if(verbose) cout << "if firstTrkPt < 10. && secondTrkPt < 10 continue" << endl;
          continue;  //filter if both cand pT are < 10GeV
        }

        //DITRK DELTA R CUT --------------------------------------------------------------------------
        float deltaEta= firstTrkEta - secondTrkEta;

        float deltaPhi = fabs(firstTrkPhi - secondTrkPhi);  //phi folding 
        if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

        deltaRK= sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
        if(deltaRK > 0.07) {
          if(verbose) cout << "deltaRK<0.07 continue" << endl;
          continue;
        }

        //OPPOSITE CHARGE - FILTER ---------------------------------------------------------------
        if(firstTrkCharge * secondTrkCharge >= 0) {
          if(verbose) cout << "firstTrkCharge * secondTrkCharge >= 0" << endl;
          continue; //choose only opposite charges
        }

        //QUADRIMOMENTUM CALCULATION ------------------------------------------------------------
        firstTrkP4  = cand1->p4(); //take quadrimomentum
        secondTrkP4 = cand2->p4();

        firstTrkPx  = firstTrkP4.px(); //take px, py, pz of the first candidate
        firstTrkPy  = firstTrkP4.py();
        firstTrkPz  = firstTrkP4.pz();
        secondTrkPx = secondTrkP4.px(); //take px, py, pz of the second candidate
        secondTrkPy = secondTrkP4.py();
        secondTrkPz = secondTrkP4.pz();

        //PIONS OR KAONS HYPOTHESIS -----------------------------------------------------------------------------------------------------------------------------------------------            
        firstTrkEnergyK   = sqrt(firstTrkPx  * firstTrkPx  + firstTrkPy  * firstTrkPy  + firstTrkPz  * firstTrkPz  + kMass  * kMass ); //Kaon hypothesis energy recalculation
        secondTrkEnergyK  = sqrt(secondTrkPx * secondTrkPx + secondTrkPy * secondTrkPy + secondTrkPz * secondTrkPz + kMass  * kMass ); //Kaon hypothesis energy recalculation
        firstTrkEnergyPi  = sqrt(firstTrkPx  * firstTrkPx  + firstTrkPy  * firstTrkPy  + firstTrkPz  * firstTrkPz  + PiMass * PiMass); //Pion hypothesis energy recalculation
        secondTrkEnergyPi = sqrt(secondTrkPx * secondTrkPx + secondTrkPy * secondTrkPy + secondTrkPz * secondTrkPz + PiMass * PiMass); //Pion hypothesis energy recalculation
        
        
        //if (verbose) {
        //  cout<<"firstTrkEnergyK   = "<<firstTrkEnergyK<<endl;
          //cout<<"firstTrkEnergyPi  = "<<firstTrkEnergyPi<<endl;
          //cout<<"secondTrkEnergyK  = "<<secondTrkEnergyK<<endl;
          //cout<<"secondTrkEnergyPi = "<<secondTrkEnergyPi<<endl;
        //}
        

        
        firstTrkP4K   = firstTrkP4.SetE(firstTrkEnergyK); //Kaon hypothesis quadrimomentum correction
        secondTrkP4K  = secondTrkP4.SetE(secondTrkEnergyK); //Kaon hypothesis quadrimomentum correction
        firstTrkP4Pi  = firstTrkP4.SetE(firstTrkEnergyPi); //Kaon hypothesis quadrimomentum correction
        secondTrkP4Pi = secondTrkP4.SetE(secondTrkEnergyPi); //Kaon hypothesis quadrimomentum correction

        pairP4K  = firstTrkP4K  + secondTrkP4K; //calculation of the couple-quadrimomentum after the correction
        PairP4Pi = firstTrkP4Pi + secondTrkP4Pi; //calculation of the couple-quadrimomentum after the correction
        
        if (verbose) {
          cout<<"KK pT = "<<pairP4K.pt()<<endl;
          cout<<"PiPi pT = "<<PairP4Pi.pt()<<endl;
        }

        //DITRK PT CUT -------------------------------------------------------------------------
        if(pairP4K.pt() < 38.) {
          if(verbose) cout<<"couplePt cut NOT passed, continue"<<endl;
          continue;
        }

        //MESON INV MASS CUT -------------------------------------------------------------------------
        isPhi = false;
        isRho = false;

        PhiMass = (pairP4K).M(); //calculate inv mass of the Phi candidate  
        if (verbose) cout<<"mKK (before the meson mass selection) =  "<<PhiMass<<endl;
        if(PhiMass > 1. && PhiMass < 1.05) isPhi = true; //filter on Phi invariant mass  

        RhoMass = (PairP4Pi).M(); //calculate inv mass of the Rho candidate  
        if (verbose) cout<<"mPiPi (before the meson mass selection) =  "<<RhoMass<<endl;
        if(RhoMass > 0.5 && RhoMass < 1.) isRho = true; //filter on Rho invariant mass   

        if (!isPhi && !isRho) {
          if (verbose) cout<<"the pair mass doesn't match any of the two mass hypothesis, continue "<<endl;
          continue; //continue if the pair mass doesn't match any of the two mass hypothesis
          }

        if (isPhi && isRho){ //if both hypothesis are true, mark it as a Phi candidate (this is done because the Phi mass window is tighter)
          isPhi = true;
          isRho = false;
        }

        //update values of quadrimomenta
        if(isPhi){
          firstTrkP4  = firstTrkP4K; 
          secondTrkP4 = secondTrkP4K;
          pairP4      = pairP4K;
        }  
        if(isRho){
          firstTrkP4  = firstTrkP4Pi;
          secondTrkP4 = secondTrkP4Pi;
          pairP4      = PairP4Pi;
        }


    K1SumPt05     = 0.;
    K1SumPt05Ch   = 0.;

    K2SumPt05     = 0.;
    K2SumPt05Ch   = 0.;

    pairSumPt05   = 0.;
    pairSumPt05Ch = 0.;

        // ISOLATION CUT-------------------------------------------------------------------------  
        for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){ //ISOLATION FORLOOP START
            if(verbose) cout << "isolation loop starts" << endl;
            if(debug){
              cout <<endl<<"ISO CALC DETAILS ---------------------"<<endl;
              cout << "pt cand_iso = "<<cand_iso->pt()<<endl;
            }
            if(cand_iso->pt() < 0.5) {
              if(verbose) cout << "cand_iso->pt() < 0.5, continue" << endl;
              continue; //do not consider tracks with pT < 500MeV
            }

          //calculate the deltaR between the track and the first candidate ---------------------------------------
          float deltaPhi_K1 = fabs(firstTrkP4.phi()-cand_iso->phi());  //phi folding 
          if (deltaPhi_K1 > M_PI) deltaPhi_K1 = 2*M_PI - deltaPhi_K1;

          float deltaRK1 = sqrt((firstTrkP4.eta()-cand_iso->eta())*(firstTrkP4.eta()-cand_iso->eta()) + deltaPhi_K1*deltaPhi_K1);
          if (debug) cout << "deltaRK1 = "<<deltaRK1<<endl;
          if(deltaRK1 < 0.0005) {
              if(verbose) cout<<"deltaRK1 < 0.0005, continue" << endl;
              continue; //remove first candidate from the sum
          }

          //calculate the deltaR between the track and the second candidate ---------------------------------------
          float deltaPhi_K2 = fabs(secondTrkP4.phi()-cand_iso->phi());  //phi folding  
          if (deltaPhi_K2 > M_PI) deltaPhi_K2 = 2*M_PI - deltaPhi_K2;

          float deltaRK2 = sqrt((secondTrkP4.eta()-cand_iso->eta())*(secondTrkP4.eta()-cand_iso->eta()) + deltaPhi_K2*deltaPhi_K2);
            if (debug) cout << "deltaRK2 = "<<deltaRK2<<endl;
            if(deltaRK2 < 0.0005) {
              if(verbose) cout << "deltaRK2 < 0.0005, continue" << endl;
            continue; //remove second candidate from the sum
            }

            //calculate the deltaR between the track and the best pair ---------------------------------------
            float deltaPhi_Couple = fabs(pairP4.phi()-cand_iso->phi());  //phi folding  
            if (deltaPhi_Couple > M_PI) deltaPhi_Couple = 2*M_PI - deltaPhi_Couple;

            float deltaR_Couple = sqrt((pairP4.eta()-cand_iso->eta())*(pairP4.eta()-cand_iso->eta()) + deltaPhi_Couple*deltaPhi_Couple);

            //sum pT of the tracks inside a cone of deltaR = 0.3 ---------------------------------------
            if(deltaRK1 <= 0.3) K1SumPt05 += cand_iso->pt();
            if(deltaRK2 <= 0.3) K2SumPt05 += cand_iso->pt();
            if(deltaR_Couple <= 0.3) pairSumPt05 += cand_iso->pt();

            //sum pT of the charged tracks inside a cone of deltaR = 0.3 ---------------------------------------
            if (debug) cout << "Charge = "<< cand_iso->charge()<<endl;
            if(cand_iso->charge() == 0) {
            if(verbose) cout << "cand_iso->charge() == 0, continue" << endl;
            continue;
            }

            if (debug) cout << "dxy = "<<fabs(cand_iso->dxy())<<" and dz = "<< fabs(cand_iso->dz())<<endl;
            if(fabs(cand_iso->dxy()) >= 0.2 || fabs(cand_iso->dz()) >= 0.5) {
            if(verbose) cout << "cand_iso->dxy()) >= 0.2 || cand_iso->dz() >= 0.5, continue" << endl;
            continue; // Requesting charged particles to come from PV
            }
          
            if(deltaRK1 <= 0.3) K1SumPt05Ch += cand_iso->pt();
            if(deltaRK2 <= 0.3) K2SumPt05Ch += cand_iso->pt();
            if (debug) cout <<"deltaR_Couple = "<<deltaR_Couple<<endl;
            if(deltaR_Couple <= 0.3){
              pairSumPt05Ch += cand_iso->pt();
              if (verbose) cout<<"Particle in the cone: SumPt = "<<pairSumPt05Ch<<endl;
            }
            if(verbose) cout << "isolations loop ends" << endl;
        } //ISOLATION FORLOOP END

        float isoCoupleCh = pairP4.pt()/(pairSumPt05Ch + pairP4.pt());
        if(verbose) cout << "pairP4.pt() = " << pairP4.pt() << " pairSumPt05Ch = " << pairSumPt05Ch << " isoCoupleCh = " << isoCoupleCh << endl;
        if(isoCoupleCh < 0.9) {
          cout<<"No isolation cut passed, continue."<<endl;
          continue; 
        }

        nEventsPairIsolationFilter++;

        //PT MAX OF THE JET - FILTER -------------------------------------------------
        if (verbose) cout<<"Current bestCoupleOfTheEvent_Pt = "<<bestCoupleOfTheJetPt<<endl;
        
        if(pairP4.pt() <= bestCoupleOfTheJetPt) {
          if(verbose) cout<<"Not passed: pT lower than the current best pair of the event. Continue"<<endl;
          continue; //choose the couple with greatest pt
        }

        //If passed, this is the pair with the largest pT of the event so far
        bestCoupleOfTheJetPt = pairP4.pt();     
        if (verbose) cout<<"pairP4.pt() = "<<bestCoupleOfTheJetPt << endl; 

        if(verbose) cout<<"This is the best pair so far!"<<endl<<"-------------------------"<<endl;
        isBestCoupleOfTheEventFound = true;


        //Save if best pair has been found  
        deltaRKChosen        = deltaRK;
        _isPhi               = isPhi;
        _isRho               = isRho;
        firstTrkCharge       = firstTrkCharge;
        secondTrkCharge      = secondTrkCharge;
        bestFirstTrkDxy      = firstTrkDxy;
        bestFirstTrkDz       = firstTrkDz;
        bestSecondTrkDxy     = secondTrkDxy;
        bestSecondTrkDz      = secondTrkDz;
        bestFirstTrkDxyErr   = firstTrkDxyErr;
        bestFirstTrkDzErr    = firstTrkDzErr;
        bestSecondTrkDxyErr  = secondTrkDxyErr;
        bestSecondTrkDzErr   = secondTrkDzErr;

        bestFirstTrkP4       = firstTrkP4; 
        bestSecondTrkP4      = secondTrkP4;
        bestPairP4           = pairP4;
        bestPairIsoCh        = isoCoupleCh;
        bestPairSumPt05Ch    = pairSumPt05Ch;

        if(verbose) cout << "2nd loop ends" << endl;
    }//SECOND TRACK LOOP ENDS
    if(verbose) cout << "1st loop ends" << endl;
}//FIRST TRACK LOOP ENDS




  }


  



  if(!isBestCoupleOfTheEventFound) {
    cout<<"No best couple detected for current event, RETURN."<<endl;
    return;
  }
  nEventsBestPairFound++;      
  if(verbose||jet_verbose) cout<<"Bool: nEventsBestPairFound: "<<nEventsBestPairFound<<endl;        

  //DATAMEMBER SAVING
  firstTrkPt   = bestFirstTrkP4.pt();
  firstTrkEta  = bestFirstTrkP4.eta();
  firstTrkPhi  = bestFirstTrkP4.phi();
  secondTrkPt  = bestSecondTrkP4.pt();       
  secondTrkEta = bestSecondTrkP4.eta();
  secondTrkPhi = bestSecondTrkP4.phi();
  bestPairPt   = bestPairP4.pt();
  bestPairEta  = bestPairP4.eta();
  bestPairPhi  = bestPairP4.phi();

  if(isJet){
  bestJetInvMass                  = slimmedJets->at(bestJet_Index).mass();
  bestJetPt                       = slimmedJets->at(bestJet_Index).pt();
  bestJetEta                      = slimmedJets->at(bestJet_Index).eta();
  bestJetPhi                      = slimmedJets->at(bestJet_Index).phi();
  bestJetnDaughters               = slimmedJets->at(bestJet_Index).numberOfDaughters();
  bestJetChargedEmEnergy          = slimmedJets->at(bestJet_Index).chargedEmEnergy();
  bestJetNeutralEmEnergy          = slimmedJets->at(bestJet_Index).neutralEmEnergy();
  bestJetChargedHadEnergy         = slimmedJets->at(bestJet_Index).chargedHadronEnergy();
  bestJetNeutralHadEnergy         = slimmedJets->at(bestJet_Index).neutralHadronEnergy();
  bestJetChargedEmEnergyFraction  = slimmedJets->at(bestJet_Index).chargedEmEnergyFraction();
  bestJetNeutralEmEnergyFraction  = slimmedJets->at(bestJet_Index).neutralEmEnergyFraction();
  bestJetChargedHadEnergyFraction = slimmedJets->at(bestJet_Index).chargedHadronEnergyFraction();
  bestJetNeutralHadEnergyFraction = slimmedJets->at(bestJet_Index).neutralHadronEnergyFraction();
  bestJetChargedHadMultiplicity   = slimmedJets->at(bestJet_Index).chargedHadronMultiplicity();
  }

  //MESON MASS CALCULATION
  mesonMass = (bestFirstTrkP4 + bestSecondTrkP4).M();

  //Z INV MASS CALCULATION
  ZMassFrom2KPhoton = (bestFirstTrkP4 + bestSecondTrkP4 + ph_p4).M(); //calculate inv mass of the Z candidate
    

  //CANDIDATES SORTING
  if(firstTrkPt < secondTrkPt){  //swap-values loop, in order to fill the tree with the candidate with max pt of the couple in firstCand branches                                 //and one with min pt in secondCand branches
    float a,b,c,d,e;
    a = firstTrkPt;
    b = firstTrkEta;
    c = firstTrkPhi;
    d = firstTrkEnergy;
    e = firstTrkCharge;
    firstTrkPt     = secondTrkPt;
    firstTrkEta    = secondTrkEta;
    firstTrkPhi    = secondTrkPhi;
    firstTrkEnergy = secondTrkEnergy;
    firstTrkCharge = secondTrkCharge;
    secondTrkPt     = a;
    secondTrkEta    = b;
    secondTrkPhi    = c;
    secondTrkEnergy = d;
    secondTrkCharge = e;
  }

  //CUTS ON CANDIDATES PT
  if(firstTrkPt < 20. || secondTrkPt < 5.) {
    cout<<"Final cut on candidates pT not passed, RETURN."<<endl;
    return;
  }
  nEventsTrkPtFilter++;

  //ISOLATION DATAMEMBER FOR TREE FILLING 
  isoK1     = firstTrkPt/(K1SumPt05 + firstTrkPt);
  isoK2     = secondTrkPt/(K2SumPt05 + secondTrkPt);
  isoPair   = bestPairPt/(pairSumPt05 + bestPairPt);
  isoK1Ch   = firstTrkPt/(K1SumPt05Ch + firstTrkPt);
  isoK2Ch   = secondTrkPt/(K2SumPt05Ch + secondTrkPt);
  isoPairCh = bestPairPt/(pairSumPt05Ch + bestPairPt);
  if(verbose || jet_verbose) cout << "bestPairPt = " << bestPairPt << " bestPairSumPt05Ch = " << bestPairSumPt05Ch << " bestPairIsoCh = " << bestPairIsoCh << endl; 

  //CUT ON PHI ISOLATION
  if(verbose){
    cout<<endl;
    cout<<"###### ISO           = "<<bestPairIsoCh<<endl;
    cout<<"###### isRho         = "<<isRho<<endl;
    cout<<"###### SUM pT        = "<<bestPairSumPt05Ch<<endl;
    cout<<"###### pT leading    = "<<firstTrkPt<<endl;
    cout<<"###### pT subleading = "<<secondTrkPt<<endl;
    cout<<"###### MesonMass     = "<<mesonMass<<endl;
    cout<<"###### ZMass         = "<<ZMassFrom2KPhoton<<endl;
  }



  //nEventsPairIsolationFilter++;



  //*************************************************************///
  //                                                             //
  //----------------------- Access MC Truth ---------------------//
  //                                                             //
  //*************************************************************//

  //In signal, identify if there's a real mu or ele from W
  is_Kplus_matched   = false;
  is_Kminus_matched  = false;
  is_Piplus_matched  = false;
  is_Piminus_matched = false;
  //is_Phi_matched     = false;
  //is_Rho_matched     = false;
  is_photon_matched  = false;
  is_meson_matched   = false;
  is_Higgs_matched   = false; 

  Kminus_phi         = -999.;
  Kplus_phi          = -999.;
  float Piminus_phi  = -999.;
  float Piplus_phi   = -999.;
  Kminus_eta         = -999.;
  Kplus_eta          = -999.;
  float Piminus_eta  = -999.;
  float Piplus_eta   = -999.;
  deltaRKplus        = -999;
  deltaR_wrong       = -999;
  deltaRKminus       = -999.;
  deltaR_Piplus      = -999.;
  deltaR_Piminus     = -999.;
  genPhoton_eT       = -999.;
  genPhoton_eta      = -999.;
  genPhoton_phi      = -999.;
  genMeson_pT        = -999.;
  genMeson_m         = -999.;
  genMeson_eta       = -999.;
  genMeson_phi       = -999.;
  KplusPt            = -999.;
  KminusPt           = -999.;
  Kplus_dxy          = -999.;
  Kplus_dz           = -999.;
  Kminus_dxy         = -999.;
  Kminus_dz          = -999.;

  theta_pol     = -10.;
  TLorentzVector mu[2];
  theta_pol_tree = 0.;

  if(!runningOnData_){
    for(auto gen = prunedGenParticles->begin(); gen != prunedGenParticles->end(); ++gen){
      if( gen->pdgId() == 321  && gen->mother()->pdgId() == 333 && gen->mother()->mother()->pdgId() == 23)  Kplus_phi   = gen->phi(), Kplus_eta   = gen->eta(), KplusPt  = gen->pt();//, Kplus_dz  = gen->dz();//(&slimmedPV->at(0))->position()
      if( gen->pdgId() == -321 && gen->mother()->pdgId() == 333 && gen->mother()->mother()->pdgId() == 23)  Kminus_phi  = gen->phi(), Kminus_eta  = gen->eta(), KminusPt = gen->pt();//, Kminus_dxy = gen->dxy(), Kminus_dz = gen->dz();
      if( gen->pdgId() == 211  && gen->mother()->pdgId() == 113 && gen->mother()->mother()->pdgId() == 23)  Piplus_phi  = gen->phi(), Piplus_eta  = gen->eta();
      if( gen->pdgId() == -211 && gen->mother()->pdgId() == 113 && gen->mother()->mother()->pdgId() == 23)  Piminus_phi = gen->phi(), Piminus_eta = gen->eta();
      if( gen->pdgId() == 333  && gen->mother()->pdgId() == 23) genMeson_pT  = gen->pt(), genMeson_phi  = gen->phi(),  genMeson_eta = gen->eta(), genMeson_m = gen->mass();
      if( gen->pdgId() == 113  && gen->mother()->pdgId() == 23) genMeson_pT  = gen->pt(), genMeson_phi  = gen->phi(),  genMeson_eta = gen->eta(), genMeson_m = gen->mass();
      if( gen->pdgId() == 22   && gen->mother()->pdgId() == 23) genPhoton_eT = gen->pt(), genPhoton_phi = gen->phi(), genPhoton_eta = gen->eta();
    }
  }


  if(!runningOnData_){
    //For the polarization reweighting
    for(auto gen = prunedGenParticles->begin(); gen != prunedGenParticles->end(); ++gen){
      if(gen->pdgId() == 23 && gen->numberOfDaughters() == 2){
        //for each daughter
        for(int i = 0; i < 2; i++){
          //if daughters are not Phi or Rho and gamma, continue
          if( !(gen->daughter(i)->pdgId() == 22 || (gen->daughter(i)->pdgId() == 333 || gen->daughter(i)->pdgId() == 113)) ) continue;
            //if daughter(i) is a Phi or a Rho
            if(gen->daughter(i)->pdgId() == 333 || gen->daughter(i)->pdgId() == 113){
              if(gen->daughter(i)->numberOfDaughters() == 2){      
                //for each Meson daughter
                for(int j = 0; j < 2; j++){
                  //if daughter(j) is a K+
                  if(gen->daughter(i)->daughter(j)->pdgId() == 321){
                    mu[1].SetPxPyPzE(gen->daughter(i)->daughter(j)->px(),gen->daughter(i)->daughter(j)->py(), gen->daughter(i)->daughter(j)->pz(),gen->daughter(i)->daughter(j)->energy());
                    mu[0].SetPxPyPzE(gen->daughter(i)->px(),gen->daughter(i)->py(), gen->daughter(i)->pz(),gen->daughter(i)->energy());
                    TVector3 trackBoost = mu[0].BoostVector();
                    mu[1].Boost(- trackBoost);
                    theta_pol = mu[0].Vect().Angle(mu[1].Vect());
                  }
                  //if daughter(j) is a pi+
                  if(gen->daughter(i)->daughter(j)->pdgId() == 211){
                    mu[1].SetPxPyPzE(gen->daughter(i)->daughter(j)->px(),gen->daughter(i)->daughter(j)->py(), gen->daughter(i)->daughter(j)->pz(),gen->daughter(i)->daughter(j)->energy());
                    mu[0].SetPxPyPzE(gen->daughter(i)->px(),gen->daughter(i)->py(), gen->daughter(i)->pz(),gen->daughter(i)->energy());
                    TVector3 trackBoost = mu[0].BoostVector();
                    mu[1].Boost(- trackBoost);
                    theta_pol = mu[0].Vect().Angle(mu[1].Vect());
                  }
               }
              }
            }
        }
      }
    }    
  }

  theta_pol_tree = theta_pol;

  
  //MC TRUTH CHECK
  if(!runningOnData_){ //ONLY FOR MC START
   
    //photon matching -----------------------------------------
    float deltaPhiPhoton = fabs(photonPhi - genPhoton_phi);
    if (deltaPhiPhoton > M_PI) deltaPhiPhoton = 2*M_PI - deltaPhiPhoton;

    float deltaR_photonGenVsReco = sqrt((photonEta - genPhoton_eta) * (photonEta - genPhoton_eta) + deltaPhiPhoton * deltaPhiPhoton);
    //float photonRelPt = 
    if (deltaR_photonGenVsReco < 0.2) is_photon_matched = true;

    //meson matching -----------------------------------------
    float deltaPhiMeson = fabs(bestPairPhi - genMeson_phi);
    if (deltaPhiMeson > M_PI) deltaPhiMeson = 2*M_PI - deltaPhiMeson;

    float deltaR_mesonGenVsReco = sqrt((bestPairEta - genMeson_eta) * (bestPairEta - genMeson_eta) + deltaPhiMeson * deltaPhiMeson);
    if (deltaR_mesonGenVsReco < 0.3) is_meson_matched = true;    

    //meson pT matching -----------------------------------------
    if (bestPairPt < genMeson_pT - 4. || bestPairPt > genMeson_pT + 4.) nEventsMesonPtNotMatched ++; 

    //Higgs matching -----------------------------------------    
    if(is_photon_matched && is_meson_matched){

      if(verbose) cout<<endl<<"**************** Z FOUND ******************"<<endl;
      if(verbose) cout<<"Z deltaR = "<<deltaR<<endl;
      nEventsZMatched++;
      is_Higgs_matched=true;
    }
    
    else 
    {
      nEventsZNotMatched++;
      if(verbose) cout<<endl<<"THAT'S NOT A Z!"<<endl;
    
    }

    //if is PhiGamma event
    if (_isPhi){

    //First cand positive and second cand negative
    if (firstTrkCharge > 0){ 
      //phi angle folding K plus
      float deltaPhi_Kplus = fabs(firstTrkPhi - Kplus_phi);
      float deltaPhi_wrong = fabs(secondTrkPhi - Kplus_phi);
      if (deltaPhi_Kplus > M_PI) deltaPhi_Kplus = 2*M_PI - deltaPhi_Kplus;
      if (deltaPhi_wrong > M_PI) deltaPhi_wrong = 2*M_PI - deltaPhi_wrong;

      //deltaR K plus
      deltaRKplus = sqrt((firstTrkEta - Kplus_eta) * (firstTrkEta - Kplus_eta) + deltaPhi_Kplus * deltaPhi_Kplus);
      deltaR_wrong = sqrt((secondTrkEta - Kplus_eta) * (secondTrkEta - Kplus_eta) + deltaPhi_wrong * deltaPhi_wrong);
      
      //phi angle folding K minus
      float deltaPhi_Kminus = fabs(secondTrkPhi - Kminus_phi);
      if (deltaPhi_Kminus > M_PI) deltaPhi_Kminus = 2*M_PI - deltaPhi_Kminus;

      //deltaR K minus
      deltaRKminus = sqrt((secondTrkEta - Kminus_eta) * (secondTrkEta - Kminus_eta) + deltaPhi_Kminus * deltaPhi_Kminus);
      cout<<endl;
      //if (firstTrkPt < 0.95*KplusPt && firstTrkPt > 1.05*KplusPt) cout<<"firstCand pT not matched"<<endl;

    }

    else{ //Second cand positive and first cand negative
      
      //phi angle folding K plus
      float deltaPhi_Kplus = fabs(secondTrkPhi - Kplus_phi);
      if (deltaPhi_Kplus > M_PI) deltaPhi_Kplus = 2*M_PI - deltaPhi_Kplus;

      //deltaR K plus
      deltaRKplus = sqrt((secondTrkEta - Kplus_eta) * (secondTrkEta - Kplus_eta) + deltaPhi_Kplus * deltaPhi_Kplus);

      cout<<endl;
      //if (firstTrkPt < 0.95*KminusPt && firstTrkPt > 1.05*KminusPt) cout<<"firstTrkPt pT not matched"<<endl;

    }

    //phi angle folding K minus
    float deltaPhi_Kminus = fabs(firstTrkPhi - Kminus_phi);
    if (deltaPhi_Kminus > M_PI) deltaPhi_Kminus = 2*M_PI - deltaPhi_Kminus;

    //deltaR K minus
    deltaRKminus = sqrt((firstTrkEta - Kminus_eta) * (firstTrkEta - Kminus_eta) + deltaPhi_Kminus * deltaPhi_Kminus);


    } //if isPhi END
  
    else{ //RhoGamma event

      //First cand positive and second cand negative
      if (firstTrkCharge > 0){ 
        //phi angle folding Pi plus
        float deltaPhi_Piplus = fabs(firstTrkPhi - Piplus_phi);
        if (deltaPhi_Piplus > M_PI) deltaPhi_Piplus = 2*M_PI - deltaPhi_Piplus;

        //deltaR Pi plus
        deltaR_Piplus = sqrt((firstTrkEta - Piplus_eta) * (firstTrkEta - Piplus_eta) + deltaPhi_Piplus * deltaPhi_Piplus);

        //phi angle folding Pi minus
        float deltaPhi_Piminus = fabs(secondTrkPhi - Piminus_phi);
        if (deltaPhi_Piminus > M_PI) deltaPhi_Piminus = 2*M_PI - deltaPhi_Piminus;

        //deltaR Pi minus
        deltaR_Piminus = sqrt((secondTrkEta - Piminus_eta) * (secondTrkEta - Piminus_eta) + deltaPhi_Piminus * deltaPhi_Piminus);
        
      }

      else{ //Second cand positive and first cand negative
      
      //phi angle folding Pi plus
      float deltaPhi_Piplus = fabs(secondTrkPhi - Piplus_phi);
      if (deltaPhi_Piplus > M_PI) deltaPhi_Piplus = 2*M_PI - deltaPhi_Piplus;

      //deltaR Pi plus
      deltaR_Piplus = sqrt((secondTrkEta - Piplus_eta) * (secondTrkEta - Piplus_eta) + deltaPhi_Piplus * deltaPhi_Piplus);
      
      
      //phi angle folding Pi minus
      float deltaPhi_Piminus = fabs(firstTrkPhi - Piminus_phi);
      if (deltaPhi_Piminus > M_PI) deltaPhi_Piminus = 2*M_PI - deltaPhi_Piminus;

      //deltaR Pi minus
      deltaR_Piminus = sqrt((firstTrkEta - Piminus_eta) * (firstTrkEta - Piminus_eta) + deltaPhi_Piminus * deltaPhi_Piminus);

      }
    } // if isRho END



    //phi angle folding
    //float deltaPhi_Kpm = fabs(Kplus_phi - Kminus_phi);
    //if (deltaPhi_Kpm > M_PI) deltaPhi_Kpm = 2*M_PI - deltaPhi_Kpm;

    //float deltaRKpm = sqrt((Kplus_eta - Kminus_eta) * (Kplus_eta - Kminus_eta) + deltaPhi_Kpm * deltaPhi_Kpm);

    //cout<<"deltaRKpm   = "<<deltaRKpm<<endl;



    //some prints
    
    if(verbose || jet_verbose){
      cout<<"Photon eT = "<<photonEt<<endl;
      //cout<<"ph_en_sigmaUP = "<< ph_en_sigmaUP<<endl;
      //cout<<"ph_en_sigmaDW = "<< ph_en_sigmaDW<<endl;
      //cout<<"ph_en_scaleUP = "<<ph_en_scaleUP<<endl;
      //cout<<"ph_en_scaleDW = "<<ph_en_scaleDW<<endl;
      cout<<"n Jets = "<<nJets25<<endl;
      if(isJet) cout<<"Jet + photon inv. mass = "<<bestJetPhotonInvMass<<endl;
      if(isJet) cout<<"n. of daughters: "<<bestJetnDaughters<<endl;
      cout<<"Best couple pT = "<<firstTrkPt + secondTrkPt<<endl;
      cout<<"firstTrkPt   = "<<firstTrkPt<<endl;
      cout<<"secondTrkPt  = "<<secondTrkPt<<endl;
      cout<<"genMeson_pT    = "<<genMeson_pT<<endl;
      cout<<"genKplus pT    = "<<KplusPt<<endl;
      cout<<"genKminus pT   = "<<KminusPt<<endl;
      cout<<"trk1 dxy       = "<<bestFirstTrkDxy<<endl;
      cout<<"trk2 dxy       = "<<bestSecondTrkDxy<<endl;
      cout<<"trk1 dxyErr    = "<<bestFirstTrkDxyErr<<endl;
      cout<<"trk2 dxyErr    = "<<bestSecondTrkDxyErr<<endl;
      cout<<"deltaDxy       = "<<abs(bestFirstTrkDxy - bestSecondTrkDxy)<<endl;
      cout<<"trk1 dz        = "<<bestFirstTrkDz<<endl;
      cout<<"trk2 dz        = "<<bestSecondTrkDz<<endl;
      cout<<"trk1 dzErr     = "<<bestFirstTrkDzErr<<endl;
      cout<<"trk2 dzErr     = "<<bestSecondTrkDzErr<<endl;
      cout<<"deltaDz        = "<<abs(bestFirstTrkDz - bestSecondTrkDz)<<endl;
      cout<<"Kplus  dxy     = "<<Kplus_dxy<<endl;
      cout<<"Kplus  dz      = "<<Kplus_dz<<endl;
      cout<<"Kminus dxy     = "<<Kminus_dxy<<endl;
      cout<<"Kminus dz      = "<<Kminus_dz<<endl;
      cout<<"Best couple DeltaR = "<<deltaRKChosen<<endl;
      cout<<"Meson candidate inv. mass  = "<<mesonMass<<endl;
      cout<<"isPhi = "<<isPhi<<" and isRho = "<<isRho<<endl;
      cout<<"Z inv. mass = "<<ZMassFrom2KPhoton<<endl;
      cout<<"--------------------------------------------------"<<endl;
      cout<<"MC Z found = "<<nEventsZMatched<<",   Z NOT matched = "<<nEventsZNotMatched<<",   mesonPt not matched = "<<nEventsMesonPtNotMatched<<endl;
      cout<<"--------------------------------------------------"<<endl<<endl;
    }
     
  }  //ONLY FOR MC END 
 
  else {//ONLY FOR DATA
 
    cout<<"CANDIDATE Z FOUND IN DATA: EVENT RECORDED!"<<endl;
    if(debug){
      cout<<"Photon eT = "<<photonEt<<endl;
    }
 }

  //cout<<endl<<"Event n = "<<eventNumber<<endl;
  mytree->Fill();

}



//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void ZMesonGamma::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isTwoProngTrigger",&isTwoProngTrigger);

//Save run number info when running on data
  if(runningOnData_){
    mytree->Branch("runNumber",&runNumber);
    mytree->Branch("eventNumber",&eventNumber);
  }
  else{
    mytree->Branch("eventNumber",&eventNumber);
  }

  mytree->Branch("nMuons10",&nMuons10);
  mytree->Branch("nMuons20",&nMuons20);
  mytree->Branch("nElectrons10",&nElectrons10);
  mytree->Branch("nElectrons20",&nElectrons20);
  mytree->Branch("nPhotons38WP80",&nPhotons38WP80);
  mytree->Branch("nPhotons20WP90",&nPhotons20WP90);
  mytree->Branch("nPhotonsChosen",&nPhotonsChosen);
  mytree->Branch("nJets30",&nJets30);
  mytree->Branch("nJets25",&nJets25);
  mytree->Branch("metPt",&metPt);
  mytree->Branch("metpuppiPt",&metpuppiPt);

  mytree->Branch("photon_eT",&photonEt);
  mytree->Branch("photon_eta",&photonEta);
  mytree->Branch("photon_etaSC",&photonEtaSC);
  mytree->Branch("photon_phi",&photonPhi);
  //mytree->Branch("photon_iso_ChargedHadron",&photonIsoChargedHadron);
  //mytree->Branch("photon_iso_NeutralHadron",&photonIsoNeutralHadron);
  //mytree->Branch("photon_iso_Photon",&photonIsoPhoton);
  //mytree->Branch("photon_iso_eArho",&photonIsoEArho);
  mytree->Branch("photonRegressionError",&photonRegressionError);

  mytree->Branch("bestJet_pT",&bestJetPt);
  mytree->Branch("bestJet_eta",&bestJetEta);
  mytree->Branch("bestJet_phi",&bestJetPhi);
  mytree->Branch("bestJet_nDaughters",&bestJetnDaughters);
  mytree->Branch("bestJet_chargedEmEnergy",&bestJetChargedEmEnergy);
  mytree->Branch("bestJet_neutralEmEnergy",&bestJetNeutralEmEnergy);
  mytree->Branch("bestJet_chargedHadEnergy",&bestJetChargedHadEnergy);
  mytree->Branch("bestJet_neutralHadEnergy",&bestJetNeutralHadEnergy);
  mytree->Branch("bestJet_chargedEmEnergyFraction",&bestJetChargedEmEnergyFraction);
  mytree->Branch("bestJet_neutralEmEnergyFraction",&bestJetNeutralEmEnergyFraction);
  mytree->Branch("bestJet_chargedHadEnergyFraction",&bestJetChargedHadEnergyFraction);
  mytree->Branch("bestJet_neutralHadEnergyFraction",&bestJetNeutralHadEnergyFraction);
  mytree->Branch("bestJet_invMass",&bestJetInvMass);
  mytree->Branch("bestJet_Photon_invMass",&bestJetPhotonInvMass);
  //mytree->Branch("bestJet_JECunc",&bestJetJECunc);
  

  mytree->Branch("firstTrkCharge",&firstTrkCharge);
  mytree->Branch("firstTrkPt",&firstTrkPt);
  mytree->Branch("firstTrkEta",&firstTrkEta);
  mytree->Branch("firstTrkPhi",&firstTrkPhi);
  mytree->Branch("firstTrkDxy",&bestFirstTrkDxy);
  mytree->Branch("firstTrkDz",&bestFirstTrkDz);
  mytree->Branch("firstTrkDxyErr",&bestFirstTrkDxyErr);
  mytree->Branch("firstTrkDzErr",&bestFirstTrkDzErr);
  mytree->Branch("secondTrkCharge",&secondTrkCharge);
  mytree->Branch("secondTrkPt",&secondTrkPt);
  mytree->Branch("secondTrkEta",&secondTrkEta);
  mytree->Branch("secondTrkPhi",&secondTrkPhi);
  mytree->Branch("secondTrkDxy",&bestSecondTrkDxy);
  mytree->Branch("secondTrkDz",&bestSecondTrkDz);
  mytree->Branch("secondTrkDxyErr",&bestSecondTrkDxyErr);
  mytree->Branch("secondTrkDzErr",&bestSecondTrkDzErr);
  mytree->Branch("bestCouplePt",&bestPairPt);
  mytree->Branch("bestCoupleEta",&bestPairEta);
  mytree->Branch("bestCouplePhi",&bestPairPhi);
  mytree->Branch("isPhi",&_isPhi);
  mytree->Branch("isRho",&_isRho);

  mytree->Branch("firstTrkEnergy",&firstTrkEnergy);
  mytree->Branch("secondTrkEnergy",&secondTrkEnergy);

  mytree->Branch("MesonMass",&mesonMass);
  mytree->Branch("ZMassFrom2KPhoton",&ZMassFrom2KPhoton);

  //mytree->Branch("K1SumPt05",&K1SumPt05);
  //mytree->Branch("K1SumPt05Ch",&K1SumPt05Ch);
  //mytree->Branch("K2SumPt05",&K2SumPt05);
  //mytree->Branch("K2SumPt05Ch",&K2SumPt05Ch);
  mytree->Branch("pairSumPt05",&pairSumPt05);
  mytree->Branch("pairSumPt05Ch",&pairSumPt05Ch);

  mytree->Branch("iso_K1",&isoK1);
  mytree->Branch("iso_K1_ch",&isoK1Ch);
  mytree->Branch("iso_K2",&isoK2);
  mytree->Branch("iso_K2_ch",&isoK2Ch);
  mytree->Branch("iso_couple",&isoPair);
  mytree->Branch("bestIso_couple_ch",&bestPairIsoCh);


  //Save MC info
  if(!runningOnData_){ //NO INFO FOR DATA
    mytree->Branch("PU_Weight",&PU_Weight);
    mytree->Branch("MC_Weight",&MC_Weight);
    mytree->Branch("minPDFWeight",&minPDFWeight);
    mytree->Branch("maxPDFWeight",&maxPDFWeight);
    mytree->Branch("minQCDWeight",&minQCDWeight);
    mytree->Branch("maxQCDWeight",&maxQCDWeight);
    mytree->Branch("isHiggsMatched",&is_Higgs_matched);
    //mytree->Branch("isKplusMatched",&is_Kplus_matched);
    //mytree->Branch("isKminusMatched",&is_Kminus_matched);
    //mytree->Branch("isPiplusMatched",&is_Piplus_matched);
    //mytree->Branch("isPiminusMatched",&is_Piminus_matched);
    //mytree->Branch("isPhiFromH",&is_Phi_fromH);
    //mytree->Branch("isRhofromH",&is_Rho_fromH);
    //mytree->Branch("isPhotonFromH",&is_Photon_fromH);
    //mytree->Branch("isPhotonTrue",&is_photon_a_photon);
    mytree->Branch("isPhotonMatched",&is_photon_matched);
    mytree->Branch("genPhoton_eT",&genPhoton_eT);
    mytree->Branch("isMesonMatched",&is_meson_matched);
    mytree->Branch("genMeson_pT",&genMeson_pT);
    mytree->Branch("genMeson_m",&genMeson_m);
    mytree->Branch("KplusPt",&KplusPt);
    mytree->Branch("KminusPt",&KminusPt);
    mytree->Branch("Kminus_eta",&Kminus_eta);
    mytree->Branch("Kplus_eta",&Kplus_eta);
    mytree->Branch("Kminus_phi",&Kminus_phi);
    mytree->Branch("Kplus_phi",&Kplus_phi);

    mytree->Branch("deltaRKplus",&deltaRKplus);
    mytree->Branch("deltaR_wrong",&deltaR_wrong);
    mytree->Branch("deltaRKminus",&deltaRKminus);
    mytree->Branch("deltaR_Piplus",&deltaR_Piplus);
    mytree->Branch("deltaR_Piminus",&deltaR_Piminus);
    mytree->Branch("theta_polarization",&theta_pol_tree);

  }

}

void ZMesonGamma::beginJob(){
  //Flag for PileUp reweighting
  if (!runningOnData_){ // PU reweighting for 2017
   //Lumiweights_ = edm::LumiReWeighting("MCpileUp_2018_25ns_UltraLegacy_PoissonOOTPU.root", "MyDataPileupHistogram.root", "pileup", "pileup");
  }
}



//*************************************************************//
//                                                             //
//------------------- Fill event loss histos ------------------//
//                                                             //
//*************************************************************//

void ZMesonGamma::endJob() {
  hEvents->Fill(0.5,nEventsProcessed);
  hEvents->Fill(1.5,nEventsTriggered);
  hEvents->Fill(2.5,nEventsIsPhoton);
  hEvents->Fill(3.5,nEventsPairIsolationFilter);
  hEvents->Fill(4.5,nEventsBestPairFound);  
  hEvents->Fill(5.5,nEventsTrkPtFilter);  
  
  hEvents->GetXaxis()->SetBinLabel(1,"processed");
  hEvents->GetXaxis()->SetBinLabel(2,"triggered");
  hEvents->GetXaxis()->SetBinLabel(3,"best photon");
  hEvents->GetXaxis()->SetBinLabel(4,"trks iso");
  hEvents->GetXaxis()->SetBinLabel(5,"best pair");
  hEvents->GetXaxis()->SetBinLabel(6,"trks pT");
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZMesonGamma);