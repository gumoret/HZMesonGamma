
//---------- class declaration----------

class ZMesonGamma : public edm::stream::EDAnalyzer<> {//////////
public:
  explicit ZMesonGamma(const edm::ParameterSet&);
  ~ZMesonGamma();

private:
  virtual void beginJob();/////// override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();//////// override;
 
 
  bool runningOnData_;
  bool verboseIdFlag_;
  const edm::InputTag packedPFCandidates_;
  const edm::InputTag slimmedMuons_; 
  const edm::InputTag prunedGenParticles_;
  const edm::InputTag genParticles_;
  const edm::InputTag slimmedPhotons_;
  const edm::InputTag slimmedElectrons_;
  const edm::InputTag slimmedJets_;
  const edm::InputTag slimmedMETs_;
  const edm::InputTag slimmedMETsPuppi_;
  const edm::InputTag pvCollection_;  
  const edm::InputTag bsCollection_;  
  const edm::InputTag PileupSrc_;
  const edm::InputTag GenInfo_;

  edm::LumiReWeighting Lumiweights_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  TH1F* hEvents;

  TH1F* hPileup;

  //debug

  bool debug;
  bool verbose;

  //Counters
  int nPV;
  int nMuons10;
  int nMuons20;
  int nElectrons10;
  int nElectrons20;
  int nPhotonsChosen;
  int nPhotons20WP90;
  int nPhotons38WP80;
  int nJets30;
  int nJets25;
  int nEventsTriggered;
  int nEventsProcessed;
  int nEventsIsTwoKaons;
  int nEventsIsPhoton;
  int nEventsZMatched; 
  int nEventsZNotMatched; 
  int nEventsZMassMatched; 
  int nEventsZMassNotMatched; 
  int nEventsMesonPtNotMatched;
  int nEventsBestPairFound;
  int nEventsTrkPtFilter;
  int nEventsPairIsolationFilter;
  
  //TTree and TTree variables
  TTree *mytree;

  int runNumber;
  int eventNumber;

  float photonEt;
  float photonEta;
  float photonEtaSC;
  float photonPhi;
  float photonIsoChargedHadron;
  float photonIsoNeutralHadron;
  float photonIsoPhoton;
  float photonIsoEArho;
  bool isPhotonWP90;
  float photonEtMax;
  float photonRegressionError;

  float firstTrkPx;
  float firstTrkPy;
  float firstTrkPz;
  float bestFirstTrkDxy;
  float bestFirstTrkDxyErr;
  float bestFirstTrkDz;
  float bestFirstTrkDzErr;
  float secondTrkPx;
  float secondTrkPy;
  float secondTrkPz;
  float bestSecondTrkDxy;
  float bestSecondTrkDxyErr;
  float bestSecondTrkDz;
  float bestSecondTrkDzErr;
  float firstTrkEnergy;
  float secondTrkEnergy;
  float firstTrkEnergyK;
  float secondTrkEnergyK;
  float firstTrkEnergyPi;
  float secondTrkEnergyPi;
  float firstTrkPt;
  float firstTrkEta;
  float firstTrkPhi;
  float firstTrkCharge;
  float secondTrkPt;
  float secondTrkEta;
  float secondTrkPhi;
  float secondTrkCharge;
  float bestPairPt;
  float bestPairEta;
  float bestPairPhi;
  float minPDFWeight;
  float maxPDFWeight;
  float minQCDWeight;
  float maxQCDWeight;

  float jetPhotonInvMass;
  float mesonMass;
  float ZMassFrom2KPhoton;
    
  float K1SumPt05;
  float K1SumPt05Ch;
  float K2SumPt05;
  float K2SumPt05Ch;
  float pairSumPt05;
  float pairSumPt05Ch;
  float isoK1;
  float isoK1Ch;
  float isoK2;
  float isoK2Ch;
  float isoPair;
  float isoPairCh;
  float bestPairSumPt05Ch;
  float bestPairIsoCh;
    
  float metPt;
  float metpuppiPt;
  
  bool isTwoProngTrigger;

  //Jet datamember
  
  float jetInvMass;
  float bestJetPt;
  float bestJetEta;
  float bestJetPhi;
  int bestJetnDaughters;
  float bestJetPtMax;
  float bestJetChargedEmEnergy;
  float bestJetNeutralEmEnergy;
  float bestJetChargedHadEnergy;
  float bestJetNeutralHadEnergy;
  float bestJetChargedEmEnergyFraction;
  float bestJetNeutralEmEnergyFraction;
  float bestJetChargedHadEnergyFraction;
  float bestJetNeutralHadEnergyFraction;
  int bestJetChargedHadMultiplicity;
  float bestJetInvMass;
  float bestJetPhotonInvMass;
  float bestJetJECunc;
  
  //MC truth
  float PU_Weight;
  float MC_Weight;
  float deltaRKplus;
  float deltaR_wrong;
  float deltaRKminus;
  float deltaR_Piplus;
  float deltaR_Piminus;
  float genPhoton_eT;
  float genPhoton_eta;
  float genPhoton_phi;
  float genMeson_pT;
  float genMeson_eta;
  float genMeson_phi;
  float genMeson_m;
  float KplusPt;
  float KminusPt;
  float Kplus_dxy;
  float Kplus_dz;
  float Kminus_dxy;
  float Kminus_dz;
  float Kminus_eta;
  float Kplus_eta;
  float Kminus_phi;
  float Kplus_phi;

  bool is_Kplus_matched;
  bool is_Kminus_matched;
  bool is_Piplus_matched;
  bool is_Piminus_matched;
  bool is_Phi_fromH;
  bool is_Rho_fromH;
  bool is_Photon_fromH;
  bool is_photon_a_photon;
  bool is_photon_matched;
  bool is_meson_matched;
  bool is_Higgs_matched;
  bool _isPhi;
  bool _isRho;
  //rho for isolation
  float rho_;

  //for VBF veto
  int nJets20;

  float theta_pol;
  float theta_pol_tree;


  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatesToken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonsToken_;   edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsPuppiToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > offlineSlimmedPrimaryVerticesToken_; 
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpotToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<LHEEventProduct> LHEEventProduct_; //LHE reader
  //edm::EDGetTokenT<std::vector<pat::PackedGenParticle> > packedGenParticlesToken_; //GENPART
  edm::EDGetTokenT<std::vector<reco::GenParticle>> prunedGenParticlesToken_;
  
  //edm::EDGetTokenT<reco::JetCorrector> jetCorrectorToken_;
  //edm::EDGetTokenT<JetCorrectionUncertainty> mJetCorrectorUnc;

//  edm::EDGetTokenT<reco::JetCorrector> jetCorrectorToken_;


  //Ele ID decisions objects
  edm::EDGetToken electronsMiniAODToken_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;

  //rho (PU energy density)
  edm::EDGetTokenT<double> rhoToken_;


  //Effective areas for isolation
  //EffectiveAreas   effectiveAreas_el_;
  //EffectiveAreas   effectiveAreas_ph_;

};