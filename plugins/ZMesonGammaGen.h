//---------- class declaration----------

class ZMesonGammaGen : public edm::stream::EDAnalyzer<> {///////// edm::EDAnalyzer {

 public:
  explicit ZMesonGammaGen(const edm::ParameterSet&);
  ~ZMesonGammaGen();

 private:
  virtual void beginJob(); /////remove override
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob(); /////remove override

  const edm::InputTag prunedGenParticles_;
  //const edm::InputTag genParticles_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ---------- member data ----------- //

  TTree *mytree;

  int genZ_ID_tree;
  float genZ_pT_tree;
  float genZ_eta_tree;
  float genZ_phi_tree;
  float genZ_E_tree;
  float genZ_mass_tree;

  int genGamma_ID_tree;
  float genGamma_pT_tree;
  float genGamma_eta_tree;
  float genGamma_phi_tree;
  float genGamma_E_tree;

  int genMeson_ID_tree;
  float genMeson_pT_tree;
  float genMeson_eta_tree;
  float genMeson_phi_tree;
  float genMeson_E_tree;
  float genMeson_mass_tree;

  int genTrackminus_ID_tree;
  float genTrackminus_pT_tree;
  float genTrackminus_eta_tree;
  float genTrackminus_phi_tree;
  float genTrackminus_E_tree;

  int genTrackplus_ID_tree;
  float genTrackplus_pT_tree;
  float genTrackplus_eta_tree;
  float genTrackplus_phi_tree;
  float genTrackplus_E_tree;

  float genTrackBig_pT_tree;
  float genTrackSmall_pT_tree;  

  float theta_pol; 
  float theta_pol_tree;

  int bosonID999;
  int nEvent;
  int particleNumber;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlesToken_; 
  //edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_; 

  //edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
};