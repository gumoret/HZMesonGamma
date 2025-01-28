//ROOT includes
#include <TH1F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>
#include <TMath.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"////////
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h" 
 
#include "ZMesonGammaGen.h"

#include <iostream>
using namespace std;  

 
// constructors and destructor
ZMesonGammaGen::ZMesonGammaGen(const edm::ParameterSet& iConfig) 
{
  prunedGenParticlesToken_ = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  //genParticlesToken_       = consumes<std::vector<reco::GenParticle> >(edm::InputTag("genParticles"));
  create_trees();
}

ZMesonGammaGen::~ZMesonGammaGen()
{
}


// ------------ method called for each event  ------------
void ZMesonGammaGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  iEvent.getByToken(prunedGenParticlesToken_, genParticles);
  //iEvent.getByToken(genParticlesToken_, genParticles);


  int genZ_ID     = -999;
  float genZ_pT   = -999.;
  float genZ_eta  = -999.;
  float genZ_phi  = -999.;
  float genZ_E    = -999.;
  float genZ_mass = -999.;

  int genGamma_ID    = -999;
  float genGamma_pT  = -999.;
  float genGamma_eta = -999.;
  float genGamma_phi = -999.;
  float genGamma_E   = -999.;

  int genMeson_ID     = -999;
  float genMeson_pT   = -999.;
  float genMeson_eta  = -999.;
  float genMeson_phi  = -999.;
  float genMeson_E    = -999.;
  float genMeson_mass = -999.;

  int genTrackminus_ID    = -999;
  float genTrackminus_pT  = -999.;
  float genTrackminus_eta = -999.;
  float genTrackminus_phi = -999.;
  float genTrackminus_E   = -999.;

  int genTrackplus_ID    = -999;
  float genTrackplus_pT  = -999.;
  float genTrackplus_eta = -999.;
  float genTrackplus_phi = -999.;
  float genTrackplus_E   = -999.;

  float genTrackBig_pT   = -999.;
  float genTrackSmall_pT= -999.;

  genZ_ID_tree   = 0;
  genZ_pT_tree   = 0.;
  genZ_eta_tree  = 0.;
  genZ_phi_tree  = 0.;
  genZ_E_tree    = 0.;
  genZ_mass_tree = 0.;

  genGamma_ID_tree  = 0;
  genGamma_pT_tree  = 0.;
  genGamma_eta_tree = 0.;
  genGamma_phi_tree = 0.;
  genGamma_E_tree   = 0.;

  genMeson_ID_tree   = 0;
  genMeson_pT_tree   = 0.;
  genMeson_eta_tree  = 0.;
  genMeson_phi_tree  = 0.;
  genMeson_E_tree    = 0.;
  genMeson_mass_tree = 0.;

  genTrackminus_ID_tree  = 0;
  genTrackminus_pT_tree  = 0.;
  genTrackminus_eta_tree = 0.;
  genTrackminus_phi_tree = 0.;
  genTrackminus_E_tree   = 0.;

  genTrackplus_ID_tree  = 0;
  genTrackplus_pT_tree  = 0.;
  genTrackplus_eta_tree = 0.;
  genTrackplus_phi_tree = 0.;
  genTrackplus_E_tree   = 0.;

  genTrackBig_pT_tree   = 0;
  genTrackSmall_pT_tree = 0.;  

  theta_pol     = -10.;
  TLorentzVector mu[2];
  theta_pol_tree = 0.;


  
  for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){

    //if it is a Z with 2 daughters (sometimes Pythia sends a Z in itself, therefore only 1 daughter)
    if(gen->pdgId() == 25 && gen->numberOfDaughters() == 2){

      //for each daughter
      for(int i = 0; i < 2; i++){

        //if daughters are not Phi or Rho and gamma, continue
        if( !(gen->daughter(i)->pdgId() == 22 || (gen->daughter(i)->pdgId() == 333 || gen->daughter(i)->pdgId() == 113)) ) continue;

        //save Z variables
        genZ_ID   = gen->pdgId();
        genZ_pT   = gen->pt();
        genZ_eta  = gen->eta();
        genZ_phi  = gen->phi();
        genZ_E    = gen->energy();
        genZ_mass = gen->p4().M();
              
        //cout << "gen->daughter(i)->pdgId() = " << gen->daughter(i)->pdgId() << endl;

        //if daughter is Gamma
        if(gen->daughter(i)->pdgId() == 22){ 
                   
          //save Gamma variables
          genGamma_ID  = gen->daughter(i)->pdgId();
          genGamma_pT  = gen->daughter(i)->pt();
          genGamma_eta = gen->daughter(i)->eta();
          genGamma_phi = gen->daughter(i)->phi();
          genGamma_E   = gen->daughter(i)->energy();
        } 

        //if daughter(i) is a Phi or a Rho
        if(gen->daughter(i)->pdgId() == 333 || gen->daughter(i)->pdgId() == 113){

          //save Phi/Rho variables
          genMeson_ID   = gen->daughter(i)->pdgId();
          genMeson_pT   = gen->daughter(i)->pt();
          genMeson_eta  = gen->daughter(i)->eta();
          genMeson_phi  = gen->daughter(i)->phi();
          genMeson_E    = gen->daughter(i)->energy();
          genMeson_mass = gen->daughter(i)->p4().M();
  
          //if Meson has two daughters
          if(gen->daughter(i)->numberOfDaughters() == 2){
            //cout << "try1" << endl;
      
            //for each Meson daughter
            for(int j = 0; j < 2; j++){

              //if daughter(j) is a K-
              if(gen->daughter(i)->daughter(j)->pdgId() == -321){

                //save K- variables
                genTrackminus_ID  = gen->daughter(i)->daughter(j)->pdgId();
                genTrackminus_pT  = gen->daughter(i)->daughter(j)->pt();
                genTrackminus_eta = gen->daughter(i)->daughter(j)->eta();
                genTrackminus_phi = gen->daughter(i)->daughter(j)->phi();
                genTrackminus_E   = gen->daughter(i)->daughter(j)->energy();
              }

              //if daughter(j) is a K+
              if(gen->daughter(i)->daughter(j)->pdgId() == 321){

                //cout << "try2" << endl;

                //save K+ variables
                genTrackplus_ID  = gen->daughter(i)->daughter(j)->pdgId();
                genTrackplus_pT  = gen->daughter(i)->daughter(j)->pt();
                genTrackplus_eta = gen->daughter(i)->daughter(j)->eta();
                genTrackplus_phi = gen->daughter(i)->daughter(j)->phi();
                genTrackplus_E   = gen->daughter(i)->daughter(j)->energy();
                
              }

       
              //if daughter(j) is a Pi-
              if(gen->daughter(i)->daughter(j)->pdgId() == -211){

                //cout << "try3" << endl;
                //save Pi- variables
                genTrackminus_ID  = gen->daughter(i)->daughter(j)->pdgId();
                genTrackminus_pT  = gen->daughter(i)->daughter(j)->pt();
                genTrackminus_eta = gen->daughter(i)->daughter(j)->eta();
                genTrackminus_phi = gen->daughter(i)->daughter(j)->phi();
                genTrackminus_E   = gen->daughter(i)->daughter(j)->energy();
              }
        
              //if daughter(j) is a Pi+
              if(gen->daughter(i)->daughter(j)->pdgId() == 211){

                //cout << "try3" << endl;
                //save Pi+ variables
                genTrackplus_ID  = gen->daughter(i)->daughter(j)->pdgId();
                genTrackplus_pT  = gen->daughter(i)->daughter(j)->pt();
                genTrackplus_eta = gen->daughter(i)->daughter(j)->eta();
                genTrackplus_phi = gen->daughter(i)->daughter(j)->phi();
                genTrackplus_E   = gen->daughter(i)->daughter(j)->energy();
              
              }

      
            }//j-forloop end

            if (genTrackplus_pT < genTrackminus_pT){ 
              genTrackBig_pT = genTrackminus_pT;
              genTrackSmall_pT = genTrackplus_pT;
            }
            else if(genTrackplus_pT > genTrackminus_pT){
              genTrackBig_pT = genTrackplus_pT;
              genTrackSmall_pT = genTrackminus_pT;
            }
            //cout << genTrackSmall_pT << endl;
            //cout << "pdgIDE" << endl;
            //cout<<genTrackplus_ID<<endl;
          }//"if Phi/Rho has two daughters" end   //quindi: il mesone è sempre un phi/rho, infatti pdgID è stampato sempre, mentre non sempre il mesone decade in due figlie,
          //cout<<"pdgID"<<endl;             // infatti pdgIDE non sempre è stampato
          //cout<<genTrackplus_ID<<endl;
          //cout<<genTrackminus_ID<<endl;
        }//"if it is a Phi/Rho" end
      }//i-forloop end
    }//"if it is a Z" end
  }//gen-forloop end
  
  genZ_ID_tree   = genZ_ID;
  genZ_pT_tree   = genZ_pT;
  genZ_eta_tree  = genZ_eta;
  genZ_phi_tree  = genZ_phi;
  genZ_E_tree    = genZ_E;
  genZ_mass_tree = genZ_mass;
  
  genGamma_ID_tree  = genGamma_ID;
  genGamma_pT_tree  = genGamma_pT;
  genGamma_eta_tree = genGamma_eta;
  genGamma_phi_tree = genGamma_phi;
  genGamma_E_tree   = genGamma_E;
  
  genMeson_ID_tree   = genMeson_ID;
  genMeson_pT_tree   = genMeson_pT;
  genMeson_eta_tree  = genMeson_eta;
  genMeson_phi_tree  = genMeson_phi;
  genMeson_E_tree    = genMeson_E;
  genMeson_mass_tree = genMeson_mass;

  genTrackminus_ID_tree  = genTrackminus_ID;
  genTrackminus_pT_tree  = genTrackminus_pT;
  genTrackminus_eta_tree = genTrackminus_eta;
  genTrackminus_phi_tree = genTrackminus_phi;
  genTrackminus_E_tree   = genTrackminus_E;

  genTrackplus_ID_tree  = genTrackplus_ID;
  genTrackplus_pT_tree  = genTrackplus_pT;
  genTrackplus_eta_tree = genTrackplus_eta;
  genTrackplus_phi_tree = genTrackplus_phi;
  genTrackplus_E_tree   = genTrackplus_E;

  genTrackBig_pT_tree   = genTrackBig_pT;
  genTrackSmall_pT_tree = genTrackSmall_pT;
  
  
  mytree->Fill();

}


//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void ZMesonGammaGen::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen info");

  mytree->Branch("genZ_ID",&genZ_ID_tree);
  mytree->Branch("genZ_pT",&genZ_pT_tree);
  mytree->Branch("genZ_eta",&genZ_eta_tree);
  mytree->Branch("genZ_phi",&genZ_phi_tree);
  mytree->Branch("genZ_E",&genZ_E_tree);
  mytree->Branch("genZ_mass",&genZ_mass_tree);

  
  mytree->Branch("genGamma_ID",&genGamma_ID_tree);
  mytree->Branch("genGamma_pT",&genGamma_pT_tree);
  mytree->Branch("genGamma_eta",&genGamma_eta_tree);
  mytree->Branch("genGamma_phi",&genGamma_phi_tree);
  mytree->Branch("genGamma_E",&genGamma_E_tree);
  
  mytree->Branch("genMeson_ID",&genMeson_ID_tree);
  mytree->Branch("genMeson_pT",&genMeson_pT_tree);
  mytree->Branch("genMeson_eta",&genMeson_eta_tree);
  mytree->Branch("genMeson_phi",&genMeson_phi_tree);
  mytree->Branch("genMeson_E",&genMeson_E_tree);
  mytree->Branch("genMeson_mass",&genMeson_mass_tree);


  mytree->Branch("genTrackminus_ID",&genTrackminus_ID_tree);
  mytree->Branch("genTrackminus_pT",&genTrackminus_pT_tree);
  mytree->Branch("genTrackminus_eta",&genTrackminus_eta_tree);
  mytree->Branch("genTrackminus_phi",&genTrackminus_phi_tree);
  mytree->Branch("genTrackminus_E",&genTrackminus_E_tree);

  mytree->Branch("genTrackplus_ID",&genTrackplus_ID_tree);
  mytree->Branch("genTrackplus_pT",&genTrackplus_pT_tree);
  mytree->Branch("genTrackplus_eta",&genTrackplus_eta_tree);
  mytree->Branch("genTrackplus_phi",&genTrackplus_phi_tree);
  mytree->Branch("genTrackplus_E",&genTrackplus_E_tree);

  mytree->Branch("genTrackBig_pT",&genTrackBig_pT_tree);
  mytree->Branch("genTrackSmall_pT",&genTrackSmall_pT_tree);

  mytree->Branch("theta_pol",&theta_pol_tree);

    
}

void ZMesonGammaGen::beginJob()
{
}

void ZMesonGammaGen::endJob() 
{
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZMesonGammaGen);
