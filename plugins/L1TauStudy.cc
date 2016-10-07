// -*- C++ -*-
//
// Package:    L1TauStudy
// Class:      L1TauStudy
// 
/**\class L1TauStudy L1TauStudy.cc SLHCUpgradeSimulations/L1TauStudy/plugins/L1TauStudy.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Thu, 01 Sep 2016 20:40:06 GMT
// $Id$
//
//


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/SLHC/interface/L1EGCrystalCluster.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/Particle/interface/ParticleTable.h"

#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTPixelTrack.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/L1TkElectronTrackMatchAlgo.h"
#include "SLHCUpgradeSimulations/L1TauStudy/plugins/helpers.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <TLorentzVector.h>
#include <memory>
#include <math.h>
#include <vector>
#include <list>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SLHCUpgradeSimulations/L1TauStudy/interface/triggerGeometryTools.hh"

//
// class declaration
//
using namespace edm;
using std::cout;
using std::endl;
using std::vector;

struct l1Object{
    double objectPt;
    double objectEta;
    double objectPhi;
    double trackPt; //Primary Track Pt
    double trackEta; //Primary Track Eta
    double trackPhi; // Primary Track Phi
    double HCALEnergy;
    double ECALEnergy;
    int tauDecayMode;
    double iso;
};

bool compareByPt (TTTrack<Ref_PixelDigi_> i,TTTrack<Ref_PixelDigi_> j) { 
  return(i.getMomentum().perp() > j.getMomentum().perp()); 
	 };

bool compareL1ObjectByPt (l1Object i,l1Object j) { 
  return(i.objectPt > j.objectPt); 
	 };


class L1TauStudy : public edm::EDAnalyzer {
  typedef std::vector<TTTrack<Ref_PixelDigi_>> L1TkTrackCollectionType;
  typedef vector<reco::GenParticle> GenParticleCollectionType;

  struct genVisTau{

    reco::Candidate::LorentzVector p4;
    int decayMode;
  };

   public:
      explicit L1TauStudy(const edm::ParameterSet&);
      ~L1TauStudy();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  void initializHcalTpgs(edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs,std::vector<TLorentzVector> &allHcalTPGs,const edm::EventSetup& es);
  void initializEcalTpgs(edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs,std::vector<TLorentzVector> &allEcalTPGs);

  void printLorentzVector(TLorentzVector input){
    std::cout<<"pt: "<<input.Pt()<<" eta: "<<input.Eta()<<" phi: "<<input.Phi();
  }

  void printL1Object(l1Object input){
    std::cout<<"pt: "<<input.objectPt<<" eta: "<<input.objectEta<<" phi: "<<input.objectPhi;
    std::cout<<"hcalEnergy "<<input.HCALEnergy<<" ecalEnergy "<<input.ECALEnergy<<" iso "<<input.iso<<std::endl;
  }

  void crossCleanThreeProngTaus(vector<l1Object> &inputL1Objects){
    //std::cout<<"BEFORE===================================="<<std::endl;
    //for(auto inputL1Object: inputL1Objects)
    //std::cout<<"inputL1Object Pt: "<<inputL1Object.objectPt<<" Eta: "<<inputL1Object.objectEta<<" Phi: "<<inputL1Object.objectPhi<<std::endl;

    //this loop must stop at the second to last element, hence -1
    for(unsigned int i = 0; i<inputL1Objects.size()-1; i++){
      //for(auto objectA : inputL1Objects){
      l1Object objectA = inputL1Objects.at(i);
      //std::cout<<"objectA Pt: "<<objectA.objectPt<<" eta: "<<objectA.objectEta<<" phi: "<<objectA.objectPhi<<std::endl;
      for(unsigned int j = i+1; j<inputL1Objects.size(); j++){
	l1Object objectB = inputL1Objects.at(j);
	//std::cout<<"objectB Pt: "<<objectB.objectPt<<" eta: "<<objectB.objectEta<<" phi: "<<objectB.objectPhi<<std::endl;
	if(objectA.objectPt  == objectB.objectPt  && 
	   objectA.objectEta == objectB.objectEta &&  
	   objectA.objectPhi == objectB.objectPhi 
	   ){
	  //if(objectA.trackPt>objectB.trackPt){
	  //std::cout<<"erase object"<<std::endl;
	  inputL1Objects.erase(inputL1Objects.begin()+j);
	  j = j-1;
	  //}else{
	  //inputL1Objects.erase(inputL1Objects.begin()+i);
	  //break;
	  //}
	}
      }
    }
    //print for sanity check
    //std::cout<<"AFTER========================================"<<std::endl;
    //for(auto inputL1Object: inputL1Objects)
    //std::cout<<"inputL1Object Pt: "<<inputL1Object.objectPt<<" Eta: "<<inputL1Object.objectEta<<" Phi: "<<inputL1Object.objectPhi<<std::endl;
  }

  float SumL1TrackPt(std::vector<TTTrack<Ref_PixelDigi_>> inputVector){
    float sum = 0;
    for(auto element : inputVector){
      sum += element.getMomentum().perp();
    }
    return(sum);
  }

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------

  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalSrc_; 
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalSrc_;

  //edm::InputTag ecalSrc_;
  //edm::InputTag hcalSrc_;
  edm::InputTag vtxLabel_;
  edm::InputTag tauSrc_;
  edm::InputTag genSrc_;
  double genMatchDeltaRcut;
  double tauHcalPt_;
  double tauEcalPt_;

  edm::InputTag L1TrackInputTag;
  edm::InputTag L1TrackPrimaryVertexTag;

  TTree* efficiencyTree;
  TTree* oneProngTree;
  TTree* threeProngTree;

  double genPt, genEta, genPhi;
  int decayMode, run, lumi, event;
  double trackPt, trackEta, trackPhi;
  double l1TauPt, l1TauEta, l1TauPhi;
  double gen1ProngPt, gen1ProngEta, gen1ProngPhi;
  double gen3ProngPt, gen3ProngEta, gen3ProngPhi;
  int gen1ProngDecayMode;
  int gen3ProngDecayMode;

  l1Object oneProngTau;
  l1Object threeProngTau;

  TH1F* nEvents;
  TH1F* track_pt;
  TH1F* track_pt_eta2p1;
  TH1F* track_pt_eta2p4;
  TH1F* l1Tau_pt;     
  TH1F* l1Tau_pt_eta2p4;     
  TH1F* l1Tau_pt_eta2p1;

  TH1F* l1SingleProngTau_pt;     
  TH1F* l1SingleProngTau_pt_eta2p4;     
  TH1F* l1SingleProngTau_pt_eta2p1;

  TH1F* l1SingleProngTauIso_pt;     
  TH1F* l1SingleProngTauIso_pt_eta2p4;     
  TH1F* l1SingleProngTauIso_pt_eta2p1;

  TH1F* l1ThreeProngTau_pt;     
  TH1F* l1ThreeProngTau_pt_eta2p4;     
  TH1F* l1ThreeProngTau_pt_eta2p1;

  TH1F* l1ThreeProngTauIso_pt;     
  TH1F* l1ThreeProngTauIso_pt_eta2p4;     
  TH1F* l1ThreeProngTauIso_pt_eta2p1;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1TauStudy::L1TauStudy(const edm::ParameterSet& cfg) :
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_((cfg.getParameter<edm::InputTag>( "vertices"    ))),
  genSrc_ ((cfg.getParameter<edm::InputTag>( "genParticles"))),
  tauHcalPt_ (cfg.getUntrackedParameter<double>("tauHcalPtCut", -1)),
  tauEcalPt_ (cfg.getUntrackedParameter<double>("tauEcalPtCut", -1)),
  genMatchDeltaRcut(cfg.getUntrackedParameter<double>("genMatchDeltaRcut", 0.1))
{
   //now do what ever initialization is needed
  L1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
  L1TrackPrimaryVertexTag = cfg.getParameter<edm::InputTag>("L1TrackPrimaryVertexTag");
  //std::cout<<"creating efficiency tree"<<std::endl;
  edm::Service<TFileService> fs;
  efficiencyTree = fs->make<TTree>("efficiencyTree", "Crystal cluster individual crystal pt values");
  efficiencyTree->Branch("run",    &run,     "run/I");
  efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
  efficiencyTree->Branch("event",  &event,   "event/I");
  
  efficiencyTree->Branch("genPt",    &genPt,   "genPt/D");
  efficiencyTree->Branch("genEta",       &genEta,   "genEta/D");
  efficiencyTree->Branch("genPhi",       &genPhi,   "genPhi/D");
  efficiencyTree->Branch("genDM",       &decayMode,   "genDM/I");

  efficiencyTree->Branch("trackPt",    &trackPt,   "trackPt/D");
  efficiencyTree->Branch("trackEta",       &trackEta,   "trackEta/D");
  efficiencyTree->Branch("trackPhi",       &trackPhi,   "trackPhi/D");

  efficiencyTree->Branch("l1TauPt",    &l1TauPt,   "l1TauPt/D");
  efficiencyTree->Branch("l1TauEta",       &l1TauEta,   "l1TauEta/D");
  efficiencyTree->Branch("l1TauPhi",       &l1TauPhi,   "l1TauPhi/D");

  oneProngTree = fs->make<TTree>("oneProngTree", "Crystal cluster individual crystal pt values");
  oneProngTree->Branch("run",    &run,     "run/I");
  oneProngTree->Branch("lumi",   &lumi,    "lumi/I");
  oneProngTree->Branch("event",  &event,   "event/I");
  
  oneProngTree->Branch("genPt",   &gen1ProngPt,   "genPt/D");
  oneProngTree->Branch("genEta",  &gen1ProngEta,   "genEta/D");
  oneProngTree->Branch("genPhi",  &gen1ProngPhi,   "genPhi/D");
  oneProngTree->Branch("genDM",   &gen1ProngDecayMode,   "genDM/I");

  oneProngTree->Branch("trackPt",  &oneProngTau.trackPt,   "trackPt/D");
  oneProngTree->Branch("trackEta", &oneProngTau.trackEta,   "trackEta/D");
  oneProngTree->Branch("trackPhi", &oneProngTau.trackPhi,   "trackPhi/D");

  oneProngTree->Branch("objectPt",  &oneProngTau.objectPt,   "objectPt/D");
  oneProngTree->Branch("objectEta", &oneProngTau.objectEta,   "objectEta/D");
  oneProngTree->Branch("objectPhi", &oneProngTau.objectPhi,   "objectPhi/D");

  oneProngTree->Branch("HCALEnergy",  &oneProngTau.HCALEnergy,   "HCALEnergy/D");
  oneProngTree->Branch("ECALEnergy",  &oneProngTau.ECALEnergy,   "ECALEnergy/D");
  
  oneProngTree->Branch("isoRaw",  &oneProngTau.iso,   "isoRaw/D");
  oneProngTree->Branch("decayMode",  &oneProngTau.tauDecayMode,   "decayMode/I");


  threeProngTree = fs->make<TTree>("threeProngTree", "Crystal cluster individual crystal pt values");
  threeProngTree->Branch("run",    &run,     "run/I");
  threeProngTree->Branch("lumi",   &lumi,    "lumi/I");
  threeProngTree->Branch("event",  &event,   "event/I");
  
  threeProngTree->Branch("genPt",   &gen3ProngPt,   "genPt/D");
  threeProngTree->Branch("genEta",  &gen3ProngEta,   "genEta/D");
  threeProngTree->Branch("genPhi",  &gen3ProngPhi,   "genPhi/D");
  threeProngTree->Branch("genDM",   &gen3ProngDecayMode,   "genDM/I");

  threeProngTree->Branch("trackPt",  &threeProngTau.trackPt,   "trackPt/D");
  threeProngTree->Branch("trackEta", &threeProngTau.trackEta,   "trackEta/D");
  threeProngTree->Branch("trackPhi", &threeProngTau.trackPhi,   "trackPhi/D");

  threeProngTree->Branch("objectPt",  &threeProngTau.objectPt,   "objectPt/D");
  threeProngTree->Branch("objectEta", &threeProngTau.objectEta,   "objectEta/D");
  threeProngTree->Branch("objectPhi", &threeProngTau.objectPhi,   "objectPhi/D");


  threeProngTree->Branch("HCALEnergy",  &threeProngTau.HCALEnergy,   "HCALEnergy/D");
  threeProngTree->Branch("ECALEnergy",  &threeProngTau.ECALEnergy,   "ECALEnergy/D");

  threeProngTree->Branch("isoRaw",     &threeProngTau.iso,   "isoRaw/D");
  threeProngTree->Branch("decayMode",  &threeProngTau.tauDecayMode,   "decayMode/I");


  nEvents     = fs->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  track_pt      = fs->make<TH1F>( "track_pt"  , "p_{t}", 200,  0., 200. );
  track_pt_eta2p4      = fs->make<TH1F>( "track_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  track_pt_eta2p1      = fs->make<TH1F>( "track_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1Tau_pt             = fs->make<TH1F>( "l1Tau_pt"  , "p_{t}", 200,  0., 200. );
  l1Tau_pt_eta2p4      = fs->make<TH1F>( "l1Tau_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1Tau_pt_eta2p1      = fs->make<TH1F>( "l1Tau_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1SingleProngTau_pt             = fs->make<TH1F>( "l1SingleProngTau"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTau_pt_eta2p4      = fs->make<TH1F>( "l1SingleProngTau_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTau_pt_eta2p1      = fs->make<TH1F>( "l1SingleProngTau_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1SingleProngTauIso_pt             = fs->make<TH1F>( "l1SingleProngTauIso"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTauIso_pt_eta2p4      = fs->make<TH1F>( "l1SingleProngTauIso_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1SingleProngTauIso_pt_eta2p1      = fs->make<TH1F>( "l1SingleProngTauIso_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1ThreeProngTau_pt             = fs->make<TH1F>( "l1ThreeProngTau_pt"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTau_pt_eta2p4      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTau_pt_eta2p1      = fs->make<TH1F>( "l1ThreeProngTau_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

  l1ThreeProngTauIso_pt             = fs->make<TH1F>( "l1ThreeProngTauIso_pt"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTauIso_pt_eta2p4      = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p4"  , "p_{t}", 200,  0., 200. );
  l1ThreeProngTauIso_pt_eta2p1      = fs->make<TH1F>( "l1ThreeProngTauIso_pt_eta2p1"  , "p_{t}", 200,  0., 200. );

}


L1TauStudy::~L1TauStudy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1TauStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   nEvents->Fill(1);

   run = iEvent.id().run();
   lumi = iEvent.id().luminosityBlock();
   event = iEvent.id().event();

   // L1 Tracks
   edm::Handle<L1TkTrackCollectionType> l1trackHandle;
   iEvent.getByLabel(L1TrackInputTag, l1trackHandle);

   // L1 Track based primary vertex
   edm::Handle<L1TkPrimaryVertexCollection> l1PrimaryVertexHandle;
   iEvent.getByLabel(L1TrackPrimaryVertexTag, l1PrimaryVertexHandle);

   //decal ecal and hcal tpgs
   edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
   edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  

   // Get genParticles
   edm::Handle<GenParticleCollectionType> genParticleHandle;
   if(!iEvent.getByLabel(genSrc_,genParticleHandle))
     std::cout<<"No gen Particles Found "<<std::endl;

   vector<TTTrack<Ref_PixelDigi_> > l1Tracks;
   l1Tracks.clear();
   //Find and sort the tracks
   for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
       edm::Ptr<TTTrack<Ref_PixelDigi_>> ptr(l1trackHandle, track_index);
       double pt = ptr->getMomentum().perp();
       double eta = ptr->getMomentum().eta();
       //std::cout<<"# "<<track_index<<": "<<pt<<std::endl;
       l1Tracks.push_back(l1trackHandle->at(track_index));       
     }
   std::sort(l1Tracks.begin(),l1Tracks.end(),compareByPt);   

   std::vector<TLorentzVector> allEcalTPGs;

   if(!iEvent.getByToken(ecalSrc_, ecalTPGs))
     std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
   else{
     initializEcalTpgs(ecalTPGs , allEcalTPGs);
   }
   
   std::vector<TLorentzVector> allHcalTPGs;

   if(!iEvent.getByToken(hcalSrc_, hcalTPGs))
     std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
   else{
     initializHcalTpgs(hcalTPGs , allHcalTPGs, iSetup);
   }

   //create high pt track distribution
   if(l1Tracks.size()>0){
     TTTrack<Ref_PixelDigi_> l1Track = l1Tracks.at(0);
     double pt = l1Track.getMomentum().perp();
     double eta = l1Track.getMomentum().eta();

     track_pt->Fill(pt);
     if(fabs(eta<2.1))
       track_pt_eta2p1->Fill(pt);
     if(fabs(eta<2.4))
       track_pt_eta2p4->Fill(pt);
   }

   vector<TTTrack<Ref_PixelDigi_> > oneProngSeeds;
   vector<TTTrack<Ref_PixelDigi_> > threeProngSeeds;
   for(auto l1Track : l1Tracks){

   if(l1Track.getMomentum().perp()>7)
     //this is used for efficiency studies
     threeProngSeeds.push_back(l1Track);

   if(l1Track.getMomentum().perp()>15)//arbitrary seed
     //this is used for efficiency studies
     oneProngSeeds.push_back(l1Track);
   }

   /* Create 1 prong Tau objects
      1) Create charged hadrons
	  -associate proper trigger tower to track
          -H/E requirement with respect to single track pt
   */
   //make seed vector

   
   vector<reco::GenParticle> genTaus;
   vector<reco::GenParticle> genParticles;
   //reco::GenParticle* genTau;
   for(unsigned int i = 0; i< genParticleHandle->size(); i++){
     edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
     genParticles.push_back(*ptr);

     if(abs(ptr->pdgId())==15){
       genTaus.push_back(*ptr);
       //const reco::GenParticle *genTau = &(*ptr);  
     }
   }



   std::vector<genVisTau> GenOneProngTaus;
   std::vector<genVisTau> GenOneProngPi0Taus;
   std::vector<genVisTau> GenThreeProngTaus;
   //Find and Sort the 1 Prong, 1 Prong + pi0 and 3 Prong Taus
   for(auto genTau: genTaus){
     std::cout<<"got tau"<<std::endl;

     std::cout<<"Tau Decay Mode "<<GetDecayMode(&genTau)<<std::endl;
     reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &genParticles);
     genVisTau Temp;
     genPt = visGenTau.pt();
     genEta = visGenTau.eta();
     genPhi = visGenTau.phi();
     decayMode = GetDecayMode(&genTau);
     Temp.p4 = visGenTau;
     Temp.decayMode = decayMode;
     std::cout<<"tau vis pt: "<<genPt<<" genEta: "<<genEta<<" genPhi: "<<genPhi<<std::endl;

     if(decayMode >21 ){
       std::cout<<"found 3 prong tau: "<<decayMode<<std::endl;
       GenThreeProngTaus.push_back(Temp);
     }

     if(decayMode == 10 ){
       GenOneProngTaus.push_back(Temp);
     }

     if(decayMode > 10 && decayMode < 20 ){
       GenOneProngPi0Taus.push_back(Temp);
     }
   }

     /*

     trackPt  = -10;
     trackEta = -10;
     trackPhi = -10;

       for(auto l1Track : l1Tracks){

	 //std::cout<<"l1track"<<std::endl;
	 double dR = deltaR( l1Track.getMomentum(), genTau.p4());
	 if( dR < 0.2){
	   std::cout<<"found matched track pt: "<<l1Track.getMomentum().perp()<<std::endl;
	   trackPt = l1Track.getMomentum().perp();
	   trackEta = l1Track.getMomentum().eta();
	   trackPhi = l1Track.getMomentum().phi();
	   
	   l1TauPt = l1Track.getMomentum().perp();
	   l1TauEta = l1Track.getMomentum().eta();
	   l1TauPhi = l1Track.getMomentum().phi();
	   
	   break;
	   
	 }
       }
     }
     ///efficiencyTree filling goes here
     
   }
     */

   std::vector<l1Object> l1OneProngTaus;
   std::vector<l1Object> l1OneProngTaus_eta2p1;
   std::vector<l1Object> l1OneProngTaus_eta2p4;

   std::vector<l1Object> l1OneProngTausIso;
   std::vector<l1Object> l1OneProngTausIso_eta2p1;
   std::vector<l1Object> l1OneProngTausIso_eta2p4;

   ///Fill 1 prong taus here
   triggerGeometryTools trigTools;
   for(auto l1Track : l1Tracks){
     if(l1Track.getMomentum().perp()>15){//arbitrary seed
       vector<TLorentzVector> isoTracks;
       vector<TLorentzVector> isoECAL;
       vector<TLorentzVector> isoHCAL;

       TLorentzVector tempTrack;
       double pt = l1Track.getMomentum().perp();
       double eta = l1Track.getMomentum().eta();
       double phi = l1Track.getMomentum().phi();
       tempTrack.SetPtEtaPhiE(pt, eta, phi, pt);

       //find closest hcalTPG
       std::vector<TLorentzVector> nearHCALTowers;
       for(auto hcalTPG : allHcalTPGs){
	 if(hcalTPG.Pt()>0)
	   if(tempTrack.DeltaR(hcalTPG)<0.087){
	     //std::cout<<"hcal tpgs ";
	     //printLorentzVector(hcalTPG);
	     //std::cout<<std::endl;
	     nearHCALTowers.push_back(hcalTPG);
	   }
       }
       float hcalEnergy = trigTools.SumTLorentzPt(nearHCALTowers);
       
       std::vector<TLorentzVector> nearECALTowers;
       for(auto ecalTPG : allEcalTPGs){
	 if(ecalTPG.Pt()>0)
	   if(tempTrack.DeltaR(ecalTPG)<0.087){
	     //std::cout<<"ecal tpgs ";
	     //printLorentzVector(ecalTPG);
	     //std::cout<<std::endl;
	     nearECALTowers.push_back(ecalTPG);
	   }
       }
       float ecalEnergy = trigTools.SumTLorentzPt(nearECALTowers);

       for(auto l1TrackB : l1Tracks){
	 TLorentzVector tempTrackB;
	 double ptB = l1TrackB.getMomentum().perp();
	 double etaB = l1TrackB.getMomentum().eta();
	 double phiB = l1TrackB.getMomentum().phi();
	 tempTrackB.SetPtEtaPhiE(ptB, etaB, phiB, ptB);
	 //check to see if it is the same track
	 if(tempTrackB == tempTrack)
	   continue;
	 if(tempTrack.DeltaR(tempTrackB)<0.4)
	   isoTracks.push_back(tempTrackB);
       }

       float trackIsoPt = trigTools.SumTLorentzPt(isoTracks);

       l1Object tempL1Object;
       tempL1Object.trackPt = pt;
       tempL1Object.trackEta = eta;
       tempL1Object.trackPhi = phi;
       tempL1Object.objectPt = pt;
       tempL1Object.objectEta = eta;
       tempL1Object.objectPhi = phi;
       tempL1Object.ECALEnergy = ecalEnergy;
       tempL1Object.HCALEnergy = hcalEnergy;
       tempL1Object.tauDecayMode = 1;
       tempL1Object.iso = trackIsoPt;

       //loose ID
       if(hcalEnergy>0&&((ecalEnergy+hcalEnergy)>(pt*0.1))){
	 l1OneProngTaus.push_back(tempL1Object);
	 if(abs(eta)<2.1)
	   l1OneProngTaus_eta2p1.push_back(tempL1Object);
	 if(abs(eta)<2.4)
	   l1OneProngTaus_eta2p4.push_back(tempL1Object);

	 if(trackIsoPt<10){
	   l1OneProngTausIso.push_back(tempL1Object);
	   if(abs(eta)<2.1)
	     l1OneProngTausIso_eta2p1.push_back(tempL1Object);
	   if(abs(eta)<2.4)
	     l1OneProngTausIso_eta2p4.push_back(tempL1Object);
	 }
       }
     }
   }

   int i = 0;
   for(auto l1ThreeProngTau: l1OneProngTaus){
     std::cout<<"i: "<<i<<std::endl;
     printL1Object(l1ThreeProngTau);
     i++;
   }


   //std::cout<<"here 1"<<std::endl;
   if(l1OneProngTaus.size()>0)
     l1SingleProngTau_pt->Fill(l1OneProngTaus.at(0).objectPt);
   
   if(l1OneProngTaus_eta2p1.size()>0)
     l1SingleProngTau_pt_eta2p1->Fill(l1OneProngTaus_eta2p1.at(0).objectPt);

   if(l1OneProngTaus_eta2p4.size()>0)
     l1SingleProngTau_pt_eta2p4->Fill(l1OneProngTaus_eta2p4.at(0).objectPt);

   if(l1OneProngTausIso.size()>0)
     l1SingleProngTauIso_pt->Fill(l1OneProngTausIso.at(0).objectPt);
   
   if(l1OneProngTausIso_eta2p1.size()>0)
     l1SingleProngTauIso_pt_eta2p1->Fill(l1OneProngTausIso_eta2p1.at(0).objectPt);

   if(l1OneProngTausIso_eta2p4.size()>0)
     l1SingleProngTauIso_pt_eta2p4->Fill(l1OneProngTausIso_eta2p4.at(0).objectPt);


   //create one prong taus when gen one prong tau is nearby
   for(auto GenOneProngTau: GenOneProngTaus){
     gen1ProngPt = GenOneProngTau.p4.Pt();
     gen1ProngEta = GenOneProngTau.p4.Eta();
     gen1ProngPhi = GenOneProngTau.p4.Phi();
     gen1ProngDecayMode = GenOneProngTau.decayMode;
     std::cout<<"One Prong Tau (decaymode "<<gen1ProngDecayMode <<") PT: "<<gen1ProngPt<<" Eta: "<<gen1ProngEta<<" Phi: "<<gen1ProngPhi<<std::endl;
     vector<TTTrack<Ref_PixelDigi_>> oneProngTaus;
     oneProngTau.HCALEnergy = 0;
     oneProngTau.ECALEnergy = 0;
     oneProngTau.trackPt  = 0;
     oneProngTau.trackEta = 0;
     oneProngTau.trackPhi = 0;
     
     triggerGeometryTools trigTools;
     for(auto l1OneProngTau : l1OneProngTaus){
       //for(auto oneProngSeed : oneProngSeeds){
       double l1ObjectPt  = l1OneProngTau.objectPt;
       double l1ObjectEta = l1OneProngTau.objectEta;
       double l1ObjectPhi = l1OneProngTau.objectPhi;

       TLorentzVector tempTrack;
       tempTrack.SetPtEtaPhiE(l1ObjectPt, l1ObjectEta, l1ObjectPhi, l1ObjectPt);
       TLorentzVector tempGen;
       tempGen.SetPtEtaPhiE(GenOneProngTau.p4.Pt(),GenOneProngTau.p4.Eta(),GenOneProngTau.p4.Phi(),GenOneProngTau.p4.Pt());
       if(tempTrack.DeltaR(tempGen)<0.4){
	 oneProngTau.trackPt  = tempTrack.Pt();
	 oneProngTau.trackEta = tempTrack.Eta();
	 oneProngTau.trackPhi = tempTrack.Phi();
	 /*	 
	 //find closest hcalTPG
	 std::vector<TLorentzVector> nearHCALTowers;
	 //TODO: check if dynamics cluster based on the track pt is OK
	 for(auto hcalTPG : allHcalTPGs){
	   if(hcalTPG.Pt()>0)
	     if(tempTrack.DeltaR(hcalTPG)<0.087){
	       //std::cout<<"hcal tpgs ";
	       //printLorentzVector(hcalTPG);
	       //std::cout<<std::endl;
	       nearHCALTowers.push_back(hcalTPG);
	     }
	     }*/
	 oneProngTau.HCALEnergy = l1OneProngTau.HCALEnergy;
	 /*
	 std::vector<TLorentzVector> nearECALTowers;
	 for(auto ecalTPG : allEcalTPGs){
	   if(ecalTPG.Pt()>0)
	     if(tempTrack.DeltaR(ecalTPG)<0.087){
	       //std::cout<<"ecal tpgs ";
	       //printLorentzVector(ecalTPG);
	       //std::cout<<std::endl;
	       nearECALTowers.push_back(ecalTPG);
	     }
	     }*/
	 oneProngTau.ECALEnergy = l1OneProngTau.ECALEnergy;

	 oneProngTau.iso = l1OneProngTau.iso;

	 oneProngTau.tauDecayMode = l1OneProngTau.tauDecayMode;
	 break;
       }
     }
     //std::cout<<"one prong tau, track pt: "<<oneProngTau.trackPt<<" eta: "<<oneProngTau.trackEta<<" phi: "<<oneProngTau.trackPhi;
     oneProngTree->Fill();
   }
   //std::cout<<"here 2"<<std::endl;
   //create the clustered 3 prong taus
   std::vector<l1Object> l1ThreeProngTaus;
   
   for(auto threeProngSeed: threeProngSeeds){
     l1Object tempL1Object;
     double pt = threeProngSeed.getMomentum().perp();
     double eta = threeProngSeed.getMomentum().eta();
     double phi = threeProngSeed.getMomentum().phi();
     tempL1Object.trackPt  = pt;
     tempL1Object.trackEta = eta;
     tempL1Object.trackPhi = phi;
     //std::cout<<"Seed pt: "<<pt<<" eta: "<<eta<<" phi: "<<phi<<std::endl;

     vector<TTTrack<Ref_PixelDigi_>> threeProngTauTracks;
     //Fill Iso Tracks if greater than 3 tracks in 0.1 cone OR if between 0.4 and 0.1
     vector<TTTrack<Ref_PixelDigi_>> threeProngTauTracks_Iso;

     threeProngTauTracks.push_back(threeProngSeed);
     int iMatched = 1;
     //loop over all tracks
     for(auto l1Track : l1Tracks){
       //don't double count
       if(l1Track.getMomentum() == threeProngSeed.getMomentum())
	 continue;

       if(deltaR(l1Track.getMomentum(),threeProngSeed.getMomentum())<0.4 &&
	  deltaR(l1Track.getMomentum(),threeProngSeed.getMomentum())>0.1){
	 threeProngTauTracks_Iso.push_back(l1Track);
       }
       else if(deltaR(l1Track.getMomentum(),threeProngSeed.getMomentum())<0.1){
	 //std::cout<<"i "<< iMatched <<" pt: "<<l1Track.getMomentum().perp()<<" eta: "<<l1Track.getMomentum().eta()<<" phi: "<<l1Track.getMomentum().phi()<<std::endl;

	 threeProngTauTracks.push_back(l1Track);
	 //incremement the number of matched tracks
	 iMatched++; 
	 
	 // when we get to 3, then create the 3 prong object and break
	 if(iMatched==3){
	   double objectPt  = 0;
	   double objectEta = 0;
	   double objectPhi = 0;
	   for(auto threeProngTrack : threeProngTauTracks){
	     double itrPt= threeProngTrack.getMomentum().perp();
	     double itrEta=threeProngTrack.getMomentum().eta();
	     double itrPhi=threeProngTrack.getMomentum().phi();
	     
	     //Get Weighted Eta Position
	     objectEta = (objectEta*objectPt + itrEta*itrPt)/(objectPt+itrPt);
	     objectPhi = (objectPhi*objectPt + itrPhi*itrPt)/(objectPt+itrPt);
	     
	     //recalculate objectPt 
	     objectPt += itrPt; 
	     //std::cout<<"Found 3 matched tracks, objectPt: "<<objectPt<<std::endl;
	   }
	   tempL1Object.objectPt  = objectPt;
	   tempL1Object.objectEta = objectEta;
	   tempL1Object.objectPhi = objectPhi;
	   tempL1Object.tauDecayMode = 10;


	   //////// now calculate local ecal and hcal
	   double hcalEnergy = 0;
	   double ecalEnergy = 0;
	   TLorentzVector tempTau;
	   double ipt =  objectPt;
	   double ieta = objectEta;
	   double iphi = objectPhi;
	   
	   tempTau.SetPtEtaPhiE(ipt, ieta, iphi, ipt);
	   std::vector<TLorentzVector> nearECALTowers;
	   nearECALTowers.clear();
	   for(auto ecalTPG : allEcalTPGs){
	     if(ecalTPG.Pt()>0)
	       if(tempTau.DeltaR(ecalTPG)<0.1){
		 nearECALTowers.push_back(ecalTPG);
	       }
	   }
	   ecalEnergy = trigTools.SumTLorentzPt(nearECALTowers);
	   //find closest hcalTPG
	   std::vector<TLorentzVector> nearHCALTowers;
	   nearHCALTowers.clear();
	   for(auto hcalTPG : allHcalTPGs){
	     if(hcalTPG.Pt()>0)
	       if(tempTau.DeltaR(hcalTPG)<0.1){
		 nearHCALTowers.push_back(hcalTPG);
	       }
	   }
	   hcalEnergy = trigTools.SumTLorentzPt(nearHCALTowers);
	   //std::cout<<"ecalEnergy "<<ecalEnergy<<" hcalEnergy "<<hcalEnergy<<std::endl;
	   tempL1Object.ECALEnergy = ecalEnergy;
	   tempL1Object.HCALEnergy = hcalEnergy;

	   ///////
	   //break; //!!
	 }
	 if(iMatched>3){
	   threeProngTauTracks_Iso.push_back(l1Track);
	 }
       }
     }
     //triggerGeometryTools trigTools;
     tempL1Object.iso = SumL1TrackPt(threeProngTauTracks_Iso);
     l1ThreeProngTaus.push_back(tempL1Object);
   }
   //std::cout<<"here 3"<<std::endl;
   //keep only the three prong tau with the highest pt primary track
   if(l1ThreeProngTaus.size()>1)
     crossCleanThreeProngTaus(l1ThreeProngTaus);

   sort(l1ThreeProngTaus.begin(), l1ThreeProngTaus.end(), compareL1ObjectByPt);


   //Getting ecal and hcal local information
   for(auto l1ThreeProngTau: l1ThreeProngTaus){

   }
   /*
   int i = 0;
   for(auto l1ThreeProngTau: l1ThreeProngTaus){
     std::cout<<"i: "<<i<<std::endl;
     printL1Object(l1ThreeProngTau);
     i++;
   }
   */
   std::vector<l1Object> l1ThreeProngTaus_eta2p1;
   std::vector<l1Object> l1ThreeProngTaus_eta2p4;

   std::vector<l1Object> l1ThreeProngTausIso;
   std::vector<l1Object> l1ThreeProngTausIso_eta2p1;
   std::vector<l1Object> l1ThreeProngTausIso_eta2p4;

   l1ThreeProngTaus_eta2p1.clear();
   l1ThreeProngTaus_eta2p4.clear();
   l1ThreeProngTausIso.clear();
   l1ThreeProngTausIso_eta2p1.clear();
   l1ThreeProngTausIso_eta2p4.clear();
   //std::cout<<"here 4"<<std::endl;
   for(auto l1ThreeProngTau : l1ThreeProngTaus){

     if(l1ThreeProngTau.objectEta<2.1)
       l1ThreeProngTaus_eta2p1.push_back(l1ThreeProngTau);

     if(l1ThreeProngTau.objectEta<2.4)
       l1ThreeProngTaus_eta2p4.push_back(l1ThreeProngTau);

     if(l1ThreeProngTau.iso<10){
       l1ThreeProngTausIso.push_back(l1ThreeProngTau);
       
       if(l1ThreeProngTau.objectEta<2.1)
	 l1ThreeProngTausIso_eta2p1.push_back(l1ThreeProngTau);

       if(l1ThreeProngTau.objectEta<2.4)
	 l1ThreeProngTausIso_eta2p4.push_back(l1ThreeProngTau);

     }//end iso
   }

     //Rates
   if(l1ThreeProngTaus.size()>0)
     l1ThreeProngTau_pt->Fill(l1ThreeProngTaus.at(0).objectPt);



   if(l1ThreeProngTaus_eta2p4.size()>0)
     l1ThreeProngTau_pt_eta2p4->Fill(l1ThreeProngTaus_eta2p4.at(0).objectPt);


   if(l1ThreeProngTaus_eta2p1.size()>0)
     l1ThreeProngTau_pt_eta2p1->Fill(l1ThreeProngTaus_eta2p1.at(0).objectPt);


   if(l1ThreeProngTausIso.size()>0)
     l1ThreeProngTauIso_pt->Fill(l1ThreeProngTausIso.at(0).objectPt);


   if(l1ThreeProngTausIso_eta2p4.size()>0)
     l1ThreeProngTauIso_pt_eta2p4->Fill(l1ThreeProngTausIso_eta2p4.at(0).objectPt);

   if(l1ThreeProngTausIso_eta2p1.size()>0)
     l1ThreeProngTauIso_pt_eta2p1->Fill(l1ThreeProngTausIso_eta2p1.at(0).objectPt);

   //std::cout<<"here 5"<<std::endl;
   //associate the threeprong taus to gen objects
   for(auto GenThreeProngTau: GenThreeProngTaus){
     gen3ProngPt = GenThreeProngTau.p4.Pt();
     gen3ProngEta = GenThreeProngTau.p4.Eta();
     gen3ProngPhi = GenThreeProngTau.p4.Phi();
     gen3ProngDecayMode = GenThreeProngTau.decayMode;

     threeProngTau.HCALEnergy = 0;
     threeProngTau.ECALEnergy = 0;
     threeProngTau.trackPt  = 0;
     threeProngTau.trackEta = 0;
     threeProngTau.trackPhi = 0;
     threeProngTau.iso = 0;
     threeProngTau.tauDecayMode = 0;

     TLorentzVector tempGen;
     tempGen.SetPtEtaPhiE(GenThreeProngTau.p4.Pt(),GenThreeProngTau.p4.Eta(),GenThreeProngTau.p4.Phi(),GenThreeProngTau.p4.Pt());
   //std::cout<<"here 6"<<std::endl;
     for(auto l1ThreeProngTau : l1ThreeProngTaus){
       if(l1ThreeProngTau.objectPt<10)
	 continue;
       TLorentzVector tempThreeProng;
       tempThreeProng.SetPtEtaPhiE(l1ThreeProngTau.objectPt, l1ThreeProngTau.objectEta, l1ThreeProngTau.objectPhi, l1ThreeProngTau.objectPt);

       if(tempThreeProng.DeltaR(tempGen)<0.3){
	 threeProngTau.trackPt  = l1ThreeProngTau.trackPt;
	 threeProngTau.trackEta = l1ThreeProngTau.trackEta;
	 threeProngTau.trackPhi = l1ThreeProngTau.trackPhi;

	 threeProngTau.objectPt  = tempThreeProng.Pt();
	 threeProngTau.objectEta = tempThreeProng.Eta();
	 threeProngTau.objectPhi = tempThreeProng.Phi();

	 threeProngTau.tauDecayMode = 10;
	 threeProngTau.iso = l1ThreeProngTau.iso;
	 //std::cout<<"gen Pt: "<<gen3ProngPt <<"eta: "<< gen3ProngEta<<" phi: "<< gen3ProngPhi<<std::endl;
	 //std::cout<<"matched L1 Pt: "<<threeProngTau.trackPt <<"eta: "<< threeProngTau.trackEta<<" phi: "<< threeProngTau.trackPhi<<std::endl;
	 break;
       }
     }
   //std::cout<<"here 7"<<std::endl;
     threeProngTree->Fill();
   }

   //std::cout<<"here 8"<<std::endl;
   
   //find 3 highest pt tracks and combine them into one object
   /*
     for(auto l1Track : l1Tracks){
       //std::cout<<"l1track"<<std::endl;
       double dR = deltaR( l1Track.getMomentum(), genTau.p4());
       if( dR < 0.5){
	 std::cout<<"found matched track pt: "<<l1Track.getMomentum().perp()<<std::endl;
	 trackPt = l1Track.getMomentum().perp();
	   trackEta = l1Track.getMomentum().eta();
	   trackPhi = l1Track.getMomentum().phi();
	   threeHighTracks.push_back(l1Track);
	   i++;
	 }
	 if(i>2){
	   float tempTauPt = 0;
	   for(auto track : threeHighTracks){
	     tempTauPt += track.getMomentum().perp(); 
	   }
	   l1TauPt = tempTauPt;
	   std::cout<<"found matched 3 prong tau pt: "<<l1TauPt<<std::endl;
	   l1TauEta = threeHighTracks.at(0).getMomentum().eta();
	   l1TauPhi = threeHighTracks.at(0).getMomentum().phi();
	   break;
	 }
       }
     }
   }
   


   //Find and Sort the Single Prong + pi0 Taus


   //for(unsigned int t = 0; t<l1Tracks.size(); t++)
   //std::cout<<"# "<<t <<": "<<l1Tracks.at(t).getMomentum().perp()<<std::endl;
   
   efficiencyTree->Fill();
   */
}

void L1TauStudy::initializEcalTpgs(edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs,std::vector<TLorentzVector> &allEcalTPGs)
{
  triggerGeometryTools trigTools;
  for (size_t i = 0; i < ecalTPGs->size(); ++i) {
    
    int cal_ieta = (*ecalTPGs)[i].id().ieta();
    int cal_iphi = (*ecalTPGs)[i].id().iphi();
    if(cal_iphi==0)
      std::cout<<"cal_phi is 0"<<std::endl;
    if(cal_ieta<-28)
      continue;
    if(cal_ieta>28)
      continue;
    int ieta = trigTools.TPGEtaRange(cal_ieta);
    short zside = (*ecalTPGs)[i].id().zside();
    // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
    // TPG ieta ideal goes from 0-55.
    double LSB = 0.5;
    double et= (*ecalTPGs)[i].compressedEt()*LSB;
    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
    }

    float eta = trigTools.getRecoEta(ieta, zside);
    float phi = trigTools.getRecoPhi(cal_iphi);    
    //if(et>0)
    //std::cout<<"et "<<et<<std::endl;
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    //if(et>5)
    //std::cout<<"Event Display tpg ecal pt() "<<temp.Pt()<< " eta " <<eta << " phi "<< phi <<std::endl;
    allEcalTPGs.push_back(temp);
  }
}  

void L1TauStudy::initializHcalTpgs(edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs,std::vector<TLorentzVector> &allHcalTPGs,const edm::EventSetup& es){

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  triggerGeometryTools trigTools;  
  for (size_t i = 0; i < hcalTPGs->size(); ++i) {
    HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    if(cal_ieta>28)continue; 
    if(cal_ieta<-28)continue; 
    int ieta = trigTools.TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside = tpg.id().zside();
    double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //if(et>0)
    //std::cout<<"HCAL ET "<<et<<std::endl;
    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
    }

    float eta = trigTools.getRecoEta(ieta, zside);
    float phi = trigTools.getRecoPhi(cal_iphi);    
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHcalTPGs.push_back(temp);
  }
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
L1TauStudy::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TauStudy::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
L1TauStudy::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
L1TauStudy::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
L1TauStudy::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
L1TauStudy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TauStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TauStudy);
