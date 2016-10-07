/*
 * \file L1TEventDisplayGenerator.cc
 *
 * \author I. Ojalvo
 * Written for Gen
 */
#include <TLorentzVector.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/deltaR.h"

//#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
//#include "L1Trigger/Stage3Ntuplizer/plugins/UCTRegionProcess.hh"
#include "SLHCUpgradeSimulations/L1TauStudy/interface/L1TEventDisplayGenerator.h"
#include "SLHCUpgradeSimulations/L1TauStudy/interface/triggerGeometryTools.hh"

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

bool compareByPtLorentz (TLorentzVector i,TLorentzVector j) { 
  return(i.Pt() > j.Pt()); 
	 };


L1TEventDisplayGenerator::L1TEventDisplayGenerator( const ParameterSet & cfg ) :
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
  crystalSrc_(consumes<l1extra::L1EmParticleCollection>(cfg.getParameter<edm::InputTag>("emCrystals"))),
  //ecalSrc_ ((cfg.getParameter<edm::InputTag>( "ecalDigis"   ))),
  //hcalSrc_ ((cfg.getParameter<edm::InputTag>( "hcalDigis"   ))),
  //vtxLabel_((cfg.getParameter<edm::InputTag>( "vertices"    ))),
  genSrc_  ((cfg.getParameter<edm::InputTag>( "genParticles"))),
  genMatchDeltaRcut(cfg.getUntrackedParameter<double>("genMatchDeltaRcut", 0.1))
  //  tauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("recoTau"))),
  //rctSource_(cfg.getParameter<edm::InputTag>("rctSource")),
  //packedPfCandsToken_(consumes<vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedPfCands"))),
  //pfCandsToken_(consumes<vector<reco::PFCandidate> >(cfg.getParameter<edm::InputTag>("pfCands"))),
  //discriminatorMu_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorMu"))),
  //discriminatorIso_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorIso"))),
  //tauSrc_(consumes<vector<reco::PFTau> >(cfg.getParameter<edm::InputTag>("recoTaus"))),
  //slimmedTauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("slimmedTaus"))),
  //jetSrc_(consumes<vector<pat::Jet> >(cfg.getParameter<edm::InputTag>("recoJets"))),
  //l1ExtraIsoTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraIsoTauSource"))),
  //l1ExtraTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraTauSource"))),
  //l1ExtraJetSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraJetSource"))),
  //regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion"))),
  //ecalCaloSrc_(consumes<vector <reco::CaloCluster> >(cfg.getParameter<edm::InputTag>("ecalCaloClusters")))
  {
    gROOT->ProcessLine(".L /cms/ojalvo/triggerPhaseII/CMSSW_6_2_0_SLHC12_patch1/src/SLHCUpgradeSimulations/L1TauStudy/test/loader.C+");
    L1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
    L1TrackPrimaryVertexTag = cfg.getParameter<edm::InputTag>("L1TrackPrimaryVertexTag");

    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    efficiencyTree = tfs_->make<TTree>("EfficiencyTree", "Efficiency Tree");

    efficiencyTree->Branch("hcalTpgs_Pt",  &hcalTpgs_Pt); 
    efficiencyTree->Branch("hcalTpgs_Eta", &hcalTpgs_Eta); 
    efficiencyTree->Branch("hcalTpgs_Phi", &hcalTpgs_Phi); 

    efficiencyTree->Branch("ecalTpgs_Pt",  &ecalTpgs_Pt); 
    efficiencyTree->Branch("ecalTpgs_Eta", &ecalTpgs_Eta); 
    efficiencyTree->Branch("ecalTpgs_Phi", &ecalTpgs_Phi); 

    efficiencyTree->Branch("sumTpgs_Pt",  &sumTpgs_Pt); 
    efficiencyTree->Branch("sumTpgs_Eta", &sumTpgs_Eta); 
    efficiencyTree->Branch("sumTpgs_Phi", &sumTpgs_Phi); 

    //putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches

    //efficiencyTree->Branch("rlxTaus", "vector<TLorentzVector>", &rlxTaus, 32000, 0); 
    //efficiencyTree->Branch("isoTaus", "vector<TLorentzVector>", &isoTaus, 32000, 0); 
    //efficiencyTree->Branch("recoTaus", "vector<TLorentzVector>", &recoTaus, 32000, 0); 
    //efficiencyTree->Branch("allRegions", "vector<TLorentzVector>", &allRegions, 32000, 0); 
    //efficiencyTree->Branch("hcalTPGs", "ROOT.std.vector(ROOT.TLorentzVector())()", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("hcalTPGs", "std::vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("ecalTPGs", "std::vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 
    //efficiencyTree->Branch("signalPFCands", "vector<TLorentzVector>", &signalPFCands, 32000, 0); 
    //efficiencyTree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0); 
    //efficiencyTree->Branch("recoJets", "vector<TLorentzVector>", &recoJets, 32000, 0); 
    //efficiencyTree->Branch("recoJetsDR", "vector<double>", &recoJetsDr, 32000, 0); 
    //efficiencyTree->Branch("caloClusters", "vector<TLorentzVector>", &caloClusters, 32000, 0); 

    efficiencyTree->Branch("l1Tracks", "vector<TLorentzVector>", &l1Tracks, 32000, 0); 
    efficiencyTree->Branch("l1EcalClusters", "vector<TLorentzVector>", &l1EcalClusters, 32000, 0); 
    efficiencyTree->Branch("l1EcalCrystals", "vector<TLorentzVector>", &l1EcalCrystals, 32000, 0); 
    efficiencyTree->Branch("genTaus", "vector<TLorentzVector>", &genHadronicTaus, 32000, 0); 

    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",   &nvtx,         "nvtx/I");

    efficiencyTree->Branch("decayMode", &decayMode,   "decayMode/I");
    
    efficiencyTree->Branch("tauEtaEcalEnt", &tauEtaEcalEnt,"tauEtaEcalEnt/D");
    efficiencyTree->Branch("tauPhiEcalEnt", &tauPhiEcalEnt,"tauPhiEcalEnt/D");

    efficiencyTree->Branch("recoPt",        &recoPt,   "recoPt/D");
    efficiencyTree->Branch("recoEta",       &recoEta,   "recoEta/D");
    efficiencyTree->Branch("recoPhi",       &recoPhi,   "recoPhi/D");
    
  }

void L1TEventDisplayGenerator::beginJob( const EventSetup & es) {
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <TLorentzVector.h>");

}

//unsigned int const L1TEventDisplayGenerator::N_TOWER_PHI = 72;
//unsigned int const L1TEventDisplayGenerator::N_TOWER_ETA = 56;

void L1TEventDisplayGenerator::analyze( const Event& evt, const EventSetup& es )
{
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <TLorentzVector.h>");
  
  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();
   // L1 Tracks
   edm::Handle<L1TkTrackCollectionType> l1trackHandle;
   evt.getByLabel(L1TrackInputTag, l1trackHandle);

   // L1 Track based primary vertex
   edm::Handle<L1TkPrimaryVertexCollection> l1PrimaryVertexHandle;
   evt.getByLabel(L1TrackPrimaryVertexTag, l1PrimaryVertexHandle);

   edm::Handle<l1extra::L1EmParticleCollection> l1EGCrystalHandle;
   evt.getByToken(crystalSrc_, l1EGCrystalHandle);

   // Get genParticles
   edm::Handle<GenParticleCollectionType> genParticleHandle;
   if(!evt.getByLabel(genSrc_,genParticleHandle))
     std::cout<<"No gen Particles Found "<<std::endl;

   vector<TTTrack<Ref_PixelDigi_> > l1TracksRef;
   
   edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
   edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  
   
   //Clear the vectors
   l1Tracks->clear(); 
   l1EcalClusters->clear(); 
   l1EcalCrystals->clear(); 
   genHadronicTaus->clear(); 

   vector<reco::GenParticle> genTaus;
   vector<reco::GenParticle> genParticles;
   reco::GenParticle* genTau;
   for(unsigned int i = 0; i< genParticleHandle->size(); i++){
     edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
     genParticles.push_back(*ptr);
     if(abs(ptr->pdgId())==15){
       genTaus.push_back(*ptr);
       //const reco::GenParticle *genTau = &(*ptr);
  
     }
   }

   for(auto genTau: genTaus){
     std::cout<<"got tau"<<std::endl;

     std::cout<<"Tau Decay Mode "<<GetDecayMode(&genTau)<<std::endl;
     decayMode = GetDecayMode(&genTau);
     //onlygetting the hadronic taus
     if(decayMode<10)
       continue;

     reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &genParticles);

     TLorentzVector genTauTemp;
     float et, eta, phi;
     et  = visGenTau.pt();
     eta = visGenTau.eta();
     phi = visGenTau.phi();
     genTauTemp.SetPtEtaPhiE(et,eta,phi,et);
     genHadronicTaus->push_back(genTauTemp);

   }


   
   //Find and sort the tracks
   for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
       edm::Ptr<TTTrack<Ref_PixelDigi_>> ptr(l1trackHandle, track_index);
       //std::cout<<"# "<<track_index<<": "<<pt<<std::endl;
       TLorentzVector l1TrackTemp;
       float et, eta, phi;
       et  = ptr->getMomentum().perp();
       eta = ptr->getMomentum().eta();
       phi = ptr->getMomentum().phi();
       l1TrackTemp.SetPtEtaPhiE(et,eta,phi,et);
       l1Tracks->push_back(l1TrackTemp);
     }
   //std::sort(l1Tracks->begin(),l1Tracks->end(),compareByPtLorentz);   

   //Find and sort the tracks
   for(size_t crystal_index=0; crystal_index< l1EGCrystalHandle->size(); ++crystal_index)
     {
       edm::Ptr<l1extra::L1EmParticle> ptr(l1EGCrystalHandle, crystal_index);
       //std::cout<<"# "<<crystal_index<<": "<<pt<<std::endl;
       TLorentzVector l1CrystalTemp;
       float et, eta, phi;
       et  = ptr->pt();
       eta = ptr->eta();
       phi = ptr->phi();
       l1CrystalTemp.SetPtEtaPhiE(et,eta,phi,et);
       l1EcalCrystals->push_back(l1CrystalTemp);
       //std::cout<<"crystal pt:"<<et<<" eta: "<<eta<<" phi: "<<phi<<std::endl;
     }
   std::sort(l1EcalCrystals->begin(),l1EcalCrystals->end(),compareByPtLorentz);   


   // electron candidate extra info from Sacha's algorithm
   //l1slhc::L1EGCrystalClusterCollection crystalClusters;
   //edm::Handle<l1slhc::L1EGCrystalClusterCollection> crystalClustersHandle;      
   //evt.getByLabel(L1CrystalClustersInputTag,crystalClustersHandle);
   //crystalClusters = (*crystalClustersHandle.product());

   // Sort clusters so we can always pick highest pt cluster matching cuts
   //std::sort(begin(crystalClusters), end(crystalClusters), [](const l1slhc::L1EGCrystalCluster& a, const l1slhc::L1EGCrystalCluster& b){return a.pt() > b.pt();});

   //if(!evt.getByToken(hcalSrc_, hcalTPGs))
   //std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
   
   // loop over jets
   
   


  if(!evt.getByToken(ecalSrc_, ecalTPGs))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < ecalTPGs->size(); ++i) {

      int cal_ieta = (*ecalTPGs)[i].id().ieta();
      int cal_iphi = (*ecalTPGs)[i].id().iphi();
      if(cal_iphi==0)
	std::cout<<"cal_phi is 0"<<std::endl;
      if(cal_ieta<-28)
	continue;
      if(cal_ieta>28)
	continue;
      int ieta = TPGEtaRange(cal_ieta);
      short zside = (*ecalTPGs)[i].id().zside();
      // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
      // TPG ieta ideal goes from 0-55.
      double LSB = 0.5;
      double et= (*ecalTPGs)[i].compressedEt()*LSB;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      triggerGeometryTools trigTools;
      float eta = trigTools.getRecoEta(ieta, zside);
      float phi = trigTools.getRecoPhi(cal_iphi);    
      //if(et>0)
      //std::cout<<"et "<<et<<std::endl;
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      //if(et>5)
      //std::cout<<"Event Display tpg ecal pt() "<<temp.Pt()<< " eta " <<eta << " phi "<< phi <<std::endl;
      allEcalTPGs->push_back(temp);
    }

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < hcalTPGs->size(); ++i) {
      HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
      int cal_ieta = tpg.id().ieta();
      int cal_iphi = tpg.id().iphi();
      if(cal_ieta>28)continue; 
      if(cal_ieta<-28)continue; 
      int ieta = TPGEtaRange(cal_ieta);
      short absieta = std::abs(tpg.id().ieta());
      short zside = tpg.id().zside();
      double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
      //if(et>0)
      //std::cout<<"HCAL ET "<<et<<std::endl;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      triggerGeometryTools trigTools;
      float eta = trigTools.getRecoEta(ieta, zside);
      float phi = trigTools.getRecoPhi(cal_iphi);    
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allHcalTPGs->push_back(temp);
    }

  efficiencyTree->Fill();
   /**/
 }




void L1TEventDisplayGenerator::endJob() {
}

L1TEventDisplayGenerator::~L1TEventDisplayGenerator(){
}

DEFINE_FWK_MODULE(L1TEventDisplayGenerator);
