// -*- C++ -*-
//
// Package:    WbDumper
// Class:      WbDumper
//
/**\class WbDumper WbDumper.cc WbAnalysis/WbDumper/src/WbDumper.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Schizzi
//         Created:  Thu Jul 18 15:17:47 CEST 2013
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
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/View.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "PhysicsTools/CandUtils/interface/CandCombiner.h"
#include "CommonTools/Utils/interface/MassRangeSelector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "Math/VectorUtil.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class WbDumper : public edm::EDAnalyzer {
   public:
      explicit WbDumper(const edm::ParameterSet&);
      ~WbDumper();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

     std::string lepton_;
     std::string pileupDT_;
     double par_;
     double par2_;
     bool usePartonFlavour_;
     bool pcut_;
     bool useDeltaR_;

  TH2F*     w_first_jet_pt;	// leading jet of any type
  TH2F*     w_first_jet_eta;
  TH2F*     w_first_jet_mass;
  TH2F*     w_second_jet_pt;
  TH2F*     w_second_jet_eta;
  TH2F*     w_second_jet_mass;

  TH2F*     w_first_jet_pt_b;	// leading jet with at least one b jet in the event
  TH2F*     w_first_jet_eta_b;
  TH2F*     w_first_jet_mass_b;
  TH2F*     w_second_jet_pt_b;
  TH2F*     w_second_jet_eta_b;
  TH2F*     w_second_jet_mass_b;

  TH2F*     w_first_jet_pt_bb;	// leading jet with at least one b jet in the event
  TH2F*     w_first_jet_eta_bb;
  TH2F*     w_first_jet_mass_bb;
  TH2F*     w_second_jet_pt_bb;
  TH2F*     w_second_jet_eta_bb;
  TH2F*     w_second_jet_mass_bb;

  TH2F*     w_first_bjet_pt;	// leading b jet
  TH2F*     w_first_bjet_eta;
  TH2F*     w_first_bjet_mass;

  TH2F*     w_single_bjet_pt;	// only 1 b jet
  TH2F*     w_single_bjet_eta;
  TH2F*     w_single_bjet_mass;

  TH2F*     w_second_bjet_pt;
  TH2F*     w_second_bjet_eta;
  TH2F*     w_second_bjet_mass;

  TH2F*     w_mt_wenu;
  TH2F*     w_mt_wmnu;
  TH2F*     w_mass_wenu_blepton;	// at least one b jet in the event
  TH2F*     w_mass_wmnu_blepton;
  TH2F*     w_mass_wenu_blepton_b;	// at least one b jet in the event
  TH2F*     w_mass_wmnu_blepton_b;
  TH2F*     w_mass_wenu_blepton_bb;	// at least one b jet in the event
  TH2F*     w_mass_wmnu_blepton_bb;
  TH2F*     w_mt_wenu_b;	// at least one b jet in the event
  TH2F*     w_mt_wmnu_b;
  TH2F*     w_mt_wenu_bb;	// at least one b jet in the event
  TH2F*     w_mt_wmnu_bb;
  TH2F*     w_delta_wenu;
  TH2F*     w_delta_wenu_b;
  TH2F*     w_delta_wenu_bb;
  TH2F*     w_delta_wenu_2b;
  TH2F*     w_delta_wmnu;
  TH2F*     w_delta_wmnu_b;
  TH2F*     w_delta_wmnu_bb;
  TH2F*     w_delta_wmnu_2b;
  TH2F*     w_deltaR_wenu;
  TH2F*     w_deltaR_wenu_b;
  TH2F*     w_deltaR_wenu_bb;
  TH2F*     w_deltaR_wenu_2b;
  TH2F*     w_deltaR_wmnu;
  TH2F*     w_deltaR_wmnu_b;
  TH2F*     w_deltaR_wmnu_bb;
  TH2F*     w_deltaR_wmnu_2b;

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
WbDumper::WbDumper(const edm::ParameterSet& iConfig) {

   lepton_           = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
   pileupDT_         = iConfig.getUntrackedParameter < std::string > ("pileupDT", "");
   par_              = iConfig.getUntrackedParameter <double> ("JEC", 0);
   par2_             = iConfig.getUntrackedParameter <double> ("JER", 0);
   pcut_             = iConfig.getUntrackedParameter <bool> ("pcut", false);   
   useDeltaR_        = iConfig.getUntrackedParameter <bool> ("useDeltaR", false);
   //now do what ever initialization is needed
   edm::Service < TFileService > fs;


  w_first_jet_pt =      fs->make < TH2F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_first_jet_eta =     fs->make < TH2F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_jet_mass =      fs->make < TH2F > ("w_first_jet_mass",    "w_first_jet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_second_jet_pt =     fs->make < TH2F > ("w_second_jet_pt",   "w_second_jet_pt;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_second_jet_eta =    fs->make < TH2F > ("w_second_jet_eta",  "w_second_jet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_jet_mass =      fs->make < TH2F > ("w_second_jet_mass",    "w_second_jet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_first_jet_pt_b =    fs->make < TH2F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_first_jet_eta_b =   fs->make < TH2F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_jet_mass_b =      fs->make < TH2F > ("w_first_jet_mass_b",    "w_first_jet_mass_b;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_second_jet_pt_b =   fs->make < TH2F > ("w_second_jet_pt_b",  "w_second_jet_pt_b;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_second_jet_eta_b =  fs->make < TH2F > ("w_second_jet_eta_b", "w_second_jet_eta_b;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_jet_mass_b =      fs->make < TH2F > ("w_second_jet_mass_b",    "w_second_jet_mass_b;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_first_jet_pt_bb =    fs->make < TH2F > ("w_first_jet_pt_bb",   "w_first_jet_pt_bb;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_first_jet_eta_bb =   fs->make < TH2F > ("w_first_jet_eta_bb",  "w_first_jet_eta_bb;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_jet_mass_bb =      fs->make < TH2F > ("w_first_jet_mass_bb",    "w_first_jet_mass_bb;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_second_jet_pt_bb =   fs->make < TH2F > ("w_second_jet_pt_bb",  "w_second_jet_pt_bb;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_second_jet_eta_bb =  fs->make < TH2F > ("w_second_jet_eta_bb", "w_second_jet_eta_bb;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_jet_mass_bb =      fs->make < TH2F > ("w_second_jet_mass_bb",    "w_second_jet_mass_bb;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_first_bjet_pt =     fs->make < TH2F > ("w_first_bjet_pt",    "w_first_bjet_pt;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_first_bjet_eta =    fs->make < TH2F > ("w_first_bjet_eta",   "w_first_bjet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_bjet_mass =      fs->make < TH2F > ("w_first_bjet_mass",    "w_first_bjet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_single_bjet_pt =    fs->make < TH2F > ("w_single_bjet_pt",    "w_single_bjet_pt;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_single_bjet_eta =   fs->make < TH2F > ("w_single_bjet_eta",   "w_single_bjet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_single_bjet_mass =      fs->make < TH2F > ("w_single_bjet_mass",    "w_single_bjet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_second_bjet_pt =    fs->make < TH2F > ("w_second_bjet_pt",   "w_second_bjet_pt;P_t [GeV]", 19, 20., 210., 19, 20., 210.);
  w_second_bjet_eta =   fs->make < TH2F > ("w_second_bjet_eta",  "w_second_bjet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_bjet_mass =      fs->make < TH2F > ("w_second_bjet_mass",    "w_second_bjet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_mt_wenu =           fs->make < TH2F > ("w_mt_wenu",         "w_mt_wenu;M_{T} [GeV]", 40, 40., 200., 40, 40., 200.);
  w_mt_wmnu =           fs->make < TH2F > ("w_mt_wmnu",         "w_mt_wmnu;M_{T} [GeV]", 40, 40., 200., 40, 40., 200.);
  w_mass_wenu_blepton =      fs->make < TH2F > ("w_mass_wenu_blepton",    "w_mass_wenu_blepton;M_{T} [GeV]", 25, 0., 250., 25, 0., 250.);
  w_mass_wmnu_blepton =      fs->make < TH2F > ("w_mass_wmnu_blepton",    "w_mass_wmnu_blepton;M_{T} [GeV]", 25, 0., 250., 25, 0., 250.);
  w_mass_wenu_blepton_b =      fs->make < TH2F > ("w_mass_wenu_blepton_b",    "w_mass_wenu_blepton_b;M_{T} [GeV]", 25, 0., 250., 25, 0., 250.);
  w_mass_wmnu_blepton_b =      fs->make < TH2F > ("w_mass_wmnu_blepton_b",    "w_mass_wmnu_blepton_b;M_{T} [GeV]", 25, 0., 250., 25, 0., 250.);
  w_mass_wenu_blepton_bb =      fs->make < TH2F > ("w_mass_wenu_blepton_bb",    "w_mass_wenu_blepton_bb;M_{T} [GeV]", 25, 0., 250., 25, 0., 250.);
  w_mass_wmnu_blepton_bb =      fs->make < TH2F > ("w_mass_wmnu_blepton_bb",    "w_mass_wmnu_blepton_bb;M_{T} [GeV]", 25, 0., 250., 25, 0., 250.);
  w_mt_wenu_b =         fs->make < TH2F > ("w_mt_wenu_b",       "w_mt_wenu_b;M_{T} [GeV]", 40, 40., 200., 40, 40., 200.);
  w_mt_wmnu_b =         fs->make < TH2F > ("w_mt_wmnu_b",       "w_mt_wmnu_b;M_{T} [GeV]", 40, 40., 200., 40, 40., 200.);
  w_mt_wenu_bb =        fs->make < TH2F > ("w_mt_wenu_bb",      "w_mt_wenu_bb;M_{T} [GeV]", 40, 40., 200., 40, 40., 200.);
  w_mt_wmnu_bb =        fs->make < TH2F > ("w_mt_wmnu_bb",      "w_mt_wmnu_bb;M_{T} [GeV]", 40, 40., 200., 40, 40., 200.);
  w_delta_wenu =          fs->make < TH2F > ("w_delta_wenu",     "w_delta_wenu",    20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wenu_b =        fs->make < TH2F > ("w_delta_wenu_b",   "w_delta_wenu_b",  20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wenu_bb =       fs->make < TH2F > ("w_delta_wenu_bb",  "w_delta_wenu_bb", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wenu_2b =       fs->make < TH2F > ("w_delta_wenu_2b",  "w_delta_wenu_2b", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu =          fs->make < TH2F > ("w_delta_wmnu",     "w_delta_wmnu",    20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu_b =        fs->make < TH2F > ("w_delta_wmnu_b",   "w_delta_wmnu_b",  20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu_bb =       fs->make < TH2F > ("w_delta_wmnu_bb",  "w_delta_wmnu_bb", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu_2b =       fs->make < TH2F > ("w_delta_wmnu_2b",  "w_delta_wmnu_2b", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_deltaR_wenu =          fs->make < TH2F > ("w_deltaR_wenu",     "w_deltaR_wenu",    25, 0., 5., 25, 0., 5.);
  w_deltaR_wenu_b =        fs->make < TH2F > ("w_deltaR_wenu_b",   "w_deltaR_wenu_b",  25, 0., 5., 25, 0., 5.);
  w_deltaR_wenu_bb =       fs->make < TH2F > ("w_deltaR_wenu_bb",  "w_deltaR_wenu_bb", 25, 0., 5., 25, 0., 5.);
  w_deltaR_wenu_2b =       fs->make < TH2F > ("w_deltaR_wenu_2b",  "w_deltaR_wenu_2b", 25, 0., 5., 25, 0., 5.);
  w_deltaR_wmnu =          fs->make < TH2F > ("w_deltaR_wmnu",     "w_deltaR_wmnu",    25, 0., 5., 25, 0., 5.);
  w_deltaR_wmnu_b =        fs->make < TH2F > ("w_deltaR_wmnu_b",   "w_deltaR_wmnu_b",  25, 0., 5., 25, 0., 5.);
  w_deltaR_wmnu_bb =       fs->make < TH2F > ("w_deltaR_wmnu_bb",  "w_deltaR_wmnu_bb", 25, 0., 5., 25, 0., 5.);
  w_deltaR_wmnu_2b =       fs->make < TH2F > ("w_deltaR_wmnu_2b",  "w_deltaR_wmnu_2b", 25, 0., 5., 25, 0., 5.);

}

WbDumper::~WbDumper() {

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void WbDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace std;

    edm::Handle <std::vector<math::XYZTLorentzVector>> electrons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> muons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> jets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> bjets;
    edm::Handle <std::vector<double>>   weight;
    edm::Handle <std::vector<double>>   bweight;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_electrons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_muons;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_jets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_jets2;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_bjets;
    edm::Handle <std::vector<math::XYZTLorentzVector>> gen_bjets2;
    edm::Handle <std::vector<double>>   gen_weight;

    string postfix = "";
    if (pileupDT_=="ee_pup" || pileupDT_=="mm_pup") postfix = "Pup";
    if (pileupDT_=="ee_pum" || pileupDT_=="mm_pum") postfix = "Pum";
    if (par_==+1) postfix = "Up";
    if (par_==-1) postfix = "Down";
    if (par2_==+1) postfix = "JerUp";
    if (par2_==-1) postfix = "JerDown";
    if (pcut_) postfix = "Pur";
    if (useDeltaR_) postfix = "DR";

   if (lepton_== "electron") {

     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myEventWeight"), weight);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myJets"), jets);

     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myBJets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoEle"+postfix,"myBJetsWeights"), bweight);

     iEvent.getByLabel (edm::InputTag("demoEleGen","myEventWeight"), gen_weight);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myElectrons"), gen_electrons);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myMuons"), gen_muons);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myJets"), gen_jets);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myJets2"), gen_jets2);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myBJets"), gen_bjets);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myBJets2"), gen_bjets2);
   }

   if (lepton_== "muon") {

     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myEventWeight"), weight);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myElectrons"), electrons);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myMuons"), muons);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myJets"), jets);

     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myBJets"), bjets);
     iEvent.getByLabel (edm::InputTag("demoMuo"+postfix,"myBJetsWeights"), bweight);

     iEvent.getByLabel (edm::InputTag("demoMuoGen","myEventWeight"), gen_weight);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myElectrons"), gen_electrons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myMuons"), gen_muons);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myJets"), gen_jets);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myJets2"), gen_jets2);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myBJets"), gen_bjets);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myBJets2"), gen_bjets2);
   }

   //match leading jet
   int k=-1;
   if (jets->size()>0) {
     double R = 0.5;
     for (unsigned int i=0; i<gen_jets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*jets)[0], (*gen_jets2)[i]) < R) {
         k=i;
         R = ROOT::Math::VectorUtil::DeltaR((*jets)[0], (*gen_jets2)[i]);
       }
     }
   }
   int k_b=-1;
   if (bjets->size()>0) {
     double R_b = 0.5;
     for (unsigned int i=0; i<gen_bjets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*bjets)[0], (*gen_bjets2)[i]) < R_b) {
         k_b=i;
         R_b = ROOT::Math::VectorUtil::DeltaR((*bjets)[0], (*gen_bjets2)[i]);
       }
     }
   }

   if (k!=-1) {
     if ((*gen_jets2)[k].pt() < 25) k=-1;
     if (fabs((*gen_jets2)[k].eta()) > 2.4) k=-1;
   }
   if (k_b!=-1) {
     if ((*gen_bjets2)[k_b].pt() < 25) k_b=-1;
     if (fabs((*gen_bjets2)[k_b].eta()) > 2.4) k_b=-1;
   }

   // match second jet
   int k2=-1;
   if (jets->size()>1) {
     double R2 = 0.5;
     for (unsigned int i=0; i<gen_jets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*jets)[1], (*gen_jets2)[i]) < R2) {
         k2=i;
         R2 = ROOT::Math::VectorUtil::DeltaR((*jets)[1], (*gen_jets2)[i]);
       }
     }
   }
   int k2_b=-1;
   if (bjets->size()>1) {
     double R2_b = 0.5;
     for (unsigned int i=0; i<gen_bjets2->size(); ++i) {
       if (ROOT::Math::VectorUtil::DeltaR((*bjets)[1], (*gen_bjets2)[i]) < R2_b) {
         k2_b=i;
         R2_b = ROOT::Math::VectorUtil::DeltaR((*bjets)[1], (*gen_bjets2)[i]);
       }
     }
   }

   if (k2!=-1) {
     if ((*gen_jets2)[k2].pt() < 25) k2=-1;
     if (fabs((*gen_jets2)[k2].eta()) > 2.4) k2=-1;
   }
   if (k2_b!=-1) {
     if ((*gen_bjets2)[k2_b].pt() < 25) k2_b=-1;
     if (fabs((*gen_bjets2)[k2_b].eta()) > 2.4) k2_b=-1;
   }

   double my_weight = weight->empty() ? ( gen_weight->empty() ? -1 : (*gen_weight)[0] ) : (*weight)[0];
   double my_bweight = my_weight * ( bweight->empty() ? 1 : (*bweight)[0] );

   if (my_weight>0) {

     bool wenu_event = lepton_ == "electron" && (electrons->size() != 0 || gen_electrons->size() != 0);
     bool wmnu_event = lepton_ == "muon" && (muons->size() != 0 || gen_muons->size() != 0);

     if (wenu_event || wmnu_event) {
       w_first_jet_pt->Fill(jets->empty() ? -1 : (*jets)[0].pt(), k<0 ? -1 : (*gen_jets2)[k].pt(), my_weight);
       w_first_jet_eta->Fill(jets->empty() ? -3 : (*jets)[0].eta(), k<0 ? -3 : (*gen_jets2)[k].eta(), my_weight);
       w_first_jet_mass->Fill(jets->empty() ? -1 : (*jets)[0].mass(), k<0 ? -1 : (*gen_jets2)[k].mass(), my_weight);
       w_first_bjet_pt->Fill(bjets->empty() ? -1 : (*bjets)[0].pt(), k_b<0 ? -1 : (*gen_bjets2)[k_b].pt(), my_bweight);
       w_first_bjet_eta->Fill(bjets->empty() ? -3 : (*bjets)[0].eta(), k_b<0 ? -3 : (*gen_bjets2)[k_b].eta(), my_bweight);
       w_first_bjet_mass->Fill(bjets->empty() ? -1 : (*bjets)[0].mass(), k_b<0 ? -1 : (*gen_bjets2)[k_b].mass(), my_bweight);
       w_second_jet_pt->Fill(jets->size()<2 ? -1 : (*jets)[1].pt(), k2<0 ? -1 : (*gen_jets2)[k2].pt(), my_weight);
       w_second_jet_eta->Fill(jets->size()<2 ? -3 : (*jets)[1].eta(), k2<0 ? -3 : (*gen_jets2)[k2].eta(), my_weight);
       w_second_jet_mass->Fill(jets->size()<2 ? -1 : (*jets)[1].mass(), k2<0 ? -1 : (*gen_jets2)[k2].mass(), my_weight);
       w_second_bjet_pt->Fill(bjets->size()<2 ? -1 : (*bjets)[1].pt(), k2_b<0 ? -1 : (*gen_bjets2)[k2_b].pt(), my_bweight);
       w_second_bjet_eta->Fill(bjets->size()<2 ? -3 : (*bjets)[1].eta(), k2_b<0 ? -3 : (*gen_bjets2)[k2_b].eta(), my_bweight);
       w_second_bjet_mass->Fill(bjets->size()<2 ? -1 : (*bjets)[1].mass(), k2_b<0 ? -1 : (*gen_bjets2)[k2_b].mass(), my_bweight);

       w_first_jet_pt_b->Fill(bjets->size()!=1 ? -1 : (*jets)[0].pt(), k<0 ? -1 : (*gen_jets2)[k].pt(), my_weight);
       w_first_jet_eta_b->Fill(bjets->size()!=1 ? -3 : (*jets)[0].eta(), k<0 ? -3 : (*gen_jets2)[k].eta(), my_weight);
       w_first_jet_mass_b->Fill(bjets->size()!=1 ? -1 : (*jets)[0].mass(), k<0 ? -1 : (*gen_jets2)[k].mass(), my_weight);
       w_second_jet_pt_b->Fill(bjets->size()!=1 ? -1 : (*jets)[1].pt(), k2<0 ? -1 : (*gen_jets2)[k2].pt(), my_weight);
       w_second_jet_eta_b->Fill(bjets->size()!=1 ? -3 : (*jets)[1].eta(), k2<0 ? -3 : (*gen_jets2)[k2].eta(), my_weight);
       w_second_jet_mass_b->Fill(bjets->size()!=1 ? -1 : (*jets)[1].mass(), k2<0 ? -1 : (*gen_jets2)[k2].mass(), my_weight);

       w_first_jet_pt_bb->Fill(bjets->size()<2 ? -1 : (*jets)[0].pt(), k<0 ? -1 : (*gen_jets2)[k].pt(), my_weight);
       w_first_jet_eta_bb->Fill(bjets->size()<2 ? -3 : (*jets)[0].eta(), k<0 ? -3 : (*gen_jets2)[k].eta(), my_weight);
       w_first_jet_mass_bb->Fill(bjets->size()<2 ? -1 : (*jets)[0].mass(), k<0 ? -1 : (*gen_jets2)[k].mass(), my_weight);
       w_second_jet_pt_bb->Fill(bjets->size()<2 ? -1 : (*jets)[1].pt(), k2<0 ? -1 : (*gen_jets2)[k2].pt(), my_weight);
       w_second_jet_eta_bb->Fill(bjets->size()<2 ? -3 : (*jets)[1].eta(), k2<0 ? -3 : (*gen_jets2)[k2].eta(), my_weight);
       w_second_jet_mass_bb->Fill(bjets->size()<2 ? -1 : (*jets)[1].mass(), k2<0 ? -1 : (*gen_jets2)[k2].mass(), my_weight);
     }

/*
     if ((ee_event || mm_event) && numB_!=1) {
       w_delta_phi_2b->Fill(delta_phi_bb->empty() ? -1 : (*delta_phi_bb)[0], gen_delta_phi_bb->empty() ? -1 : (*gen_delta_phi_bb)[0], my_bweight);       
     }
     
     if ((ee_event || mm_event) && numB_==2) {
       w_DR_bb->Fill(DR_bb->empty() ? -1 : (*DR_bb)[0], gen_DR_bb->empty() ? -1 : (*gen_DR_bb)[0], my_bweight);
     }
 
     if (ee_event || mm_event) {
       w_Ht->Fill(Ht->empty() ? -1 : (*Ht)[0], gen_Ht->empty() ? -1 : (*gen_Ht)[0], my_weight);
       w_Ht_b->Fill(Ht_b->empty() ? -1 : (*Ht_b)[0], gen_Ht_b->empty() ? -1 : (*gen_Ht_b)[0], my_bweight);
     }

//     if (ee_event) {
//       w_delta_ee->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi->empty() ? -1 : (*gen_delta_phi)[0], my_weight);
//       w_delta_ee_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi->empty() ? -1 : (*gen_bdelta_phi)[0], my_bweight);
//     }

     if (ee_event) {
       math::XYZTLorentzVector gen_z;
       if (gen_electrons->size() != 0) gen_z = (*gen_electrons)[0] + (*gen_electrons)[1];
       double gen_delta_phi = (gen_electrons->empty() || k<0) ? -1 : fabs(gen_z.phi() - (*gen_jets2)[k].phi());
       double gen_bdelta_phi = (gen_electrons->empty() || k_b<0) ? -1 : fabs(gen_z.phi() - (*gen_bjets2)[k_b].phi());
       if (gen_delta_phi > acos (-1)) gen_delta_phi = 2 * acos (-1) - gen_delta_phi;
       if (gen_bdelta_phi > acos (-1)) gen_bdelta_phi = 2 * acos (-1) - gen_bdelta_phi;
       w_delta_ee->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi, my_weight);
       w_delta_ee_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi, my_bweight);
     }

     if (mm_event) {
       math::XYZTLorentzVector gen_z;
       if (gen_muons->size() != 0) gen_z = (*gen_muons)[0] + (*gen_muons)[1];
       double gen_delta_phi = (gen_muons->empty() || k<0) ? -1 : fabs(gen_z.phi() - (*gen_jets2)[k].phi());
       double gen_bdelta_phi = (gen_muons->empty() || k_b<0) ? -1 : fabs(gen_z.phi() - (*gen_bjets2)[k_b].phi());
       if (gen_delta_phi > acos (-1)) gen_delta_phi = 2 * acos (-1) - gen_delta_phi;
       if (gen_bdelta_phi > acos (-1)) gen_bdelta_phi = 2 * acos (-1) - gen_bdelta_phi;
       w_delta_mm->Fill(delta_phi->empty() ? -1 : (*delta_phi)[0], gen_delta_phi, my_weight);
       w_delta_mm_b->Fill(bdelta_phi->empty() ? -1 : (*bdelta_phi)[0], gen_bdelta_phi, my_bweight);
     }
*/

   }
}


// ------------ method called once each job just before starting event loop  ------------
void WbDumper::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void WbDumper::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void WbDumper::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void WbDumper::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void WbDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void WbDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(WbDumper);
