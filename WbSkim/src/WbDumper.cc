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
  TH2F*     w_dijet_pt;
  TH2F*     w_dijet_eta;
  TH2F*     w_dijet_mass;

  TH2F*     w_first_jet_pt_b;	// leading jet with at least one b jet in the event
  TH2F*     w_first_jet_eta_b;
  TH2F*     w_first_jet_mass_b;
  TH2F*     w_second_jet_pt_b;
  TH2F*     w_second_jet_eta_b;
  TH2F*     w_second_jet_mass_b;
  TH2F*     w_dijet_pt_b;
  TH2F*     w_dijet_eta_b;
  TH2F*     w_dijet_mass_b;

  TH2F*     w_first_jet_pt_bb;	// leading jet with at least one b jet in the event
  TH2F*     w_first_jet_eta_bb;
  TH2F*     w_first_jet_mass_bb;
  TH2F*     w_second_jet_pt_bb;
  TH2F*     w_second_jet_eta_bb;
  TH2F*     w_second_jet_mass_bb;
  TH2F*     w_dijet_pt_bb;
  TH2F*     w_dijet_eta_bb;
  TH2F*     w_dijet_mass_bb;

  TH2F*     w_first_bjet_pt;	// leading b jet
  TH2F*     w_first_bjet_eta;
  TH2F*     w_first_bjet_mass;

  TH2F*     w_second_bjet_pt;
  TH2F*     w_second_bjet_eta;
  TH2F*     w_second_bjet_mass;

  TH2F*     w_mt_wenu;
  TH2F*     w_mt_wmnu;
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

  TH2F*     w_pt_W_wenu;
  TH2F*     w_pt_W_wmnu;
  TH2F*     w_pt_W_wenu_b;
  TH2F*     w_pt_W_wmnu_b;
  TH2F*     w_pt_W_wenu_bb;
  TH2F*     w_pt_W_wmnu_bb;
  TH2F*     w_eta_W_wenu;
  TH2F*     w_eta_W_wmnu;
  TH2F*     w_eta_W_wenu_b;
  TH2F*     w_eta_W_wmnu_b;
  TH2F*     w_eta_W_wenu_bb;
  TH2F*     w_eta_W_wmnu_bb;

  TH2F*     w_Ht;
  TH2F*     w_Ht_b;
  TH2F*     w_Ht_bb;
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


  w_first_jet_pt =      fs->make < TH2F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_first_jet_eta =     fs->make < TH2F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_jet_mass =      fs->make < TH2F > ("w_first_jet_mass",    "w_first_jet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_second_jet_pt =     fs->make < TH2F > ("w_second_jet_pt",   "w_second_jet_pt;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_second_jet_eta =    fs->make < TH2F > ("w_second_jet_eta",  "w_second_jet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_jet_mass =      fs->make < TH2F > ("w_second_jet_mass",    "w_second_jet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_dijet_pt =     fs->make < TH2F > ("w_dijet_pt",   "w_dijet_pt;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_dijet_eta =    fs->make < TH2F > ("w_dijet_eta",  "w_dijet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_dijet_mass =      fs->make < TH2F > ("w_dijet_mass",    "w_dijet_mass;Mass [GeV]", 20, 20., 260., 20, 20., 260.);

  w_first_jet_pt_b =    fs->make < TH2F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_first_jet_eta_b =   fs->make < TH2F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_jet_mass_b =      fs->make < TH2F > ("w_first_jet_mass_b",    "w_first_jet_mass_b;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_second_jet_pt_b =   fs->make < TH2F > ("w_second_jet_pt_b",  "w_second_jet_pt_b;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_second_jet_eta_b =  fs->make < TH2F > ("w_second_jet_eta_b", "w_second_jet_eta_b;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_jet_mass_b =      fs->make < TH2F > ("w_second_jet_mass_b",    "w_second_jet_mass_b;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_dijet_pt_b =   fs->make < TH2F > ("w_dijet_pt_b",  "w_dijet_pt_b;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_dijet_eta_b =  fs->make < TH2F > ("w_dijet_eta_b", "w_dijet_eta_b;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_dijet_mass_b =      fs->make < TH2F > ("w_dijet_mass_b",    "w_dijet_mass_b;Mass [GeV]", 20, 20., 260., 20, 20., 260.);

  w_first_jet_pt_bb =    fs->make < TH2F > ("w_first_jet_pt_bb",   "w_first_jet_pt_bb;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_first_jet_eta_bb =   fs->make < TH2F > ("w_first_jet_eta_bb",  "w_first_jet_eta_bb;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_jet_mass_bb =      fs->make < TH2F > ("w_first_jet_mass_bb",    "w_first_jet_mass_bb;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_second_jet_pt_bb =   fs->make < TH2F > ("w_second_jet_pt_bb",  "w_second_jet_pt_bb;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_second_jet_eta_bb =  fs->make < TH2F > ("w_second_jet_eta_bb", "w_second_jet_eta_bb;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_jet_mass_bb =      fs->make < TH2F > ("w_second_jet_mass_bb",    "w_second_jet_mass_bb;Mass [GeV]", 18, 0., 36., 18, 0., 36.);
  w_dijet_pt_bb =   fs->make < TH2F > ("w_dijet_pt_bb",  "w_dijet_pt_bb;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_dijet_eta_bb =  fs->make < TH2F > ("w_dijet_eta_bb", "w_dijet_eta_bb;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_dijet_mass_bb =      fs->make < TH2F > ("w_dijet_mass_bb",    "w_dijet_mass_bb;Mass [GeV]", 20, 20., 260., 20, 20., 260.);

  w_first_bjet_pt =     fs->make < TH2F > ("w_first_bjet_pt",    "w_first_bjet_pt;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_first_bjet_eta =    fs->make < TH2F > ("w_first_bjet_eta",   "w_first_bjet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_first_bjet_mass =      fs->make < TH2F > ("w_first_bjet_mass",    "w_first_bjet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_second_bjet_pt =    fs->make < TH2F > ("w_second_bjet_pt",   "w_second_bjet_pt;P_t [GeV]", 20, 20., 220., 20, 20., 220.);
  w_second_bjet_eta =   fs->make < TH2F > ("w_second_bjet_eta",  "w_second_bjet_eta;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_second_bjet_mass =      fs->make < TH2F > ("w_second_bjet_mass",    "w_second_bjet_mass;Mass [GeV]", 18, 0., 36., 18, 0., 36.);

  w_mt_wenu =           fs->make < TH2F > ("w_mt_wenu",         "w_mt_wenu;M_{T} [GeV]", 20, 45., 205., 20, 45., 205.);
  w_mt_wmnu =           fs->make < TH2F > ("w_mt_wmnu",         "w_mt_wmnu;M_{T} [GeV]", 20, 45., 205., 20, 45., 205.);
  w_mt_wenu_b =         fs->make < TH2F > ("w_mt_wenu_b",       "w_mt_wenu_b;M_{T} [GeV]", 20, 45., 205., 20, 45., 205.);
  w_mt_wmnu_b =         fs->make < TH2F > ("w_mt_wmnu_b",       "w_mt_wmnu_b;M_{T} [GeV]", 20, 45., 205., 20, 45., 205.);
  w_mt_wenu_bb =        fs->make < TH2F > ("w_mt_wenu_bb",      "w_mt_wenu_bb;M_{T} [GeV]", 20, 45., 205., 20, 45., 205.);
  w_mt_wmnu_bb =        fs->make < TH2F > ("w_mt_wmnu_bb",      "w_mt_wmnu_bb;M_{T} [GeV]", 20, 45., 205., 20, 45., 205.);
  w_delta_wenu =          fs->make < TH2F > ("w_delta_wenu",     "w_delta_wenu",    20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wenu_b =        fs->make < TH2F > ("w_delta_wenu_b",   "w_delta_wenu_b",  20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wenu_bb =       fs->make < TH2F > ("w_delta_wenu_bb",  "w_delta_wenu_bb", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wenu_2b =       fs->make < TH2F > ("w_delta_wenu_2b",  "w_delta_wenu_2b", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu =          fs->make < TH2F > ("w_delta_wmnu",     "w_delta_wmnu",    20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu_b =        fs->make < TH2F > ("w_delta_wmnu_b",   "w_delta_wmnu_b",  20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu_bb =       fs->make < TH2F > ("w_delta_wmnu_bb",  "w_delta_wmnu_bb", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_delta_wmnu_2b =       fs->make < TH2F > ("w_delta_wmnu_2b",  "w_delta_wmnu_2b", 20, 0., TMath::Pi (),  20, 0., TMath::Pi ());
  w_deltaR_wenu =          fs->make < TH2F > ("w_deltaR_wenu",     "w_deltaR_wenu",    24, 0., 4.8, 24, 0., 4.8);
  w_deltaR_wenu_b =        fs->make < TH2F > ("w_deltaR_wenu_b",   "w_deltaR_wenu_b",  24, 0., 4.8, 24, 0., 4.8);
  w_deltaR_wenu_bb =       fs->make < TH2F > ("w_deltaR_wenu_bb",  "w_deltaR_wenu_bb", 24, 0., 4.8, 24, 0., 4.8);
  w_deltaR_wenu_2b =       fs->make < TH2F > ("w_deltaR_wenu_2b",  "w_deltaR_wenu_2b", 24, 0., 4.8, 24, 0., 4.8);
  w_deltaR_wmnu =          fs->make < TH2F > ("w_deltaR_wmnu",     "w_deltaR_wmnu",    24, 0., 4.8, 24, 0., 4.8);
  w_deltaR_wmnu_b =        fs->make < TH2F > ("w_deltaR_wmnu_b",   "w_deltaR_wmnu_b",  24, 0., 4.8, 24, 0., 4.8);
  w_deltaR_wmnu_bb =       fs->make < TH2F > ("w_deltaR_wmnu_bb",  "w_deltaR_wmnu_bb", 24, 0., 4.8, 24, 0., 4.8);
  w_deltaR_wmnu_2b =       fs->make < TH2F > ("w_deltaR_wmnu_2b",  "w_deltaR_wmnu_2b", 24, 0., 4.8, 24, 0., 4.8);

  w_pt_W_wenu =           fs->make < TH2F > ("w_pt_W_wenu",         "w_pt_W_wenu;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_pt_W_wmnu =           fs->make < TH2F > ("w_pt_W_wmnu",         "w_pt_W_wmnu;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_pt_W_wenu_b =         fs->make < TH2F > ("w_pt_W_wenu_b",       "w_pt_W_wenu_b;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_pt_W_wmnu_b =         fs->make < TH2F > ("w_pt_W_wmnu_b",       "w_pt_W_wmnu_b;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_pt_W_wenu_bb =        fs->make < TH2F > ("w_pt_W_wenu_bb",      "w_pt_W_wenu_bb;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_pt_W_wmnu_bb =        fs->make < TH2F > ("w_pt_W_wmnu_bb",      "w_pt_W_wmnu_bb;P_t [GeV]", 20, 0., 200., 20, 0., 200.);
  w_eta_W_wenu =           fs->make < TH2F > ("w_eta_W_wenu",         "w_eta_W_wenu;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_eta_W_wmnu =           fs->make < TH2F > ("w_eta_W_wmnu",         "w_eta_W_wmnu;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_eta_W_wenu_b =         fs->make < TH2F > ("w_eta_W_wenu_b",       "w_eta_W_wenu_b;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_eta_W_wmnu_b =         fs->make < TH2F > ("w_eta_W_wmnu_b",       "w_eta_W_wmnu_b;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_eta_W_wenu_bb =        fs->make < TH2F > ("w_eta_W_wenu_bb",      "w_eta_W_wenu_bb;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);
  w_eta_W_wmnu_bb =        fs->make < TH2F > ("w_eta_W_wmnu_bb",      "w_eta_W_wmnu_bb;Eta", 20, -2.4, 2.4, 20, -2.4, 2.4);

  w_Ht =                fs->make < TH2F > ("w_Ht",              "w_Ht [GeV]", 20, 20., 220., 20, 20., 220.);
  w_Ht_b =              fs->make < TH2F > ("w_Ht_b",            "w_Ht_b [GeV]", 20, 20., 220., 20, 20., 220.);
  w_Ht_bb =             fs->make < TH2F > ("w_Ht_bb",           "w_Ht_bb [GeV]", 20, 20., 220., 20, 20., 220.);
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

    edm::Handle <std::vector<double>>   deltaphi_ej;
    edm::Handle <std::vector<double>>   deltaphi_ebj;
    edm::Handle <std::vector<double>>   deltaphi_ebjbj;
    edm::Handle <std::vector<double>>   deltaphi_mj;
    edm::Handle <std::vector<double>>   deltaphi_mbj;
    edm::Handle <std::vector<double>>   deltaphi_mbjbj;
    edm::Handle <std::vector<double>>   deltaR_ej;
    edm::Handle <std::vector<double>>   deltaR_ebj;
    edm::Handle <std::vector<double>>   deltaR_ebjbj;
    edm::Handle <std::vector<double>>   deltaR_mj;
    edm::Handle <std::vector<double>>   deltaR_mbj;
    edm::Handle <std::vector<double>>   deltaR_mbjbj;
    edm::Handle <std::vector<double>>   pt_W_wenu;
    edm::Handle <std::vector<double>>   pt_W_wmnu;
    edm::Handle <std::vector<double>>   eta_W_wenu;
    edm::Handle <std::vector<double>>   eta_W_wmnu;

    edm::Handle <std::vector<double>>   gen_deltaphi_ej;
    edm::Handle <std::vector<double>>   gen_deltaphi_ebj;
    edm::Handle <std::vector<double>>   gen_deltaphi_ebjbj;
    edm::Handle <std::vector<double>>   gen_deltaphi_mj;
    edm::Handle <std::vector<double>>   gen_deltaphi_mbj;
    edm::Handle <std::vector<double>>   gen_deltaphi_mbjbj;
    edm::Handle <std::vector<double>>   gen_deltaR_ej;
    edm::Handle <std::vector<double>>   gen_deltaR_ebj;
    edm::Handle <std::vector<double>>   gen_deltaR_ebjbj;
    edm::Handle <std::vector<double>>   gen_deltaR_mj;
    edm::Handle <std::vector<double>>   gen_deltaR_mbj;
    edm::Handle <std::vector<double>>   gen_deltaR_mbjbj;
    edm::Handle <std::vector<double>>   gen_pt_W_wenu;
    edm::Handle <std::vector<double>>   gen_pt_W_wmnu;
    edm::Handle <std::vector<double>>   gen_eta_W_wenu;
    edm::Handle <std::vector<double>>   gen_eta_W_wmnu;

    edm::Handle <std::vector<double>>   Ht;
    edm::Handle <std::vector<double>>   gen_Ht;

    edm::Handle <std::vector<double>>   dijet_pt;
    edm::Handle <std::vector<double>>   dijet_eta;
    edm::Handle <std::vector<double>>   dijet_mass;
    edm::Handle <std::vector<double>>   gen_dijet_pt;
    edm::Handle <std::vector<double>>   gen_dijet_eta;
    edm::Handle <std::vector<double>>   gen_dijet_mass;

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

     iEvent.getByLabel (edm::InputTag("demoEle","myDeltaPhiEJ"), deltaphi_ej);
     iEvent.getByLabel (edm::InputTag("demoEle","myDeltaPhiEBJ"), deltaphi_ebj);
     iEvent.getByLabel (edm::InputTag("demoEle","myDeltaPhiEBJBJ"), deltaphi_ebjbj);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDeltaPhiEJ"), gen_deltaphi_ej);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDeltaPhiEBJ"), gen_deltaphi_ebj);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDeltaPhiEBJBJ"), gen_deltaphi_ebjbj);
     iEvent.getByLabel (edm::InputTag("demoEle","myDeltaREJ"), deltaR_ej);
     iEvent.getByLabel (edm::InputTag("demoEle","myDeltaREBJ"), deltaR_ebj);
     iEvent.getByLabel (edm::InputTag("demoEle","myDeltaREBJBJ"), deltaR_ebjbj);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDeltaREJ"), gen_deltaR_ej);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDeltaREBJ"), gen_deltaR_ebj);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDeltaREBJBJ"), gen_deltaR_ebjbj);

     iEvent.getByLabel (edm::InputTag("demoEle","myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myHt"), gen_Ht);

     iEvent.getByLabel (edm::InputTag("demoEle","myWenuPt"), pt_W_wenu);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myWenuPt"), gen_pt_W_wenu);
     iEvent.getByLabel (edm::InputTag("demoEle","myWenuEta"), eta_W_wenu);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myWenuEta"), gen_eta_W_wenu);

     iEvent.getByLabel (edm::InputTag("demoEle","myDijetPt"), dijet_pt);
     iEvent.getByLabel (edm::InputTag("demoEle","myDijetEta"), dijet_eta);
     iEvent.getByLabel (edm::InputTag("demoEle","myDijetMass"), dijet_mass);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDijetPt"), gen_dijet_pt);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDijetEta"), gen_dijet_eta);
     iEvent.getByLabel (edm::InputTag("demoEleGen","myDijetMass"), gen_dijet_mass);
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

     iEvent.getByLabel (edm::InputTag("demoMuo","myDeltaPhiMJ"), deltaphi_mj);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDeltaPhiMBJ"), deltaphi_mbj);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDeltaPhiMBJBJ"), deltaphi_mbjbj);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDeltaPhiMJ"), gen_deltaphi_mj);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDeltaPhiMBJ"), gen_deltaphi_mbj);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDeltaPhiMBJBJ"), gen_deltaphi_mbjbj);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDeltaRMJ"), deltaR_mj);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDeltaRMBJ"), deltaR_mbj);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDeltaRMBJBJ"), deltaR_mbjbj);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDeltaRMJ"), gen_deltaR_mj);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDeltaRMBJ"), gen_deltaR_mbj);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDeltaRMBJBJ"), gen_deltaR_mbjbj);

     iEvent.getByLabel (edm::InputTag("demoMuo","myHt"), Ht);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myHt"), gen_Ht);

     iEvent.getByLabel (edm::InputTag("demoMuo","myWmnuPt"), pt_W_wmnu);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myWmnuPt"), gen_pt_W_wmnu);
     iEvent.getByLabel (edm::InputTag("demoMuo","myWmnuEta"), eta_W_wmnu);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myWmnuEta"), gen_eta_W_wmnu);

     iEvent.getByLabel (edm::InputTag("demoMuo","myDijetPt"), dijet_pt);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDijetEta"), dijet_eta);
     iEvent.getByLabel (edm::InputTag("demoMuo","myDijetMass"), dijet_mass);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDijetPt"), gen_dijet_pt);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDijetEta"), gen_dijet_eta);
     iEvent.getByLabel (edm::InputTag("demoMuoGen","myDijetMass"), gen_dijet_mass);
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

   if (my_weight>0) {

     bool wenu_event = lepton_ == "electron" && (electrons->size() != 0 || gen_electrons->size() != 0);
     bool wmnu_event = lepton_ == "muon" && (muons->size() != 0 || gen_muons->size() != 0);

     if (wenu_event || wmnu_event) {
       w_first_jet_pt->Fill(jets->empty() ? -1 : (*jets)[0].pt(), k<0 ? -1 : (*gen_jets2)[k].pt(), my_weight);
       w_first_jet_eta->Fill(jets->empty() ? -3 : (*jets)[0].eta(), k<0 ? -3 : (*gen_jets2)[k].eta(), my_weight);
       w_first_jet_mass->Fill(jets->empty() ? -1 : (*jets)[0].mass(), k<0 ? -1 : (*gen_jets2)[k].mass(), my_weight);
       w_first_bjet_pt->Fill(bjets->empty() ? -1 : (*bjets)[0].pt(), k_b<0 ? -1 : (*gen_bjets2)[k_b].pt(), my_weight);
       w_first_bjet_eta->Fill(bjets->empty() ? -3 : (*bjets)[0].eta(), k_b<0 ? -3 : (*gen_bjets2)[k_b].eta(), my_weight);
       w_first_bjet_mass->Fill(bjets->empty() ? -1 : (*bjets)[0].mass(), k_b<0 ? -1 : (*gen_bjets2)[k_b].mass(), my_weight);
       w_second_jet_pt->Fill(jets->size()<2 ? -1 : (*jets)[1].pt(), k2<0 ? -1 : (*gen_jets2)[k2].pt(), my_weight);
       w_second_jet_eta->Fill(jets->size()<2 ? -3 : (*jets)[1].eta(), k2<0 ? -3 : (*gen_jets2)[k2].eta(), my_weight);
       w_second_jet_mass->Fill(jets->size()<2 ? -1 : (*jets)[1].mass(), k2<0 ? -1 : (*gen_jets2)[k2].mass(), my_weight);
       w_second_bjet_pt->Fill(bjets->size()<2 ? -1 : (*bjets)[1].pt(), k2_b<0 ? -1 : (*gen_bjets2)[k2_b].pt(), my_weight);
       w_second_bjet_eta->Fill(bjets->size()<2 ? -3 : (*bjets)[1].eta(), k2_b<0 ? -3 : (*gen_bjets2)[k2_b].eta(), my_weight);
       w_second_bjet_mass->Fill(bjets->size()<2 ? -1 : (*bjets)[1].mass(), k2_b<0 ? -1 : (*gen_bjets2)[k2_b].mass(), my_weight);

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

     if (wenu_event || wmnu_event) {
       w_Ht->Fill(Ht->empty() ? -1 : (*Ht)[0], gen_Ht->empty() ? -1 : (*gen_Ht)[0], my_weight);
       w_Ht_b->Fill(Ht->empty() ? -1 : (*Ht)[0], gen_Ht->empty() ? -1 : (*gen_Ht)[0], my_weight);
       w_Ht_bb->Fill(Ht->empty() ? -1 : (*Ht)[0], gen_Ht->empty() ? -1 : (*gen_Ht)[0], my_weight);
       w_dijet_pt->Fill(dijet_pt->empty() ? -1 : (*dijet_pt)[0], gen_dijet_pt->empty() ? -1 : (*gen_dijet_pt)[0], my_weight);
       w_dijet_eta->Fill(dijet_eta->empty() ? -1 : (*dijet_eta)[0], gen_dijet_eta->empty() ? -1 : (*gen_dijet_eta)[0], my_weight);
       w_dijet_mass->Fill(dijet_mass->empty() ? -1 : (*dijet_mass)[0], gen_dijet_mass->empty() ? -1 : (*gen_dijet_mass)[0], my_weight);
       w_dijet_pt_b->Fill(dijet_pt->empty() ? -1 : (*dijet_pt)[0], gen_dijet_pt->empty() ? -1 : (*gen_dijet_pt)[0], my_weight);
       w_dijet_eta_b->Fill(dijet_eta->empty() ? -1 : (*dijet_eta)[0], gen_dijet_eta->empty() ? -1 : (*gen_dijet_eta)[0], my_weight);
       w_dijet_mass_b->Fill(dijet_mass->empty() ? -1 : (*dijet_mass)[0], gen_dijet_mass->empty() ? -1 : (*gen_dijet_mass)[0], my_weight);
       w_dijet_pt_bb->Fill(dijet_pt->empty() ? -1 : (*dijet_pt)[0], gen_dijet_pt->empty() ? -1 : (*gen_dijet_pt)[0], my_weight);
       w_dijet_eta_bb->Fill(dijet_eta->empty() ? -1 : (*dijet_eta)[0], gen_dijet_eta->empty() ? -1 : (*gen_dijet_eta)[0], my_weight);
       w_dijet_mass_bb->Fill(dijet_mass->empty() ? -1 : (*dijet_mass)[0], gen_dijet_mass->empty() ? -1 : (*gen_dijet_mass)[0], my_weight);
     }

     if (wenu_event) {
       w_delta_wenu->Fill(deltaphi_ej->empty() ? -1 : (*deltaphi_ej)[0], gen_deltaphi_ej->empty() ? -1 : (*gen_deltaphi_ej)[0], my_weight);
       w_delta_wenu_b->Fill(deltaphi_ebj->empty() ? -1 : (*deltaphi_ebj)[0], gen_deltaphi_ebj->empty() ? -1 : (*gen_deltaphi_ebj)[0], my_weight);
       w_delta_wenu_bb->Fill(deltaphi_ebj->empty() ? -1 : (*deltaphi_ebj)[0], gen_deltaphi_ebj->empty() ? -1 : (*gen_deltaphi_ebj)[0], my_weight);
       w_delta_wenu_2b->Fill(deltaphi_ebjbj->empty() ? -1 : (*deltaphi_ebjbj)[0], gen_deltaphi_ebjbj->empty() ? -1 : (*gen_deltaphi_ebjbj)[0], my_weight);
       w_deltaR_wenu->Fill(deltaR_ej->empty() ? -1 : (*deltaR_ej)[0], gen_deltaR_ej->empty() ? -1 : (*gen_deltaR_ej)[0], my_weight);
       w_deltaR_wenu_b->Fill(deltaR_ebj->empty() ? -1 : (*deltaR_ebj)[0], gen_deltaR_ebj->empty() ? -1 : (*gen_deltaR_ebj)[0], my_weight);
       w_deltaR_wenu_bb->Fill(deltaR_ebj->empty() ? -1 : (*deltaR_ebj)[0], gen_deltaR_ebj->empty() ? -1 : (*gen_deltaR_ebj)[0], my_weight);
       w_deltaR_wenu_2b->Fill(deltaR_ebjbj->empty() ? -1 : (*deltaR_ebjbj)[0], gen_deltaR_ebjbj->empty() ? -1 : (*gen_deltaR_ebjbj)[0], my_weight);
       w_pt_W_wenu->Fill(pt_W_wenu->empty() ? -1 : (*pt_W_wenu)[0], gen_pt_W_wenu->empty() ? -1 : (*gen_pt_W_wenu)[0], my_weight);
       w_pt_W_wenu_b->Fill(pt_W_wenu->empty() ? -1 : (*pt_W_wenu)[0], gen_pt_W_wenu->empty() ? -1 : (*gen_pt_W_wenu)[0], my_weight);
       w_pt_W_wenu_bb->Fill(pt_W_wenu->empty() ? -1 : (*pt_W_wenu)[0], gen_pt_W_wenu->empty() ? -1 : (*gen_pt_W_wenu)[0], my_weight);
       w_eta_W_wenu->Fill(eta_W_wenu->empty() ? -1 : (*eta_W_wenu)[0], gen_eta_W_wenu->empty() ? -1 : (*gen_eta_W_wenu)[0], my_weight);
       w_eta_W_wenu_b->Fill(eta_W_wenu->empty() ? -1 : (*eta_W_wenu)[0], gen_eta_W_wenu->empty() ? -1 : (*gen_eta_W_wenu)[0], my_weight);
       w_eta_W_wenu_bb->Fill(eta_W_wenu->empty() ? -1 : (*eta_W_wenu)[0], gen_eta_W_wenu->empty() ? -1 : (*gen_eta_W_wenu)[0], my_weight);
     }

     if (wmnu_event) {
       w_delta_wmnu->Fill(deltaphi_mj->empty() ? -1 : (*deltaphi_mj)[0], gen_deltaphi_mj->empty() ? -1 : (*gen_deltaphi_mj)[0], my_weight);
       w_delta_wmnu_b->Fill(deltaphi_mbj->empty() ? -1 : (*deltaphi_mbj)[0], gen_deltaphi_mbj->empty() ? -1 : (*gen_deltaphi_mbj)[0], my_weight);
       w_delta_wmnu_bb->Fill(deltaphi_mbj->empty() ? -1 : (*deltaphi_mbj)[0], gen_deltaphi_mbj->empty() ? -1 : (*gen_deltaphi_mbj)[0], my_weight);
       w_delta_wmnu_2b->Fill(deltaphi_mbjbj->empty() ? -1 : (*deltaphi_mbjbj)[0], gen_deltaphi_mbjbj->empty() ? -1 : (*gen_deltaphi_mbjbj)[0], my_weight);
       w_deltaR_wmnu->Fill(deltaR_mj->empty() ? -1 : (*deltaR_mj)[0], gen_deltaR_mj->empty() ? -1 : (*gen_deltaR_mj)[0], my_weight);
       w_deltaR_wmnu_b->Fill(deltaR_mbj->empty() ? -1 : (*deltaR_mbj)[0], gen_deltaR_mbj->empty() ? -1 : (*gen_deltaR_mbj)[0], my_weight);
       w_deltaR_wmnu_bb->Fill(deltaR_mbj->empty() ? -1 : (*deltaR_mbj)[0], gen_deltaR_mbj->empty() ? -1 : (*gen_deltaR_mbj)[0], my_weight);
       w_deltaR_wmnu_2b->Fill(deltaR_mbjbj->empty() ? -1 : (*deltaR_mbjbj)[0], gen_deltaR_mbjbj->empty() ? -1 : (*gen_deltaR_mbjbj)[0], my_weight);
       w_pt_W_wmnu->Fill(pt_W_wmnu->empty() ? -1 : (*pt_W_wmnu)[0], gen_pt_W_wmnu->empty() ? -1 : (*gen_pt_W_wmnu)[0], my_weight);
       w_pt_W_wmnu_b->Fill(pt_W_wmnu->empty() ? -1 : (*pt_W_wmnu)[0], gen_pt_W_wmnu->empty() ? -1 : (*gen_pt_W_wmnu)[0], my_weight);
       w_pt_W_wmnu_bb->Fill(pt_W_wmnu->empty() ? -1 : (*pt_W_wmnu)[0], gen_pt_W_wmnu->empty() ? -1 : (*gen_pt_W_wmnu)[0], my_weight);
       w_eta_W_wmnu->Fill(eta_W_wmnu->empty() ? -1 : (*eta_W_wmnu)[0], gen_eta_W_wmnu->empty() ? -1 : (*gen_eta_W_wmnu)[0], my_weight);
       w_eta_W_wmnu_b->Fill(eta_W_wmnu->empty() ? -1 : (*eta_W_wmnu)[0], gen_eta_W_wmnu->empty() ? -1 : (*gen_eta_W_wmnu)[0], my_weight);
       w_eta_W_wmnu_bb->Fill(eta_W_wmnu->empty() ? -1 : (*eta_W_wmnu)[0], gen_eta_W_wmnu->empty() ? -1 : (*gen_eta_W_wmnu)[0], my_weight);
     }

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
