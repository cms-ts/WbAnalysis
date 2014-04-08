// -*- C++ -*-
//
// Package: GenbAnalyzer
// Class: GenbAnalyzer
//
/**\class GenbAnalyzer GenWbAnalyzer.cc WbAnalysis/WbSkim/src/GenWbAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]

*/
//
// Original Author: Andrea Schizzi
// Created: Thu Jan 10 15:57:03 CET 2013
// $Id: GenWbAnalyzer.cc,v 1.37 2013/07/22 10:27:10 dellaric Exp $
//
//

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stddef.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
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
#include "Math/VectorUtil.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "PhysicsTools/CandUtils/interface/CandCombiner.h"
#include "CommonTools/Utils/interface/MassRangeSelector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "Rivet/Projections/FastJets.hh"

//
// class declaration
//

class GenWbAnalyzer:public  edm::EDProducer {

public:

  explicit GenWbAnalyzer (const edm::ParameterSet &);
  ~GenWbAnalyzer ();

private:

  virtual void beginJob ();
  virtual void produce (edm::Event &, const edm::EventSetup &);
  virtual void endJob ();

  virtual void beginRun (edm::Run &, edm::EventSetup const &);
  virtual void endRun (edm::Run  &, edm::EventSetup const &);
  virtual void beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &);

  struct order_leptons { bool operator() (const TLorentzVector &p1, const TLorentzVector &p2) const {
      return (p1.Pt() > p2.Pt());
    }
  };
  struct order_jets { bool operator() (const fastjet::PseudoJet &j1, const fastjet::PseudoJet &j2) const {
      return (j1.pt() > j2.pt());
    }
  };

  void fill(TH1F* histogram, double value, double weight=1.0) {
    TAxis* axis = histogram->GetXaxis();
    Int_t nx = histogram->GetNbinsX();
    if (axis->FindBin(value) <= 0) {
      histogram->Fill(histogram->GetBinCenter(1), weight);
    } else if (axis->FindBin(value) >= nx+1) {
      histogram->Fill(histogram->GetBinCenter(nx), weight);
    } else {
      histogram->Fill(value, weight);
    }
  };

  bool isB (reco::GenParticle gp) {
    int pid = gp.pdgId();
    if ((abs(pid)/100)%10 == 5 || (abs(pid)/1000)%10 == 5) {
      return true;
    } else {
      return false;
    }
  }

  const reco::GenParticle* getBAncestors (const reco::GenParticle* gp) {
    if (gp->status()<1||gp->status()>3) return NULL;
    if (isB(*gp)) return gp;
    if (gp->numberOfMothers()==0) return NULL;
    for (unsigned int im=0; im<gp->numberOfMothers(); im++) {
      const reco::GenParticle *p = gp->motherRef(im).get();
      const reco::GenParticle* mom = getBAncestors(p);
      if (mom != NULL) return mom; 
    }
    return NULL;
  }

  // ----------member data ---------------------------

  std::string pileupMC_;
  std::string pileupDT_;
  std::string lepton_;
  std::string path_;
  double numB_;
  bool rivet_;

  edm::LumiReWeighting LumiWeights_;

  int nprup;

  TH1F*     h_gen_weights;

  TH1F*     h_eventYields;

  TH1F* h_nmult0;
  TH1F* h_nmult1;

  TH1F*     w_jetmultiplicity;
  TH1F*     b_jetmultiplicity;
  TH1F*     c_jetmultiplicity;
  TH1F*     t_jetmultiplicity;

  TH1F*     w_first_jet_pt;	// leading jet of any type
  TH1F*     w_first_jet_eta;
  TH1F*     w_first_jet_mass;
  TH1F*     w_second_jet_pt;
  TH1F*     w_second_jet_eta;
  TH1F*     w_second_jet_mass;
  TH1F*     w_third_jet_pt;
  TH1F*     w_third_jet_eta;
  TH1F*     w_third_jet_mass;

  TH1F*     w_first_jet_pt_b;	// leading jet with at least one b jet in the event
  TH1F*     w_first_jet_eta_b;
  TH1F*     w_first_jet_mass_b;
  TH1F*     w_second_jet_pt_b;
  TH1F*     w_second_jet_eta_b;
  TH1F*     w_second_jet_mass_b;
  TH1F*     w_third_jet_pt_b;
  TH1F*     w_third_jet_eta_b;
  TH1F*     w_third_jet_mass_b;

  TH1F*     w_bjetmultiplicity;

  TH1F*     w_first_bjet_pt;	// leading b jet
  TH1F*     w_first_bjet_eta;
  TH1F*     w_first_bjet_mass;

  TH1F*     w_single_bjet_pt;	// only 1 b jet
  TH1F*     w_single_bjet_eta;
  TH1F*     w_single_bjet_mass;

  TH1F*     w_second_bjet_pt;
  TH1F*     w_second_bjet_eta;
  TH1F*     w_second_bjet_mass;

  TH1F*     w_third_bjet_pt;
  TH1F*     w_third_bjet_eta;
  TH1F*     w_third_bjet_mass;

  TH1F*     w_mt_wenu;
  TH1F*     w_mt_wmnu;
  TH1F*     w_mass_ee_b_wide;	// at least one b jet in the event
  TH1F*     w_mass_mm_b_wide;
  TH1F*     w_mt_wenu_b_wide;	// at least one b jet in the event
  TH1F*     w_mt_wmnu_b_wide;
  TH1F*     w_mt_wenu_bb_wide;	// at least one b jet in the event
  TH1F*     w_mt_wmnu_bb_wide;
  TH1F*     w_mass_wenu_blepton;	// at least one b jet in the event
  TH1F*     w_mass_wmnu_blepton;
  TH1F*     w_mass_wenu_blepton_b;	// at least one b jet in the event
  TH1F*     w_mass_wmnu_blepton_b;
  TH1F*     w_mass_wenu_blepton_bb;	// at least one b jet in the event
  TH1F*     w_mass_wmnu_blepton_bb;
  TH1F*     w_mt_wenu_b;	// at least one b jet in the event
  TH1F*     w_mt_wmnu_b;
  TH1F*     w_mt_wenu_bb;	// at least one b jet in the event
  TH1F*     w_mt_wmnu_bb;
  TH1F*     w_delta_wenu;
  TH1F*     w_delta_wenu_b;
  TH1F*     w_delta_wenu_bb;
  TH1F*     w_delta_wenu_2b;
  TH1F*     w_delta_wmnu;
  TH1F*     w_delta_wmnu_b;
  TH1F*     w_delta_wmnu_bb;
  TH1F*     w_delta_wmnu_2b;
  TH1F*     w_deltaR_wenu;
  TH1F*     w_deltaR_wenu_b;
  TH1F*     w_deltaR_wenu_bb;
  TH1F*     w_deltaR_wenu_2b;
  TH1F*     w_deltaR_wmnu;
  TH1F*     w_deltaR_wmnu_b;
  TH1F*     w_deltaR_wmnu_bb;
  TH1F*     w_deltaR_wmnu_2b;
  TH1F*     w_single_delta_wenu_b;
  TH1F*     w_single_delta_wmnu_b;
  TH1F*     w_single_deltaR_wenu_b;
  TH1F*     w_single_deltaR_wmnu_b;

};

using namespace  pat;

//
// constants, enums and typedefs
//


//
// static data member definitions
//


//
// constructors and destructor
//
GenWbAnalyzer::GenWbAnalyzer (const edm::ParameterSet & iConfig) {

  pileupMC_ = iConfig.getUntrackedParameter < std::string > ("pileupMC", "S10");
  pileupDT_ = iConfig.getUntrackedParameter < std::string > ("pileupDT", "");
  lepton_ = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  path_ =   iConfig.getUntrackedParameter < std::string > ("path", "/gpfs/cms/users/candelis/work/ZbSkim/test");
  numB_ =  iConfig.getUntrackedParameter <double> ("numB", 0);

  rivet_  = iConfig.getUntrackedParameter < bool > ("rivet", false);

  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_gen_weights     =   fs->make < TH1F > ("h_gen_weights",      "h_gen_weights", 2, 0, 2);

  h_eventYields =   fs->make < TH1F > ("h_eventYields", "h_eventYields;selection", 8, 0.5, 8.5);

  h_nmult0 =            fs->make < TH1F > ("h_nmult0", "h_nmult0", 8, -0.5, 7.5);
  h_nmult1 =            fs->make < TH1F > ("h_nmult1", "h_nmult1", 8, -0.5, 7.5);


  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity", "w_jetmultiplicity;N_jets", 8, 0.5, 8.5);

  w_first_jet_pt =      fs->make < TH1F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 70, 0., 350.);
  w_first_jet_eta =     fs->make < TH1F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 50, -2.5, 2.5);
  w_first_jet_mass =      fs->make < TH1F > ("w_first_jet_mass",    "w_first_jet_mass;Mass [GeV]", 50, 0., 75.);
  w_second_jet_pt =     fs->make < TH1F > ("w_second_jet_pt",   "w_second_jet_pt;P_t [GeV]", 50, 0., 250.);
  w_second_jet_eta =    fs->make < TH1F > ("w_second_jet_eta",  "w_second_jet_eta;Eta", 50, -2.5, 2.5);
  w_second_jet_mass =      fs->make < TH1F > ("w_second_jet_mass",    "w_second_jet_mass;Mass [GeV]", 50, 0., 75.);
  w_third_jet_pt =      fs->make < TH1F > ("w_third_jet_pt",    "w_third_jet_pt;P_t [GeV]", 50, 0., 150.);
  w_third_jet_eta =     fs->make < TH1F > ("w_third_jet_eta",   "w_third_jet_eta;Eta", 50, -2.5, 2.5);
  w_third_jet_mass =      fs->make < TH1F > ("w_third_jet_mass",    "w_third_jet_mass;Mass [GeV]", 50, 0., 75.);

  w_first_jet_pt_b =    fs->make < TH1F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 70, 0., 350.);
  w_first_jet_eta_b =   fs->make < TH1F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 50, -2.5, 2.5);
  w_first_jet_mass_b =      fs->make < TH1F > ("w_first_jet_mass_b",    "w_first_jet_mass_b;Mass [GeV]", 50, 0., 75.);
  w_second_jet_pt_b =   fs->make < TH1F > ("w_second_jet_pt_b",  "w_second_jet_pt_b;P_t [GeV]", 50, 0., 250.);
  w_second_jet_eta_b =  fs->make < TH1F > ("w_second_jet_eta_b", "w_second_jet_eta_b;Eta", 50, -2.5, 2.5);
  w_second_jet_mass_b =      fs->make < TH1F > ("w_second_jet_mass_b",    "w_second_jet_mass_b;Mass [GeV]", 50, 0., 75.);
  w_third_jet_pt_b =    fs->make < TH1F > ("w_third_jet_pt_b",   "w_third_jet_pt_b;P_t [GeV]", 50, 0., 150.);
  w_third_jet_eta_b =   fs->make < TH1F > ("w_third_jet_eta_b",  "w_third_jet_eta_b;Eta", 50, -2.5, 2.5);
  w_third_jet_mass_b =      fs->make < TH1F > ("w_third_jet_mass_b",    "w_third_jet_mass_b;Mass [GeV]", 50, 0., 75.);

  w_bjetmultiplicity =  fs->make < TH1F > ("w_bjetmultiplicity", "w_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);

  w_first_bjet_pt =     fs->make < TH1F > ("w_first_bjet_pt",    "w_first_bjet_pt;P_t [GeV]", 70, 0., 350.);
  w_first_bjet_eta =    fs->make < TH1F > ("w_first_bjet_eta",   "w_first_bjet_eta;Eta", 50, -2.5, 2.5);
  w_first_bjet_mass =      fs->make < TH1F > ("w_first_bjet_mass",    "w_first_bjet_mass;Mass [GeV]", 50, 0., 75.);

  w_single_bjet_pt =    fs->make < TH1F > ("w_single_bjet_pt",    "w_single_bjet_pt;P_t [GeV]", 70, 0., 350.);
  w_single_bjet_eta =   fs->make < TH1F > ("w_single_bjet_eta",   "w_single_bjet_eta;Eta", 50, -2.5, 2.5);
  w_single_bjet_mass =      fs->make < TH1F > ("w_single_bjet_mass",    "w_single_bjet_mass;Mass [GeV]", 50, 0., 75.);

  w_second_bjet_pt =    fs->make < TH1F > ("w_second_bjet_pt",   "w_second_bjet_pt;P_t [GeV]", 50, 0., 250.);
  w_second_bjet_eta =   fs->make < TH1F > ("w_second_bjet_eta",  "w_second_bjet_eta;Eta", 50, -2.5, 2.5);
  w_second_bjet_mass =      fs->make < TH1F > ("w_second_bjet_mass",    "w_second_bjet_mass;Mass [GeV]", 50, 0., 75.);

  w_third_bjet_pt =     fs->make < TH1F > ("w_third_bjet_pt",    "w_third_bjet_pt;P_t [GeV]", 50, 0., 150.);
  w_third_bjet_eta =    fs->make < TH1F > ("w_third_bjet_eta",   "w_third_bjet_eta;Eta", 50, -2.5, 2.5);
  w_third_bjet_mass =      fs->make < TH1F > ("w_third_bjet_mass",    "w_third_bjet_mass;Mass [GeV]", 50, 0., 75.);

  w_mt_wenu =           fs->make < TH1F > ("w_mt_wenu",         "w_mt_wenu;M_{T} [GeV]", 40, 40., 200.);
  w_mt_wmnu =           fs->make < TH1F > ("w_mt_wmnu",         "w_mt_wmnu;M_{T} [GeV]", 40, 40., 200.);
  w_mt_wenu_b_wide =    fs->make < TH1F > ("w_mt_wenu_b_wide",  "w_mt_wenu_b_wide;M_{T} [GeV]", 40, 0., 160.);
  w_mt_wmnu_b_wide =    fs->make < TH1F > ("w_mt_wmnu_b_wide",  "w_mt_wmnu_b_wide;M_{T} [GeV]", 40, 0., 160.);
  w_mt_wenu_bb_wide =   fs->make < TH1F > ("w_mt_wenu_bb_wide", "w_mt_wenu_bb_wide;M_{T} [GeV]", 40, 0., 160.);
  w_mt_wmnu_bb_wide =   fs->make < TH1F > ("w_mt_wmnu_bb_wide", "w_mt_wmnu_bb_wide;M_{T} [GeV]", 40, 0., 160.);
  w_mass_wenu_blepton =      fs->make < TH1F > ("w_mass_wenu_blepton",    "w_mass_wenu_blepton;M_{T} [GeV]", 50, 0., 250.);
  w_mass_wmnu_blepton =      fs->make < TH1F > ("w_mass_wmnu_blepton",    "w_mass_wmnu_blepton;M_{T} [GeV]", 50, 0., 250.);
  w_mass_wenu_blepton_b =      fs->make < TH1F > ("w_mass_wenu_blepton_b",    "w_mass_wenu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  w_mass_wmnu_blepton_b =      fs->make < TH1F > ("w_mass_wmnu_blepton_b",    "w_mass_wmnu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  w_mass_wenu_blepton_bb =      fs->make < TH1F > ("w_mass_wenu_blepton_bb",    "w_mass_wenu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  w_mass_wmnu_blepton_bb =      fs->make < TH1F > ("w_mass_wmnu_blepton_bb",    "w_mass_wmnu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  w_mt_wenu_b =         fs->make < TH1F > ("w_mt_wenu_b",       "w_mt_wenu_b;M_{T} [GeV]", 40, 40., 200.);
  w_mt_wmnu_b =         fs->make < TH1F > ("w_mt_wmnu_b",       "w_mt_wmnu_b;M_{T} [GeV]", 40, 40., 200.);
  w_mt_wenu_bb =        fs->make < TH1F > ("w_mt_wenu_bb",      "w_mt_wenu_bb;M_{T} [GeV]", 40, 40., 200.);
  w_mt_wmnu_bb =        fs->make < TH1F > ("w_mt_wmnu_bb",      "w_mt_wmnu_bb;M_{T} [GeV]", 40, 40., 200.);
  w_delta_wenu =          fs->make < TH1F > ("w_delta_wenu",     "w_delta_wenu",    30, 0, TMath::Pi ());
  w_delta_wenu_b =        fs->make < TH1F > ("w_delta_wenu_b",   "w_delta_wenu_b",  30, 0, TMath::Pi ());
  w_delta_wenu_bb =       fs->make < TH1F > ("w_delta_wenu_bb",  "w_delta_wenu_bb", 30, 0, TMath::Pi ());
  w_delta_wenu_2b =       fs->make < TH1F > ("w_delta_wenu_2b",  "w_delta_wenu_2b", 30, 0, TMath::Pi ());
  w_delta_wmnu =          fs->make < TH1F > ("w_delta_wmnu",     "w_delta_wmnu",    30, 0, TMath::Pi ());
  w_delta_wmnu_b =        fs->make < TH1F > ("w_delta_wmnu_b",   "w_delta_wmnu_b",  30, 0, TMath::Pi ());
  w_delta_wmnu_bb =       fs->make < TH1F > ("w_delta_wmnu_bb",  "w_delta_wmnu_bb", 30, 0, TMath::Pi ());
  w_delta_wmnu_2b =       fs->make < TH1F > ("w_delta_wmnu_2b",  "w_delta_wmnu_2b", 30, 0, TMath::Pi ());
  w_deltaR_wenu =          fs->make < TH1F > ("w_deltaR_wenu",     "w_deltaR_wenu",    30, 0, 3.);
  w_deltaR_wenu_b =        fs->make < TH1F > ("w_deltaR_wenu_b",   "w_deltaR_wenu_b",  30, 0, 3.);
  w_deltaR_wenu_bb =       fs->make < TH1F > ("w_deltaR_wenu_bb",  "w_deltaR_wenu_bb", 30, 0, 3.);
  w_deltaR_wenu_2b =       fs->make < TH1F > ("w_deltaR_wenu_2b",  "w_deltaR_wenu_2b", 30, 0, 3.);
  w_deltaR_wmnu =          fs->make < TH1F > ("w_deltaR_wmnu",     "w_deltaR_wmnu",    30, 0, 3.);
  w_deltaR_wmnu_b =        fs->make < TH1F > ("w_deltaR_wmnu_b",   "w_deltaR_wmnu_b",  30, 0, 3.);
  w_deltaR_wmnu_bb =       fs->make < TH1F > ("w_deltaR_wmnu_bb",  "w_deltaR_wmnu_bb", 30, 0, 3.);
  w_deltaR_wmnu_2b =       fs->make < TH1F > ("w_deltaR_wmnu_2b",  "w_deltaR_wmnu_2b", 30, 0, 3.);
  w_single_delta_wenu_b =        fs->make < TH1F > ("w_single_delta_wenu_b",  "w_single_delta_wenu_b", 12, 0, TMath::Pi ());
  w_single_delta_wmnu_b =        fs->make < TH1F > ("w_single_delta_wmnu_b",  "w_single_delta_wmnu_b", 12, 0, TMath::Pi ());
  w_single_deltaR_wenu_b =        fs->make < TH1F > ("w_single_deltaR_wenu_b",  "w_single_deltaR_wenu_b", 12, 0, TMath::Pi ());
  w_single_deltaR_wmnu_b =        fs->make < TH1F > ("w_single_deltaR_wmnu_b",  "w_single_deltaR_wmnu_b", 12, 0, TMath::Pi ());

  produces<std::vector<double>>("myEventWeight");

  produces<std::vector<math::XYZTLorentzVector>>("myElectrons");
  produces<std::vector<math::XYZTLorentzVector>>("myMuons");

  produces<std::vector<double>>("myPtZ");
  produces<std::vector<double>>("myPtZb");
  
  produces<std::vector<double>>("myYZ");
  produces<std::vector<double>>("myYZb");

  produces<std::vector<double>>("myMassZj");
  produces<std::vector<double>>("myMassZb");

  produces<std::vector<math::XYZTLorentzVector>>("myJets");
  produces<std::vector<math::XYZTLorentzVector>>("myJets2");
  produces<std::vector<double>>("myDeltaPhi");

  produces<std::vector<double>>("myHt");
  produces<std::vector<double>>("myHtb");

  produces<std::vector<math::XYZTLorentzVector>>("myBJets");
  produces<std::vector<math::XYZTLorentzVector>>("myBJets2");
  produces<std::vector<double>>("myBDeltaPhi");
 
}

GenWbAnalyzer::~GenWbAnalyzer () {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void GenWbAnalyzer::produce (edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace edm;
  using namespace std;

  edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);

  std::auto_ptr<std::vector<double>> myEventWeight( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myElectrons( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myMuons( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myPtZ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myPtZb( new std::vector<double> );
  
  std::auto_ptr<std::vector<double>> myYZ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myYZb( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myMassZj( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myMassZb( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets2( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<double>> myDeltaPhi( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myHt( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myHtb( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets2( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<double>> myBDeltaPhi( new std::vector<double> );

  bool wenu_event = false;
  bool wmnu_event = false;
  bool ist = false;

  int Nj = 0;
  int Nj2 = 0;
  int Nb = 0;
  int Nb2 = 0;

  double MyWeight = 1;

  double Ht = 0;

  struct pt_and_particles {
    TLorentzVector p_part;
    vector<unsigned int> lepton_photon;
  };

  struct pt_and_particles ele_dres;
  vector<unsigned int> ele_photons;
  vector<unsigned int> ele_photons_canc;

  struct pt_and_particles mu_dres;
  vector<unsigned int> mu_photons;
  vector<unsigned int> mu_photons_canc;

  double lepton1_eta = -9999;
  double lepton1_phi = -9999;

  // ++++++ Pile-Up

  MyWeight = 1.0;

  Handle < vector < PileupSummaryInfo > > PupInfo;

  if (iEvent.getByLabel (edm::InputTag ("addPileupInfo"), PupInfo))  {

    float Tnpv = -1;

    for (vector < PileupSummaryInfo >::const_iterator PVI = PupInfo->begin (); PVI != PupInfo->end (); ++PVI) {
      int BX = PVI->getBunchCrossing ();
      if (BX == 0) {
        Tnpv = PVI->getTrueNumInteractions ();
        continue;
      }
    }

    MyWeight = LumiWeights_.weight (Tnpv);

  }

  if (rivet_) MyWeight = 1.0;

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;

  if (iEvent.getByLabel ("generator", genEventInfoHandle)) {

    double mcWeight = genEventInfoHandle->weight();

    h_gen_weights->Fill(0.5, 1.0);
    h_gen_weights->Fill(1.5, mcWeight);

    MyWeight = MyWeight*mcWeight;

  }

  Handle<LHEEventProduct> evt;

  if (iEvent.getByLabel ("source", evt)) {

    const lhef::HEPEUP hepeup = evt->hepeup();
    const int nup = hepeup.NUP;
    const std::vector<int> idup = hepeup.IDUP;

    if (nup>=5) {
      if (abs(idup[2])==24 && abs(idup[3])>=11 && abs(idup[3])<=16) {

        double lheWeight = 1.0;
        int nmult = nup-5;
        h_nmult0->Fill(nmult);

        if (nprup==10) {
          if (nmult!=0) lheWeight = 0.0;
        }

        h_nmult1->Fill(nmult,lheWeight);
        MyWeight = MyWeight*lheWeight;

      }
    }

  }

  // +++++++++ MET

  TLorentzVector met;
  met.SetPtEtaPhiM(0.,0.,0.,0.);

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
    if ((fabs(itgen->pdgId())==12 || fabs(itgen->pdgId())==14 || fabs(itgen->pdgId())==16) && itgen->status()==1) { // loop over gen electrons
      TLorentzVector temp_neutrino;
      temp_neutrino.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());
      met += temp_neutrino;
    }
  }

  // +++++++++ ELECTRONS

  //vector < unsigned int > lepton_photon;
  vector <reco::GenParticle> part_vect;

  unsigned int index_ele=0;
  unsigned int index_ele2=0;
  unsigned int index_goodele=0;

  vector < TLorentzVector > vect_ele;
  ele_photons_canc.clear();

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
    if (fabs(itgen->pdgId())==11 && itgen->status()==1) { // loop over gen electrons
      ele_photons.clear();
      TLorentzVector ele;
      ele.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());

      ele_photons.push_back(index_ele);
      if (ele.Pt()>10. && fabs(ele.Eta())<2.4) ele_photons_canc.push_back(index_ele);

      unsigned int index_gamma=0;
      // Loop over photons: FSR dressing for electrons
      for (vector<reco::GenParticle>::const_iterator itgen2=genPart->begin(); itgen2!=genPart->end(); itgen2++) {
        if (fabs(itgen2->pdgId())==22 && itgen2->status()==1) { // loop over primary gen photon
       	  TLorentzVector gam;
          gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
	  double deltaPhi_eg = fabs(gam.Phi()- itgen->phi());
	  if (deltaPhi_eg > acos(-1)) deltaPhi_eg= 2*acos(-1) - deltaPhi_eg;
	  double deltaR_eg = sqrt( deltaPhi_eg*deltaPhi_eg + pow(fabs(gam.Eta()-itgen->eta()),2));
	  if (deltaR_eg< 0.1) {
	    ele_photons.push_back(index_gamma);
	    ele += gam;
	    if (ele.Pt()>10. && fabs(ele.Eta())<2.4) ele_photons_canc.push_back(index_gamma);
	  }
	}
	index_gamma++;
      }

      if (ele.Pt()>30. && fabs(ele.Eta())<2.1) {
	index_goodele++;
	if (ele_dres.lepton_photon.empty()) {
	  ele_dres.p_part = ele;
	  ele_dres.lepton_photon = ele_photons;
	} else {
	  if (ele.Pt()>ele_dres.p_part.Pt()) {
	    ele_dres.p_part = ele;
	    ele_dres.lepton_photon = ele_photons;
	  }
	}
      } else {
	if (ele.Pt()>30. && fabs(ele.Eta())<2.1) index_ele2++;
      }
    }
    index_ele++;
  }

  // Computing Mt:
  double op_met = met.Et();
  double elept = 0.;
  double deltaPhiMetEle = 0.;
  double mt_wenu = 0.;
  bool mt_cut_wenu = false;
  if ( !ele_dres.lepton_photon.empty() ) {
    vect_ele.push_back(ele_dres.p_part);
    elept = vect_ele[0].Pt();
    if (op_met>0. && elept>0.) deltaPhiMetEle = fabs(vect_ele[0].Phi() - met.Phi());
    if (deltaPhiMetEle > acos(-1)) deltaPhiMetEle = 2*acos(-1) - deltaPhiMetEle;
    mt_wenu = sqrt(2*elept*op_met*(1-TMath::Cos(deltaPhiMetEle)));
    mt_cut_wenu = mt_wenu > 45.;
  }

  // +++++++++ MUONS

  vector < TLorentzVector > vect_muon;
  unsigned int index_mu = 0;
  unsigned int index_mu2 = 0;
  unsigned int index_goodmu = 0;

  mu_photons_canc.clear();

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
    if (fabs(itgen->pdgId())==13 && itgen->status()==1) { // loop over gen muons
      mu_photons.clear();
      TLorentzVector muon;
      muon.SetPtEtaPhiM(itgen->pt(),itgen->eta(),itgen->phi(),itgen->mass());

      mu_photons.push_back(index_mu);
      if (muon.Pt()>10. && fabs(muon.Eta())<2.4) mu_photons_canc.push_back(index_mu);

      // Loop over photons: FSR dressing for muons
      unsigned int index_gammamu = 0;
      for (vector<reco::GenParticle>::const_iterator itgen2=genPart->begin(); itgen2!=genPart->end(); itgen2++) {
        if (fabs(itgen2->pdgId())==22 && itgen2->status()==1) { // loop over primary gen photon
	  TLorentzVector gam;
	  gam.SetPtEtaPhiM(itgen2->pt(),itgen2->eta(),itgen2->phi(),itgen2->mass());
	  double deltaPhi_mg = fabs(gam.Phi()- itgen->phi());
	  if (deltaPhi_mg > acos(-1)) deltaPhi_mg= 2*acos(-1) - deltaPhi_mg;
	  double deltaR_mg = sqrt( deltaPhi_mg*deltaPhi_mg + pow(fabs(gam.Eta()-itgen->eta()),2));
	  if (deltaR_mg< 0.1) {
	    mu_photons.push_back(index_gammamu);
	    muon += gam;
	    if (muon.Pt()>10. && fabs(muon.Eta())<2.4) mu_photons_canc.push_back(index_gammamu);
	  }
        }
	index_gammamu++;
      }
      if (muon.Pt()>30. && fabs(muon.Eta())<2.1) {
	index_goodmu++;
	if (mu_dres.lepton_photon.empty()) {
	  mu_dres.p_part = muon;
	  mu_dres.lepton_photon = mu_photons;
	} else {
	  if (muon.Pt()>mu_dres.p_part.Pt()) {
	    mu_dres.p_part = muon;
	    mu_dres.lepton_photon = mu_photons;
	  }
	}
      } else {
	if (muon.Pt()>10. && fabs(muon.Eta())<2.1) index_mu2++;
      }
    }
    index_mu++;
  }

  // Computing Mt:
  op_met = met.Et();
  double muopt = 0.;
  double deltaPhiMetMuo = 0.;
  double mt_wmnu = 0.;
  bool mt_cut_wmnu = false;
  if ( !mu_dres.lepton_photon.empty() ) {
    vect_muon.push_back(mu_dres.p_part);
    muopt = vect_muon[0].Pt();
    if (op_met>0. && muopt>0.) deltaPhiMetMuo = fabs(vect_muon[0].Phi() - met.Phi());
    if (deltaPhiMetMuo > acos(-1)) deltaPhiMetMuo = 2*acos(-1) - deltaPhiMetMuo;
    mt_wmnu = sqrt(2*muopt*op_met*(1-TMath::Cos(deltaPhiMetMuo)));
    mt_cut_wmnu = mt_wmnu > 45.;
  }

  // +++++++++ Decisions:

  wenu_event = (lepton_ == "electron") && index_goodele==1 && index_ele2==0 && index_goodmu==0 && index_mu2==0;
  wmnu_event = (lepton_ == "muon") && index_goodmu==1 && index_mu2==0 && index_goodele==0 && index_ele2==0;

  if (wenu_event) {
    lepton1_eta = ele_dres.p_part.Eta();
    lepton1_phi = ele_dres.p_part.Phi();
   }

  if (wmnu_event) {
    lepton1_eta = mu_dres.p_part.Eta();
    lepton1_phi = mu_dres.p_part.Phi();
   }

  vector<reco::GenParticle> part_jets;
  vector<reco::GenParticle> part_jets_st1;

  unsigned int index_nu = 0;
  unsigned index_ch = 0;

  vector<unsigned int> neutrino;
  vector<unsigned int> ch_part;

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
	if (itgen->status()==1 && (fabs(itgen->pdgId())==12 || fabs(itgen->pdgId())==14 || fabs(itgen->pdgId())==16)) {
		neutrino.push_back(index_nu);
	}
	index_nu++;
/*
	if (itgen->charge()!=0 && itgen->pt()<0.25 && itgen->status()==1) {
		ch_part.push_back(index_ch);
	}
*/
	index_ch++;
  }

  for (vector<reco::GenParticle>::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
	part_jets.push_back(*itgen);
  }

  vector<unsigned int> canc_part;

  if (!mu_photons_canc.empty()) {
    //    canc_part = ele_dres.lepton_photon;
    for (unsigned int i =0; i < mu_photons_canc.size(); i++) {
      canc_part.push_back(mu_photons_canc.at(i));
    }
  }

  if (!ele_photons_canc.empty()) {
    //    canc_part = mu_dres.lepton_photon;
    for (unsigned int l = 0; l < ele_photons_canc.size(); l++) {
      canc_part.push_back(ele_photons_canc.at(l));
    }
  }

  if (!neutrino.empty()) {
	for (unsigned int j=0; j < neutrino.size(); j++ )
		canc_part.push_back(neutrino.at(j));
  }

  if (!ch_part.empty()) {
	for (unsigned int p = 0; p < ch_part.size(); p++ )
		canc_part.push_back(ch_part.at(p));
  }

  //ordino il vettore in modo crescente
  stable_sort(canc_part.begin(), canc_part.end());

  //scorro il vettore dal basso e per ogni elemento cancello il corrispondente elemento da part_jets
  for (int i = (canc_part.size() - 1); i >= 0; i--)
	part_jets.erase(part_jets.begin() + canc_part.at(i));

  for (unsigned int j=0; j<part_jets.size(); j++) {
	if (part_jets.at(j).status()==1)
		part_jets_st1.push_back(part_jets.at(j));
  }

  // part_jets_st1 contiene ora tutte le particelle su cui fare la riclusterizzazione
  std::vector<fastjet::PseudoJet> vecs;

  for (unsigned int l = 0; l < part_jets_st1.size(); l++) {
	TLorentzVector part;
	part.SetPtEtaPhiM(part_jets_st1.at(l).pt(),part_jets_st1.at(l).eta(),part_jets_st1.at(l).phi(),part_jets_st1.at(l).mass());

	fastjet::PseudoJet pseudoJet(part_jets_st1.at(l).px(), part_jets_st1.at(l).py(), part_jets_st1.at(l).pz(), part.E());
	pseudoJet.set_user_index(l);
	vecs.push_back(pseudoJet);
  }


  // ++++++++ JETS

  vector<fastjet::PseudoJet> vect_jets;
  vector<fastjet::PseudoJet> vect_jets2;
  fastjet::ClusterSequence cseq(vecs, fastjet::JetDefinition(fastjet:: antikt_algorithm, 0.5));
  vector<fastjet::PseudoJet> jets = sorted_by_pt(cseq.inclusive_jets(0.0));
  for (unsigned int i = 0; i < jets.size(); i++) {
	double etaj = jets[i].eta();
	double phij = jets[i].phi();
	double ptj = jets[i].perp();
	       
	if (fabs(etaj) < 2.4 && ptj > 25) {
          double delta_eta1 = lepton1_eta - etaj;
          double delta_phi1 = fabs(lepton1_phi - phij);
          if (delta_phi1 > acos(-1)) delta_phi1 = 2*acos(-1) - delta_phi1;

          double deltaR_jl1 = sqrt(pow(delta_eta1,2) + pow(delta_phi1,2));
	  
	  if (deltaR_jl1 > 0.5 && (wenu_event || wmnu_event)) { 
	    Nj++;
            vect_jets.push_back(jets[i]);
            Ht = Ht + jets[i].perp();
          }
        } else {
	  if (fabs(etaj) > 2.4 && fabs(etaj) < 5.0 && ptj > 25) {
	    Nj2++;
	    vect_jets2.push_back(jets[i]);
	  }
	}
  }

  // ++++++++ BJETS

  /*loop over gen particles, find the b*/

  vector<fastjet::PseudoJet> vect_bjets;
  vector<fastjet::PseudoJet> vect_bjets2;

  bool Bjet_found = false;
  bool Bjet2_found = false;

  for(unsigned int k = 0; k < vect_jets.size() ; k++) {
    Bjet_found = false;
    double B_eta = 9999;
    double B_phi = 9999;
    double deltaR_Bb = 9999;

    vector<fastjet::PseudoJet> constituents = cseq.constituents(vect_jets[k]);

    for (unsigned int c = 0; c < constituents.size() && !Bjet_found; c++) {
      int index = constituents.at(c).user_index();
      reco::GenParticle gp = part_jets_st1.at(index);
      const reco::GenParticle* bpart = getBAncestors(&gp);
      if (bpart!=NULL) {
        Bjet_found = true;
        B_eta = bpart->eta();
        B_phi = bpart->phi();
     }
    }

    if (Bjet_found) {
      double eta_bj = vect_jets[k].eta();
      double phi_bj = vect_jets[k].phi();

      double deltaEta_Bb = eta_bj - B_eta; 
      double deltaPhi_Bb = fabs(phi_bj - B_phi);
      if (deltaPhi_Bb > acos(-1)) deltaPhi_Bb = 2*acos(-1) - deltaPhi_Bb;

      deltaR_Bb = sqrt(pow(deltaEta_Bb,2) + pow(deltaPhi_Bb,2)); 

      if (deltaR_Bb < 0.5) {
        vect_bjets.push_back(vect_jets[k]);
        Nb++;
      }
    }  
  }  

  for(unsigned int k = 0; k < vect_jets2.size() ; k++) {
    Bjet2_found = false;
    double B_eta2 = 9999;
    double B_phi2 = 9999;
    double deltaR_Bb2 = 9999;

    vector<fastjet::PseudoJet> constituents = cseq.constituents(vect_jets2[k]);

    for (unsigned int c = 0; c < constituents.size() && !Bjet2_found; c++) {
      int index = constituents.at(c).user_index();
      reco::GenParticle gp = part_jets_st1.at(index);
      const reco::GenParticle * bpart2 = getBAncestors(&gp);
      if (bpart2!=NULL) {
        Bjet2_found = true;
        B_eta2 = bpart2->eta();
        B_phi2 = bpart2->phi();
      }
    }
    if (Bjet2_found) {
      double eta_bj = vect_jets2[k].eta();
      double phi_bj = vect_jets2[k].phi();

      double deltaEta_Bb = eta_bj - B_eta2;
      double deltaPhi_Bb = fabs(phi_bj - B_phi2);
      if (deltaPhi_Bb > acos(-1)) deltaPhi_Bb = 2*acos(-1) - deltaPhi_Bb;

      deltaR_Bb2 = sqrt(pow(deltaEta_Bb,2) + pow(deltaPhi_Bb,2));

      if (deltaR_Bb2 < 0.5) {
        vect_bjets2.push_back(vect_jets2[k]);
        Nb2++;
      }
    }
  }

  // Sort b-jets in pT
  std::sort( vect_bjets.begin(), vect_bjets.end(), order_jets() );
  std::sort( vect_bjets2.begin(), vect_bjets2.end(), order_jets() );

  for (std::vector <reco::GenParticle>::const_iterator thepart = genPart->begin(); thepart != genPart->end(); thepart++) {
    if (thepart->pdgId()==24) {
      for (UInt_t i=0; i<thepart->numberOfDaughters(); i++){
        if (abs(thepart->daughter(i)->pdgId())==15 && thepart->daughter(i)->status()==3){
          ist = true;
	}
      }
    }
  }

  wenu_event = wenu_event && !ist && Nj>1 && Nb>0 && Nj2==0;
  wmnu_event = wmnu_event && !ist && Nj>1 && Nb>0 && Nj2==0;


  // ++++++++ EVENT YIELDS:

  h_eventYields->Fill(1);
  if (lepton_ == "muon" && index_goodmu==1) h_eventYields->Fill(2);
  if (lepton_ == "muon" && index_goodmu==1 && index_goodele==0) h_eventYields->Fill(3);
  if (lepton_ == "muon" && index_goodmu==1 && index_goodele==0 && mt_cut_wmnu) h_eventYields->Fill(4);
  if (lepton_ == "muon" && index_goodmu==1 && index_goodele==0 && mt_cut_wmnu && !ist) h_eventYields->Fill(5);
  if (lepton_ == "muon" && index_goodmu==1 && index_goodele==0 && mt_cut_wmnu && !ist && Nj>1) h_eventYields->Fill(6);
  if (lepton_ == "muon" && index_goodmu==1 && index_goodele==0 && mt_cut_wmnu && !ist && Nj>1 && Nb>0) h_eventYields->Fill(7);
  if (wmnu_event && mt_cut_wmnu) h_eventYields->Fill(8);


  // ++++++++ WeNu PLOTS

  if (wenu_event && mt_cut_wenu) {
    w_mt_wenu->Fill (mt_wenu, MyWeight);
    double delta_phi_ej = fabs(vect_ele[0].Phi() - vect_jets[0].phi());
    double delta_eta_ej = fabs(vect_ele[0].Eta() - vect_jets[0].eta());
    if (delta_phi_ej > acos (-1)) delta_phi_ej = 2 * acos (-1) - delta_phi_ej;
    double DR_ej = TMath::Sqrt(delta_phi_ej*delta_phi_ej + delta_eta_ej*delta_eta_ej);
    w_delta_wenu->Fill (delta_phi_ej, MyWeight);
    w_deltaR_wenu->Fill (DR_ej, MyWeight);
    double delta_phi_ebj = fabs(vect_ele[0].Phi() - vect_bjets[0].phi());
    double delta_eta_ebj = fabs(vect_ele[0].Eta() - vect_bjets[0].eta());
    if (delta_phi_ebj > acos (-1)) delta_phi_ebj = 2 * acos (-1) - delta_phi_ebj;
    double DR_ebj = TMath::Sqrt(delta_phi_ebj*delta_phi_ebj + delta_eta_ebj*delta_eta_ebj);
    w_delta_wenu_b->Fill (delta_phi_ebj, MyWeight);
    w_deltaR_wenu_b->Fill (DR_ebj, MyWeight);
    TLorentzVector belectron(vect_ele[0].Px(),vect_ele[0].Py(),vect_ele[0].Pz(),vect_ele[0].E());
    TLorentzVector belectron2(vect_bjets[0].px(),vect_bjets[0].py(),vect_bjets[0].pz(),vect_bjets[0].E());
    belectron += belectron2;
    w_mass_wenu_blepton->Fill(belectron.M(), MyWeight);
    if (Nb == 1) {
      w_mt_wenu_b->Fill (mt_wenu, MyWeight);
      w_single_delta_wenu_b->Fill (delta_phi_ebj, MyWeight);
      w_single_deltaR_wenu_b->Fill (DR_ebj, MyWeight);
      //      w_mass_wenu_blepton_b->Fill(belectron.mass(), MyWeight);
    }
    if (Nb > 1) {
      w_mt_wenu_bb->Fill (mt_wenu, MyWeight);
      w_delta_wenu_bb->Fill (delta_phi_ebj, MyWeight);
      w_deltaR_wenu_bb->Fill (DR_ebj, MyWeight);
      double delta_phi_ebjbj = fabs(vect_bjets[0].phi() - vect_bjets[1].phi());
      double delta_eta_ebjbj = fabs(vect_bjets[0].eta() - vect_bjets[1].eta());
      if (delta_phi_ebjbj > acos (-1)) delta_phi_ebjbj = 2 * acos (-1) - delta_phi_ebjbj;
      double DR_ebjbj = TMath::Sqrt(delta_phi_ebjbj*delta_phi_ebjbj + delta_eta_ebjbj*delta_eta_ebjbj);
      w_delta_wenu_2b->Fill (delta_phi_ebjbj, MyWeight);
      w_deltaR_wenu_2b->Fill (DR_ebjbj, MyWeight);
      //      w_mass_wenu_blepton_bb->Fill(belectron.mass(), MyWeight);
    }
  }


  // ++++++++ WmNu PLOTS

  if (wmnu_event && mt_cut_wmnu) {
    w_mt_wmnu->Fill (mt_wmnu, MyWeight);
    double delta_phi_mj = fabs(vect_muon[0].Phi() - vect_jets[0].phi());
    double delta_eta_mj = fabs(vect_muon[0].Eta() - vect_jets[0].eta());
    if (delta_phi_mj > acos (-1)) delta_phi_mj = 2 * acos (-1) - delta_phi_mj;
    double DR_mj = TMath::Sqrt(delta_phi_mj*delta_phi_mj + delta_eta_mj*delta_eta_mj);
    w_delta_wmnu->Fill (delta_phi_mj, MyWeight);
    w_deltaR_wmnu->Fill (DR_mj, MyWeight);
    double delta_phi_mbj = fabs(vect_muon[0].Phi() - vect_bjets[0].phi());
    double delta_eta_mbj = fabs(vect_muon[0].Eta() - vect_bjets[0].eta());
    if (delta_phi_mbj > acos (-1)) delta_phi_mbj = 2 * acos (-1) - delta_phi_mbj;
    double DR_mbj = TMath::Sqrt(delta_phi_mbj*delta_phi_mbj + delta_eta_mbj*delta_eta_mbj);
    w_delta_wmnu_b->Fill (delta_phi_mbj, MyWeight);
    w_deltaR_wmnu_b->Fill (DR_mbj, MyWeight);
    TLorentzVector bmuon(vect_muon[0].Px(),vect_muon[0].Py(),vect_muon[0].Pz(),vect_muon[0].E());
    TLorentzVector bmuon2(vect_bjets[0].px(),vect_bjets[0].py(),vect_bjets[0].pz(),vect_bjets[0].E());
    bmuon += bmuon2;
    w_mass_wmnu_blepton->Fill(bmuon.M(), MyWeight);
    if (Nb == 1) {
      w_mt_wmnu_b->Fill (mt_wmnu, MyWeight);
      w_single_delta_wmnu_b->Fill (delta_phi_mbj, MyWeight);
      w_single_deltaR_wmnu_b->Fill (DR_mbj, MyWeight);
      //      w_mass_wmnu_blepton_b->Fill(bmuon.mass(), MyWeight);
    }
    if (Nb > 1) {
      w_mt_wmnu_bb->Fill (mt_wmnu, MyWeight);
      w_delta_wmnu_bb->Fill (delta_phi_mbj, MyWeight);
      w_deltaR_wmnu_bb->Fill (DR_mbj, MyWeight);
      double delta_phi_mbjbj = fabs(vect_bjets[0].phi() - vect_bjets[1].phi());
      double delta_eta_mbjbj = fabs(vect_bjets[0].eta() - vect_bjets[1].eta());
      if (delta_phi_mbjbj > acos (-1)) delta_phi_mbjbj = 2 * acos (-1) - delta_phi_mbjbj;
      double DR_mbjbj = TMath::Sqrt(delta_phi_mbjbj*delta_phi_mbjbj + delta_eta_mbjbj*delta_eta_mbjbj);
      w_delta_wmnu_2b->Fill (delta_phi_mbjbj, MyWeight);
      w_deltaR_wmnu_2b->Fill (DR_mbjbj, MyWeight);
      //      w_mass_wmnu_blepton_bb->Fill(bmuon.mass(), MyWeight);
    }
  }


  // ++++++++ JETS PLOTS

  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) {
    w_jetmultiplicity->Fill (Nj, MyWeight);
    w_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight);
    //    w_first_jet_mass->Fill (vect_jets[0].mass(), MyWeight);
  }

  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) {
    w_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight);
    w_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight);
    //    w_second_jet_mass->Fill (vect_jets[1].mass(), MyWeight);
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nj > 2) {
    w_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight);
    w_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight);
    //    w_third_jet_mass->Fill (vect_jets[2].mass(), MyWeight);
  }

  // ++++++++ B JETS PLOTS

  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) {
    w_bjetmultiplicity->Fill (Nb, MyWeight);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight);
    //    w_first_jet_mass_b->Fill (vect_jets[0].mass(), MyWeight);
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    //    w_first_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight);
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 1) {
    w_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight);
    w_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight);
    //    w_second_jet_mass_b->Fill (vect_jets[1].mass(), MyWeight);
    w_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight);
    w_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight);
    //    w_second_bjet_mass->Fill (vect_bjets[1].mass(), MyWeight);
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 2) {
    w_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight);
    w_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight);
    //    w_third_jet_mass_b->Fill (vect_jets[2].mass(), MyWeight);
    w_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight);
    w_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight);
    //    w_third_bjet_mass->Fill (vect_bjets[2].mass(), MyWeight);
  }

  // ++++++++ SINGLE BJET

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb == 1) {
    w_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight);
    w_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight);
    //    w_single_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight);
  }


//  // ++++++++ OUTPUT COLLECTIONS
//
//  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) {
//     myEventWeight->push_back(MyWeight);
//  }
//
//  if (wenu_event && mt_cut_wenu) {
//    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[0].Px(),vect_ele[0].Py(),vect_ele[0].Pz(),vect_ele[0].E()));
//  }
//
//  if (mm_event && Nj > 0) {
//    myMuons->push_back(math::XYZTLorentzVector(vect_muon[0].Px(),vect_muon[0].Py(),vect_muon[0].Pz(),vect_muon[0].E()));
//  }
//
//  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) {
//    for (unsigned int i=0; i<vect_jets.size(); ++i) {
//      myJets->push_back(math::XYZTLorentzVector(vect_jets[i].px(),vect_jets[i].py(),vect_jets[i].pz(),vect_jets[i].e()));
//    }
//    for (unsigned int i=0; i<vect_bjets.size(); ++i) {
//      myBJets->push_back(math::XYZTLorentzVector(vect_bjets[i].px(),vect_bjets[i].py(),vect_bjets[i].pz(),vect_bjets[i].e()));
//    }
//  }
//
//  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) {
//    myHt->push_back(Ht);
//  }
//
//  iEvent.put( myEventWeight, "myEventWeight" );
//
//  iEvent.put( myElectrons, "myElectrons" );
//  iEvent.put( myMuons, "myMuons" );
//
//  iEvent.put( myJets, "myJets" );
//
//  iEvent.put( myHt, "myHt" );
//
//  iEvent.put( myBJets, "myBJets" );

}

// ------------ method called once each job just before starting event loop ------------
void GenWbAnalyzer::beginJob () {

  LumiWeights_ = edm::LumiReWeighting(path_ + "/" + "pileup_" + pileupMC_ + ".root", path_ + "/" + "pileup_2012_" + pileupDT_ + ".root", "pileup", "pileup");

}

// ------------ method called once each job just after ending the event loop ------------
void GenWbAnalyzer::endJob () {
}

// ------------ method called when starting to processes a run ------------
void GenWbAnalyzer::beginRun (edm::Run & iRun, edm::EventSetup const & iSetup) {

  nprup = 0;

  edm::Handle<LHERunInfoProduct> run;
  if (iRun.getByLabel ("source", run)) {
    const lhef::HEPRUP heprup = run->heprup();
    nprup = heprup.NPRUP;
  }

}

// ------------ method called when ending the processing of a run ------------
void GenWbAnalyzer::endRun (edm::Run &, edm::EventSetup const &) {
}

// ------------ method called when starting to processes a luminosity block ------------
void GenWbAnalyzer::beginLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a luminosity block ------------
void GenWbAnalyzer::endLuminosityBlock (edm::LuminosityBlock const &, edm::EventSetup const &) {
}

// define this as a plug-in
DEFINE_FWK_MODULE (GenWbAnalyzer);

