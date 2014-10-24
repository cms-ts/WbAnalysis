// -*- C++ -*-
//
// Package: WbAnalyzer
// Class: WbAnalyzer
//
/**\class WbAnalyzer WbAnalyzer.cc WbAnalysis/WbAnalyzer/src/WbAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]

*/
//
// Original Author: Andrea Schizzi
// Created: Thu Jan 10 15:57:03 CET 2013
// $Id: WbAnalyzer.cc,v 1.111 2013/07/20 07:23:57 dellaric Exp $
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
#include "Math/VectorUtil.h"
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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
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

#include "table.h"

//
// class declaration
//

class WbAnalyzer : public edm::EDProducer {

public:

  explicit WbAnalyzer (const edm::ParameterSet &);
  ~WbAnalyzer ();

private:

  virtual void beginJob ();
  virtual void produce (edm::Event &, const edm::EventSetup &);
  virtual void endJob ();

  virtual void beginRun (edm::Run &, edm::EventSetup const &);
  virtual void endRun (edm::Run &, edm::EventSetup const &);
  virtual void beginLuminosityBlock (edm::LuminosityBlock &, edm::EventSetup const &);
  virtual void endLuminosityBlock (edm::LuminosityBlock &, edm::EventSetup const &);

  double btagSF(bool isMC, std::vector < pat::Jet >& jets, int k) {
    if (isMC == false) return 1.0;

    double w0n=1.0;
    double w1n=0.0;

    for (unsigned int i=0;i<jets.size();i++) {
      if ((fabs(jets[i].partonFlavour()) == 5 || fabs(jets[i].partonFlavour()) == 4)) {
        w0n = w0n * (1.0 - BtSF_->Val(jets[i].pt(), jets[i].eta()));
      } else {
        w0n = w0n * (1.0 - LtSF_->Val(jets[i].pt(), jets[i].eta()));
      }
      double w=1.0;
      for (unsigned int j=0;j<jets.size();j++) {
        if (i!=j) {
          if ((fabs(jets[j].partonFlavour()) == 5 || fabs(jets[j].partonFlavour()) == 4)) {
            w = w * (1.0 - BtSF_->Val(jets[j].pt(), jets[j].eta()));
          } else {
            w = w * (1.0 - LtSF_->Val(jets[j].pt(), jets[j].eta()));
          }
        }
      }
      if ((fabs(jets[i].partonFlavour()) == 5 || fabs(jets[i].partonFlavour()) == 4)) {
        w = w * BtSF_->Val(jets[i].pt(), jets[i].eta());
      } else {
        w = w * LtSF_->Val(jets[i].pt(), jets[i].eta());
      }
      w1n = w1n + w;
    }

    if (k>=1) return (1.0-w0n);     // >= 1 b tagged jet
    if (k==2) return (1.0-w0n-w1n); // >= 2 b tagged jets
    if (k==3) return (1.0-w0n-w1n); // >= 3 b tagged jets // FIXME //
    return (0);

  };

  double jetResolutionCorrection(double jetEta, double jetPt, double jetPtGen, int syst) {

    if (jetPt <= 0 || jetPtGen <= 0) return jetPt;

    double correctionFactor[7]     = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
    double correctionFactorUp[7]   = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
    double correctionFactorDown[7] = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};

    int index = 0;

    if (                      fabs(jetEta) <= 0.5) index = 0;
    if (fabs(jetEta) > 0.5 && fabs(jetEta) <= 1.1) index = 1;
    if (fabs(jetEta) > 1.1 && fabs(jetEta) <= 1.7) index = 2;
    if (fabs(jetEta) > 1.7 && fabs(jetEta) <= 2.3) index = 3;
    if (fabs(jetEta) > 2.3 && fabs(jetEta) <= 2.8) index = 4;
    if (fabs(jetEta) > 2.8 && fabs(jetEta) <= 3.2) index = 5;
    if (fabs(jetEta) > 3.2 && fabs(jetEta) <= 5.0) index = 6;

    double jetPtNew = jetPt;

    if (syst ==  0) jetPtNew = jetPtGen + correctionFactor[index] * (jetPt-jetPtGen);
    if (syst == +1) jetPtNew = jetPtGen + correctionFactorUp[index] * (jetPt-jetPtGen);
    if (syst == -1) jetPtNew = jetPtGen + correctionFactorDown[index] * (jetPt-jetPtGen);

    return fmax(0.0, jetPtNew/jetPt);

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

  // ----------member data ---------------------------

  std::string pileupMC_;
  std::string pileupDT_;
  std::string lepton_;
  double par_;
  double par2_;
  bool usePartonFlavour_;
  bool pcut_;
  bool useDeltaR_;
  std::string path_;
  unsigned int icut_;

  JetCorrectionUncertainty *jetCorrectionUncertainty_;
  edm::LumiReWeighting LumiWeights_;

  int nprup;

  table* ElSF_;
  table* ElSF2_;
  table* MuSF_;
  table* MuSF2_;
  table* MuSF3_;
  table* BtSF_;
  table* LtSF_;

  /****************** LEGEND *************************

   w_plotname_b => at least a b in the event, Nb > 0
   b_plotname   => b quark fraction, isb
   c_plotname   => c quark fraction, isc

   ***************************************************/

  TH1F* h_nmult0;
  TH1F* h_nmult1;

  TH1F*     h_jetmultiplicity;

  TH1F*     h_eventYields;
  TH1F*     w_eventYields;

  TH1F*     h_trgMatchEle;
  TH1F*     h_trgMatchMuo;

  TH1F*     h_pu_weights;

  TH1F*     h_tracks;
  TH1F*     w_tracks;
  TH1F*     h_recoVTX;
  TH1F*     w_recoVTX;

  TH1F*     w_jetmultiplicity;
  TH1F*     b_jetmultiplicity;
  TH1F*     c_jetmultiplicity;
  TH1F*     t_jetmultiplicity;

  TH1F*     h_first_jet_pt;	// leading jet of any type
  TH1F*     w_first_jet_pt;
  TH1F*     b_first_jet_pt;
  TH1F*     c_first_jet_pt;
  TH1F*     t_first_jet_pt;
  TH1F*     w_first_jet_eta;
  TH1F*     b_first_jet_eta;
  TH1F*     c_first_jet_eta;
  TH1F*     t_first_jet_eta;
  TH1F*     w_first_jet_mass;
  TH1F*     b_first_jet_mass;
  TH1F*     c_first_jet_mass;
  TH1F*     t_first_jet_mass;
  TH1F*     h_second_jet_pt;
  TH1F*     w_second_jet_pt;
  TH1F*     b_second_jet_pt;
  TH1F*     c_second_jet_pt;
  TH1F*     t_second_jet_pt;
  TH1F*     w_second_jet_eta;
  TH1F*     b_second_jet_eta;
  TH1F*     c_second_jet_eta;
  TH1F*     t_second_jet_eta;
  TH1F*     w_second_jet_mass;
  TH1F*     b_second_jet_mass;
  TH1F*     c_second_jet_mass;
  TH1F*     t_second_jet_mass;
  TH1F*     h_dijet_pt;
  TH1F*     w_dijet_pt;
  TH1F*     b_dijet_pt;
  TH1F*     c_dijet_pt;
  TH1F*     t_dijet_pt;
  TH1F*     w_dijet_eta;
  TH1F*     b_dijet_eta;
  TH1F*     c_dijet_eta;
  TH1F*     t_dijet_eta;
  TH1F*     w_dijet_mass;
  TH1F*     b_dijet_mass;
  TH1F*     c_dijet_mass;
  TH1F*     t_dijet_mass;

  TH1F*     h_first_jet_pt_b;	// leading jet with at least one b jet in the event
  TH1F*     w_first_jet_pt_b;
  TH1F*     b_first_jet_pt_b;
  TH1F*     c_first_jet_pt_b;
  TH1F*     t_first_jet_pt_b;
  TH1F*     w_first_jet_eta_b;
  TH1F*     b_first_jet_eta_b;
  TH1F*     c_first_jet_eta_b;
  TH1F*     t_first_jet_eta_b;
  TH1F*     w_first_jet_mass_b;
  TH1F*     b_first_jet_mass_b;
  TH1F*     c_first_jet_mass_b;
  TH1F*     t_first_jet_mass_b;
  TH1F*     h_second_jet_pt_b;
  TH1F*     w_second_jet_pt_b;
  TH1F*     b_second_jet_pt_b;
  TH1F*     c_second_jet_pt_b;
  TH1F*     t_second_jet_pt_b;
  TH1F*     w_second_jet_eta_b;
  TH1F*     b_second_jet_eta_b;
  TH1F*     c_second_jet_eta_b;
  TH1F*     t_second_jet_eta_b;
  TH1F*     w_second_jet_mass_b;
  TH1F*     b_second_jet_mass_b;
  TH1F*     c_second_jet_mass_b;
  TH1F*     t_second_jet_mass_b;
  TH1F*     h_dijet_pt_b;
  TH1F*     w_dijet_pt_b;
  TH1F*     b_dijet_pt_b;
  TH1F*     c_dijet_pt_b;
  TH1F*     t_dijet_pt_b;
  TH1F*     w_dijet_eta_b;
  TH1F*     b_dijet_eta_b;
  TH1F*     c_dijet_eta_b;
  TH1F*     t_dijet_eta_b;
  TH1F*     w_dijet_mass_b;
  TH1F*     b_dijet_mass_b;
  TH1F*     c_dijet_mass_b;
  TH1F*     t_dijet_mass_b;

  TH1F*     h_first_jet_pt_bb;	// leading jet with at least one b jet in the event
  TH1F*     w_first_jet_pt_bb;
  TH1F*     b_first_jet_pt_bb;
  TH1F*     c_first_jet_pt_bb;
  TH1F*     t_first_jet_pt_bb;
  TH1F*     w_first_jet_eta_bb;
  TH1F*     b_first_jet_eta_bb;
  TH1F*     c_first_jet_eta_bb;
  TH1F*     t_first_jet_eta_bb;
  TH1F*     w_first_jet_mass_bb;
  TH1F*     b_first_jet_mass_bb;
  TH1F*     c_first_jet_mass_bb;
  TH1F*     t_first_jet_mass_bb;
  TH1F*     h_second_jet_pt_bb;
  TH1F*     w_second_jet_pt_bb;
  TH1F*     b_second_jet_pt_bb;
  TH1F*     c_second_jet_pt_bb;
  TH1F*     t_second_jet_pt_bb;
  TH1F*     w_second_jet_eta_bb;
  TH1F*     b_second_jet_eta_bb;
  TH1F*     c_second_jet_eta_bb;
  TH1F*     t_second_jet_eta_bb;
  TH1F*     w_second_jet_mass_bb;
  TH1F*     b_second_jet_mass_bb;
  TH1F*     c_second_jet_mass_bb;
  TH1F*     t_second_jet_mass_bb;
  TH1F*     h_dijet_pt_bb;
  TH1F*     w_dijet_pt_bb;
  TH1F*     b_dijet_pt_bb;
  TH1F*     c_dijet_pt_bb;
  TH1F*     t_dijet_pt_bb;
  TH1F*     w_dijet_eta_bb;
  TH1F*     b_dijet_eta_bb;
  TH1F*     c_dijet_eta_bb;
  TH1F*     t_dijet_eta_bb;
  TH1F*     w_dijet_mass_bb;
  TH1F*     b_dijet_mass_bb;
  TH1F*     c_dijet_mass_bb;
  TH1F*     t_dijet_mass_bb;

  TH1F*     w_bjetmultiplicity;
  TH1F*     b_bjetmultiplicity;
  TH1F*     c_bjetmultiplicity;
  TH1F*     t_bjetmultiplicity;

  TH1F*     h_first_bjet_pt;	// leading b jet
  TH1F*     w_first_bjet_pt;
  TH1F*     b_first_bjet_pt;
  TH1F*     c_first_bjet_pt;
  TH1F*     t_first_bjet_pt;
  TH1F*     w_first_bjet_eta;
  TH1F*     b_first_bjet_eta;
  TH1F*     c_first_bjet_eta;
  TH1F*     t_first_bjet_eta;
  TH1F*     w_first_bjet_mass;
  TH1F*     b_first_bjet_mass;
  TH1F*     c_first_bjet_mass;
  TH1F*     t_first_bjet_mass;

  TH1F*     w_single_bjet_pt;	// only 1 b jet
  TH1F*     b_single_bjet_pt;
  TH1F*     c_single_bjet_pt;
  TH1F*     t_single_bjet_pt;
  TH1F*     w_single_bjet_eta;
  TH1F*     b_single_bjet_eta;
  TH1F*     c_single_bjet_eta;
  TH1F*     t_single_bjet_eta;
  TH1F*     w_single_bjet_mass;
  TH1F*     b_single_bjet_mass;
  TH1F*     c_single_bjet_mass;
  TH1F*     t_single_bjet_mass;

  TH1F*     h_second_bjet_pt;
  TH1F*     w_second_bjet_pt;
  TH1F*     b_second_bjet_pt;
  TH1F*     c_second_bjet_pt;
  TH1F*     t_second_bjet_pt;
  TH1F*     w_second_bjet_eta;
  TH1F*     b_second_bjet_eta;
  TH1F*     c_second_bjet_eta;
  TH1F*     t_second_bjet_eta;
  TH1F*     w_second_bjet_mass;
  TH1F*     b_second_bjet_mass;
  TH1F*     c_second_bjet_mass;
  TH1F*     t_second_bjet_mass;

  TH1F*     w_first_ele_pt;
  TH1F*     w_first_ele_pt_b;
  TH1F*     w_first_ele_pt_bb;
  TH1F*     h_first_ele_pt;
  TH1F*     h_first_ele_pt_b;
  TH1F*     h_first_ele_pt_bb;
  TH1F*     b_first_ele_pt;
  TH1F*     c_first_ele_pt;
  TH1F*     t_first_ele_pt;
  TH1F*     w_second_ele_pt;
  TH1F*     b_second_ele_pt;
  TH1F*     c_second_ele_pt;
  TH1F*     t_second_ele_pt;
  TH1F*     w_first_muon_pt;
  TH1F*     w_first_muon_pt_b;
  TH1F*     w_first_muon_pt_bb;
  TH1F*     h_first_muon_pt;
  TH1F*     h_first_muon_pt_b;
  TH1F*     h_first_muon_pt_bb;
  TH1F*     b_first_muon_pt;
  TH1F*     c_first_muon_pt;
  TH1F*     t_first_muon_pt;
  TH1F*     w_second_muon_pt;
  TH1F*     b_second_muon_pt;
  TH1F*     c_second_muon_pt;
  TH1F*     t_second_muon_pt;
  TH1F*     w_first_ele_eta;
  TH1F*     b_first_ele_eta;
  TH1F*     c_first_ele_eta;
  TH1F*     t_first_ele_eta;
  TH1F*     w_second_ele_eta;
  TH1F*     b_second_ele_eta;
  TH1F*     c_second_ele_eta;
  TH1F*     t_second_ele_eta;
  TH1F*     w_first_ele_iso;
  TH1F*     b_first_ele_iso;
  TH1F*     c_first_ele_iso;
  TH1F*     t_first_ele_iso;
  TH1F*     w_first_muon_eta;
  TH1F*     b_first_muon_eta;
  TH1F*     c_first_muon_eta;
  TH1F*     t_first_muon_eta;
  TH1F*     w_second_muon_eta;
  TH1F*     b_second_muon_eta;
  TH1F*     c_second_muon_eta;
  TH1F*     t_second_muon_eta;
  TH1F*     w_first_muon_iso;
  TH1F*     b_first_muon_iso;
  TH1F*     c_first_muon_iso;
  TH1F*     t_first_muon_iso;

  TH1F*     w_mass_ee_wide;
  TH1F*     b_mass_ee_wide;
  TH1F*     c_mass_ee_wide;
  TH1F*     t_mass_ee_wide;

  TH1F*     w_mass_mm_wide;
  TH1F*     b_mass_mm_wide;
  TH1F*     c_mass_mm_wide;
  TH1F*     t_mass_mm_wide;

  TH1F*     h_mt_wenu_wide;
  TH1F*     w_mt_wenu_wide;
  TH1F*     b_mt_wenu_wide;
  TH1F*     c_mt_wenu_wide;
  TH1F*     t_mt_wenu_wide;

  TH1F*     h_mt_wmnu_wide;
  TH1F*     w_mt_wmnu_wide;
  TH1F*     b_mt_wmnu_wide;
  TH1F*     c_mt_wmnu_wide;
  TH1F*     t_mt_wmnu_wide;

  TH1F*     h_mass_ee;
  TH1F*     w_mass_ee;
  TH1F*     b_mass_ee;
  TH1F*     c_mass_ee;
  TH1F*     t_mass_ee;

  TH1F*     h_mass_mm;
  TH1F*     w_mass_mm;
  TH1F*     b_mass_mm;
  TH1F*     c_mass_mm;
  TH1F*     t_mass_mm;

  TH1F*     h_mt_wenu;
  TH1F*     w_mt_wenu;
  TH1F*     b_mt_wenu;
  TH1F*     c_mt_wenu;
  TH1F*     t_mt_wenu;
  TH1F*     h_mt_wmnu;
  TH1F*     w_mt_wmnu;
  TH1F*     b_mt_wmnu;
  TH1F*     c_mt_wmnu;
  TH1F*     t_mt_wmnu;

  TH1F*     w_pt_W_wenu;
  TH1F*     b_pt_W_wenu;
  TH1F*     c_pt_W_wenu;
  TH1F*     t_pt_W_wenu;

  TH1F*     w_pt_W_wmnu;
  TH1F*     b_pt_W_wmnu;
  TH1F*     c_pt_W_wmnu;
  TH1F*     t_pt_W_wmnu;

  TH1F*     w_pt_W_wenu_b;
  TH1F*     b_pt_W_wenu_b;
  TH1F*     c_pt_W_wenu_b;
  TH1F*     t_pt_W_wenu_b;

  TH1F*     w_pt_W_wmnu_b;
  TH1F*     b_pt_W_wmnu_b;
  TH1F*     c_pt_W_wmnu_b;
  TH1F*     t_pt_W_wmnu_b;

  TH1F*     w_pt_W_wenu_bb;
  TH1F*     b_pt_W_wenu_bb;
  TH1F*     c_pt_W_wenu_bb;
  TH1F*     t_pt_W_wenu_bb;

  TH1F*     w_pt_W_wmnu_bb;
  TH1F*     b_pt_W_wmnu_bb;
  TH1F*     c_pt_W_wmnu_bb;
  TH1F*     t_pt_W_wmnu_bb;

  TH1F*     w_eta_W_wenu;
  TH1F*     b_eta_W_wenu;
  TH1F*     c_eta_W_wenu;
  TH1F*     t_eta_W_wenu;

  TH1F*     w_eta_W_wmnu;
  TH1F*     b_eta_W_wmnu;
  TH1F*     c_eta_W_wmnu;
  TH1F*     t_eta_W_wmnu;

  TH1F*     w_eta_W_wenu_b;
  TH1F*     b_eta_W_wenu_b;
  TH1F*     c_eta_W_wenu_b;
  TH1F*     t_eta_W_wenu_b;

  TH1F*     w_eta_W_wmnu_b;
  TH1F*     b_eta_W_wmnu_b;
  TH1F*     c_eta_W_wmnu_b;
  TH1F*     t_eta_W_wmnu_b;

  TH1F*     w_eta_W_wenu_bb;
  TH1F*     b_eta_W_wenu_bb;
  TH1F*     c_eta_W_wenu_bb;
  TH1F*     t_eta_W_wenu_bb;

  TH1F*     w_eta_W_wmnu_bb;
  TH1F*     b_eta_W_wmnu_bb;
  TH1F*     c_eta_W_wmnu_bb;
  TH1F*     t_eta_W_wmnu_bb;

  TH1F*     w_pt_Z_ee;
  TH1F*     b_pt_Z_ee;
  TH1F*     c_pt_Z_ee;
  TH1F*     t_pt_Z_ee;

  TH1F*     w_pt_Z_mm;
  TH1F*     b_pt_Z_mm;
  TH1F*     c_pt_Z_mm;
  TH1F*     t_pt_Z_mm;

  TH1F*     w_single_pt_Z_ee_b;
  TH1F*     b_single_pt_Z_ee_b;
  TH1F*     c_single_pt_Z_ee_b;
  TH1F*     t_single_pt_Z_ee_b;

  TH1F*     w_single_pt_Z_mm_b;
  TH1F*     b_single_pt_Z_mm_b;
  TH1F*     c_single_pt_Z_mm_b;
  TH1F*     t_single_pt_Z_mm_b;

  TH1F*     w_mass_ee_b_wide;	// at least one b jet in the event
  TH1F*     b_mass_ee_b_wide;
  TH1F*     c_mass_ee_b_wide;
  TH1F*     t_mass_ee_b_wide;

  TH1F*     w_mass_mm_b_wide;
  TH1F*     b_mass_mm_b_wide;
  TH1F*     c_mass_mm_b_wide;
  TH1F*     t_mass_mm_b_wide;

  TH1F*     h_mt_wenu_b_wide;	// at least one b jet in the event
  TH1F*     w_mt_wenu_b_wide;
  TH1F*     b_mt_wenu_b_wide;
  TH1F*     c_mt_wenu_b_wide;
  TH1F*     t_mt_wenu_b_wide;

  TH1F*     h_mt_wmnu_b_wide;
  TH1F*     w_mt_wmnu_b_wide;
  TH1F*     b_mt_wmnu_b_wide;
  TH1F*     c_mt_wmnu_b_wide;
  TH1F*     t_mt_wmnu_b_wide;

  TH1F*     h_mt_wenu_bb_wide;	// at least one b jet in the event
  TH1F*     w_mt_wenu_bb_wide;
  TH1F*     b_mt_wenu_bb_wide;
  TH1F*     c_mt_wenu_bb_wide;
  TH1F*     t_mt_wenu_bb_wide;

  TH1F*     h_mt_wmnu_bb_wide;
  TH1F*     w_mt_wmnu_bb_wide;
  TH1F*     b_mt_wmnu_bb_wide;
  TH1F*     c_mt_wmnu_bb_wide;
  TH1F*     t_mt_wmnu_bb_wide;

  TH1F*     w_mass_wenu_blepton;	// at least one b jet in the event
  TH1F*     b_mass_wenu_blepton;
  TH1F*     c_mass_wenu_blepton;
  TH1F*     t_mass_wenu_blepton;

  TH1F*     w_mass_wmnu_blepton;
  TH1F*     b_mass_wmnu_blepton;
  TH1F*     c_mass_wmnu_blepton;
  TH1F*     t_mass_wmnu_blepton;

  TH1F*     w_mass_wenu_blepton_b;	// at least one b jet in the event
  TH1F*     b_mass_wenu_blepton_b;
  TH1F*     c_mass_wenu_blepton_b;
  TH1F*     t_mass_wenu_blepton_b;

  TH1F*     w_mass_wmnu_blepton_b;
  TH1F*     b_mass_wmnu_blepton_b;
  TH1F*     c_mass_wmnu_blepton_b;
  TH1F*     t_mass_wmnu_blepton_b;

  TH1F*     w_mass_wenu_blepton_bb;	// at least one b jet in the event
  TH1F*     b_mass_wenu_blepton_bb;
  TH1F*     c_mass_wenu_blepton_bb;
  TH1F*     t_mass_wenu_blepton_bb;

  TH1F*     w_mass_wmnu_blepton_bb;
  TH1F*     b_mass_wmnu_blepton_bb;
  TH1F*     c_mass_wmnu_blepton_bb;
  TH1F*     t_mass_wmnu_blepton_bb;

  TH1F*     w_mass_ee_b;	// at least one b jet in the event
  TH1F*     b_mass_ee_b;
  TH1F*     c_mass_ee_b;
  TH1F*     t_mass_ee_b;

  TH1F*     w_mass_mm_b;
  TH1F*     b_mass_mm_b;
  TH1F*     c_mass_mm_b;
  TH1F*     t_mass_mm_b;

  TH1F*     w_mass_ee_bb;	// at least one b jet in the event
  TH1F*     b_mass_ee_bb;
  TH1F*     c_mass_ee_bb;
  TH1F*     t_mass_ee_bb;

  TH1F*     w_mass_mm_bb;
  TH1F*     b_mass_mm_bb;
  TH1F*     c_mass_mm_bb;
  TH1F*     t_mass_mm_bb;

  TH1F*     h_mt_wenu_b;	// at least one b jet in the event
  TH1F*     w_mt_wenu_b;	// at least one b jet in the event
  TH1F*     b_mt_wenu_b;
  TH1F*     c_mt_wenu_b;
  TH1F*     t_mt_wenu_b;

  TH1F*     h_mt_wmnu_b;
  TH1F*     w_mt_wmnu_b;
  TH1F*     b_mt_wmnu_b;
  TH1F*     c_mt_wmnu_b;
  TH1F*     t_mt_wmnu_b;

  TH1F*     h_mt_wenu_bb;	// at least one b jet in the event
  TH1F*     w_mt_wenu_bb;	// at least one b jet in the event
  TH1F*     b_mt_wenu_bb;
  TH1F*     c_mt_wenu_bb;
  TH1F*     t_mt_wenu_bb;

  TH1F*     h_mt_wmnu_bb;
  TH1F*     w_mt_wmnu_bb;
  TH1F*     b_mt_wmnu_bb;
  TH1F*     c_mt_wmnu_bb;
  TH1F*     t_mt_wmnu_bb;

  TH1F*     w_mass_Zj_ee;
  TH1F*     b_mass_Zj_ee;
  TH1F*     c_mass_Zj_ee;
  TH1F*     t_mass_Zj_ee;

  TH1F*     w_mass_Zj_mm;
  TH1F*     b_mass_Zj_mm;
  TH1F*     c_mass_Zj_mm;
  TH1F*     t_mass_Zj_mm;

  TH1F*     w_mass_Zj_ee_b;
  TH1F*     b_mass_Zj_ee_b;
  TH1F*     c_mass_Zj_ee_b;
  TH1F*     t_mass_Zj_ee_b;

  TH1F*     w_mass_Zj_mm_b;
  TH1F*     b_mass_Zj_mm_b;
  TH1F*     c_mass_Zj_mm_b;
  TH1F*     t_mass_Zj_mm_b;

  TH1F*     w_mass_Zj_ee_bb;
  TH1F*     b_mass_Zj_ee_bb;
  TH1F*     c_mass_Zj_ee_bb;
  TH1F*     t_mass_Zj_ee_bb;

  TH1F*     w_mass_Zj_mm_bb;
  TH1F*     b_mass_Zj_mm_bb;
  TH1F*     c_mass_Zj_mm_bb;
  TH1F*     t_mass_Zj_mm_bb;

  TH1F*     w_pt_Z_ee_b;
  TH1F*     b_pt_Z_ee_b;
  TH1F*     c_pt_Z_ee_b;
  TH1F*     t_pt_Z_ee_b;

  TH1F*     w_pt_Z_mm_b;
  TH1F*     b_pt_Z_mm_b;
  TH1F*     c_pt_Z_mm_b;
  TH1F*     t_pt_Z_mm_b;

  TH1F*     w_pt_Z_ee_bb;
  TH1F*     b_pt_Z_ee_bb;
  TH1F*     c_pt_Z_ee_bb;
  TH1F*     t_pt_Z_ee_bb;

  TH1F*     w_pt_Z_mm_bb;
  TH1F*     b_pt_Z_mm_bb;
  TH1F*     c_pt_Z_mm_bb;
  TH1F*     t_pt_Z_mm_bb;

  TH1F*     w_delta_ee;
  TH1F*     b_delta_ee;
  TH1F*     c_delta_ee;
  TH1F*     t_delta_ee;
  TH1F*     w_delta_ee_b;
  TH1F*     b_delta_ee_b;
  TH1F*     c_delta_ee_b;
  TH1F*     t_delta_ee_b;
  TH1F*     w_delta_ee_bb;
  TH1F*     b_delta_ee_bb;
  TH1F*     c_delta_ee_bb;
  TH1F*     t_delta_ee_bb;

  TH1F*     w_delta_mm;
  TH1F*     b_delta_mm;
  TH1F*     c_delta_mm;
  TH1F*     t_delta_mm;
  TH1F*     w_delta_mm_b;
  TH1F*     b_delta_mm_b;
  TH1F*     c_delta_mm_b;
  TH1F*     t_delta_mm_b;
  TH1F*     w_delta_mm_bb;
  TH1F*     b_delta_mm_bb;
  TH1F*     c_delta_mm_bb;
  TH1F*     t_delta_mm_bb;

  TH1F*     w_delta_wenu;
  TH1F*     b_delta_wenu;
  TH1F*     c_delta_wenu;
  TH1F*     t_delta_wenu;
  TH1F*     w_delta_wenu_b;
  TH1F*     b_delta_wenu_b;
  TH1F*     c_delta_wenu_b;
  TH1F*     t_delta_wenu_b;
  TH1F*     w_delta_wenu_bb;
  TH1F*     b_delta_wenu_bb;
  TH1F*     c_delta_wenu_bb;
  TH1F*     t_delta_wenu_bb;
  TH1F*     w_delta_wenu_2b;
  TH1F*     b_delta_wenu_2b;
  TH1F*     c_delta_wenu_2b;
  TH1F*     t_delta_wenu_2b;

  TH1F*     w_delta_wmnu;
  TH1F*     b_delta_wmnu;
  TH1F*     c_delta_wmnu;
  TH1F*     t_delta_wmnu;
  TH1F*     w_delta_wmnu_b;
  TH1F*     b_delta_wmnu_b;
  TH1F*     c_delta_wmnu_b;
  TH1F*     t_delta_wmnu_b;
  TH1F*     w_delta_wmnu_bb;
  TH1F*     b_delta_wmnu_bb;
  TH1F*     c_delta_wmnu_bb;
  TH1F*     t_delta_wmnu_bb;
  TH1F*     w_delta_wmnu_2b;
  TH1F*     b_delta_wmnu_2b;
  TH1F*     c_delta_wmnu_2b;
  TH1F*     t_delta_wmnu_2b;

  TH1F*     w_deltaR_wenu;
  TH1F*     b_deltaR_wenu;
  TH1F*     c_deltaR_wenu;
  TH1F*     t_deltaR_wenu;
  TH1F*     w_deltaR_wenu_b;
  TH1F*     b_deltaR_wenu_b;
  TH1F*     c_deltaR_wenu_b;
  TH1F*     t_deltaR_wenu_b;
  TH1F*     w_deltaR_wenu_bb;
  TH1F*     b_deltaR_wenu_bb;
  TH1F*     c_deltaR_wenu_bb;
  TH1F*     t_deltaR_wenu_bb;
  TH1F*     w_deltaR_wenu_2b;
  TH1F*     b_deltaR_wenu_2b;
  TH1F*     c_deltaR_wenu_2b;
  TH1F*     t_deltaR_wenu_2b;

  TH1F*     w_deltaR_wmnu;
  TH1F*     b_deltaR_wmnu;
  TH1F*     c_deltaR_wmnu;
  TH1F*     t_deltaR_wmnu;
  TH1F*     w_deltaR_wmnu_b;
  TH1F*     b_deltaR_wmnu_b;
  TH1F*     c_deltaR_wmnu_b;
  TH1F*     t_deltaR_wmnu_b;
  TH1F*     w_deltaR_wmnu_bb;
  TH1F*     b_deltaR_wmnu_bb;
  TH1F*     c_deltaR_wmnu_bb;
  TH1F*     t_deltaR_wmnu_bb;
  TH1F*     w_deltaR_wmnu_2b;
  TH1F*     b_deltaR_wmnu_2b;
  TH1F*     c_deltaR_wmnu_2b;
  TH1F*     t_deltaR_wmnu_2b;

  TH1F*     w_single_delta_ee_b;
  TH1F*     b_single_delta_ee_b;
  TH1F*     c_single_delta_ee_b;
  TH1F*     t_single_delta_ee_b;

  TH1F*     w_single_delta_mm_b;
  TH1F*     b_single_delta_mm_b;
  TH1F*     c_single_delta_mm_b;
  TH1F*     t_single_delta_mm_b;

  TH1F*     w_single_delta_wenu_b;
  TH1F*     b_single_delta_wenu_b;
  TH1F*     c_single_delta_wenu_b;
  TH1F*     t_single_delta_wenu_b;

  TH1F*     w_single_delta_wmnu_b;
  TH1F*     b_single_delta_wmnu_b;
  TH1F*     c_single_delta_wmnu_b;
  TH1F*     t_single_delta_wmnu_b;

  TH1F*     w_single_deltaR_wenu_b;
  TH1F*     b_single_deltaR_wenu_b;
  TH1F*     c_single_deltaR_wenu_b;
  TH1F*     t_single_deltaR_wenu_b;

  TH1F*     w_single_deltaR_wmnu_b;
  TH1F*     b_single_deltaR_wmnu_b;
  TH1F*     c_single_deltaR_wmnu_b;
  TH1F*     t_single_deltaR_wmnu_b;

  TH1F*     h_secondvtx_N;
  TH1F*     w_secondvtx_N;
  TH1F*     b_secondvtx_N;
  TH1F*     c_secondvtx_N;
  TH1F*     t_secondvtx_N;

  TH1F*     w_secondvtx_N_zoom;
  TH1F*     b_secondvtx_N_zoom;
  TH1F*     c_secondvtx_N_zoom;
  TH1F*     t_secondvtx_N_zoom;

  TH1F*     w_secondvtx_N_mass;
  TH1F*     b_secondvtx_N_mass;
  TH1F*     c_secondvtx_N_mass;
  TH1F*     t_secondvtx_N_mass;

  TH1F*     w_secondvtx_N_nomass;
  TH1F*     b_secondvtx_N_nomass;
  TH1F*     c_secondvtx_N_nomass;
  TH1F*     t_secondvtx_N_nomass;

  TH1F*     w_SVTX_mass_jet;
  TH1F*     b_SVTX_mass_jet;
  TH1F*     c_SVTX_mass_jet;
  TH1F*     t_SVTX_mass_jet;

  TH1F*	    w_SVTX_mass_trk;
  TH1F*     b_SVTX_mass_trk;
  TH1F*     c_SVTX_mass_trk;
  TH1F*     t_SVTX_mass_trk;

  TH1F*     w_SVTX_mass_jet_b;
  TH1F*     b_SVTX_mass_jet_b;
  TH1F*     c_SVTX_mass_jet_b;
  TH1F*     t_SVTX_mass_jet_b;

  TH1F*	    w_SVTX_mass_trk_b;
  TH1F*     b_SVTX_mass_trk_b;
  TH1F*     c_SVTX_mass_trk_b;
  TH1F*     t_SVTX_mass_trk_b;

  TH1F*     w_SVTX_mass_jet_bb;
  TH1F*     b_SVTX_mass_jet_bb;
  TH1F*     c_SVTX_mass_jet_bb;
  TH1F*     t_SVTX_mass_jet_bb;

  TH1F*	    w_SVTX_mass_trk_bb;
  TH1F*     b_SVTX_mass_trk_bb;
  TH1F*     c_SVTX_mass_trk_bb;
  TH1F*     t_SVTX_mass_trk_bb;

  TH1F*     w_SVTX_mass;
  TH1F*     b_SVTX_mass;
  TH1F*     c_SVTX_mass;
  TH1F*     t_SVTX_mass;

  TH1F*     w_SVTX_mass_b;
  TH1F*     b_SVTX_mass_b;
  TH1F*     c_SVTX_mass_b;
  TH1F*     t_SVTX_mass_b;

  TH1F*     w_SVTX_mass_bb;
  TH1F*     b_SVTX_mass_bb;
  TH1F*     c_SVTX_mass_bb;
  TH1F*     t_SVTX_mass_bb;

  TH1F*     w_BJP;
  TH1F*     b_BJP;
  TH1F*     c_BJP;
  TH1F*     t_BJP;

  TH1F*     w_BJP_b;
  TH1F*     b_BJP_b;
  TH1F*     c_BJP_b;
  TH1F*     t_BJP_b;

  TH1F*     w_BJP_bb;
  TH1F*     b_BJP_bb;
  TH1F*     c_BJP_bb;
  TH1F*     t_BJP_bb;

  TH1F*     w_JBP;
  TH1F*     b_JBP;
  TH1F*     c_JBP;
  TH1F*     t_JBP;

  TH1F*     w_JBP_b;
  TH1F*     b_JBP_b;
  TH1F*     c_JBP_b;
  TH1F*     t_JBP_b;

  TH1F*     w_JBP_bb;
  TH1F*     b_JBP_bb;
  TH1F*     c_JBP_bb;
  TH1F*     t_JBP_bb;

  TH1F*     w_BJP0;
  TH1F*     b_BJP0;
  TH1F*     c_BJP0;
  TH1F*     t_BJP0;

  TH1F*     w_BJP1;
  TH1F*     b_BJP1;
  TH1F*     c_BJP1;
  TH1F*     t_BJP1;

  TH1F*     w_BJP2;
  TH1F*     b_BJP2;
  TH1F*     c_BJP2;
  TH1F*     t_BJP2;

  TH1F*     w_BJP_mass;
  TH1F*     b_BJP_mass;
  TH1F*     c_BJP_mass;
  TH1F*     t_BJP_mass;

  TH1F*     w_JBP_mass;
  TH1F*     b_JBP_mass;
  TH1F*     c_JBP_mass;
  TH1F*     t_JBP_mass;

  TH1F*     w_BJP_nomass;
  TH1F*     b_BJP_nomass;
  TH1F*     c_BJP_nomass;
  TH1F*     t_BJP_nomass;

  TH1F*     w_JBP_nomass;
  TH1F*     b_JBP_nomass;
  TH1F*     c_JBP_nomass;
  TH1F*     t_JBP_nomass;

  TH1F*     w_BJP_mass_b;
  TH1F*     b_BJP_mass_b;
  TH1F*     c_BJP_mass_b;
  TH1F*     t_BJP_mass_b;

  TH1F*     w_JBP_mass_b;
  TH1F*     b_JBP_mass_b;
  TH1F*     c_JBP_mass_b;
  TH1F*     t_JBP_mass_b;

  TH1F*     w_BJP_nomass_b;
  TH1F*     b_BJP_nomass_b;
  TH1F*     c_BJP_nomass_b;
  TH1F*     t_BJP_nomass_b;

  TH1F*     w_JBP_nomass_b;
  TH1F*     b_JBP_nomass_b;
  TH1F*     c_JBP_nomass_b;
  TH1F*     t_JBP_nomass_b;

  TH1F*     w_BJP_mass_bb;
  TH1F*     b_BJP_mass_bb;
  TH1F*     c_BJP_mass_bb;
  TH1F*     t_BJP_mass_bb;

  TH1F*     w_JBP_mass_bb;
  TH1F*     b_JBP_mass_bb;
  TH1F*     c_JBP_mass_bb;
  TH1F*     t_JBP_mass_bb;

  TH1F*     w_BJP_nomass_bb;
  TH1F*     b_BJP_nomass_bb;
  TH1F*     c_BJP_nomass_bb;
  TH1F*     t_BJP_nomass_bb;

  TH1F*     w_JBP_nomass_bb;
  TH1F*     b_JBP_nomass_bb;
  TH1F*     c_JBP_nomass_bb;
  TH1F*     t_JBP_nomass_bb;

  TH1F*     w_Ht;
  TH1F*     b_Ht;
  TH1F*     c_Ht;
  TH1F*     t_Ht;

  TH1F*     w_single_Ht_b;
  TH1F*     b_single_Ht_b;
  TH1F*     c_single_Ht_b;
  TH1F*     t_single_Ht_b;

  TH1F*     w_Ht_b; // at least one b jet in the event
  TH1F*     b_Ht_b;
  TH1F*     c_Ht_b;
  TH1F*     t_Ht_b;

  TH1F*     w_Ht_bb; // at least one b jet in the event
  TH1F*     b_Ht_bb;
  TH1F*     c_Ht_bb;
  TH1F*     t_Ht_bb;

  TH1F*     h_MET;
  TH1F*     w_MET;
  TH1F*     b_MET;
  TH1F*     c_MET;
  TH1F*     t_MET;
  TH1F*     h_MET_phi;
  TH1F*     w_MET_phi;
  TH1F*     b_MET_phi;
  TH1F*     c_MET_phi;
  TH1F*     t_MET_phi;
  TH1F*     w_MET_sign;
  TH1F*     b_MET_sign;
  TH1F*     c_MET_sign;
  TH1F*     t_MET_sign;

  TH1F*     h_MET_b;
  TH1F*     w_MET_b;
  TH1F*     b_MET_b;
  TH1F*     c_MET_b;
  TH1F*     t_MET_b;
  TH1F*     h_MET_phi_b;
  TH1F*     w_MET_phi_b;
  TH1F*     b_MET_phi_b;
  TH1F*     c_MET_phi_b;
  TH1F*     t_MET_phi_b;
  TH1F*     w_MET_sign_b;
  TH1F*     b_MET_sign_b;
  TH1F*     c_MET_sign_b;
  TH1F*     t_MET_sign_b;

  TH1F*     h_MET_bb;
  TH1F*     w_MET_bb;
  TH1F*     b_MET_bb;
  TH1F*     c_MET_bb;
  TH1F*     t_MET_bb;
  TH1F*     h_MET_phi_bb;
  TH1F*     w_MET_phi_bb;
  TH1F*     b_MET_phi_bb;
  TH1F*     c_MET_phi_bb;
  TH1F*     t_MET_phi_bb;
  TH1F*     w_MET_sign_bb;
  TH1F*     b_MET_sign_bb;
  TH1F*     c_MET_sign_bb;
  TH1F*     t_MET_sign_bb;

  TH1F*     w_Afb;

  TH1F*     h_scaleFactor_first_ele;
  TH1F*     b_scaleFactor_first_ele;
  TH1F*     h_scaleFactor_first_muon;
  TH1F*     b_scaleFactor_first_muon;
  TH1F*     h_scaleFactor_second_ele;
  TH1F*     b_scaleFactor_second_ele;
  TH1F*     h_scaleFactor_second_muon;
  TH1F*     b_scaleFactor_second_muon;

  TH1F*     h_JEC_uncert;

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
WbAnalyzer::WbAnalyzer (const edm::ParameterSet & iConfig) {

  pileupMC_         = iConfig.getUntrackedParameter < std::string > ("pileupMC", "S10");
  pileupDT_         = iConfig.getUntrackedParameter < std::string > ("pileupDT", "");
  lepton_           = iConfig.getUntrackedParameter < std::string > ("lepton", "electron");
  par_              = iConfig.getUntrackedParameter <double> ("JEC", 0);
  par2_             = iConfig.getUntrackedParameter <double> ("JER", 0);
  path_             = iConfig.getUntrackedParameter < std::string > ("path", "/gpfs/cms/users/schizzi/Wbb2012/test/");
  icut_             = iConfig.getUntrackedParameter <unsigned int> ("icut", 0);
  usePartonFlavour_ = iConfig.getUntrackedParameter <bool> ("usePartonFlavour", false);
  pcut_             = iConfig.getUntrackedParameter <bool> ("pcut", false);
  useDeltaR_        = iConfig.getUntrackedParameter <bool> ("useDeltaR", false);
  // now do what ever initialization is needed
  edm::Service < TFileService > fs;

  h_nmult0 =            fs->make < TH1F > ("h_nmult0", "h_nmult0", 8, -0.5, 7.5);
  h_nmult1 =            fs->make < TH1F > ("h_nmult1", "h_nmult1", 8, -0.5, 7.5);

  h_jetmultiplicity =   fs->make < TH1F > ("h_jetmultiplicity", "h_jetmultiplicity;N_jets", 8, 0.5, 8.5);

  h_eventYields =   fs->make < TH1F > ("h_eventYields", "h_eventYields;selection", 8, 0.5, 8.5);
  w_eventYields =   fs->make < TH1F > ("w_eventYields", "w_eventYields;selection", 8, 0.5, 8.5);

  h_trgMatchEle =   fs->make < TH1F > ("h_trgMatchEle", "h_trgMatchEle;selection", 8, -0.5, 7.5);
  h_trgMatchMuo =   fs->make < TH1F > ("h_trgMatchMuo", "h_trgMatchMuo;selection", 8, -0.5, 7.5);

  h_pu_weights =        fs->make < TH1F > ("h_pu_weights",      "h_pu_weights;PU weight", 10, 0, 10);

  h_tracks =            fs->make < TH1F > ("h_tracks",          "h_tracks;N_tracks", 100, 0, 2500);
  w_tracks = 	        fs->make < TH1F > ("w_tracks",          "w_tracks;N_tracks", 100, 0, 2500);
  h_recoVTX =           fs->make < TH1F > ("h_recoVTX",         "h_recoVTX;N_vtx", 40, 0., 40.);
  w_recoVTX =           fs->make < TH1F > ("w_recoVTX",         "w_recoVTX;N_vtx", 40, 0., 40.);

  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity", "w_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  b_jetmultiplicity =   fs->make < TH1F > ("b_jetmultiplicity", "b_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  c_jetmultiplicity =   fs->make < TH1F > ("c_jetmultiplicity", "c_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  t_jetmultiplicity =   fs->make < TH1F > ("t_jetmultiplicity", "t_jetmultiplicity;N_jets", 8, 0.5, 8.5);

  w_first_jet_pt =      fs->make < TH1F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 20, 20., 220.);
  h_first_jet_pt =      fs->make < TH1F > ("h_first_jet_pt",    "h_first_jet_pt;P_t [GeV]", 20, 20., 220.);
  b_first_jet_pt =      fs->make < TH1F > ("b_first_jet_pt",    "b_first_jet_pt;P_t [GeV]", 20, 20., 220.);
  c_first_jet_pt =      fs->make < TH1F > ("c_first_jet_pt",    "c_first_jet_pt;P_t [GeV]", 20, 20., 220.);
  t_first_jet_pt =      fs->make < TH1F > ("t_first_jet_pt",    "t_first_jet_pt;P_t [GeV]", 20, 20., 220.);
  w_first_jet_eta =     fs->make < TH1F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 20, -2.4, 2.4);
  b_first_jet_eta =     fs->make < TH1F > ("b_first_jet_eta",   "b_first_jet_eta;Eta", 20, -2.4, 2.4);
  c_first_jet_eta =     fs->make < TH1F > ("c_first_jet_eta",   "c_first_jet_eta;Eta", 20, -2.4, 2.4);
  t_first_jet_eta =     fs->make < TH1F > ("t_first_jet_eta",   "t_first_jet_eta;Eta", 20, -2.4, 2.4);
  w_first_jet_mass =      fs->make < TH1F > ("w_first_jet_mass",    "w_first_jet_mass;Mass [GeV]", 18, 0., 36.);
  b_first_jet_mass =      fs->make < TH1F > ("b_first_jet_mass",    "b_first_jet_mass;Mass [GeV]", 18, 0., 36.);
  c_first_jet_mass =      fs->make < TH1F > ("c_first_jet_mass",    "c_first_jet_mass;Mass [GeV]", 18, 0., 36.);
  t_first_jet_mass =      fs->make < TH1F > ("t_first_jet_mass",    "t_first_jet_mass;Mass [GeV]", 18, 0., 36.);
  h_second_jet_pt =     fs->make < TH1F > ("h_second_jet_pt",   "h_second_jet_pt;P_t [GeV]", 20, 20., 220.);
  w_second_jet_pt =     fs->make < TH1F > ("w_second_jet_pt",   "w_second_jet_pt;P_t [GeV]", 20, 20., 220.);
  b_second_jet_pt =     fs->make < TH1F > ("b_second_jet_pt",   "b_second_jet_pt;P_t [GeV]", 20, 20., 220.);
  c_second_jet_pt =     fs->make < TH1F > ("c_second_jet_pt",   "c_second_jet_pt;P_t [GeV]", 20, 20., 220.);
  t_second_jet_pt =     fs->make < TH1F > ("t_second_jet_pt",   "t_second_jet_pt;P_t [GeV]", 20, 20., 220.);
  w_second_jet_eta =    fs->make < TH1F > ("w_second_jet_eta",  "w_second_jet_eta;Eta", 20, -2.4, 2.4);
  b_second_jet_eta =    fs->make < TH1F > ("b_second_jet_eta",  "b_second_jet_eta;Eta", 20, -2.4, 2.4);
  c_second_jet_eta =    fs->make < TH1F > ("c_second_jet_eta",  "c_second_jet_eta;Eta", 20, -2.4, 2.4);
  t_second_jet_eta =    fs->make < TH1F > ("t_second_jet_eta",  "t_second_jet_eta;Eta", 20, -2.4, 2.4);
  w_second_jet_mass =      fs->make < TH1F > ("w_second_jet_mass",    "w_second_jet_mass;Mass [GeV]", 18, 0., 36.);
  b_second_jet_mass =      fs->make < TH1F > ("b_second_jet_mass",    "b_second_jet_mass;Mass [GeV]", 18, 0., 36.);
  c_second_jet_mass =      fs->make < TH1F > ("c_second_jet_mass",    "c_second_jet_mass;Mass [GeV]", 18, 0., 36.);
  t_second_jet_mass =      fs->make < TH1F > ("t_second_jet_mass",    "t_second_jet_mass;Mass [GeV]", 18, 0., 36.);
  w_dijet_pt =      fs->make < TH1F > ("w_dijet_pt",    "w_dijet_pt;P_t [GeV]", 20, 0., 200.);
  h_dijet_pt =      fs->make < TH1F > ("h_dijet_pt",    "h_dijet_pt;P_t [GeV]", 20, 0., 200.);
  b_dijet_pt =      fs->make < TH1F > ("b_dijet_pt",    "b_dijet_pt;P_t [GeV]", 20, 0., 200.);
  c_dijet_pt =      fs->make < TH1F > ("c_dijet_pt",    "c_dijet_pt;P_t [GeV]", 20, 0., 200.);
  t_dijet_pt =      fs->make < TH1F > ("t_dijet_pt",    "t_dijet_pt;P_t [GeV]", 20, 0., 200.);
  w_dijet_eta =     fs->make < TH1F > ("w_dijet_eta",   "w_dijet_eta;Eta", 20, -2.4, 2.4);
  b_dijet_eta =     fs->make < TH1F > ("b_dijet_eta",   "b_dijet_eta;Eta", 20, -2.4, 2.4);
  c_dijet_eta =     fs->make < TH1F > ("c_dijet_eta",   "c_dijet_eta;Eta", 20, -2.4, 2.4);
  t_dijet_eta =     fs->make < TH1F > ("t_dijet_eta",   "t_dijet_eta;Eta", 20, -2.4, 2.4);
  w_dijet_mass =      fs->make < TH1F > ("w_dijet_mass",    "w_dijet_mass;Mass [GeV]", 20, 20., 260.);
  b_dijet_mass =      fs->make < TH1F > ("b_dijet_mass",    "b_dijet_mass;Mass [GeV]", 20, 20., 260.);
  c_dijet_mass =      fs->make < TH1F > ("c_dijet_mass",    "c_dijet_mass;Mass [GeV]", 20, 20., 260.);
  t_dijet_mass =      fs->make < TH1F > ("t_dijet_mass",    "t_dijet_mass;Mass [GeV]", 20, 20., 260.);

  h_first_jet_pt_b =    fs->make < TH1F > ("h_first_jet_pt_b",   "h_first_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  w_first_jet_pt_b =    fs->make < TH1F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  b_first_jet_pt_b =    fs->make < TH1F > ("b_first_jet_pt_b",   "b_first_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  c_first_jet_pt_b =    fs->make < TH1F > ("c_first_jet_pt_b",   "c_first_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  t_first_jet_pt_b =    fs->make < TH1F > ("t_first_jet_pt_b",   "t_first_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  w_first_jet_eta_b =   fs->make < TH1F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 20, -2.4, 2.4);
  b_first_jet_eta_b =   fs->make < TH1F > ("b_first_jet_eta_b",  "b_first_jet_eta_b;Eta", 20, -2.4, 2.4);
  c_first_jet_eta_b =   fs->make < TH1F > ("c_first_jet_eta_b",  "c_first_jet_eta_b;Eta", 20, -2.4, 2.4);
  t_first_jet_eta_b =   fs->make < TH1F > ("t_first_jet_eta_b",  "t_first_jet_eta_b;Eta", 20, -2.4, 2.4);
  w_first_jet_mass_b =      fs->make < TH1F > ("w_first_jet_mass_b",    "w_first_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  b_first_jet_mass_b =      fs->make < TH1F > ("b_first_jet_mass_b",    "b_first_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  c_first_jet_mass_b =      fs->make < TH1F > ("c_first_jet_mass_b",    "c_first_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  t_first_jet_mass_b =      fs->make < TH1F > ("t_first_jet_mass_b",    "t_first_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  h_second_jet_pt_b =   fs->make < TH1F > ("h_second_jet_pt_b",  "h_second_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  w_second_jet_pt_b =   fs->make < TH1F > ("w_second_jet_pt_b",  "w_second_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  b_second_jet_pt_b =   fs->make < TH1F > ("b_second_jet_pt_b",  "b_second_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  c_second_jet_pt_b =   fs->make < TH1F > ("c_second_jet_pt_b",  "c_second_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  t_second_jet_pt_b =   fs->make < TH1F > ("t_second_jet_pt_b",  "t_second_jet_pt_b;P_t [GeV]", 20, 20., 220.);
  w_second_jet_eta_b =  fs->make < TH1F > ("w_second_jet_eta_b", "w_second_jet_eta_b;Eta", 20, -2.4, 2.4);
  b_second_jet_eta_b =  fs->make < TH1F > ("b_second_jet_eta_b", "b_second_jet_eta_b;Eta", 20, -2.4, 2.4);
  c_second_jet_eta_b =  fs->make < TH1F > ("c_second_jet_eta_b", "c_second_jet_eta_b;Eta", 20, -2.4, 2.4);
  t_second_jet_eta_b =  fs->make < TH1F > ("t_second_jet_eta_b", "t_second_jet_eta_b;Eta", 20, -2.4, 2.4);
  w_second_jet_mass_b =      fs->make < TH1F > ("w_second_jet_mass_b",    "w_second_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  b_second_jet_mass_b =      fs->make < TH1F > ("b_second_jet_mass_b",    "b_second_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  c_second_jet_mass_b =      fs->make < TH1F > ("c_second_jet_mass_b",    "c_second_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  t_second_jet_mass_b =      fs->make < TH1F > ("t_second_jet_mass_b",    "t_second_jet_mass_b;Mass [GeV]", 18, 0., 36.);
  h_dijet_pt_b =    fs->make < TH1F > ("h_dijet_pt_b",   "h_dijet_pt_b;P_t [GeV]", 20, 0., 200.);
  w_dijet_pt_b =    fs->make < TH1F > ("w_dijet_pt_b",   "w_dijet_pt_b;P_t [GeV]", 20, 0., 200.);
  b_dijet_pt_b =    fs->make < TH1F > ("b_dijet_pt_b",   "b_dijet_pt_b;P_t [GeV]", 20, 0., 200.);
  c_dijet_pt_b =    fs->make < TH1F > ("c_dijet_pt_b",   "c_dijet_pt_b;P_t [GeV]", 20, 0., 200.);
  t_dijet_pt_b =    fs->make < TH1F > ("t_dijet_pt_b",   "t_dijet_pt_b;P_t [GeV]", 20, 0., 200.);
  w_dijet_eta_b =   fs->make < TH1F > ("w_dijet_eta_b",  "w_dijet_eta_b;Eta", 20, -2.4, 2.4);
  b_dijet_eta_b =   fs->make < TH1F > ("b_dijet_eta_b",  "b_dijet_eta_b;Eta", 20, -2.4, 2.4);
  c_dijet_eta_b =   fs->make < TH1F > ("c_dijet_eta_b",  "c_dijet_eta_b;Eta", 20, -2.4, 2.4);
  t_dijet_eta_b =   fs->make < TH1F > ("t_dijet_eta_b",  "t_dijet_eta_b;Eta", 20, -2.4, 2.4);
  w_dijet_mass_b =      fs->make < TH1F > ("w_dijet_mass_b",    "w_dijet_mass_b;Mass [GeV]", 20, 20., 260.);
  b_dijet_mass_b =      fs->make < TH1F > ("b_dijet_mass_b",    "b_dijet_mass_b;Mass [GeV]", 20, 20., 260.);
  c_dijet_mass_b =      fs->make < TH1F > ("c_dijet_mass_b",    "c_dijet_mass_b;Mass [GeV]", 20, 20., 260.);
  t_dijet_mass_b =      fs->make < TH1F > ("t_dijet_mass_b",    "t_dijet_mass_b;Mass [GeV]", 20, 20., 260.);

  h_first_jet_pt_bb =    fs->make < TH1F > ("h_first_jet_pt_bb",   "h_first_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  w_first_jet_pt_bb =    fs->make < TH1F > ("w_first_jet_pt_bb",   "w_first_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  b_first_jet_pt_bb =    fs->make < TH1F > ("b_first_jet_pt_bb",   "b_first_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  c_first_jet_pt_bb =    fs->make < TH1F > ("c_first_jet_pt_bb",   "c_first_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  t_first_jet_pt_bb =    fs->make < TH1F > ("t_first_jet_pt_bb",   "t_first_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  w_first_jet_eta_bb =   fs->make < TH1F > ("w_first_jet_eta_bb",  "w_first_jet_eta_bb;Eta", 20, -2.4, 2.4);
  b_first_jet_eta_bb =   fs->make < TH1F > ("b_first_jet_eta_bb",  "b_first_jet_eta_bb;Eta", 20, -2.4, 2.4);
  c_first_jet_eta_bb =   fs->make < TH1F > ("c_first_jet_eta_bb",  "c_first_jet_eta_bb;Eta", 20, -2.4, 2.4);
  t_first_jet_eta_bb =   fs->make < TH1F > ("t_first_jet_eta_bb",  "t_first_jet_eta_bb;Eta", 20, -2.4, 2.4);
  w_first_jet_mass_bb =      fs->make < TH1F > ("w_first_jet_mass_bb",    "w_first_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  b_first_jet_mass_bb =      fs->make < TH1F > ("b_first_jet_mass_bb",    "b_first_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  c_first_jet_mass_bb =      fs->make < TH1F > ("c_first_jet_mass_bb",    "c_first_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  t_first_jet_mass_bb =      fs->make < TH1F > ("t_first_jet_mass_bb",    "t_first_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  h_second_jet_pt_bb =   fs->make < TH1F > ("h_second_jet_pt_bb",  "h_second_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  w_second_jet_pt_bb =   fs->make < TH1F > ("w_second_jet_pt_bb",  "w_second_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  b_second_jet_pt_bb =   fs->make < TH1F > ("b_second_jet_pt_bb",  "b_second_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  c_second_jet_pt_bb =   fs->make < TH1F > ("c_second_jet_pt_bb",  "c_second_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  t_second_jet_pt_bb =   fs->make < TH1F > ("t_second_jet_pt_bb",  "t_second_jet_pt_bb;P_t [GeV]", 20, 20., 220.);
  w_second_jet_eta_bb =  fs->make < TH1F > ("w_second_jet_eta_bb", "w_second_jet_eta_bb;Eta", 20, -2.4, 2.4);
  b_second_jet_eta_bb =  fs->make < TH1F > ("b_second_jet_eta_bb", "b_second_jet_eta_bb;Eta", 20, -2.4, 2.4);
  c_second_jet_eta_bb =  fs->make < TH1F > ("c_second_jet_eta_bb", "c_second_jet_eta_bb;Eta", 20, -2.4, 2.4);
  t_second_jet_eta_bb =  fs->make < TH1F > ("t_second_jet_eta_bb", "t_second_jet_eta_bb;Eta", 20, -2.4, 2.4);
  w_second_jet_mass_bb =      fs->make < TH1F > ("w_second_jet_mass_bb",    "w_second_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  b_second_jet_mass_bb =      fs->make < TH1F > ("b_second_jet_mass_bb",    "b_second_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  c_second_jet_mass_bb =      fs->make < TH1F > ("c_second_jet_mass_bb",    "c_second_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  t_second_jet_mass_bb =      fs->make < TH1F > ("t_second_jet_mass_bb",    "t_second_jet_mass_bb;Mass [GeV]", 18, 0., 36.);
  h_dijet_pt_bb =    fs->make < TH1F > ("h_dijet_pt_bb",   "h_dijet_pt_bb;P_t [GeV]", 20, 0., 200.);
  w_dijet_pt_bb =    fs->make < TH1F > ("w_dijet_pt_bb",   "w_dijet_pt_bb;P_t [GeV]", 20, 0., 200.);
  b_dijet_pt_bb =    fs->make < TH1F > ("b_dijet_pt_bb",   "b_dijet_pt_bb;P_t [GeV]", 20, 0., 200.);
  c_dijet_pt_bb =    fs->make < TH1F > ("c_dijet_pt_bb",   "c_dijet_pt_bb;P_t [GeV]", 20, 0., 200.);
  t_dijet_pt_bb =    fs->make < TH1F > ("t_dijet_pt_bb",   "t_dijet_pt_bb;P_t [GeV]", 20, 0., 200.);
  w_dijet_eta_bb =   fs->make < TH1F > ("w_dijet_eta_bb",  "w_dijet_eta_bb;Eta", 20, -2.4, 2.4);
  b_dijet_eta_bb =   fs->make < TH1F > ("b_dijet_eta_bb",  "b_dijet_eta_bb;Eta", 20, -2.4, 2.4);
  c_dijet_eta_bb =   fs->make < TH1F > ("c_dijet_eta_bb",  "c_dijet_eta_bb;Eta", 20, -2.4, 2.4);
  t_dijet_eta_bb =   fs->make < TH1F > ("t_dijet_eta_bb",  "t_dijet_eta_bb;Eta", 20, -2.4, 2.4);
  w_dijet_mass_bb =      fs->make < TH1F > ("w_dijet_mass_bb",    "w_dijet_mass_bb;Mass [GeV]", 20, 20., 260.);
  b_dijet_mass_bb =      fs->make < TH1F > ("b_dijet_mass_bb",    "b_dijet_mass_bb;Mass [GeV]", 20, 20., 260.);
  c_dijet_mass_bb =      fs->make < TH1F > ("c_dijet_mass_bb",    "c_dijet_mass_bb;Mass [GeV]", 20, 20., 260.);
  t_dijet_mass_bb =      fs->make < TH1F > ("t_dijet_mass_bb",    "t_dijet_mass_bb;Mass [GeV]", 20, 20., 260.);

  w_bjetmultiplicity =  fs->make < TH1F > ("w_bjetmultiplicity", "w_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  b_bjetmultiplicity =  fs->make < TH1F > ("b_bjetmultiplicity", "b_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  c_bjetmultiplicity =  fs->make < TH1F > ("c_bjetmultiplicity", "c_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  t_bjetmultiplicity =  fs->make < TH1F > ("t_bjetmultiplicity", "t_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);

  h_first_bjet_pt =     fs->make < TH1F > ("h_first_bjet_pt",    "h_first_bjet_pt;P_t [GeV]", 20, 20., 220.);
  w_first_bjet_pt =     fs->make < TH1F > ("w_first_bjet_pt",    "w_first_bjet_pt;P_t [GeV]", 20, 20., 220.);
  b_first_bjet_pt =     fs->make < TH1F > ("b_first_bjet_pt",    "b_first_bjet_pt;P_t [GeV]", 20, 20., 220.);
  c_first_bjet_pt =     fs->make < TH1F > ("c_first_bjet_pt",    "c_first_bjet_pt;P_t [GeV]", 20, 20., 220.);
  t_first_bjet_pt =     fs->make < TH1F > ("t_first_bjet_pt",    "t_first_bjet_pt;P_t [GeV]", 20, 20., 220.);
  w_first_bjet_eta =    fs->make < TH1F > ("w_first_bjet_eta",   "w_first_bjet_eta;Eta", 20, -2.4, 2.4);
  b_first_bjet_eta =    fs->make < TH1F > ("b_first_bjet_eta",   "b_first_bjet_eta;Eta", 20, -2.4, 2.4);
  c_first_bjet_eta =    fs->make < TH1F > ("c_first_bjet_eta",   "c_first_bjet_eta;Eta", 20, -2.4, 2.4);
  t_first_bjet_eta =    fs->make < TH1F > ("t_first_bjet_eta",   "t_first_bjet_eta;Eta", 20, -2.4, 2.4);
  w_first_bjet_mass =      fs->make < TH1F > ("w_first_bjet_mass",    "w_first_bjet_mass;Mass [GeV]", 18, 0., 36.);
  b_first_bjet_mass =      fs->make < TH1F > ("b_first_bjet_mass",    "b_first_bjet_mass;Mass [GeV]", 18, 0., 36.);
  c_first_bjet_mass =      fs->make < TH1F > ("c_first_bjet_mass",    "c_first_bjet_mass;Mass [GeV]", 18, 0., 36.);
  t_first_bjet_mass =      fs->make < TH1F > ("t_first_bjet_mass",    "t_first_bjet_mass;Mass [GeV]", 18, 0., 36.);

  w_single_bjet_pt =    fs->make < TH1F > ("w_single_bjet_pt",    "w_single_bjet_pt;P_t [GeV]", 20, 20., 220.);
  b_single_bjet_pt =    fs->make < TH1F > ("b_single_bjet_pt",    "b_single_bjet_pt;P_t [GeV]", 20, 20., 220.);
  c_single_bjet_pt =    fs->make < TH1F > ("c_single_bjet_pt",    "c_single_bjet_pt;P_t [GeV]", 20, 20., 220.);
  t_single_bjet_pt =    fs->make < TH1F > ("t_single_bjet_pt",    "t_single_bjet_pt;P_t [GeV]", 20, 20., 220.);
  w_single_bjet_eta =   fs->make < TH1F > ("w_single_bjet_eta",   "w_single_bjet_eta;Eta", 20, -2.4, 2.4);
  b_single_bjet_eta =   fs->make < TH1F > ("b_single_bjet_eta",   "b_single_bjet_eta;Eta", 20, -2.4, 2.4);
  c_single_bjet_eta =   fs->make < TH1F > ("c_single_bjet_eta",   "c_single_bjet_eta;Eta", 20, -2.4, 2.4);
  t_single_bjet_eta =   fs->make < TH1F > ("t_single_bjet_eta",   "t_single_bjet_eta;Eta", 20, -2.4, 2.4);
  w_single_bjet_mass =      fs->make < TH1F > ("w_single_bjet_mass",    "w_single_bjet_mass;Mass [GeV]", 18, 0., 36.);
  b_single_bjet_mass =      fs->make < TH1F > ("b_single_bjet_mass",    "b_single_bjet_mass;Mass [GeV]", 18, 0., 36.);
  c_single_bjet_mass =      fs->make < TH1F > ("c_single_bjet_mass",    "c_single_bjet_mass;Mass [GeV]", 18, 0., 36.);
  t_single_bjet_mass =      fs->make < TH1F > ("t_single_bjet_mass",    "t_single_bjet_mass;Mass [GeV]", 18, 0., 36.);

  h_second_bjet_pt =    fs->make < TH1F > ("h_second_bjet_pt",   "h_second_bjet_pt;P_t [GeV]", 20, 20., 220.);
  w_second_bjet_pt =    fs->make < TH1F > ("w_second_bjet_pt",   "w_second_bjet_pt;P_t [GeV]", 20, 20., 220.);
  b_second_bjet_pt =    fs->make < TH1F > ("b_second_bjet_pt",   "b_second_bjet_pt;P_t [GeV]", 20, 20., 220.);
  c_second_bjet_pt =    fs->make < TH1F > ("c_second_bjet_pt",   "c_second_bjet_pt;P_t [GeV]", 20, 20., 220.);
  t_second_bjet_pt =    fs->make < TH1F > ("t_second_bjet_pt",   "t_second_bjet_pt;P_t [GeV]", 20, 20., 220.);
  w_second_bjet_eta =   fs->make < TH1F > ("w_second_bjet_eta",  "w_second_bjet_eta;Eta", 20, -2.4, 2.4);
  b_second_bjet_eta =   fs->make < TH1F > ("b_second_bjet_eta",  "b_second_bjet_eta;Eta", 20, -2.4, 2.4);
  c_second_bjet_eta =   fs->make < TH1F > ("c_second_bjet_eta",  "c_second_bjet_eta;Eta", 20, -2.4, 2.4);
  t_second_bjet_eta =   fs->make < TH1F > ("t_second_bjet_eta",  "t_second_bjet_eta;Eta", 20, -2.4, 2.4);
  w_second_bjet_mass =      fs->make < TH1F > ("w_second_bjet_mass",    "w_second_bjet_mass;Mass [GeV]", 18, 0., 36.);
  b_second_bjet_mass =      fs->make < TH1F > ("b_second_bjet_mass",    "b_second_bjet_mass;Mass [GeV]", 18, 0., 36.);
  c_second_bjet_mass =      fs->make < TH1F > ("c_second_bjet_mass",    "c_second_bjet_mass;Mass [GeV]", 18, 0., 36.);
  t_second_bjet_mass =      fs->make < TH1F > ("t_second_bjet_mass",    "t_second_bjet_mass;Mass [GeV]", 18, 0., 36.);

  w_first_ele_pt =      fs->make < TH1F > ("w_first_ele_pt",     "w_first_ele_pt;P_t [GeV]", 20, 20., 220.);
  w_first_ele_pt_b =    fs->make < TH1F > ("w_first_ele_pt_b",   "w_first_ele_pt_b;P_t [GeV]", 20, 20., 220.);
  w_first_ele_pt_bb =   fs->make < TH1F > ("w_first_ele_pt_bb",  "w_first_ele_pt_bb;P_t [GeV]", 20, 20., 220.);
  h_first_ele_pt =      fs->make < TH1F > ("h_first_ele_pt",     "h_first_ele_pt;P_t [GeV]", 20, 20., 220.);
  h_first_ele_pt_b =    fs->make < TH1F > ("h_first_ele_pt_b",   "h_first_ele_pt_b;P_t [GeV]", 20, 20., 220.);
  h_first_ele_pt_bb =   fs->make < TH1F > ("h_first_ele_pt_bb",  "h_first_ele_pt_bb;P_t [GeV]", 20, 20., 220.);
  b_first_ele_pt =      fs->make < TH1F > ("b_first_ele_pt",     "b_first_ele_pt;P_t [GeV]", 20, 20., 220.);
  c_first_ele_pt =      fs->make < TH1F > ("c_first_ele_pt",     "c_first_ele_pt;P_t [GeV]", 20, 20., 220.);
  t_first_ele_pt =      fs->make < TH1F > ("t_first_ele_pt",     "t_first_ele_pt;P_t [GeV]", 20, 20., 220.);
  w_second_ele_pt =     fs->make < TH1F > ("w_second_ele_pt",    "w_second_ele_pt;P_t [GeV]", 20, 20., 220.);
  b_second_ele_pt =     fs->make < TH1F > ("b_second_ele_pt",    "b_second_ele_pt;P_t [GeV]", 20, 20., 220.);
  c_second_ele_pt =     fs->make < TH1F > ("c_second_ele_pt",    "c_second_ele_pt;P_t [GeV]", 20, 20., 220.);
  t_second_ele_pt =     fs->make < TH1F > ("t_second_ele_pt",    "t_second_ele_pt;P_t [GeV]", 20, 20., 220.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",    "w_first_muon_pt;P_t [GeV]", 20, 20., 220.);
  w_first_muon_pt_b =   fs->make < TH1F > ("w_first_muon_pt_b",  "w_first_muon_pt_b [GeV]", 20, 20., 220.);
  w_first_muon_pt_bb =  fs->make < TH1F > ("w_first_muon_pt_bb",  "w_first_muon_pt_bb [GeV]", 20, 20., 220.);
  h_first_muon_pt =     fs->make < TH1F > ("h_first_muon_pt",    "h_first_muon_pt;P_t [GeV]", 20, 20., 220.);
  h_first_muon_pt_b =   fs->make < TH1F > ("h_first_muon_pt_b",  "h_first_muon_pt_b [GeV]", 20, 20., 220.);
  h_first_muon_pt_bb =  fs->make < TH1F > ("h_first_muon_pt_bb", "h_first_muon_pt_bb [GeV]", 20, 20., 220.);
  b_first_muon_pt =     fs->make < TH1F > ("b_first_muon_pt",    "b_first_muon_pt;P_t [GeV]", 20, 20., 220.);
  c_first_muon_pt =     fs->make < TH1F > ("c_first_muon_pt",    "c_first_muon_pt;P_t [GeV]", 20, 20., 220.);
  t_first_muon_pt =     fs->make < TH1F > ("t_first_muon_pt",    "t_first_muon_pt;P_t [GeV]", 20, 20., 220.);
  w_second_muon_pt =    fs->make < TH1F > ("w_second_muon_pt",   "w_second_muon_pt;P_t [GeV]", 20, 20., 220.);
  b_second_muon_pt =    fs->make < TH1F > ("b_second_muon_pt",   "b_second_muon_pt;P_t [GeV]", 20, 20., 220.);
  c_second_muon_pt =    fs->make < TH1F > ("c_second_muon_pt",   "c_second_muon_pt;P_t [GeV]", 20, 20., 220.);
  t_second_muon_pt =    fs->make < TH1F > ("t_second_muon_pt",   "t_second_muon_pt;P_t [GeV]", 20, 20., 220.);
  w_first_ele_eta =     fs->make < TH1F > ("w_first_ele_eta",    "w_first_ele_eta;Eta", 20, -2.4, 2.4);
  b_first_ele_eta =     fs->make < TH1F > ("b_first_ele_eta",    "b_first_ele_eta;Eta", 20, -2.4, 2.4);
  c_first_ele_eta =     fs->make < TH1F > ("c_first_ele_eta",    "c_first_ele_eta;Eta", 20, -2.4, 2.4);
  t_first_ele_eta =     fs->make < TH1F > ("t_first_ele_eta",    "t_first_ele_eta;Eta", 20, -2.4, 2.4);
  w_second_ele_eta =    fs->make < TH1F > ("w_second_ele_eta",   "w_second_ele_eta;Eta", 20, -2.4, 2.4);
  b_second_ele_eta =    fs->make < TH1F > ("b_second_ele_eta",   "b_second_ele_eta;Eta", 20, -2.4, 2.4);
  c_second_ele_eta =    fs->make < TH1F > ("c_second_ele_eta",   "c_second_ele_eta;Eta", 20, -2.4, 2.4);
  t_second_ele_eta =    fs->make < TH1F > ("t_second_ele_eta",   "t_second_ele_eta;Eta", 20, -2.4, 2.4);
  w_first_ele_iso =     fs->make < TH1F > ("w_first_ele_iso",    "w_first_ele_iso;Iso", 25, 0., 0.25);
  b_first_ele_iso =     fs->make < TH1F > ("b_first_ele_iso",    "b_first_ele_iso;Iso", 25, 0., 0.25);
  c_first_ele_iso =     fs->make < TH1F > ("c_first_ele_iso",    "c_first_ele_iso;Iso", 25, 0., 0.25);
  t_first_ele_iso =     fs->make < TH1F > ("t_first_ele_iso",    "t_first_ele_iso;Iso", 25, 0., 0.25);
  w_first_muon_eta =    fs->make < TH1F > ("w_first_muon_eta",   "w_first_muon_eta;Eta", 20, -2.4, 2.4);
  b_first_muon_eta =    fs->make < TH1F > ("b_first_muon_eta",   "b_first_muon_eta;Eta", 20, -2.4, 2.4);
  c_first_muon_eta =    fs->make < TH1F > ("c_first_muon_eta",   "c_first_muon_eta;Eta", 20, -2.4, 2.4);
  t_first_muon_eta =    fs->make < TH1F > ("t_first_muon_eta",   "t_first_muon_eta;Eta", 20, -2.4, 2.4);
  w_second_muon_eta =   fs->make < TH1F > ("w_second_muon_eta",  "w_second_muon_eta;Eta", 20, -2.4, 2.4);
  b_second_muon_eta =   fs->make < TH1F > ("b_second_muon_eta",  "b_second_muon_eta;Eta", 20, -2.4, 2.4);
  c_second_muon_eta =   fs->make < TH1F > ("c_second_muon_eta",  "c_second_muon_eta;Eta", 20, -2.4, 2.4);
  t_second_muon_eta =   fs->make < TH1F > ("t_second_muon_eta",  "t_second_muon_eta;Eta", 20, -2.4, 2.4);
  w_first_muon_iso =    fs->make < TH1F > ("w_first_muon_iso",   "w_first_muon_iso;Iso", 25, 0., 0.25);
  b_first_muon_iso =    fs->make < TH1F > ("b_first_muon_iso",   "b_first_muon_iso;Iso", 25, 0., 0.25);
  c_first_muon_iso =    fs->make < TH1F > ("c_first_muon_iso",   "c_first_muon_iso;Iso", 25, 0., 0.25);
  t_first_muon_iso =    fs->make < TH1F > ("t_first_muon_iso",   "t_first_muon_iso;Iso", 25, 0., 0.25);

  w_mass_ee_wide =      fs->make < TH1F > ("w_mass_ee_wide",    "w_mass_ee_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_ee_wide =      fs->make < TH1F > ("b_mass_ee_wide",    "b_mass_ee_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_ee_wide =      fs->make < TH1F > ("c_mass_ee_wide",    "c_mass_ee_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_ee_wide =      fs->make < TH1F > ("t_mass_ee_wide",    "t_mass_ee_wide;Mass [GeV]", 40, 50., 200.);

  w_mass_mm_wide =      fs->make < TH1F > ("w_mass_mm_wide",    "w_mass_mm_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_mm_wide =      fs->make < TH1F > ("b_mass_mm_wide",    "b_mass_mm_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_mm_wide =      fs->make < TH1F > ("c_mass_mm_wide",    "c_mass_mm_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_mm_wide =      fs->make < TH1F > ("t_mass_mm_wide",    "t_mass_mm_wide;Mass [GeV]", 40, 50., 200.);

  h_mt_wenu_wide =      fs->make < TH1F > ("h_mt_wenu_wide",    "h_mt_wenu_wide;M_{T} [GeV]", 50, 0., 200.);
  w_mt_wenu_wide =      fs->make < TH1F > ("w_mt_wenu_wide",    "w_mt_wenu_wide;M_{T} [GeV]", 50, 0., 200.);
  b_mt_wenu_wide =      fs->make < TH1F > ("b_mt_wenu_wide",    "b_mt_wenu_wide;M_{T} [GeV]", 50, 0., 200.);
  c_mt_wenu_wide =      fs->make < TH1F > ("c_mt_wenu_wide",    "c_mt_wenu_wide;M_{T} [GeV]", 50, 0., 200.);
  t_mt_wenu_wide =      fs->make < TH1F > ("t_mt_wenu_wide",    "t_mt_wenu_wide;M_{T} [GeV]", 50, 0., 200.);

  h_mt_wmnu_wide =      fs->make < TH1F > ("h_mt_wmnu_wide",    "h_mt_wmnu_wide;M_{T} [GeV]", 50, 0., 200.);
  w_mt_wmnu_wide =      fs->make < TH1F > ("w_mt_wmnu_wide",    "w_mt_wmnu_wide;M_{T} [GeV]", 50, 0., 200.);
  b_mt_wmnu_wide =      fs->make < TH1F > ("b_mt_wmnu_wide",    "b_mt_wmnu_wide;M_{T} [GeV]", 50, 0., 200.);
  c_mt_wmnu_wide =      fs->make < TH1F > ("c_mt_wmnu_wide",    "c_mt_wmnu_wide;M_{T} [GeV]", 50, 0., 200.);
  t_mt_wmnu_wide =      fs->make < TH1F > ("t_mt_wmnu_wide",    "t_mt_wmnu_wide;M_{T} [GeV]", 50, 0., 200.);

  h_mass_ee =           fs->make < TH1F > ("h_mass_ee",         "h_mass_ee;Mass [GeV]", 80, 71., 111.);
  w_mass_ee =           fs->make < TH1F > ("w_mass_ee",         "w_mass_ee;Mass [GeV]", 80, 71., 111.);
  b_mass_ee =           fs->make < TH1F > ("b_mass_ee",         "b_mass_ee;Mass [GeV]", 80, 71., 111.);
  c_mass_ee =           fs->make < TH1F > ("c_mass_ee",         "c_mass_ee;Mass [GeV]", 80, 71., 111.);
  t_mass_ee =           fs->make < TH1F > ("t_mass_ee",         "t_mass_ee;Mass [GeV]", 80, 71., 111.);

  h_mass_mm =           fs->make < TH1F > ("h_mass_mm",         "h_mass_mm;Mass [GeV]", 80, 71., 111.);
  w_mass_mm =           fs->make < TH1F > ("w_mass_mm",         "w_mass_mm;Mass [GeV]", 80, 71., 111.);
  b_mass_mm =           fs->make < TH1F > ("b_mass_mm",         "b_mass_mm;Mass [GeV]", 80, 71., 111.);
  c_mass_mm =           fs->make < TH1F > ("c_mass_mm",         "c_mass_mm;Mass [GeV]", 80, 71., 111.);
  t_mass_mm =           fs->make < TH1F > ("t_mass_mm",         "t_mass_mm;Mass [GeV]", 80, 71., 111.);

  h_mt_wenu =           fs->make < TH1F > ("h_mt_wenu",         "h_mt_wenu;M_{T} [GeV]", 50, 0., 200.);
  w_mt_wenu =           fs->make < TH1F > ("w_mt_wenu",         "w_mt_wenu;M_{T} [GeV]", 20, 45., 205.);
  b_mt_wenu =           fs->make < TH1F > ("b_mt_wenu",         "b_mt_wenu;M_{T} [GeV]", 20, 45., 205.);
  c_mt_wenu =           fs->make < TH1F > ("c_mt_wenu",         "c_mt_wenu;M_{T} [GeV]", 20, 45., 205.);
  t_mt_wenu =           fs->make < TH1F > ("t_mt_wenu",         "t_mt_wenu;M_{T} [GeV]", 20, 45., 205.);
  h_mt_wmnu =           fs->make < TH1F > ("h_mt_wmnu",         "h_mt_wmnu;M_{T} [GeV]", 20, 45., 205.);
  w_mt_wmnu =           fs->make < TH1F > ("w_mt_wmnu",         "w_mt_wmnu;M_{T} [GeV]", 20, 45., 205.);
  b_mt_wmnu =           fs->make < TH1F > ("b_mt_wmnu",         "b_mt_wmnu;M_{T} [GeV]", 20, 45., 205.);
  c_mt_wmnu =           fs->make < TH1F > ("c_mt_wmnu",         "c_mt_wmnu;M_{T} [GeV]", 20, 45., 205.);
  t_mt_wmnu =           fs->make < TH1F > ("t_mt_wmnu",         "t_mt_wmnu;M_{T} [GeV]", 20, 45., 205.);

  w_pt_W_wenu =           fs->make < TH1F > ("w_pt_W_wenu",         "w_pt_W_wenu;P_t [GeV]", 20, 0., 200.);
  b_pt_W_wenu =           fs->make < TH1F > ("b_pt_W_wenu",         "b_pt_W_wenu;P_t [GeV]", 20, 0., 200.);
  c_pt_W_wenu =           fs->make < TH1F > ("c_pt_W_wenu",         "c_pt_W_wenu;P_t [GeV]", 20, 0., 200.);
  t_pt_W_wenu =           fs->make < TH1F > ("t_pt_W_wenu",         "t_pt_W_wenu;P_t [GeV]", 20, 0., 200.);

  w_pt_W_wmnu =           fs->make < TH1F > ("w_pt_W_wmnu",         "w_pt_W_wmnu;P_t [GeV]", 20, 0., 200.);
  b_pt_W_wmnu =           fs->make < TH1F > ("b_pt_W_wmnu",         "b_pt_W_wmnu;P_t [GeV]", 20, 0., 200.);
  c_pt_W_wmnu =           fs->make < TH1F > ("c_pt_W_wmnu",         "c_pt_W_wmnu;P_t [GeV]", 20, 0., 200.);
  t_pt_W_wmnu =           fs->make < TH1F > ("t_pt_W_wmnu",         "t_pt_W_wmnu;P_t [GeV]", 20, 0., 200.);

  w_pt_W_wenu_b =           fs->make < TH1F > ("w_pt_W_wenu_b",         "w_pt_W_wenu_b;P_t [GeV]", 20, 0., 200.);
  b_pt_W_wenu_b =           fs->make < TH1F > ("b_pt_W_wenu_b",         "b_pt_W_wenu_b;P_t [GeV]", 20, 0., 200.);
  c_pt_W_wenu_b =           fs->make < TH1F > ("c_pt_W_wenu_b",         "c_pt_W_wenu_b;P_t [GeV]", 20, 0., 200.);
  t_pt_W_wenu_b =           fs->make < TH1F > ("t_pt_W_wenu_b",         "t_pt_W_wenu_b;P_t [GeV]", 20, 0., 200.);

  w_pt_W_wmnu_b =           fs->make < TH1F > ("w_pt_W_wmnu_b",         "w_pt_W_wmnu_b;P_t [GeV]", 20, 0., 200.);
  b_pt_W_wmnu_b =           fs->make < TH1F > ("b_pt_W_wmnu_b",         "b_pt_W_wmnu_b;P_t [GeV]", 20, 0., 200.);
  c_pt_W_wmnu_b =           fs->make < TH1F > ("c_pt_W_wmnu_b",         "c_pt_W_wmnu_b;P_t [GeV]", 20, 0., 200.);
  t_pt_W_wmnu_b =           fs->make < TH1F > ("t_pt_W_wmnu_b",         "t_pt_W_wmnu_b;P_t [GeV]", 20, 0., 200.);

  w_pt_W_wenu_bb =           fs->make < TH1F > ("w_pt_W_wenu_bb",         "w_pt_W_wenu_bb;P_t [GeV]", 20, 0., 200.);
  b_pt_W_wenu_bb =           fs->make < TH1F > ("b_pt_W_wenu_bb",         "b_pt_W_wenu_bb;P_t [GeV]", 20, 0., 200.);
  c_pt_W_wenu_bb =           fs->make < TH1F > ("c_pt_W_wenu_bb",         "c_pt_W_wenu_bb;P_t [GeV]", 20, 0., 200.);
  t_pt_W_wenu_bb =           fs->make < TH1F > ("t_pt_W_wenu_bb",         "t_pt_W_wenu_bb;P_t [GeV]", 20, 0., 200.);

  w_pt_W_wmnu_bb =           fs->make < TH1F > ("w_pt_W_wmnu_bb",         "w_pt_W_wmnu_bb;P_t [GeV]", 20, 0., 200.);
  b_pt_W_wmnu_bb =           fs->make < TH1F > ("b_pt_W_wmnu_bb",         "b_pt_W_wmnu_bb;P_t [GeV]", 20, 0., 200.);
  c_pt_W_wmnu_bb =           fs->make < TH1F > ("c_pt_W_wmnu_bb",         "c_pt_W_wmnu_bb;P_t [GeV]", 20, 0., 200.);
  t_pt_W_wmnu_bb =           fs->make < TH1F > ("t_pt_W_wmnu_bb",         "t_pt_W_wmnu_bb;P_t [GeV]", 20, 0., 200.);

  w_eta_W_wenu =           fs->make < TH1F > ("w_eta_W_wenu",         "w_eta_W_wenu;Eta", 20, -2.4, 2.4);
  b_eta_W_wenu =           fs->make < TH1F > ("b_eta_W_wenu",         "b_eta_W_wenu;Eta", 20, -2.4, 2.4);
  c_eta_W_wenu =           fs->make < TH1F > ("c_eta_W_wenu",         "c_eta_W_wenu;Eta", 20, -2.4, 2.4);
  t_eta_W_wenu =           fs->make < TH1F > ("t_eta_W_wenu",         "t_eta_W_wenu;Eta", 20, -2.4, 2.4);

  w_eta_W_wmnu =           fs->make < TH1F > ("w_eta_W_wmnu",         "w_eta_W_wmnu;Eta", 20, -2.4, 2.4);
  b_eta_W_wmnu =           fs->make < TH1F > ("b_eta_W_wmnu",         "b_eta_W_wmnu;Eta", 20, -2.4, 2.4);
  c_eta_W_wmnu =           fs->make < TH1F > ("c_eta_W_wmnu",         "c_eta_W_wmnu;Eta", 20, -2.4, 2.4);
  t_eta_W_wmnu =           fs->make < TH1F > ("t_eta_W_wmnu",         "t_eta_W_wmnu;Eta", 20, -2.4, 2.4);

  w_eta_W_wenu_b =           fs->make < TH1F > ("w_eta_W_wenu_b",         "w_eta_W_wenu_b;Eta", 20, -2.4, 2.4);
  b_eta_W_wenu_b =           fs->make < TH1F > ("b_eta_W_wenu_b",         "b_eta_W_wenu_b;Eta", 20, -2.4, 2.4);
  c_eta_W_wenu_b =           fs->make < TH1F > ("c_eta_W_wenu_b",         "c_eta_W_wenu_b;Eta", 20, -2.4, 2.4);
  t_eta_W_wenu_b =           fs->make < TH1F > ("t_eta_W_wenu_b",         "t_eta_W_wenu_b;Eta", 20, -2.4, 2.4);

  w_eta_W_wmnu_b =           fs->make < TH1F > ("w_eta_W_wmnu_b",         "w_eta_W_wmnu_b;Eta", 20, -2.4, 2.4);
  b_eta_W_wmnu_b =           fs->make < TH1F > ("b_eta_W_wmnu_b",         "b_eta_W_wmnu_b;Eta", 20, -2.4, 2.4);
  c_eta_W_wmnu_b =           fs->make < TH1F > ("c_eta_W_wmnu_b",         "c_eta_W_wmnu_b;Eta", 20, -2.4, 2.4);
  t_eta_W_wmnu_b =           fs->make < TH1F > ("t_eta_W_wmnu_b",         "t_eta_W_wmnu_b;Eta", 20, -2.4, 2.4);

  w_eta_W_wenu_bb =           fs->make < TH1F > ("w_eta_W_wenu_bb",         "w_eta_W_wenu_bb;Eta", 20, -2.4, 2.4);
  b_eta_W_wenu_bb =           fs->make < TH1F > ("b_eta_W_wenu_bb",         "b_eta_W_wenu_bb;Eta", 20, -2.4, 2.4);
  c_eta_W_wenu_bb =           fs->make < TH1F > ("c_eta_W_wenu_bb",         "c_eta_W_wenu_bb;Eta", 20, -2.4, 2.4);
  t_eta_W_wenu_bb =           fs->make < TH1F > ("t_eta_W_wenu_bb",         "t_eta_W_wenu_bb;Eta", 20, -2.4, 2.4);

  w_eta_W_wmnu_bb =           fs->make < TH1F > ("w_eta_W_wmnu_bb",         "w_eta_W_wmnu_bb;Eta", 20, -2.4, 2.4);
  b_eta_W_wmnu_bb =           fs->make < TH1F > ("b_eta_W_wmnu_bb",         "b_eta_W_wmnu_bb;Eta", 20, -2.4, 2.4);
  c_eta_W_wmnu_bb =           fs->make < TH1F > ("c_eta_W_wmnu_bb",         "c_eta_W_wmnu_bb;Eta", 20, -2.4, 2.4);
  t_eta_W_wmnu_bb =           fs->make < TH1F > ("t_eta_W_wmnu_bb",         "t_eta_W_wmnu_bb;Eta", 20, -2.4, 2.4);

  w_pt_Z_ee =           fs->make < TH1F > ("w_pt_Z_ee",         "w_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_ee =           fs->make < TH1F > ("b_pt_Z_ee",         "b_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_ee =           fs->make < TH1F > ("c_pt_Z_ee",         "c_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_ee =           fs->make < TH1F > ("t_pt_Z_ee",         "t_pt_Z_ee;P_t [GeV]", 40, 0., 400.);

  w_pt_Z_mm =           fs->make < TH1F > ("w_pt_Z_mm",         "w_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_mm =           fs->make < TH1F > ("b_pt_Z_mm",         "b_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_mm =           fs->make < TH1F > ("c_pt_Z_mm",         "c_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_mm =           fs->make < TH1F > ("t_pt_Z_mm",         "t_pt_Z_mm;P_t [GeV]", 40, 0., 400.);

  w_single_pt_Z_ee_b =  fs->make < TH1F > ("w_single_pt_Z_ee_b",       "w_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  b_single_pt_Z_ee_b =  fs->make < TH1F > ("b_single_pt_Z_ee_b",       "b_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  c_single_pt_Z_ee_b =  fs->make < TH1F > ("c_single_pt_Z_ee_b",       "c_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  t_single_pt_Z_ee_b =  fs->make < TH1F > ("t_single_pt_Z_ee_b",       "t_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);

  w_single_pt_Z_mm_b =  fs->make < TH1F > ("w_single_pt_Z_mm_b",       "w_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  b_single_pt_Z_mm_b =  fs->make < TH1F > ("b_single_pt_Z_mm_b",       "b_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  c_single_pt_Z_mm_b =  fs->make < TH1F > ("c_single_pt_Z_mm_b",       "c_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  t_single_pt_Z_mm_b =  fs->make < TH1F > ("t_single_pt_Z_mm_b",       "t_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);

  w_mass_ee_b_wide =    fs->make < TH1F > ("w_mass_ee_b_wide",  "w_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_ee_b_wide =    fs->make < TH1F > ("b_mass_ee_b_wide",  "b_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_ee_b_wide =    fs->make < TH1F > ("c_mass_ee_b_wide",  "c_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_ee_b_wide =    fs->make < TH1F > ("t_mass_ee_b_wide",  "t_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);

  w_mass_mm_b_wide =    fs->make < TH1F > ("w_mass_mm_b_wide",  "w_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_mm_b_wide =    fs->make < TH1F > ("b_mass_mm_b_wide",  "b_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_mm_b_wide =    fs->make < TH1F > ("c_mass_mm_b_wide",  "c_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_mm_b_wide =    fs->make < TH1F > ("t_mass_mm_b_wide",  "t_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);

  h_mt_wenu_b_wide =    fs->make < TH1F > ("h_mt_wenu_b_wide",  "h_mt_wenu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  w_mt_wenu_b_wide =    fs->make < TH1F > ("w_mt_wenu_b_wide",  "w_mt_wenu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  b_mt_wenu_b_wide =    fs->make < TH1F > ("b_mt_wenu_b_wide",  "b_mt_wenu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  c_mt_wenu_b_wide =    fs->make < TH1F > ("c_mt_wenu_b_wide",  "c_mt_wenu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  t_mt_wenu_b_wide =    fs->make < TH1F > ("t_mt_wenu_b_wide",  "t_mt_wenu_b_wide;M_{T} [GeV]", 50, 0., 200.);

  h_mt_wmnu_b_wide =    fs->make < TH1F > ("h_mt_wmnu_b_wide",  "h_mt_wmnu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  w_mt_wmnu_b_wide =    fs->make < TH1F > ("w_mt_wmnu_b_wide",  "w_mt_wmnu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  b_mt_wmnu_b_wide =    fs->make < TH1F > ("b_mt_wmnu_b_wide",  "b_mt_wmnu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  c_mt_wmnu_b_wide =    fs->make < TH1F > ("c_mt_wmnu_b_wide",  "c_mt_wmnu_b_wide;M_{T} [GeV]", 50, 0., 200.);
  t_mt_wmnu_b_wide =    fs->make < TH1F > ("t_mt_wmnu_b_wide",  "t_mt_wmnu_b_wide;M_{T} [GeV]", 50, 0., 200.);

  h_mt_wenu_bb_wide =   fs->make < TH1F > ("h_mt_wenu_bb_wide", "h_mt_wenu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  w_mt_wenu_bb_wide =   fs->make < TH1F > ("w_mt_wenu_bb_wide", "w_mt_wenu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  b_mt_wenu_bb_wide =   fs->make < TH1F > ("b_mt_wenu_bb_wide", "b_mt_wenu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  c_mt_wenu_bb_wide =   fs->make < TH1F > ("c_mt_wenu_bb_wide", "c_mt_wenu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  t_mt_wenu_bb_wide =   fs->make < TH1F > ("t_mt_wenu_bb_wide", "t_mt_wenu_bb_wide;M_{T} [GeV]", 50, 0., 200.);

  h_mt_wmnu_bb_wide =   fs->make < TH1F > ("h_mt_wmnu_bb_wide", "h_mt_wmnu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  w_mt_wmnu_bb_wide =   fs->make < TH1F > ("w_mt_wmnu_bb_wide", "w_mt_wmnu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  b_mt_wmnu_bb_wide =   fs->make < TH1F > ("b_mt_wmnu_bb_wide", "b_mt_wmnu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  c_mt_wmnu_bb_wide =   fs->make < TH1F > ("c_mt_wmnu_bb_wide", "c_mt_wmnu_bb_wide;M_{T} [GeV]", 50, 0., 200.);
  t_mt_wmnu_bb_wide =   fs->make < TH1F > ("t_mt_wmnu_bb_wide", "t_mt_wmnu_bb_wide;M_{T} [GeV]", 50, 0., 200.);

  w_mass_wenu_blepton =      fs->make < TH1F > ("w_mass_wenu_blepton",    "w_mass_wenu_blepton;M_{T} [GeV]", 50, 0., 250.);
  b_mass_wenu_blepton =      fs->make < TH1F > ("b_mass_wenu_blepton",    "b_mass_wenu_blepton;M_{T} [GeV]", 50, 0., 250.);
  c_mass_wenu_blepton =      fs->make < TH1F > ("c_mass_wenu_blepton",    "c_mass_wenu_blepton;M_{T} [GeV]", 50, 0., 250.);
  t_mass_wenu_blepton =      fs->make < TH1F > ("t_mass_wenu_blepton",    "t_mass_wenu_blepton;M_{T} [GeV]", 50, 0., 250.);

  w_mass_wmnu_blepton =      fs->make < TH1F > ("w_mass_wmnu_blepton",    "w_mass_wmnu_blepton;M_{T} [GeV]", 50, 0., 250.);
  b_mass_wmnu_blepton =      fs->make < TH1F > ("b_mass_wmnu_blepton",    "b_mass_wmnu_blepton;M_{T} [GeV]", 50, 0., 250.);
  c_mass_wmnu_blepton =      fs->make < TH1F > ("c_mass_wmnu_blepton",    "c_mass_wmnu_blepton;M_{T} [GeV]", 50, 0., 250.);
  t_mass_wmnu_blepton =      fs->make < TH1F > ("t_mass_wmnu_blepton",    "t_mass_wmnu_blepton;M_{T} [GeV]", 50, 0., 250.);

  w_mass_wenu_blepton_b =      fs->make < TH1F > ("w_mass_wenu_blepton_b",    "w_mass_wenu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  b_mass_wenu_blepton_b =      fs->make < TH1F > ("b_mass_wenu_blepton_b",    "b_mass_wenu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  c_mass_wenu_blepton_b =      fs->make < TH1F > ("c_mass_wenu_blepton_b",    "c_mass_wenu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  t_mass_wenu_blepton_b =      fs->make < TH1F > ("t_mass_wenu_blepton_b",    "t_mass_wenu_blepton_b;M_{T} [GeV]", 50, 0., 250.);

  w_mass_wmnu_blepton_b =      fs->make < TH1F > ("w_mass_wmnu_blepton_b",    "w_mass_wmnu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  b_mass_wmnu_blepton_b =      fs->make < TH1F > ("b_mass_wmnu_blepton_b",    "b_mass_wmnu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  c_mass_wmnu_blepton_b =      fs->make < TH1F > ("c_mass_wmnu_blepton_b",    "c_mass_wmnu_blepton_b;M_{T} [GeV]", 50, 0., 250.);
  t_mass_wmnu_blepton_b =      fs->make < TH1F > ("t_mass_wmnu_blepton_b",    "t_mass_wmnu_blepton_b;M_{T} [GeV]", 50, 0., 250.);

  w_mass_wenu_blepton_bb =      fs->make < TH1F > ("w_mass_wenu_blepton_bb",    "w_mass_wenu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  b_mass_wenu_blepton_bb =      fs->make < TH1F > ("b_mass_wenu_blepton_bb",    "b_mass_wenu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  c_mass_wenu_blepton_bb =      fs->make < TH1F > ("c_mass_wenu_blepton_bb",    "c_mass_wenu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  t_mass_wenu_blepton_bb =      fs->make < TH1F > ("t_mass_wenu_blepton_bb",    "t_mass_wenu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);

  w_mass_wmnu_blepton_bb =      fs->make < TH1F > ("w_mass_wmnu_blepton_bb",    "w_mass_wmnu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  b_mass_wmnu_blepton_bb =      fs->make < TH1F > ("b_mass_wmnu_blepton_bb",    "b_mass_wmnu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  c_mass_wmnu_blepton_bb =      fs->make < TH1F > ("c_mass_wmnu_blepton_bb",    "c_mass_wmnu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);
  t_mass_wmnu_blepton_bb =      fs->make < TH1F > ("t_mass_wmnu_blepton_bb",    "t_mass_wmnu_blepton_bb;M_{T} [GeV]", 50, 0., 250.);

  w_mass_ee_b =         fs->make < TH1F > ("w_mass_ee_b",       "w_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  b_mass_ee_b =         fs->make < TH1F > ("b_mass_ee_b",       "b_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  c_mass_ee_b =         fs->make < TH1F > ("c_mass_ee_b",       "c_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  t_mass_ee_b =         fs->make < TH1F > ("t_mass_ee_b",       "t_mass_mm_b;Mass [GeV]", 80, 71., 111.);

  w_mass_mm_b =         fs->make < TH1F > ("w_mass_mm_b",       "w_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  b_mass_mm_b =         fs->make < TH1F > ("b_mass_mm_b",       "b_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  c_mass_mm_b =         fs->make < TH1F > ("c_mass_mm_b",       "c_mass_mm_b;Mass [GeV]", 80, 71., 111.);
  t_mass_mm_b =         fs->make < TH1F > ("t_mass_mm_b",       "t_mass_mm_b;Mass [GeV]", 80, 71., 111.);

  w_mass_ee_bb =         fs->make < TH1F > ("w_mass_ee_bb",       "w_mass_mm_bb;Mass [GeV]", 80, 71., 111.);
  b_mass_ee_bb =         fs->make < TH1F > ("b_mass_ee_bb",       "b_mass_mm_bb;Mass [GeV]", 80, 71., 111.);
  c_mass_ee_bb =         fs->make < TH1F > ("c_mass_ee_bb",       "c_mass_mm_bb;Mass [GeV]", 80, 71., 111.);
  t_mass_ee_bb =         fs->make < TH1F > ("t_mass_ee_bb",       "t_mass_mm_bb;Mass [GeV]", 80, 71., 111.);

  w_mass_mm_bb =         fs->make < TH1F > ("w_mass_mm_bb",       "w_mass_mm_bb;Mass [GeV]", 80, 71., 111.);
  b_mass_mm_bb =         fs->make < TH1F > ("b_mass_mm_bb",       "b_mass_mm_bb;Mass [GeV]", 80, 71., 111.);
  c_mass_mm_bb =         fs->make < TH1F > ("c_mass_mm_bb",       "c_mass_mm_bb;Mass [GeV]", 80, 71., 111.);
  t_mass_mm_bb =         fs->make < TH1F > ("t_mass_mm_bb",       "t_mass_mm_bb;Mass [GeV]", 80, 71., 111.);

  h_mt_wenu_b =         fs->make < TH1F > ("h_mt_wenu_b",       "h_mt_wenu_b;M_{T} [GeV]", 20, 45., 205.);
  w_mt_wenu_b =         fs->make < TH1F > ("w_mt_wenu_b",       "w_mt_wenu_b;M_{T} [GeV]", 20, 45., 205.);
  b_mt_wenu_b =         fs->make < TH1F > ("b_mt_wenu_b",       "b_mt_wenu_b;M_{T} [GeV]", 20, 45., 205.);
  c_mt_wenu_b =         fs->make < TH1F > ("c_mt_wenu_b",       "c_mt_wenu_b;M_{T} [GeV]", 20, 45., 205.);
  t_mt_wenu_b =         fs->make < TH1F > ("t_mt_wenu_b",       "t_mt_wenu_b;M_{T} [GeV]", 20, 45., 205.);

  h_mt_wmnu_b =         fs->make < TH1F > ("h_mt_wmnu_b",       "h_mt_wmnu_b;M_{T} [GeV]", 20, 45., 205.);
  w_mt_wmnu_b =         fs->make < TH1F > ("w_mt_wmnu_b",       "w_mt_wmnu_b;M_{T} [GeV]", 20, 45., 205.);
  b_mt_wmnu_b =         fs->make < TH1F > ("b_mt_wmnu_b",       "b_mt_wmnu_b;M_{T} [GeV]", 20, 45., 205.);
  c_mt_wmnu_b =         fs->make < TH1F > ("c_mt_wmnu_b",       "c_mt_wmnu_b;M_{T} [GeV]", 20, 45., 205.);
  t_mt_wmnu_b =         fs->make < TH1F > ("t_mt_wmnu_b",       "t_mt_wmnu_b;M_{T} [GeV]", 20, 45., 205.);

  h_mt_wenu_bb =        fs->make < TH1F > ("h_mt_wenu_bb",      "h_mt_wenu_bb;M_{T} [GeV]", 20, 45., 205.);
  w_mt_wenu_bb =        fs->make < TH1F > ("w_mt_wenu_bb",      "w_mt_wenu_bb;M_{T} [GeV]", 20, 45., 205.);
  b_mt_wenu_bb =        fs->make < TH1F > ("b_mt_wenu_bb",      "b_mt_wenu_bb;M_{T} [GeV]", 20, 45., 205.);
  c_mt_wenu_bb =        fs->make < TH1F > ("c_mt_wenu_bb",      "c_mt_wenu_bb;M_{T} [GeV]", 20, 45., 205.);
  t_mt_wenu_bb =        fs->make < TH1F > ("t_mt_wenu_bb",      "t_mt_wenu_bb;M_{T} [GeV]", 20, 45., 205.);

  h_mt_wmnu_bb =        fs->make < TH1F > ("h_mt_wmnu_bb",      "h_mt_wmnu_bb;M_{T} [GeV]", 20, 45., 205.);
  w_mt_wmnu_bb =        fs->make < TH1F > ("w_mt_wmnu_bb",      "w_mt_wmnu_bb;M_{T} [GeV]", 20, 45., 205.);
  b_mt_wmnu_bb =        fs->make < TH1F > ("b_mt_wmnu_bb",      "b_mt_wmnu_bb;M_{T} [GeV]", 20, 45., 205.);
  c_mt_wmnu_bb =        fs->make < TH1F > ("c_mt_wmnu_bb",      "c_mt_wmnu_bb;M_{T} [GeV]", 20, 45., 205.);
  t_mt_wmnu_bb =        fs->make < TH1F > ("t_mt_wmnu_bb",      "t_mt_wmnu_bb;M_{T} [GeV]", 20, 45., 205.);

  w_mass_Zj_ee =        fs->make < TH1F > ("w_mass_Zj_ee",      "w_mass_Zj_ee", 60, 100., 330.);
  b_mass_Zj_ee =        fs->make < TH1F > ("b_mass_Zj_ee",      "b_mass_Zj_ee", 60, 100., 330.);
  c_mass_Zj_ee =        fs->make < TH1F > ("c_mass_Zj_ee",      "c_mass_Zj_ee", 60, 100., 330.);
  t_mass_Zj_ee =        fs->make < TH1F > ("t_mass_Zj_ee",      "t_mass_Zj_ee", 60, 100., 330.);

  w_mass_Zj_mm =        fs->make < TH1F > ("w_mass_Zj_mm",      "w_mass_Zj_mm", 60, 100., 330.);
  b_mass_Zj_mm =        fs->make < TH1F > ("b_mass_Zj_mm",      "b_mass_Zj_mm", 60, 100., 330.);
  c_mass_Zj_mm =        fs->make < TH1F > ("c_mass_Zj_mm",      "c_mass_Zj_mm", 60, 100., 330.);
  t_mass_Zj_mm =        fs->make < TH1F > ("t_mass_Zj_mm",      "t_mass_Zj_mm", 60, 100., 330.);

  w_mass_Zj_ee_b =      fs->make < TH1F > ("w_mass_Zj_ee_b",    "w_mass_Zj_ee_b", 60, 100., 330.);
  b_mass_Zj_ee_b =      fs->make < TH1F > ("b_mass_Zj_ee_b",    "b_mass_Zj_ee_b", 60, 100., 330.);
  c_mass_Zj_ee_b =      fs->make < TH1F > ("c_mass_Zj_ee_b",    "c_mass_Zj_ee_b", 60, 100., 330.);
  t_mass_Zj_ee_b =      fs->make < TH1F > ("t_mass_Zj_ee_b",    "t_mass_Zj_ee_b", 60, 100., 330.);

  w_mass_Zj_mm_b =      fs->make < TH1F > ("w_mass_Zj_mm_b",    "w_mass_Zj_mm_b", 60, 100., 330.);
  b_mass_Zj_mm_b =      fs->make < TH1F > ("b_mass_Zj_mm_b",    "b_mass_Zj_mm_b", 60, 100., 330.);
  c_mass_Zj_mm_b =      fs->make < TH1F > ("c_mass_Zj_mm_b",    "c_mass_Zj_mm_b", 60, 100., 330.);
  t_mass_Zj_mm_b =      fs->make < TH1F > ("t_mass_Zj_mm_b",    "t_mass_Zj_mm_b", 60, 100., 330.);

  w_mass_Zj_ee_bb =      fs->make < TH1F > ("w_mass_Zj_ee_bb",    "w_mass_Zj_ee_bb", 60, 100., 330.);
  b_mass_Zj_ee_bb =      fs->make < TH1F > ("b_mass_Zj_ee_bb",    "b_mass_Zj_ee_bb", 60, 100., 330.);
  c_mass_Zj_ee_bb =      fs->make < TH1F > ("c_mass_Zj_ee_bb",    "c_mass_Zj_ee_bb", 60, 100., 330.);
  t_mass_Zj_ee_bb =      fs->make < TH1F > ("t_mass_Zj_ee_bb",    "t_mass_Zj_ee_bb", 60, 100., 330.);

  w_mass_Zj_mm_bb =      fs->make < TH1F > ("w_mass_Zj_mm_bb",    "w_mass_Zj_mm_bb", 60, 100., 330.);
  b_mass_Zj_mm_bb =      fs->make < TH1F > ("b_mass_Zj_mm_bb",    "b_mass_Zj_mm_bb", 60, 100., 330.);
  c_mass_Zj_mm_bb =      fs->make < TH1F > ("c_mass_Zj_mm_bb",    "c_mass_Zj_mm_bb", 60, 100., 330.);
  t_mass_Zj_mm_bb =      fs->make < TH1F > ("t_mass_Zj_mm_bb",    "t_mass_Zj_mm_bb", 60, 100., 330.);

  w_pt_Z_ee_b =         fs->make < TH1F > ("w_pt_Z_ee_b",       "w_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_ee_b =         fs->make < TH1F > ("b_pt_Z_ee_b",       "b_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_ee_b =         fs->make < TH1F > ("c_pt_Z_ee_b",       "c_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_ee_b =         fs->make < TH1F > ("t_pt_Z_ee_b",       "t_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);

  w_pt_Z_mm_b =         fs->make < TH1F > ("w_pt_Z_mm_b",       "w_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_mm_b =         fs->make < TH1F > ("b_pt_Z_mm_b",       "b_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_mm_b =         fs->make < TH1F > ("c_pt_Z_mm_b",       "c_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_mm_b =         fs->make < TH1F > ("t_pt_Z_mm_b",       "t_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);

  w_pt_Z_ee_bb =         fs->make < TH1F > ("w_pt_Z_ee_bb",       "w_pt_Z_ee_bb;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_ee_bb =         fs->make < TH1F > ("b_pt_Z_ee_bb",       "b_pt_Z_ee_bb;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_ee_bb =         fs->make < TH1F > ("c_pt_Z_ee_bb",       "c_pt_Z_ee_bb;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_ee_bb =         fs->make < TH1F > ("t_pt_Z_ee_bb",       "t_pt_Z_ee_bb;P_t [GeV]", 40, 0., 400.);

  w_pt_Z_mm_bb =         fs->make < TH1F > ("w_pt_Z_mm_bb",       "w_pt_Z_mm_bb;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_mm_bb =         fs->make < TH1F > ("b_pt_Z_mm_bb",       "b_pt_Z_mm_bb;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_mm_bb =         fs->make < TH1F > ("c_pt_Z_mm_bb",       "c_pt_Z_mm_bb;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_mm_bb =         fs->make < TH1F > ("t_pt_Z_mm_bb",       "t_pt_Z_mm_bb;P_t [GeV]", 40, 0., 400.);

  w_delta_ee =          fs->make < TH1F > ("w_delta_ee",    "w_delta_ee",   20, 0., TMath::Pi ());
  b_delta_ee =          fs->make < TH1F > ("b_delta_ee",    "b_delta_ee",   20, 0., TMath::Pi ());
  c_delta_ee =          fs->make < TH1F > ("c_delta_ee",    "c_delta_ee",   20, 0., TMath::Pi ());
  t_delta_ee =          fs->make < TH1F > ("t_delta_ee",    "t_delta_ee",   20, 0., TMath::Pi ());
  w_delta_ee_b =        fs->make < TH1F > ("w_delta_ee_b",  "w_delta_ee_b", 20, 0., TMath::Pi ());
  b_delta_ee_b =        fs->make < TH1F > ("b_delta_ee_b",  "b_delta_ee_b", 20, 0., TMath::Pi ());
  c_delta_ee_b =        fs->make < TH1F > ("c_delta_ee_b",  "c_delta_ee_b", 20, 0., TMath::Pi ());
  t_delta_ee_b =        fs->make < TH1F > ("t_delta_ee_b",  "t_delta_ee_b", 20, 0., TMath::Pi ());
  w_delta_ee_bb =        fs->make < TH1F > ("w_delta_ee_bb",  "w_delta_ee_bb", 20, 0., TMath::Pi ());
  b_delta_ee_bb =        fs->make < TH1F > ("b_delta_ee_bb",  "b_delta_ee_bb", 20, 0., TMath::Pi ());
  c_delta_ee_bb =        fs->make < TH1F > ("c_delta_ee_bb",  "c_delta_ee_bb", 20, 0., TMath::Pi ());
  t_delta_ee_bb =        fs->make < TH1F > ("t_delta_ee_bb",  "t_delta_ee_bb", 20, 0., TMath::Pi ());

  w_delta_mm =          fs->make < TH1F > ("w_delta_mm",    "w_delta_mm",   20, 0., TMath::Pi ());
  b_delta_mm =          fs->make < TH1F > ("b_delta_mm",    "b_delta_mm",   20, 0., TMath::Pi ());
  c_delta_mm =          fs->make < TH1F > ("c_delta_mm",    "c_delta_mm",   20, 0., TMath::Pi ());
  t_delta_mm =          fs->make < TH1F > ("t_delta_mm",    "t_delta_mm",   20, 0., TMath::Pi ());
  w_delta_mm_b =        fs->make < TH1F > ("w_delta_mm_b",  "w_delta_mm_b", 20, 0., TMath::Pi ());
  b_delta_mm_b =        fs->make < TH1F > ("b_delta_mm_b",  "b_delta_mm_b", 20, 0., TMath::Pi ());
  c_delta_mm_b =        fs->make < TH1F > ("c_delta_mm_b",  "c_delta_mm_b", 20, 0., TMath::Pi ());
  t_delta_mm_b =        fs->make < TH1F > ("t_delta_mm_b",  "t_delta_mm_b", 20, 0., TMath::Pi ());
  w_delta_mm_bb =        fs->make < TH1F > ("w_delta_mm_bb",  "w_delta_mm_bb", 20, 0., TMath::Pi ());
  b_delta_mm_bb =        fs->make < TH1F > ("b_delta_mm_bb",  "b_delta_mm_bb", 20, 0., TMath::Pi ());
  c_delta_mm_bb =        fs->make < TH1F > ("c_delta_mm_bb",  "c_delta_mm_bb", 20, 0., TMath::Pi ());
  t_delta_mm_bb =        fs->make < TH1F > ("t_delta_mm_bb",  "t_delta_mm_bb", 20, 0., TMath::Pi ());

  w_delta_wenu =          fs->make < TH1F > ("w_delta_wenu",     "w_delta_wenu",    20, 0., TMath::Pi ());
  b_delta_wenu =          fs->make < TH1F > ("b_delta_wenu",     "b_delta_wenu",    20, 0., TMath::Pi ());
  c_delta_wenu =          fs->make < TH1F > ("c_delta_wenu",     "c_delta_wenu",    20, 0., TMath::Pi ());
  t_delta_wenu =          fs->make < TH1F > ("t_delta_wenu",     "t_delta_wenu",    20, 0., TMath::Pi ());
  w_delta_wenu_b =        fs->make < TH1F > ("w_delta_wenu_b",   "w_delta_wenu_b",  20, 0., TMath::Pi ());
  b_delta_wenu_b =        fs->make < TH1F > ("b_delta_wenu_b",   "b_delta_wenu_b",  20, 0., TMath::Pi ());
  c_delta_wenu_b =        fs->make < TH1F > ("c_delta_wenu_b",   "c_delta_wenu_b",  20, 0., TMath::Pi ());
  t_delta_wenu_b =        fs->make < TH1F > ("t_delta_wenu_b",   "t_delta_wenu_b",  20, 0., TMath::Pi ());
  w_delta_wenu_bb =       fs->make < TH1F > ("w_delta_wenu_bb",  "w_delta_wenu_bb", 20, 0., TMath::Pi ());
  b_delta_wenu_bb =       fs->make < TH1F > ("b_delta_wenu_bb",  "b_delta_wenu_bb", 20, 0., TMath::Pi ());
  c_delta_wenu_bb =       fs->make < TH1F > ("c_delta_wenu_bb",  "c_delta_wenu_bb", 20, 0., TMath::Pi ());
  t_delta_wenu_bb =       fs->make < TH1F > ("t_delta_wenu_bb",  "t_delta_wenu_bb", 20, 0., TMath::Pi ());
  w_delta_wenu_2b =       fs->make < TH1F > ("w_delta_wenu_2b",  "w_delta_wenu_2b", 20, 0., TMath::Pi ());
  b_delta_wenu_2b =       fs->make < TH1F > ("b_delta_wenu_2b",  "b_delta_wenu_2b", 20, 0., TMath::Pi ());
  c_delta_wenu_2b =       fs->make < TH1F > ("c_delta_wenu_2b",  "c_delta_wenu_2b", 20, 0., TMath::Pi ());
  t_delta_wenu_2b =       fs->make < TH1F > ("t_delta_wenu_2b",  "t_delta_wenu_2b", 20, 0., TMath::Pi ());

  w_delta_wmnu =          fs->make < TH1F > ("w_delta_wmnu",     "w_delta_wmnu",    20, 0., TMath::Pi ());
  b_delta_wmnu =          fs->make < TH1F > ("b_delta_wmnu",     "b_delta_wmnu",    20, 0., TMath::Pi ());
  c_delta_wmnu =          fs->make < TH1F > ("c_delta_wmnu",     "c_delta_wmnu",    20, 0., TMath::Pi ());
  t_delta_wmnu =          fs->make < TH1F > ("t_delta_wmnu",     "t_delta_wmnu",    20, 0., TMath::Pi ());
  w_delta_wmnu_b =        fs->make < TH1F > ("w_delta_wmnu_b",   "w_delta_wmnu_b",  20, 0., TMath::Pi ());
  b_delta_wmnu_b =        fs->make < TH1F > ("b_delta_wmnu_b",   "b_delta_wmnu_b",  20, 0., TMath::Pi ());
  c_delta_wmnu_b =        fs->make < TH1F > ("c_delta_wmnu_b",   "c_delta_wmnu_b",  20, 0., TMath::Pi ());
  t_delta_wmnu_b =        fs->make < TH1F > ("t_delta_wmnu_b",   "t_delta_wmnu_b",  20, 0., TMath::Pi ());
  w_delta_wmnu_bb =       fs->make < TH1F > ("w_delta_wmnu_bb",  "w_delta_wmnu_bb", 20, 0., TMath::Pi ());
  b_delta_wmnu_bb =       fs->make < TH1F > ("b_delta_wmnu_bb",  "b_delta_wmnu_bb", 20, 0., TMath::Pi ());
  c_delta_wmnu_bb =       fs->make < TH1F > ("c_delta_wmnu_bb",  "c_delta_wmnu_bb", 20, 0., TMath::Pi ());
  t_delta_wmnu_bb =       fs->make < TH1F > ("t_delta_wmnu_bb",  "t_delta_wmnu_bb", 20, 0., TMath::Pi ());
  w_delta_wmnu_2b =       fs->make < TH1F > ("w_delta_wmnu_2b",  "w_delta_wmnu_2b", 20, 0., TMath::Pi ());
  b_delta_wmnu_2b =       fs->make < TH1F > ("b_delta_wmnu_2b",  "b_delta_wmnu_2b", 20, 0., TMath::Pi ());
  c_delta_wmnu_2b =       fs->make < TH1F > ("c_delta_wmnu_2b",  "c_delta_wmnu_2b", 20, 0., TMath::Pi ());
  t_delta_wmnu_2b =       fs->make < TH1F > ("t_delta_wmnu_2b",  "t_delta_wmnu_2b", 20, 0., TMath::Pi ());

  w_deltaR_wenu =          fs->make < TH1F > ("w_deltaR_wenu",     "w_deltaR_wenu",    24, 0., 4.8);
  b_deltaR_wenu =          fs->make < TH1F > ("b_deltaR_wenu",     "b_deltaR_wenu",    24, 0., 4.8);
  c_deltaR_wenu =          fs->make < TH1F > ("c_deltaR_wenu",     "c_deltaR_wenu",    24, 0., 4.8);
  t_deltaR_wenu =          fs->make < TH1F > ("t_deltaR_wenu",     "t_deltaR_wenu",    24, 0., 4.8);
  w_deltaR_wenu_b =        fs->make < TH1F > ("w_deltaR_wenu_b",   "w_deltaR_wenu_b",  24, 0., 4.8);
  b_deltaR_wenu_b =        fs->make < TH1F > ("b_deltaR_wenu_b",   "b_deltaR_wenu_b",  24, 0., 4.8);
  c_deltaR_wenu_b =        fs->make < TH1F > ("c_deltaR_wenu_b",   "c_deltaR_wenu_b",  24, 0., 4.8);
  t_deltaR_wenu_b =        fs->make < TH1F > ("t_deltaR_wenu_b",   "t_deltaR_wenu_b",  24, 0., 4.8);
  w_deltaR_wenu_bb =       fs->make < TH1F > ("w_deltaR_wenu_bb",  "w_deltaR_wenu_bb", 24, 0., 4.8);
  b_deltaR_wenu_bb =       fs->make < TH1F > ("b_deltaR_wenu_bb",  "b_deltaR_wenu_bb", 24, 0., 4.8);
  c_deltaR_wenu_bb =       fs->make < TH1F > ("c_deltaR_wenu_bb",  "c_deltaR_wenu_bb", 24, 0., 4.8);
  t_deltaR_wenu_bb =       fs->make < TH1F > ("t_deltaR_wenu_bb",  "t_deltaR_wenu_bb", 24, 0., 4.8);
  w_deltaR_wenu_2b =       fs->make < TH1F > ("w_deltaR_wenu_2b",  "w_deltaR_wenu_2b", 24, 0., 4.8);
  b_deltaR_wenu_2b =       fs->make < TH1F > ("b_deltaR_wenu_2b",  "b_deltaR_wenu_2b", 24, 0., 4.8);
  c_deltaR_wenu_2b =       fs->make < TH1F > ("c_deltaR_wenu_2b",  "c_deltaR_wenu_2b", 24, 0., 4.8);
  t_deltaR_wenu_2b =       fs->make < TH1F > ("t_deltaR_wenu_2b",  "t_deltaR_wenu_2b", 24, 0., 4.8);

  w_deltaR_wmnu =          fs->make < TH1F > ("w_deltaR_wmnu",     "w_deltaR_wmnu",    24, 0., 4.8);
  b_deltaR_wmnu =          fs->make < TH1F > ("b_deltaR_wmnu",     "b_deltaR_wmnu",    24, 0., 4.8);
  c_deltaR_wmnu =          fs->make < TH1F > ("c_deltaR_wmnu",     "c_deltaR_wmnu",    24, 0., 4.8);
  t_deltaR_wmnu =          fs->make < TH1F > ("t_deltaR_wmnu",     "t_deltaR_wmnu",    24, 0., 4.8);
  w_deltaR_wmnu_b =        fs->make < TH1F > ("w_deltaR_wmnu_b",   "w_deltaR_wmnu_b",  24, 0., 4.8);
  b_deltaR_wmnu_b =        fs->make < TH1F > ("b_deltaR_wmnu_b",   "b_deltaR_wmnu_b",  24, 0., 4.8);
  c_deltaR_wmnu_b =        fs->make < TH1F > ("c_deltaR_wmnu_b",   "c_deltaR_wmnu_b",  24, 0., 4.8);
  t_deltaR_wmnu_b =        fs->make < TH1F > ("t_deltaR_wmnu_b",   "t_deltaR_wmnu_b",  24, 0., 4.8);
  w_deltaR_wmnu_bb =       fs->make < TH1F > ("w_deltaR_wmnu_bb",  "w_deltaR_wmnu_bb", 24, 0., 4.8);
  b_deltaR_wmnu_bb =       fs->make < TH1F > ("b_deltaR_wmnu_bb",  "b_deltaR_wmnu_bb", 24, 0., 4.8);
  c_deltaR_wmnu_bb =       fs->make < TH1F > ("c_deltaR_wmnu_bb",  "c_deltaR_wmnu_bb", 24, 0., 4.8);
  t_deltaR_wmnu_bb =       fs->make < TH1F > ("t_deltaR_wmnu_bb",  "t_deltaR_wmnu_bb", 24, 0., 4.8);
  w_deltaR_wmnu_2b =       fs->make < TH1F > ("w_deltaR_wmnu_2b",  "w_deltaR_wmnu_2b", 24, 0., 4.8);
  b_deltaR_wmnu_2b =       fs->make < TH1F > ("b_deltaR_wmnu_2b",  "b_deltaR_wmnu_2b", 24, 0., 4.8);
  c_deltaR_wmnu_2b =       fs->make < TH1F > ("c_deltaR_wmnu_2b",  "c_deltaR_wmnu_2b", 24, 0., 4.8);
  t_deltaR_wmnu_2b =       fs->make < TH1F > ("t_deltaR_wmnu_2b",  "t_deltaR_wmnu_2b", 24, 0., 4.8);

  w_single_delta_ee_b =      fs->make < TH1F > ("w_single_delta_ee_b", "w_single_delta_ee_b", 20, 0., TMath::Pi ());
  b_single_delta_ee_b =      fs->make < TH1F > ("b_single_delta_ee_b", "b_single_delta_ee_b", 20, 0., TMath::Pi ());
  c_single_delta_ee_b =      fs->make < TH1F > ("c_single_delta_ee_b", "c_single_delta_ee_b", 20, 0., TMath::Pi ());
  t_single_delta_ee_b =      fs->make < TH1F > ("t_single_delta_ee_b", "t_single_delta_ee_b", 20, 0., TMath::Pi ());

  w_single_delta_mm_b =      fs->make < TH1F > ("w_single_delta_mm_b", "w_single_delta_mm_b", 20, 0., TMath::Pi ());
  b_single_delta_mm_b =      fs->make < TH1F > ("b_single_delta_mm_b", "b_single_delta_mm_b", 20, 0., TMath::Pi ());
  c_single_delta_mm_b =      fs->make < TH1F > ("c_single_delta_mm_b", "c_single_delta_mm_b", 20, 0., TMath::Pi ());
  t_single_delta_mm_b =      fs->make < TH1F > ("t_single_delta_mm_b", "t_single_delta_mm_b", 20, 0., TMath::Pi ());

  w_single_delta_wenu_b =        fs->make < TH1F > ("w_single_delta_wenu_b",  "w_single_delta_wenu_b", 20, 0., TMath::Pi ());
  b_single_delta_wenu_b =        fs->make < TH1F > ("b_single_delta_wenu_b",  "b_single_delta_wenu_b", 20, 0., TMath::Pi ());
  c_single_delta_wenu_b =        fs->make < TH1F > ("c_single_delta_wenu_b",  "c_single_delta_wenu_b", 20, 0., TMath::Pi ());
  t_single_delta_wenu_b =        fs->make < TH1F > ("t_single_delta_wenu_b",  "t_single_delta_wenu_b", 20, 0., TMath::Pi ());

  w_single_delta_wmnu_b =        fs->make < TH1F > ("w_single_delta_wmnu_b",  "w_single_delta_wmnu_b", 20, 0., TMath::Pi ());
  b_single_delta_wmnu_b =        fs->make < TH1F > ("b_single_delta_wmnu_b",  "b_single_delta_wmnu_b", 20, 0., TMath::Pi ());
  c_single_delta_wmnu_b =        fs->make < TH1F > ("c_single_delta_wmnu_b",  "c_single_delta_wmnu_b", 20, 0., TMath::Pi ());
  t_single_delta_wmnu_b =        fs->make < TH1F > ("t_single_delta_wmnu_b",  "t_single_delta_wmnu_b", 20, 0., TMath::Pi ());

  w_single_deltaR_wenu_b =        fs->make < TH1F > ("w_single_deltaR_wenu_b",  "w_single_deltaR_wenu_b", 20, 0., TMath::Pi ());
  b_single_deltaR_wenu_b =        fs->make < TH1F > ("b_single_deltaR_wenu_b",  "b_single_deltaR_wenu_b", 20, 0., TMath::Pi ());
  c_single_deltaR_wenu_b =        fs->make < TH1F > ("c_single_deltaR_wenu_b",  "c_single_deltaR_wenu_b", 20, 0., TMath::Pi ());
  t_single_deltaR_wenu_b =        fs->make < TH1F > ("t_single_deltaR_wenu_b",  "t_single_deltaR_wenu_b", 20, 0., TMath::Pi ());

  w_single_deltaR_wmnu_b =        fs->make < TH1F > ("w_single_deltaR_wmnu_b",  "w_single_deltaR_wmnu_b", 20, 0., TMath::Pi ());
  b_single_deltaR_wmnu_b =        fs->make < TH1F > ("b_single_deltaR_wmnu_b",  "b_single_deltaR_wmnu_b", 20, 0., TMath::Pi ());
  c_single_deltaR_wmnu_b =        fs->make < TH1F > ("c_single_deltaR_wmnu_b",  "c_single_deltaR_wmnu_b", 20, 0., TMath::Pi ());
  t_single_deltaR_wmnu_b =        fs->make < TH1F > ("t_single_deltaR_wmnu_b",  "t_single_deltaR_wmnu_b", 20, 0., TMath::Pi ());

  h_secondvtx_N =         fs->make < TH1F > ("h_secondvtx_N",        "h_secondvtx_N", 50, 0, 1);
  w_secondvtx_N =         fs->make < TH1F > ("w_secondvtx_N",        "w_secondvtx_N", 50, 0, 1);
  w_secondvtx_N_zoom =    fs->make < TH1F > ("w_secondvtx_N_zoom",   "w_secondvtx_N_zoom", 20, 0.898, 1);
  w_secondvtx_N_mass =    fs->make < TH1F > ("w_secondvtx_N_mass",   "w_secondvtx_N_mass", 20, 0.898, 1);
  w_secondvtx_N_nomass =  fs->make < TH1F > ("w_secondvtx_N_nomass", "w_secondvtx_N_nomass", 20, 0.898, 1);

  b_secondvtx_N =         fs->make < TH1F > ("b_secondvtx_N",        "b_secondvtx_N", 50, 0, 1);
  b_secondvtx_N_zoom =    fs->make < TH1F > ("b_secondvtx_N_zoom",   "b_secondvtx_N_zoom", 20, 0.898, 1);
  b_secondvtx_N_mass =    fs->make < TH1F > ("b_secondvtx_N_mass",   "b_secondvtx_N_mass", 20, 0.898, 1);
  b_secondvtx_N_nomass =  fs->make < TH1F > ("b_secondvtx_N_nomass", "b_secondvtx_N_nomass", 20, 0.898, 1);

  c_secondvtx_N =         fs->make < TH1F > ("c_secondvtx_N",        "c_secondvtx_N", 50, 0, 1);
  c_secondvtx_N_zoom =    fs->make < TH1F > ("c_secondvtx_N_zoom",   "c_secondvtx_N_zoom", 20, 0.898, 1);
  c_secondvtx_N_mass =    fs->make < TH1F > ("c_secondvtx_N_mass",   "c_secondvtx_N_mass", 20, 0.898, 1);
  c_secondvtx_N_nomass =  fs->make < TH1F > ("c_secondvtx_N_nomass", "c_secondvtx_N_nomass", 20, 0.898, 1);

  t_secondvtx_N =         fs->make < TH1F > ("t_secondvtx_N",        "t_secondvtx_N", 50, 0, 1);
  t_secondvtx_N_zoom =    fs->make < TH1F > ("t_secondvtx_N_zoom",   "t_secondvtx_N_zoom", 20, 0.898, 1);
  t_secondvtx_N_mass =    fs->make < TH1F > ("t_secondvtx_N_mass",   "t_secondvtx_N_mass", 20, 0.898, 1);
  t_secondvtx_N_nomass =  fs->make < TH1F > ("t_secondvtx_N_nomass", "t_secondvtx_N_nomass", 20, 0.898, 1);

  w_SVTX_mass_jet =     fs->make < TH1F > ("w_SVTX_mass_jet",   "w_SVTX_mass_jet;Mass [GeV]", 50, 0, 6);
  b_SVTX_mass_jet =     fs->make < TH1F > ("b_SVTX_mass_jet",   "b_SVTX_mass_jet;Mass [GeV]", 50, 0, 6);
  c_SVTX_mass_jet =     fs->make < TH1F > ("c_SVTX_mass_jet",   "c_SVTX_mass_jet;Mass [GeV]", 50, 0, 6);
  t_SVTX_mass_jet =     fs->make < TH1F > ("t_SVTX_mass_jet",   "t_SVTX_mass_jet;Mass [GeV]", 50, 0, 6);

  w_SVTX_mass_trk =     fs->make < TH1F > ("w_SVTX_mass_trk",   "w_SVTX_mass_trk;Mass [GeV]", 50, 0, 50);
  b_SVTX_mass_trk =     fs->make < TH1F > ("b_SVTX_mass_trk",   "b_SVTX_mass_trk;Mass [GeV]", 50, 0, 50);
  c_SVTX_mass_trk =     fs->make < TH1F > ("c_SVTX_mass_trk",   "c_SVTX_mass_trk;Mass [GeV]", 50, 0, 50);
  t_SVTX_mass_trk =     fs->make < TH1F > ("t_SVTX_mass_trk",   "t_SVTX_mass_trk;Mass [GeV]", 50, 0, 50);

  w_SVTX_mass_jet_b =     fs->make < TH1F > ("w_SVTX_mass_jet_b",   "w_SVTX_mass_jet_b;Mass [GeV]", 50, 0, 6);
  b_SVTX_mass_jet_b =     fs->make < TH1F > ("b_SVTX_mass_jet_b",   "b_SVTX_mass_jet_b;Mass [GeV]", 50, 0, 6);
  c_SVTX_mass_jet_b =     fs->make < TH1F > ("c_SVTX_mass_jet_b",   "c_SVTX_mass_jet_b;Mass [GeV]", 50, 0, 6);
  t_SVTX_mass_jet_b =     fs->make < TH1F > ("t_SVTX_mass_jet_b",   "t_SVTX_mass_jet_b;Mass [GeV]", 50, 0, 6);

  w_SVTX_mass_trk_b =     fs->make < TH1F > ("w_SVTX_mass_trk_b",   "w_SVTX_mass_trk_b;Mass [GeV]", 50, 0, 50);
  b_SVTX_mass_trk_b =     fs->make < TH1F > ("b_SVTX_mass_trk_b",   "b_SVTX_mass_trk_b;Mass [GeV]", 50, 0, 50);
  c_SVTX_mass_trk_b =     fs->make < TH1F > ("c_SVTX_mass_trk_b",   "c_SVTX_mass_trk_b;Mass [GeV]", 50, 0, 50);
  t_SVTX_mass_trk_b =     fs->make < TH1F > ("t_SVTX_mass_trk_b",   "t_SVTX_mass_trk_b;Mass [GeV]", 50, 0, 50);

  w_SVTX_mass_jet_bb =     fs->make < TH1F > ("w_SVTX_mass_jet_bb",   "w_SVTX_mass_jet_bb;Mass [GeV]", 50, 0, 6);
  b_SVTX_mass_jet_bb =     fs->make < TH1F > ("b_SVTX_mass_jet_bb",   "b_SVTX_mass_jet_bb;Mass [GeV]", 50, 0, 6);
  c_SVTX_mass_jet_bb =     fs->make < TH1F > ("c_SVTX_mass_jet_bb",   "c_SVTX_mass_jet_bb;Mass [GeV]", 50, 0, 6);
  t_SVTX_mass_jet_bb =     fs->make < TH1F > ("t_SVTX_mass_jet_bb",   "t_SVTX_mass_jet_bb;Mass [GeV]", 50, 0, 6);

  w_SVTX_mass_trk_bb =     fs->make < TH1F > ("w_SVTX_mass_trk_bb",   "w_SVTX_mass_trk_bb;Mass [GeV]", 50, 0, 50);
  b_SVTX_mass_trk_bb =     fs->make < TH1F > ("b_SVTX_mass_trk_bb",   "b_SVTX_mass_trk_bb;Mass [GeV]", 50, 0, 50);
  c_SVTX_mass_trk_bb =     fs->make < TH1F > ("c_SVTX_mass_trk_bb",   "c_SVTX_mass_trk_bb;Mass [GeV]", 50, 0, 50);
  t_SVTX_mass_trk_bb =     fs->make < TH1F > ("t_SVTX_mass_trk_bb",   "t_SVTX_mass_trk_bb;Mass [GeV]", 50, 0, 50);

  w_SVTX_mass     =     fs->make < TH1F > ("w_SVTX_mass",       "w_SVTX_mass;Mass [GeV]", 50, 0, 6);
  b_SVTX_mass     =     fs->make < TH1F > ("b_SVTX_mass",       "b_SVTX_mass;Mass [GeV]", 50, 0, 6);
  c_SVTX_mass     =     fs->make < TH1F > ("c_SVTX_mass",       "c_SVTX_mass;Mass [GeV]", 50, 0, 6);
  t_SVTX_mass     =     fs->make < TH1F > ("t_SVTX_mass",       "t_SVTX_mass;Mass [GeV]", 50, 0, 6);

  w_SVTX_mass_b     =     fs->make < TH1F > ("w_SVTX_mass_b",       "w_SVTX_mass_b;Mass [GeV]", 50, 0, 6);
  b_SVTX_mass_b     =     fs->make < TH1F > ("b_SVTX_mass_b",       "b_SVTX_mass_b;Mass [GeV]", 50, 0, 6);
  c_SVTX_mass_b     =     fs->make < TH1F > ("c_SVTX_mass_b",       "c_SVTX_mass_b;Mass [GeV]", 50, 0, 6);
  t_SVTX_mass_b     =     fs->make < TH1F > ("t_SVTX_mass_b",       "t_SVTX_mass_b;Mass [GeV]", 50, 0, 6);

  w_SVTX_mass_bb     =     fs->make < TH1F > ("w_SVTX_mass_bb",       "w_SVTX_mass_bb;Mass [GeV]", 50, 0, 6);
  b_SVTX_mass_bb     =     fs->make < TH1F > ("b_SVTX_mass_bb",       "b_SVTX_mass_bb;Mass [GeV]", 50, 0, 6);
  c_SVTX_mass_bb     =     fs->make < TH1F > ("c_SVTX_mass_bb",       "c_SVTX_mass_bb;Mass [GeV]", 50, 0, 6);
  t_SVTX_mass_bb     =     fs->make < TH1F > ("t_SVTX_mass_bb",       "t_SVTX_mass_bb;Mass [GeV]", 50, 0, 6);

  w_BJP       =     fs->make < TH1F > ("w_BJP",   "w_BJP", 50, 0, 10);
  b_BJP       =     fs->make < TH1F > ("b_BJP",   "b_BJP", 50, 0, 10);
  c_BJP       =     fs->make < TH1F > ("c_BJP",   "c_BJP", 50, 0, 10);
  t_BJP       =     fs->make < TH1F > ("t_BJP",   "t_BJP", 50, 0, 10);

  w_BJP_b       =     fs->make < TH1F > ("w_BJP_b",   "w_BJP_b", 50, 0, 10);
  b_BJP_b       =     fs->make < TH1F > ("b_BJP_b",   "b_BJP_b", 50, 0, 10);
  c_BJP_b       =     fs->make < TH1F > ("c_BJP_b",   "c_BJP_b", 50, 0, 10);
  t_BJP_b       =     fs->make < TH1F > ("t_BJP_b",   "t_BJP_b", 50, 0, 10);

  w_BJP_bb       =     fs->make < TH1F > ("w_BJP_bb",   "w_BJP_bb", 50, 0, 10);
  b_BJP_bb       =     fs->make < TH1F > ("b_BJP_bb",   "b_BJP_bb", 50, 0, 10);
  c_BJP_bb       =     fs->make < TH1F > ("c_BJP_bb",   "c_BJP_bb", 50, 0, 10);
  t_BJP_bb       =     fs->make < TH1F > ("t_BJP_bb",   "t_BJP_bb", 50, 0, 10);

  w_JBP       =     fs->make < TH1F > ("w_JBP",   "w_JBP", 50, 0, 3);
  b_JBP       =     fs->make < TH1F > ("b_JBP",   "b_JBP", 50, 0, 3);
  c_JBP       =     fs->make < TH1F > ("c_JBP",   "c_JBP", 50, 0, 3);
  t_JBP       =     fs->make < TH1F > ("t_JBP",   "t_JBP", 50, 0, 3);

  w_JBP_b       =     fs->make < TH1F > ("w_JBP_b",   "w_JBP_b", 50, 0, 3);
  b_JBP_b       =     fs->make < TH1F > ("b_JBP_b",   "b_JBP_b", 50, 0, 3);
  c_JBP_b       =     fs->make < TH1F > ("c_JBP_b",   "c_JBP_b", 50, 0, 3);
  t_JBP_b       =     fs->make < TH1F > ("t_JBP_b",   "t_JBP_b", 50, 0, 3);

  w_JBP_bb       =     fs->make < TH1F > ("w_JBP_bb",   "w_JBP_bb", 50, 0, 3);
  b_JBP_bb       =     fs->make < TH1F > ("b_JBP_bb",   "b_JBP_bb", 50, 0, 3);
  c_JBP_bb       =     fs->make < TH1F > ("c_JBP_bb",   "c_JBP_bb", 50, 0, 3);
  t_JBP_bb       =     fs->make < TH1F > ("t_JBP_bb",   "t_JBP_bb", 50, 0, 3);

  w_BJP0       =     fs->make < TH1F > ("w_BJP0",   "w_BJP0", 50, 0, 10);
  b_BJP0       =     fs->make < TH1F > ("b_BJP0",   "b_BJP0", 50, 0, 10);
  c_BJP0       =     fs->make < TH1F > ("c_BJP0",   "c_BJP0", 50, 0, 10);
  t_BJP0       =     fs->make < TH1F > ("t_BJP0",   "t_BJP0", 50, 0, 10);

  w_BJP1       =     fs->make < TH1F > ("w_BJP1",   "w_BJP1", 50, 0, 10);
  b_BJP1       =     fs->make < TH1F > ("b_BJP1",   "b_BJP1", 50, 0, 10);
  c_BJP1       =     fs->make < TH1F > ("c_BJP1",   "c_BJP1", 50, 0, 10);
  t_BJP1       =     fs->make < TH1F > ("t_BJP1",   "t_BJP1", 50, 0, 10);

  w_BJP2       =     fs->make < TH1F > ("w_BJP2",   "w_BJP2", 50, 0, 10);
  b_BJP2       =     fs->make < TH1F > ("b_BJP2",   "b_BJP2", 50, 0, 10);
  c_BJP2       =     fs->make < TH1F > ("c_BJP2",   "c_BJP2", 50, 0, 10);
  t_BJP2       =     fs->make < TH1F > ("t_BJP2",   "t_BJP2", 50, 0, 10);

  w_BJP_mass  =     fs->make < TH1F > ("w_BJP_mass",   "w_BJP_mass", 50, 0, 10);
  b_BJP_mass  =     fs->make < TH1F > ("b_BJP_mass",   "b_BJP_mass", 50, 0, 10);
  c_BJP_mass  =     fs->make < TH1F > ("c_BJP_mass",   "c_BJP_mass", 50, 0, 10);
  t_BJP_mass  =     fs->make < TH1F > ("t_BJP_mass",   "t_BJP_mass", 50, 0, 10);

  w_JBP_mass  =     fs->make < TH1F > ("w_JBP_mass",   "w_JBP_mass", 50, 0, 3);
  b_JBP_mass  =     fs->make < TH1F > ("b_JBP_mass",   "b_JBP_mass", 50, 0, 3);
  c_JBP_mass  =     fs->make < TH1F > ("c_JBP_mass",   "c_JBP_mass", 50, 0, 3);
  t_JBP_mass  =     fs->make < TH1F > ("t_JBP_mass",   "t_JBP_mass", 50, 0, 3);

  w_BJP_nomass  =     fs->make < TH1F > ("w_BJP_nomass",   "w_BJP_nomass", 50, 0, 10);
  b_BJP_nomass  =     fs->make < TH1F > ("b_BJP_nomass",   "b_BJP_nomass", 50, 0, 10);
  c_BJP_nomass  =     fs->make < TH1F > ("c_BJP_nomass",   "c_BJP_nomass", 50, 0, 10);
  t_BJP_nomass  =     fs->make < TH1F > ("t_BJP_nomass",   "t_BJP_nomass", 50, 0, 10);

  w_JBP_nomass  =     fs->make < TH1F > ("w_JBP_nomass",   "w_JBP_nomass", 50, 0, 3);
  b_JBP_nomass  =     fs->make < TH1F > ("b_JBP_nomass",   "b_JBP_nomass", 50, 0, 3);
  c_JBP_nomass  =     fs->make < TH1F > ("c_JBP_nomass",   "c_JBP_nomass", 50, 0, 3);
  t_JBP_nomass  =     fs->make < TH1F > ("t_JBP_nomass",   "t_JBP_nomass", 50, 0, 3);

  w_BJP_mass_b  =     fs->make < TH1F > ("w_BJP_mass_b",   "w_BJP_mass_b", 50, 0, 10);
  b_BJP_mass_b  =     fs->make < TH1F > ("b_BJP_mass_b",   "b_BJP_mass_b", 50, 0, 10);
  c_BJP_mass_b  =     fs->make < TH1F > ("c_BJP_mass_b",   "c_BJP_mass_b", 50, 0, 10);
  t_BJP_mass_b  =     fs->make < TH1F > ("t_BJP_mass_b",   "t_BJP_mass_b", 50, 0, 10);

  w_JBP_mass_b  =     fs->make < TH1F > ("w_JBP_mass_b",   "w_JBP_mass_b", 50, 0, 3);
  b_JBP_mass_b  =     fs->make < TH1F > ("b_JBP_mass_b",   "b_JBP_mass_b", 50, 0, 3);
  c_JBP_mass_b  =     fs->make < TH1F > ("c_JBP_mass_b",   "c_JBP_mass_b", 50, 0, 3);
  t_JBP_mass_b  =     fs->make < TH1F > ("t_JBP_mass_b",   "t_JBP_mass_b", 50, 0, 3);

  w_BJP_nomass_b  =     fs->make < TH1F > ("w_BJP_nomass_b",   "w_BJP_nomass_b", 50, 0, 10);
  b_BJP_nomass_b  =     fs->make < TH1F > ("b_BJP_nomass_b",   "b_BJP_nomass_b", 50, 0, 10);
  c_BJP_nomass_b  =     fs->make < TH1F > ("c_BJP_nomass_b",   "c_BJP_nomass_b", 50, 0, 10);
  t_BJP_nomass_b  =     fs->make < TH1F > ("t_BJP_nomass_b",   "t_BJP_nomass_b", 50, 0, 10);

  w_JBP_nomass_b  =     fs->make < TH1F > ("w_JBP_nomass_b",   "w_JBP_nomass_b", 50, 0, 3);
  b_JBP_nomass_b  =     fs->make < TH1F > ("b_JBP_nomass_b",   "b_JBP_nomass_b", 50, 0, 3);
  c_JBP_nomass_b  =     fs->make < TH1F > ("c_JBP_nomass_b",   "c_JBP_nomass_b", 50, 0, 3);
  t_JBP_nomass_b  =     fs->make < TH1F > ("t_JBP_nomass_b",   "t_JBP_nomass_b", 50, 0, 3);

  w_BJP_mass_bb  =     fs->make < TH1F > ("w_BJP_mass_bb",   "w_BJP_mass_bb", 50, 0, 10);
  b_BJP_mass_bb  =     fs->make < TH1F > ("b_BJP_mass_bb",   "b_BJP_mass_bb", 50, 0, 10);
  c_BJP_mass_bb  =     fs->make < TH1F > ("c_BJP_mass_bb",   "c_BJP_mass_bb", 50, 0, 10);
  t_BJP_mass_bb  =     fs->make < TH1F > ("t_BJP_mass_bb",   "t_BJP_mass_bb", 50, 0, 10);

  w_JBP_mass_bb  =     fs->make < TH1F > ("w_JBP_mass_bb",   "w_JBP_mass_bb", 50, 0, 3);
  b_JBP_mass_bb  =     fs->make < TH1F > ("b_JBP_mass_bb",   "b_JBP_mass_bb", 50, 0, 3);
  c_JBP_mass_bb  =     fs->make < TH1F > ("c_JBP_mass_bb",   "c_JBP_mass_bb", 50, 0, 3);
  t_JBP_mass_bb  =     fs->make < TH1F > ("t_JBP_mass_bb",   "t_JBP_mass_bb", 50, 0, 3);

  w_BJP_nomass_bb  =     fs->make < TH1F > ("w_BJP_nomass_bb",   "w_BJP_nomass_bb", 50, 0, 10);
  b_BJP_nomass_bb  =     fs->make < TH1F > ("b_BJP_nomass_bb",   "b_BJP_nomass_bb", 50, 0, 10);
  c_BJP_nomass_bb  =     fs->make < TH1F > ("c_BJP_nomass_bb",   "c_BJP_nomass_bb", 50, 0, 10);
  t_BJP_nomass_bb  =     fs->make < TH1F > ("t_BJP_nomass_bb",   "t_BJP_nomass_bb", 50, 0, 10);

  w_JBP_nomass_bb  =     fs->make < TH1F > ("w_JBP_nomass_bb",   "w_JBP_nomass_bb", 50, 0, 3);
  b_JBP_nomass_bb  =     fs->make < TH1F > ("b_JBP_nomass_bb",   "b_JBP_nomass_bb", 50, 0, 3);
  c_JBP_nomass_bb  =     fs->make < TH1F > ("c_JBP_nomass_bb",   "c_JBP_nomass_bb", 50, 0, 3);
  t_JBP_nomass_bb  =     fs->make < TH1F > ("t_JBP_nomass_bb",   "t_JBP_nomass_bb", 50, 0, 3);

  w_Ht =                fs->make < TH1F > ("w_Ht",              "w_Ht [GeV]", 20, 20., 220.);
  b_Ht =                fs->make < TH1F > ("b_Ht",              "b_Ht [GeV]", 20, 20., 220.);
  c_Ht =                fs->make < TH1F > ("c_Ht",              "c_Ht [GeV]", 20, 20., 220.);
  t_Ht =                fs->make < TH1F > ("t_Ht",              "t_Ht [GeV]", 20, 20., 220.);

  w_Ht_b =              fs->make < TH1F > ("w_Ht_b",            "w_Ht_b [GeV]", 20, 20., 220.);
  b_Ht_b =              fs->make < TH1F > ("b_Ht_b",            "b_Ht_b [GeV]", 20, 20., 220.);
  c_Ht_b =              fs->make < TH1F > ("c_Ht_b",            "c_Ht_b [GeV]", 20, 20., 220.);
  t_Ht_b =              fs->make < TH1F > ("t_Ht_b",            "t_Ht_b [GeV]", 20, 20., 220.);

  w_Ht_bb =             fs->make < TH1F > ("w_Ht_bb",           "w_Ht_bb [GeV]", 20, 20., 220.);
  b_Ht_bb =             fs->make < TH1F > ("b_Ht_bb",           "b_Ht_bb [GeV]", 20, 20., 220.);
  c_Ht_bb =             fs->make < TH1F > ("c_Ht_bb",           "c_Ht_bb [GeV]", 20, 20., 220.);
  t_Ht_bb =             fs->make < TH1F > ("t_Ht_bb",           "t_Ht_bb [GeV]", 20, 20., 220.);

  w_single_Ht_b =       fs->make < TH1F > ("w_single_Ht_b",            "w_single_Ht_b [GeV]", 20, 20., 220.);
  b_single_Ht_b =       fs->make < TH1F > ("b_single_Ht_b",            "b_single_Ht_b [GeV]", 20, 20., 220.);
  c_single_Ht_b =       fs->make < TH1F > ("c_single_Ht_b",            "c_single_Ht_b [GeV]", 20, 20., 220.);
  t_single_Ht_b =       fs->make < TH1F > ("t_single_Ht_b",            "t_single_Ht_b [GeV]", 20, 20., 220.);

  h_MET =               fs->make < TH1F > ("h_MET",             "h_MET;MET [GeV]", 50, 0., 200.);
  w_MET =               fs->make < TH1F > ("w_MET",             "w_MET;MET [GeV]", 50, 0., 200.);
  b_MET =               fs->make < TH1F > ("b_MET",             "b_MET;MET [GeV]", 50, 0., 200.);
  c_MET =               fs->make < TH1F > ("c_MET",             "c_MET;MET [GeV]", 50, 0., 200.);
  t_MET =               fs->make < TH1F > ("t_MET",             "t_MET;MET [GeV]", 50, 0., 200.);
  h_MET_phi =               fs->make < TH1F > ("h_MET_phi",             "h_MET_phi;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  w_MET_phi =               fs->make < TH1F > ("w_MET_phi",             "w_MET_phi;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  b_MET_phi =               fs->make < TH1F > ("b_MET_phi",             "b_MET_phi;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  c_MET_phi =               fs->make < TH1F > ("c_MET_phi",             "c_MET_phi;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  t_MET_phi =               fs->make < TH1F > ("t_MET_phi",             "t_MET_phi;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  w_MET_sign = 	        fs->make < TH1F > ("w_MET_sign",        "w_MET_sign;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign = 	        fs->make < TH1F > ("b_MET_sign",        "b_MET_sign;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign = 	        fs->make < TH1F > ("c_MET_sign",        "c_MET_sign;MET significance [GeV]", 50, 0., 100.);
  t_MET_sign = 	        fs->make < TH1F > ("t_MET_sign",        "t_MET_sign;MET significance [GeV]", 50, 0., 100.);

  h_MET_b =             fs->make < TH1F > ("h_MET_b",         "h_MET_b;MET [GeV]", 50, 0., 200.);
  w_MET_b =             fs->make < TH1F > ("w_MET_b",         "w_MET_b;MET [GeV]", 50, 0., 200.);
  b_MET_b =             fs->make < TH1F > ("b_MET_b",         "b_MET_b;MET [GeV]", 50, 0., 200.);
  c_MET_b =             fs->make < TH1F > ("c_MET_b",         "c_MET_b;MET [GeV]", 50, 0., 200.);
  t_MET_b =             fs->make < TH1F > ("t_MET_b",         "t_MET_b;MET [GeV]", 50, 0., 200.);
  h_MET_phi_b =             fs->make < TH1F > ("h_MET_phi_b",         "h_MET_phi_b;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  w_MET_phi_b =             fs->make < TH1F > ("w_MET_phi_b",         "w_MET_phi_b;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  b_MET_phi_b =             fs->make < TH1F > ("b_MET_phi_b",         "b_MET_phi_b;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  c_MET_phi_b =             fs->make < TH1F > ("c_MET_phi_b",         "c_MET_phi_b;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  t_MET_phi_b =             fs->make < TH1F > ("t_MET_phi_b",         "t_MET_phi_b;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  w_MET_sign_b = 	fs->make < TH1F > ("w_MET_sign_b",    "w_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign_b = 	fs->make < TH1F > ("b_MET_sign_b",    "b_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign_b = 	fs->make < TH1F > ("c_MET_sign_b",    "c_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  t_MET_sign_b = 	fs->make < TH1F > ("t_MET_sign_b",    "t_MET_sign_b;MET significance [GeV]", 50, 0., 100.);

  h_MET_bb =            fs->make < TH1F > ("h_MET_bb",         "h_MET_bb;MET [GeV]", 50, 0., 200.);
  w_MET_bb =            fs->make < TH1F > ("w_MET_bb",         "w_MET_bb;MET [GeV]", 50, 0., 200.);
  b_MET_bb =            fs->make < TH1F > ("b_MET_bb",         "b_MET_bb;MET [GeV]", 50, 0., 200.);
  c_MET_bb =            fs->make < TH1F > ("c_MET_bb",         "c_MET_bb;MET [GeV]", 50, 0., 200.);
  t_MET_bb =            fs->make < TH1F > ("t_MET_bb",         "t_MET_bb;MET [GeV]", 50, 0., 200.);
  h_MET_phi_bb =            fs->make < TH1F > ("h_MET_phi_bb",         "h_MET_phi_bb;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  w_MET_phi_bb =            fs->make < TH1F > ("w_MET_phi_bb",         "w_MET_phi_bb;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  b_MET_phi_bb =            fs->make < TH1F > ("b_MET_phi_bb",         "b_MET_phi_bb;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  c_MET_phi_bb =            fs->make < TH1F > ("c_MET_phi_bb",         "c_MET_phi_bb;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  t_MET_phi_bb =            fs->make < TH1F > ("t_MET_phi_bb",         "t_MET_phi_bb;MET phi", 24, -TMath::Pi (), TMath::Pi ());
  w_MET_sign_bb = 	fs->make < TH1F > ("w_MET_sign_bb",    "w_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign_bb = 	fs->make < TH1F > ("b_MET_sign_bb",    "b_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign_bb = 	fs->make < TH1F > ("c_MET_sign_bb",    "c_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);
  t_MET_sign_bb = 	fs->make < TH1F > ("t_MET_sign_bb",    "t_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);

  w_Afb =               fs->make < TH1F > ("b_asymmetry",       "b_asymmetry", 10, -1, 1);

  h_scaleFactor_first_ele =   fs->make < TH1F > ("h_scaleFactor_first_ele",   "h_scaleFactor_first_ele", 50, 0.95, 1.05);
  b_scaleFactor_first_ele =   fs->make < TH1F > ("b_scaleFactor_first_ele",   "b_scaleFactor_first_ele", 50, 0.95, 1.05);
  h_scaleFactor_first_muon =  fs->make < TH1F > ("h_scaleFactor_first_muon",  "h_scaleFactor_first_muon", 50, 0.95, 1.05);
  b_scaleFactor_first_muon =  fs->make < TH1F > ("b_scaleFactor_first_muon",  "b_scaleFactor_first_muon", 50, 0.95, 1.05);
  h_scaleFactor_second_ele =  fs->make < TH1F > ("h_scaleFactor_second_ele",  "h_scaleFactor_second_ele", 50, 0.95, 1.05);
  b_scaleFactor_second_ele =  fs->make < TH1F > ("b_scaleFactor_second_ele",  "b_scaleFactor_second_ele", 50, 0.95, 1.05);
  h_scaleFactor_second_muon = fs->make < TH1F > ("h_scaleFactor_second_muon", "h_scaleFactor_second_muon", 50, 0.95, 1.05);
  b_scaleFactor_second_muon = fs->make < TH1F > ("b_scaleFactor_second_muon", "b_scaleFactor_second_muon", 50, 0.95, 1.05);

  h_JEC_uncert =        fs->make < TH1F > ("JEC uncert", "JEC uncert", 10, -0.5, 0.5);

  produces<std::vector<double>>("myEventWeight");

  produces<std::vector<math::XYZTLorentzVector>>("myElectrons");
  produces<std::vector<math::XYZTLorentzVector>>("myMuons");

  produces<std::vector<math::XYZTLorentzVector>>("myJets");

  produces<std::vector<double>>("myHt");

  produces<std::vector<double>>("myWenuPt");
  produces<std::vector<double>>("myWenuEta");
  produces<std::vector<double>>("myWmnuPt");
  produces<std::vector<double>>("myWmnuEta");

  produces<std::vector<double>>("myDijetPt");
  produces<std::vector<double>>("myDijetEta");
  produces<std::vector<double>>("myDijetMass");

  produces<std::vector<double>>("myBJetsWeights");

  produces<std::vector<math::XYZTLorentzVector>>("myBJets");

  produces<std::vector<double>>("myDeltaPhiEJ");
  produces<std::vector<double>>("myDeltaPhiEBJ");
  produces<std::vector<double>>("myDeltaPhiEBJBJ");
  produces<std::vector<double>>("myDeltaPhiMJ");
  produces<std::vector<double>>("myDeltaPhiMBJ");
  produces<std::vector<double>>("myDeltaPhiMBJBJ");

  produces<std::vector<double>>("myDeltaREJ");
  produces<std::vector<double>>("myDeltaREBJ");
  produces<std::vector<double>>("myDeltaREBJBJ");
  produces<std::vector<double>>("myDeltaRMJ");
  produces<std::vector<double>>("myDeltaRMBJ");
  produces<std::vector<double>>("myDeltaRMBJBJ");

}

WbAnalyzer::~WbAnalyzer () {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void WbAnalyzer::produce (edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace edm;
  using namespace std;

  bool debug = false;
  if (debug) cout << "Processing new event..." << endl;

  // Get electron collection
  edm::Handle < pat::ElectronCollection > electrons;
  iEvent.getByLabel ("matchedElectrons", electrons);

  if (lepton_ == "electronQCD") iEvent.getByLabel ("matchedElectronsQCD", electrons);

  // Get muon collection
  edm::Handle < pat::MuonCollection > muons;
  iEvent.getByLabel ("matchedMuons", muons);

  if (lepton_ == "muonQCD") iEvent.getByLabel ("matchedMuonsQCD", muons);

  // Get jet collection
  edm::Handle < vector < pat::Jet > > jets;
  iEvent.getByLabel ("goodJets", jets);

  // Get tracks
  edm::Handle < vector < reco::Track > > tracks;
  iEvent.getByLabel ("generalTracks", tracks);

  // Get METs
  edm::Handle < vector < reco::PFMET > > mets;
  iEvent.getByLabel (edm::InputTag ("pfType1CorrectedMet"), mets);

  if (debug && mets->empty()) cout << "Warning: empty MET collection." << endl;
  if (debug && mets->size()>1) cout << "Warning: MET collection size > 1." << endl;

  edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);

  edm::Handle<vector<reco::GenJet> > gJets;
  iEvent.getByLabel(edm::InputTag("goodJets","genJets"), gJets);

  std::auto_ptr<std::vector<double>> myEventWeight( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myElectrons( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myMuons( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myHt( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myWenuPt( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myWenuEta( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myWmnuPt( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myWmnuEta( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myDijetPt( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDijetEta( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDijetMass( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myBJetsWeights( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myDeltaPhiEJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaPhiEBJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaPhiEBJBJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaPhiMJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaPhiMBJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaPhiMBJBJ( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myDeltaREJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaREBJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaREBJBJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaRMJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaRMBJ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myDeltaRMBJBJ( new std::vector<double> );

  bool wenu_event = false;
  bool wmnu_event = false;
  bool ee_event = false;
  bool mm_event = false;

  int Nj = 0;
  int Nb = 0;

  double diele_mass = 0;
  double diele_phi = 0;
  double diele_pt = 0;

  double dimuon_mass = 0;
  double dimuon_phi = 0;
  double dimuon_pt = 0;

  double Ht = 0;

  double MyWeight = 1;

  double scalFac_first_e = 1;
  double scalFac_second_e = 1;
  double scalFac_first_m = 1;
  double scalFac_second_m = 1;
  double scalFac_b = 1;

  // ++++++ Pile-Up

  bool isMC = false;

  MyWeight = 1.0;

  Handle < vector < PileupSummaryInfo > > PupInfo;

  if (iEvent.getByLabel (edm::InputTag ("addPileupInfo"), PupInfo))  {

    isMC = true;

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

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;

  if (iEvent.getByLabel ("generator", genEventInfoHandle)) {

    double mcWeight = genEventInfoHandle->weight();

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

  // +++++++++ ELECTRONS

  vector < pat::Electron > vect_ele;
  vector < pat::Electron > vect_ele2;
  int ntrgMatchesEle = 0;

  for (pat::ElectronCollection::const_iterator ele = electrons->begin (); ele != electrons->end (); ++ele) {

    if ( lepton_ != "electronQCD" ) {
      if (ele->pt()>30. &&
	  ((fabs(ele->superCluster()->eta())<=1.442 &&
	    fabs(ele->deltaEtaSuperClusterTrackAtVtx())<0.004 &&
	    fabs(ele->deltaPhiSuperClusterTrackAtVtx())<0.03 &&
	    ele->sigmaIetaIeta()<0.01 &&
	    ele->hadronicOverEm()<0.12) ||
	   (fabs(ele->superCluster()->eta())>=1.566 && fabs(ele->eta())<2.1 &&
	    fabs(ele->deltaEtaSuperClusterTrackAtVtx())<0.005 &&
	    fabs(ele->deltaPhiSuperClusterTrackAtVtx())<0.02 &&
	    ele->sigmaIetaIeta()<0.03 &&
	    ele->hadronicOverEm()<0.10)) &&
	  fabs(ele->dB())<0.02 &&
	  fabs(1./ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy())<0.05 &&
	  (ele->chargedHadronIso() + fmax(ele->neutralHadronIso() + ele->photonIso() - 0.5*ele->puChargedHadronIso(),0))/ele->et() < 0.10 &&
	  ele->passConversionVeto() &&
	  ele->gsfTrack()->trackerExpectedHitsInner().numberOfHits() < 1 &&
	  ele->triggerObjectMatches().size()>0) {
	vect_ele.push_back (*ele);
      } else {
	vect_ele2.push_back (*ele);
      }
    } else {
      if (ele->pt()>30 &&
	  ((fabs(ele->superCluster()->eta())<=1.442 &&
	    fabs(ele->deltaEtaSuperClusterTrackAtVtx())<0.004 &&
	    fabs(ele->deltaPhiSuperClusterTrackAtVtx())<0.03 &&
	    ele->sigmaIetaIeta()<0.01 &&
	    ele->hadronicOverEm()<0.12) ||
	   (fabs(ele->superCluster()->eta())>=1.566 && fabs(ele->eta())<2.1 &&
	    fabs(ele->deltaEtaSuperClusterTrackAtVtx())<0.005 &&
	    fabs(ele->deltaPhiSuperClusterTrackAtVtx())<0.02 &&
	    ele->sigmaIetaIeta()<0.03 &&
	    ele->hadronicOverEm()<0.10)) &&
	  fabs(ele->dB())<0.02 &&
	  fabs(1./ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy())<0.05 &&
	  (ele->chargedHadronIso() + fmax(ele->neutralHadronIso() + ele->photonIso() - 0.5*ele->puChargedHadronIso(),0))/ele->et() >= 0.15 &&
	  ele->passConversionVeto() &&
	  ele->gsfTrack()->trackerExpectedHitsInner().numberOfHits() < 1) {
        vect_ele.push_back (*ele);
      } else {
	vect_ele2.push_back (*ele);
      }
    }
    if (ele->triggerObjectMatches().size()>0) ntrgMatchesEle++;

  }

  // Computing Mt:
  double op_met = mets->empty() ? 0. : (*mets)[0].et();
  double elept = vect_ele.empty() ? 0. : vect_ele[0].pt();
  double deltaPhiMetEle = 0.;
  if (op_met>0. && elept>0.) deltaPhiMetEle = fabs(vect_ele[0].phi() - (*mets)[0].phi());
  if (deltaPhiMetEle > acos(-1)) deltaPhiMetEle = 2*acos(-1) - deltaPhiMetEle;
  double mt_wenu = sqrt(2*elept*op_met*(1-TMath::Cos(deltaPhiMetEle)));
  bool mt_cut_wenu = mt_wenu > 45.;

  // Computing m_ee
  int iele0=0;
  int iele1=-1;

  for (unsigned int i=1; i<vect_ele.size(); ++i) {
    if (vect_ele[i].charge()*vect_ele[iele0].charge()<0 && iele1==-1) iele1=i;
  }

  math::XYZTLorentzVector z_ee;

  if (iele1!=-1) {
    z_ee = vect_ele[iele0].p4() + vect_ele[iele1].p4();
    diele_mass = z_ee.mass();
    diele_pt = z_ee.pt();
    diele_phi = z_ee.phi();
    if (diele_mass>71 && diele_mass<111) ee_event = true;
  }

  // computing W pt:

  math::XYZTLorentzVector w_enu;

  double wenu_pt = 0.;
  double wenu_eta = 0.;

  if (!vect_ele.empty() && !mets->empty()) {
    w_enu = vect_ele[0].p4() + (*mets)[0].p4();
    wenu_pt = w_enu.pt();
    wenu_eta = w_enu.eta();
  }

  // +++++++++ MUONS

  vector < pat::Muon > vect_muon;
  vector < pat::Muon > vect_muon2;
  int ntrgMatchesMuo=0;

  for (pat::MuonCollection::const_iterator muon = muons->begin (); muon != muons->end (); ++muon) {

    if ( lepton_ != "muonQCD" ) {
      if (muon->pt()>30 && fabs(muon->eta())<2.1 &&
	  muon->isGlobalMuon() && muon->isPFMuon() &&
	  muon->globalTrack()->normalizedChi2() < 10 &&
	  muon->track()->hitPattern().trackerLayersWithMeasurement() > 5 &&
	  muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
	  muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
	  fabs(muon->dB()) < 0.2 &&
	  muon->numberOfMatchedStations() > 1 &&
	  (muon->chargedHadronIso() + fmax(muon->neutralHadronIso() + muon->photonIso() - 0.5*muon->puChargedHadronIso(),0))/muon->pt() < 0.12 &&
	  muon->triggerObjectMatches().size()>0) {
	vect_muon.push_back (*muon);
      } else {
	vect_muon2.push_back (*muon);
      }
    } else {
      if (muon->pt()>30 && fabs(muon->eta())<2.1 &&
	  muon->isGlobalMuon() && muon->isPFMuon() &&
	  muon->globalTrack()->normalizedChi2() < 10 &&
	  muon->track()->hitPattern().trackerLayersWithMeasurement() > 5 &&
	  muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
	  muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
	  fabs(muon->dB()) < 0.2 &&
	  muon->numberOfMatchedStations() > 1 &&
	  (muon->chargedHadronIso() + fmax(muon->neutralHadronIso() + muon->photonIso() - 0.5*muon->puChargedHadronIso(),0))/muon->pt() >= 0.20) {
        vect_muon.push_back (*muon);
      } else {
	vect_muon2.push_back (*muon);
      }
    }
    if (muon->triggerObjectMatches().size()>0) ntrgMatchesMuo++;

  }

  // Computing Mt:
  double muopt = vect_muon.empty() ? 0. : vect_muon[0].pt();
  double deltaPhiMetMuo = 0.;
  if (op_met>0. && muopt>0.) deltaPhiMetMuo = fabs(vect_muon[0].phi() - (*mets)[0].phi());
  if (deltaPhiMetMuo > acos(-1)) deltaPhiMetMuo = 2*acos(-1) - deltaPhiMetMuo;
  double mt_wmnu = sqrt(2*muopt*op_met*(1-TMath::Cos(deltaPhiMetMuo)));
  bool mt_cut_wmnu = mt_wmnu > 45.;

  // Computing m_ee
  int imuon0=0;
  int imuon1=-1;

  for (unsigned int i=1; i<vect_muon.size(); ++i) {
    if (vect_muon[i].charge()*vect_muon[imuon0].charge()<0 && imuon1==-1) imuon1=i;
  }

  math::XYZTLorentzVector z_mm;

  if (imuon1!=-1) {
    z_mm = vect_muon[imuon0].p4() + vect_muon[imuon1].p4();
    dimuon_mass = z_mm.mass();
    dimuon_pt = z_mm.pt();
    dimuon_phi = z_mm.phi();
    if (dimuon_mass>71 && dimuon_mass<111) mm_event = true;
  }

  // computing W pt:

  math::XYZTLorentzVector w_mnu;

  double wmnu_pt = 0.;
  double wmnu_eta = 0.;

  if (!vect_muon.empty() && !mets->empty()) {
    w_mnu = vect_muon[0].p4() + (*mets)[0].p4();
    wmnu_pt = w_mnu.pt();
    wmnu_eta = w_mnu.eta();
  }


  // +++++++++ Decisions:

  if (debug) cout << "Decisions..." << endl;

  if (vect_ele.size()==0 && vect_muon.size()==0) {if (debug) cout << "No isolated leptons!!" << endl;}

  if (lepton_ == "electron" || lepton_ == "muon") {
    if (vect_ele.size()==1 && vect_ele2.size()==0 && vect_muon.size()==0 && vect_muon2.size()==0) wenu_event = true;
    if (vect_muon.size()==1 && vect_muon2.size()==0 && vect_ele.size()==0 && vect_ele2.size()==0) wmnu_event = true;
  }
  if (lepton_ == "electronQCD" || lepton_ == "muonQCD") {
    if (vect_ele.size()==1 && vect_ele2.size()==0 && vect_muon.size()==0 && vect_muon2.size()==0) wenu_event = true;
    if (vect_muon.size()==1 && vect_muon2.size()==0 && vect_ele.size()==0 && vect_ele2.size()==0) wmnu_event = true;
  }
  if (lepton_ == "electronFWD" || lepton_ == "muonFWD") {
    if (vect_ele.size()==1 && vect_ele2.size()==0 && vect_muon.size()==0 && vect_muon2.size()==0) wenu_event = true;
    if (vect_muon.size()==1 && vect_muon2.size()==0 && vect_ele.size()==0 && vect_ele2.size()==0) wmnu_event = true;
  }
  if (lepton_ == "electronTOP" || lepton_ == "muonTOP") {
    if (vect_ele.size()==1 && vect_ele2.size()==0 && vect_muon.size()==1 && vect_muon2.size()==0) wenu_event = true;
    if (vect_muon.size()==1 && vect_muon2.size()==0 && vect_ele.size()==1 && vect_ele2.size()==0) wmnu_event = true;
  }

  ee_event = ee_event && (lepton_ == "electron" || lepton_ == "electronQCD" || lepton_ == "electronFWD");
  mm_event = mm_event && (lepton_ == "muon" || lepton_ == "muonQCD" || lepton_ == "muonFWD");

  // +++++++++ SCALE FACTORS:

  if (isMC) {
    if (wenu_event) {
      scalFac_first_e = ElSF_->Val (vect_ele[0].pt(), vect_ele[0].eta()) * ElSF2_->Val (vect_ele[0].pt(), vect_ele[0].eta());
      MyWeight = MyWeight * scalFac_first_e;
    }
    if (wmnu_event) {
      scalFac_first_m = MuSF_->Val (vect_muon[0].pt(), vect_muon[0].eta()) * MuSF2_->Val (vect_muon[0].pt(), vect_muon[0].eta()) * MuSF3_->Val (vect_muon[0].pt(), vect_muon[0].eta());
      MyWeight = MyWeight * scalFac_first_m;
    }
  }

  if (isMC) {
    if (ee_event) {
      scalFac_first_e  = ElSF_->Val (vect_ele[iele0].pt(), vect_ele[iele0].eta()) * ElSF2_->Val (vect_ele[iele0].pt(), vect_ele[iele0].eta());
      scalFac_second_e = ElSF_->Val (vect_ele[iele1].pt(), vect_ele[iele1].eta()) * ElSF2_->Val (vect_ele[iele1].pt(), vect_ele[iele1].eta());
      MyWeight = MyWeight * scalFac_first_e * scalFac_second_e;
    }
    if (mm_event) {
      scalFac_first_m  = MuSF_->Val (vect_muon[imuon0].pt(), vect_muon[imuon0].eta()) * MuSF2_->Val (vect_muon[imuon0].pt(), vect_muon[imuon0].eta()) * MuSF3_->Val (vect_muon[imuon0].pt(), vect_muon[imuon0].eta());
      scalFac_second_m = MuSF_->Val (vect_muon[imuon1].pt(), vect_muon[imuon1].eta()) * MuSF2_->Val (vect_muon[imuon1].pt(), vect_muon[imuon1].eta()) * MuSF3_->Val (vect_muon[imuon1].pt(), vect_muon[imuon1].eta());
      MyWeight = MyWeight * scalFac_first_m * scalFac_second_m;
    }
  }

  // ++++++++ VERTICES

  bool vtx_cut = true;

  edm::Handle < vector < reco::Vertex > > vertices;
  iEvent.getByLabel (edm::InputTag ("goodOfflinePrimaryVertices"), vertices);

  if (vertices->size() > 0) {
    const reco::Vertex* theVertex = &(vertices->front());
    if (theVertex->ndof() < 5) vtx_cut = false;
    if (fabs(theVertex->z()) > 24.0) vtx_cut = false;
    if (fabs(theVertex->position().rho()) > 2.0) vtx_cut = false;
  } else {
    vtx_cut = false;
  }

  int NVtx = 0;

  if (vtx_cut) {
    for (vector < reco::Vertex >::const_iterator itv = vertices->begin (); itv != vertices->end (); ++itv) {
      if (itv->ndof() < 5) continue;
      if (fabs(itv->z()) > 50.0) continue;
      if (fabs(itv->position().rho()) > 2.0) continue;
      ++NVtx;
    }
  }

  // ++++++++ LOOP OVER GEN PARTICLES

  //  bool isb = false;
  //  bool isc = false;
  bool ist = false;

  if (isMC) {
    for (std::vector <reco::GenParticle>::const_iterator thepart = genPart->begin(); thepart != genPart->end(); thepart++) {
//      if ((int) (abs(thepart->pdgId() / 100) % 10) == 5 || (int) (abs(thepart->pdgId() / 1000) % 10) == 5) {
//        isb = true;
//      }
//      if ((int) (abs(thepart->pdgId() / 100) % 10 ) == 4 || (int) (abs(thepart->pdgId() / 1000) % 10) == 4) {
//        isc = true;
//      }
      if ((int) abs(thepart->pdgId()) == 24) {
	for (UInt_t i=0; i<thepart->numberOfDaughters(); i++){
	  if (abs(thepart->daughter(i)->pdgId()) == 15 && thepart->daughter(i)->status()==3){
	    ist = true;
  	  }
	}
      }
    }
  }

  // ++++++++ MET CUT

   bool met_cut = mets->empty() ? true : (*mets)[0].significance() < 30.;

  // ++++++++ JETS
  if (debug) cout << "Jets..." << endl;

  vector < pat::Jet > vect_jets;
  vector < pat::Jet > vect_jets2;
  vector < pat::Jet > vect_bjets;

  for (vector < pat::Jet >::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {

// check for no neutrinos
//    if (isMC && jet->genJet()) {
//      vector <const reco::GenParticle*> listGenP = jet->genJet()->getGenConstituents();
//      for (unsigned int i=0; i<listGenP.size(); i++) {
//        cout << i << " " << listGenP.size() << " : " << listGenP[i]->pdgId();
//        if (fabs(listGenP[i]->pdgId())==12 || fabs(listGenP[i]->pdgId())==14 || fabs(listGenP[i]->pdgId())==16) cout << " +++ Found neutrino in GenJet +++ ";
//        cout << endl;
//      }
//    }

    // JEC uncertainty

    double jecUnc = 0.0;
    if (!isMC) {
      jetCorrectionUncertainty_->setJetPt(jet->pt());
      jetCorrectionUncertainty_->setJetEta(jet->eta());
      jecUnc = jetCorrectionUncertainty_->getUncertainty(true);
    }
    h_JEC_uncert->Fill (jecUnc);
    //cout<< "JEC syst =" << unc << endl;

    // JER corrections

    double jerCor = 1.0;
    if (isMC && jet->genJet()) jerCor = jetResolutionCorrection(jet->eta(), jet->pt(), jet->genJet()->pt(), par2_);

    pat::Jet jetNew = (*jet);
    math::XYZTLorentzVector jetNew_p4 = jetNew.p4();

    jetNew_p4 = jetNew_p4 * (1.0 + jecUnc * par_) * jerCor;

    jetNew.setP4(jetNew_p4);

    // Jet-Lepton DR

    double etaj = jetNew.eta();
    double phij = jetNew.phi();

    double delta_eta1 = -999.;
    double delta_phi1 = -999.;
    double delta_eta2 = -999.;
    double delta_phi2 = -999.;

    if (wenu_event) {
      delta_eta1 = vect_ele[0].eta() - etaj;
      delta_phi1 = fabs(vect_ele[0].phi() - phij);
    }
    if (wmnu_event) {
      delta_eta1 = vect_muon[0].eta() - etaj;
      delta_phi1 = fabs(vect_muon[0].phi() - phij);
    }

    if (wenu_event && wmnu_event) {
      delta_eta1 = vect_ele[0].eta() - etaj;
      delta_phi1 = fabs(vect_ele[0].phi() - phij);
      delta_eta2 = vect_muon[0].eta() - etaj;
      delta_phi2 = fabs(vect_muon[0].phi() - phij);
    }
   
    if (ee_event) {
      delta_eta1 = vect_ele[0].eta() - etaj;
      delta_phi1 = fabs(vect_ele[0].phi() - phij);
      delta_eta2 = vect_ele[1].eta() - etaj;
      delta_phi2 = fabs(vect_ele[1].phi() - phij);
    }
    if (mm_event) {
      delta_eta1 = vect_muon[0].eta() - etaj;
      delta_phi1 = fabs(vect_muon[0].phi() - phij);
      delta_eta2 = vect_muon[1].eta() - etaj;
      delta_phi2 = fabs(vect_muon[1].phi() - phij);
    }

    if (delta_phi1 > acos(-1)) delta_phi1 = 2*acos(-1) - delta_phi1;
    if (delta_phi2 > acos(-1)) delta_phi2 = 2*acos(-1) - delta_phi2;
    
    double deltaR_jl1 = sqrt(pow(delta_eta1,2) + pow(delta_phi1,2));
    double deltaR_jl2 = sqrt(pow(delta_eta2,2) + pow(delta_phi2,2));    

    if (fabs(jetNew.eta()) < 2.4 && jetNew.pt() > 25 && deltaR_jl1 > 0.5 && deltaR_jl2 > 0.5) {

      ++Nj;

      Ht += jetNew.pt();

      vect_jets.push_back (jetNew);

      double discrCSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
      //cout << discrCSV << endl;

      if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
        h_secondvtx_N->Fill (discrCSV);
        w_secondvtx_N->Fill (discrCSV, MyWeight);
	if (ist) {
	  t_secondvtx_N->Fill (discrCSV, MyWeight);
	}
	if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
	  b_secondvtx_N->Fill (discrCSV, MyWeight);
	}
	if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
	  c_secondvtx_N->Fill (discrCSV, MyWeight);
	}
      }

      if (discrCSV > 0.898 || (usePartonFlavour_ && isMC && fabs(vect_jets[0].partonFlavour()) == 5)) {

	++Nb;
	//cout << Nb << endl;
        vect_bjets.push_back (jetNew);
      }
    }

    // Fill vector with forward jets (most likely from single-top):
    if (fabs(jetNew.eta()) >= 2.4 && fabs(jetNew.eta()) < 5.0 && jetNew.pt() > 25) {

      vect_jets2.push_back (jetNew);

    }
  }

  bool iflag_ee=false;
  bool iflag_mm=false;
  if (icut_==0 || Nb==0) {
    iflag_ee=true;
    iflag_mm=true;
  }
  if (icut_==1 && Nb>0 && diele_pt>0   && diele_pt<30)   iflag_ee=true;
  if (icut_==2 && Nb>0 && diele_pt>30  && diele_pt<50)   iflag_ee=true;
  if (icut_==3 && Nb>0 && diele_pt>50  && diele_pt<80)   iflag_ee=true;
  if (icut_==4 && Nb>0 && diele_pt>80  && diele_pt<120)  iflag_ee=true;
  if (icut_==5 && Nb>0 && diele_pt>120 && diele_pt<400)  iflag_ee=true;

  if (icut_==1 && Nb>0 && dimuon_pt>0   && dimuon_pt<30)   iflag_mm=true;
  if (icut_==2 && Nb>0 && dimuon_pt>30  && dimuon_pt<50)   iflag_mm=true;
  if (icut_==3 && Nb>0 && dimuon_pt>50  && dimuon_pt<80)   iflag_mm=true;
  if (icut_==4 && Nb>0 && dimuon_pt>80  && dimuon_pt<120)  iflag_mm=true;
  if (icut_==5 && Nb>0 && dimuon_pt>120 && dimuon_pt<400)  iflag_mm=true;

  if (icut_==6 && Nb>0 && vect_bjets[0].eta()> -2.5   && vect_bjets[0].eta()< -1.5)  iflag_ee=true;
  if (icut_==7 && Nb>0 && vect_bjets[0].eta()> -1.5   && vect_bjets[0].eta()< -1.0)  iflag_ee=true;
  if (icut_==8 && Nb>0 && vect_bjets[0].eta()> -1.0   && vect_bjets[0].eta()< -0.5)  iflag_ee=true;
  if (icut_==9 && Nb>0 && vect_bjets[0].eta()> -0.5   && vect_bjets[0].eta()<  0.0)  iflag_ee=true;
  if (icut_==10 && Nb>0 && vect_bjets[0].eta()>  0.0   && vect_bjets[0].eta()<  0.5) iflag_ee=true;
  if (icut_==11 && Nb>0 && vect_bjets[0].eta()>  0.5   && vect_bjets[0].eta()<  1.0) iflag_ee=true;
  if (icut_==12 && Nb>0 && vect_bjets[0].eta()>  1.0   && vect_bjets[0].eta()<  1.5) iflag_ee=true;
  if (icut_==13 && Nb>0 && vect_bjets[0].eta()>  1.5   && vect_bjets[0].eta()<  2.5) iflag_ee=true;

  if (icut_==6  && Nb>0 && vect_bjets[0].eta()> -2.5   && vect_bjets[0].eta()< -1.5)  iflag_mm=true;
  if (icut_==7  && Nb>0 && vect_bjets[0].eta()> -1.5   && vect_bjets[0].eta()< -1.0)  iflag_mm=true;
  if (icut_==8  && Nb>0 && vect_bjets[0].eta()> -1.0   && vect_bjets[0].eta()< -0.5)  iflag_mm=true;
  if (icut_==9  && Nb>0 && vect_bjets[0].eta()> -0.5   && vect_bjets[0].eta()<  0.0)  iflag_mm=true;
  if (icut_==10 && Nb>0 && vect_bjets[0].eta()>  0.0   && vect_bjets[0].eta()<  0.5)  iflag_mm=true;
  if (icut_==11 && Nb>0 && vect_bjets[0].eta()>  0.5   && vect_bjets[0].eta()<  1.0)  iflag_mm=true;
  if (icut_==12 && Nb>0 && vect_bjets[0].eta()>  1.0   && vect_bjets[0].eta()<  1.5)  iflag_mm=true;
  if (icut_==13 && Nb>0 && vect_bjets[0].eta()>  1.5   && vect_bjets[0].eta()<  2.5)  iflag_mm=true;

  if (icut_==14 && Nb>0 && vect_bjets[0].pt()> 30    && vect_bjets[0].pt()< 35)  iflag_ee=true;
  if (icut_==15 && Nb>0 && vect_bjets[0].pt()> 35    && vect_bjets[0].pt()< 40)  iflag_ee=true;
  if (icut_==16 && Nb>0 && vect_bjets[0].pt()> 40    && vect_bjets[0].pt()< 45)  iflag_ee=true;
  if (icut_==17 && Nb>0 && vect_bjets[0].pt()> 45    && vect_bjets[0].pt()< 50)  iflag_ee=true;
  if (icut_==18 && Nb>0 && vect_bjets[0].pt()> 50    && vect_bjets[0].pt()< 60)  iflag_ee=true;
  if (icut_==19 && Nb>0 && vect_bjets[0].pt()> 60    && vect_bjets[0].pt()< 80)  iflag_ee=true;
  if (icut_==20 && Nb>0 && vect_bjets[0].pt()> 80    && vect_bjets[0].pt()< 350) iflag_ee=true;

  if (icut_==14 && Nb>0 && vect_bjets[0].pt()> 30    && vect_bjets[0].pt()< 35)  iflag_mm=true;
  if (icut_==15 && Nb>0 && vect_bjets[0].pt()> 35    && vect_bjets[0].pt()< 40)  iflag_mm=true;
  if (icut_==16 && Nb>0 && vect_bjets[0].pt()> 40    && vect_bjets[0].pt()< 45)  iflag_mm=true;
  if (icut_==17 && Nb>0 && vect_bjets[0].pt()> 45    && vect_bjets[0].pt()< 50)  iflag_mm=true;
  if (icut_==18 && Nb>0 && vect_bjets[0].pt()> 50    && vect_bjets[0].pt()< 60)  iflag_mm=true;
  if (icut_==19 && Nb>0 && vect_bjets[0].pt()> 60    && vect_bjets[0].pt()< 80)  iflag_mm=true;
  if (icut_==20 && Nb>0 && vect_bjets[0].pt()> 80    && vect_bjets[0].pt()< 350) iflag_mm=true;

  double R = 0.5;
  double DR05_ej = 0;
  double DR05_mj = 0;
  double DEta_ej = 0;
  double DPhi_ej = 0;
  double DEta_mj = 0;
  double DPhi_mj = 0;

  if (useDeltaR_) {
    for (unsigned int i=0; i<vect_jets.size(); ++i) {
      for (unsigned int j=0; j<vect_ele.size(); ++j) {
        DEta_ej = fabs(vect_ele[j].eta() - vect_jets[i].eta());
        DPhi_ej = fabs(vect_ele[j].phi() - vect_jets[i].phi());
        if (DPhi_ej > TMath::ACos(-1)) DPhi_ej= 2*TMath::ACos(-1) - DPhi_ej;
        DR05_ej   = sqrt(DEta_ej*DEta_ej + DPhi_ej*DPhi_ej);
        if (DR05_ej < R) iflag_ee = false;
      }
      for (unsigned int k=0; k<vect_muon.size(); ++k) {
        DEta_mj = fabs(vect_muon[k].eta() - vect_jets[i].eta());
        DPhi_mj = fabs(vect_muon[k].phi() - vect_jets[i].phi());
        if (DPhi_mj > TMath::ACos(-1)) DPhi_mj= 2*TMath::ACos(-1) - DPhi_mj;
        DR05_mj   = sqrt(DEta_mj*DEta_mj + DPhi_mj*DPhi_mj);
        if (DR05_mj < R) iflag_mm = false;
      }
    }
  }

  if (Nb > 0 && pcut_) {
    double discrBJP = vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags");
    if (discrBJP <= 5.) {
      Nb = 0;
      vect_bjets.clear();
    }
  }

  ee_event = ee_event && Nb>0 && Nj>1 && iflag_ee;
  mm_event = mm_event && Nb>0 && Nj>1 && iflag_mm;

  if (lepton_=="electronFWD" && lepton_=="muonFWD") {
    wenu_event = wenu_event && Nb>0 && Nj==1;
    wmnu_event = wmnu_event && Nb>0 && Nj==1;
  } else {
    wenu_event = wenu_event && Nb>0 && Nj==2;
    wmnu_event = wmnu_event && Nb>0 && Nj==2;
  }

  wenu_event = wenu_event && ((vect_jets2.size()==0 && (lepton_=="electron" || lepton_=="electronQCD" || lepton_=="electronTOP"))
  			      || (vect_jets2.size()==1 && Nb==1 && lepton_=="electronFWD"));
  wmnu_event = wmnu_event && ((vect_jets2.size()==0 && (lepton_=="muon" || lepton_=="muonQCD" || lepton_=="muonTOP"))
  			      || (vect_jets2.size()==1 && Nb==1 && lepton_=="muonFWD"));

  if (debug && Nj<1) cout << "Warning: 0 Jets in the event!" << endl;
  if (debug && Nb<1) cout << "Warning: 0 b-Jets in the event!" << endl;

  // ++++++++ EVENT YIELDS:
  if (debug) cout << "Start filling plots..." << endl;

  h_eventYields->Fill(1);
  if (wenu_event || wmnu_event) h_eventYields->Fill(2);
  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) h_eventYields->Fill(3);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) h_eventYields->Fill(4);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) h_eventYields->Fill(5);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb>1) h_eventYields->Fill(6);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb==1) h_eventYields->Fill(7);

  scalFac_b = btagSF(isMC, vect_bjets, 1);
  w_eventYields->Fill(1, MyWeight*scalFac_b);
  if (wenu_event || wmnu_event) w_eventYields->Fill(2, MyWeight*scalFac_b);
  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) w_eventYields->Fill(3, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) w_eventYields->Fill(4, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) w_eventYields->Fill(5, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb>1) w_eventYields->Fill(6, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb==1) w_eventYields->Fill(7, MyWeight*scalFac_b);

  if (wenu_event && mt_cut_wenu && vtx_cut) h_trgMatchEle->Fill(ntrgMatchesEle);
  if (wmnu_event && mt_cut_wmnu && vtx_cut) h_trgMatchMuo->Fill(ntrgMatchesMuo);

  // ++++++++ MET PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_MET->Fill (mets->empty() ? 0 : (*mets)[0].et());
    w_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
    h_MET_phi->Fill (mets->empty() ? 0 : (*mets)[0].phi());
    w_MET_phi->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
    w_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    if (ist) {
      t_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      t_MET_phi->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
      t_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      b_MET_phi->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
      b_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      c_MET_phi->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
      c_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      h_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et());
      w_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      h_MET_phi_b->Fill (mets->empty() ? 0 : (*mets)[0].phi());
      w_MET_phi_b->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
      w_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      if (ist) {
	t_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
	t_MET_phi_b->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
	t_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        b_MET_phi_b->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
        b_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        c_MET_phi_b->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
        c_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      h_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et());
      w_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      h_MET_phi_bb->Fill (mets->empty() ? 0 : (*mets)[0].phi());
      w_MET_phi_bb->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
      w_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      if (ist) {
	t_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
	t_MET_phi_bb->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
	t_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        b_MET_phi_bb->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
        b_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        c_MET_phi_bb->Fill (mets->empty() ? 0 : (*mets)[0].phi(), MyWeight*scalFac_b);
        c_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ HT PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_Ht->Fill (Ht, MyWeight*scalFac_b);
    if (ist) {
      t_Ht->Fill (Ht, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_Ht->Fill (Ht, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_Ht->Fill (Ht, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (Nj == 1) w_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (ist) {
        t_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) t_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) b_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) c_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      w_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      if (ist) {
        t_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ WeNu PLOTS

  if (wenu_event && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_mt_wenu_wide->Fill (mt_wenu);
    w_mt_wenu_wide->Fill (mt_wenu, MyWeight*scalFac_b);
    if (ist) {
      t_mt_wenu_wide->Fill (mt_wenu, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mt_wenu_wide->Fill (mt_wenu, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mt_wenu_wide->Fill (mt_wenu, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      h_mt_wenu_b_wide->Fill (mt_wenu);
      w_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      h_mt_wenu_bb_wide->Fill (mt_wenu);
      w_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
    }
  }

  double delta_phi_ej=0;
  double delta_eta_ej=0;
  double DR_ej=0;
  double delta_phi_ebj=0;
  double delta_eta_ebj=0;
  double DR_ebj=0;
  double delta_phi_ebjbj=0;
  double delta_eta_ebjbj=0;
  double DR_ebjbj=0;

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_mt_wenu->Fill (mt_wenu);
    w_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
    w_pt_W_wenu->Fill (wenu_pt, MyWeight*scalFac_b);
    w_eta_W_wenu->Fill (wenu_eta, MyWeight*scalFac_b);
    delta_phi_ej = fabs(vect_ele[0].phi() - vect_jets[0].phi());
    delta_eta_ej = fabs(vect_ele[0].eta() - vect_jets[0].eta());
    if (delta_phi_ej > acos (-1)) delta_phi_ej = 2 * acos (-1) - delta_phi_ej;
    DR_ej = TMath::Sqrt(delta_phi_ej*delta_phi_ej + delta_eta_ej*delta_eta_ej);
    delta_phi_ebj = fabs(vect_ele[0].phi() - vect_bjets[0].phi());
    delta_eta_ebj = fabs(vect_ele[0].eta() - vect_bjets[0].eta());
    if (delta_phi_ebj > acos (-1)) delta_phi_ebj = 2 * acos (-1) - delta_phi_ebj;
    DR_ebj = TMath::Sqrt(delta_phi_ebj*delta_phi_ebj + delta_eta_ebj*delta_eta_ebj);
    w_delta_wenu->Fill (delta_phi_ebj, MyWeight*scalFac_b);
    w_deltaR_wenu->Fill (DR_ebj, MyWeight*scalFac_b);
    math::XYZTLorentzVector belectron;
    belectron = vect_ele[0].p4() + vect_bjets[0].p4();
    w_mass_wenu_blepton->Fill(belectron.mass(), MyWeight*scalFac_b);
    if (ist) {
      t_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
      t_pt_W_wenu->Fill (wenu_pt, MyWeight*scalFac_b);
      t_eta_W_wenu->Fill (wenu_eta, MyWeight*scalFac_b);
      t_delta_wenu->Fill (delta_phi_ebj, MyWeight*scalFac_b);
      t_deltaR_wenu->Fill (DR_ebj, MyWeight*scalFac_b);
      t_mass_wenu_blepton->Fill(belectron.mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
      b_pt_W_wenu->Fill (wenu_pt, MyWeight*scalFac_b);
      b_eta_W_wenu->Fill (wenu_eta, MyWeight*scalFac_b);
      b_delta_wenu->Fill (delta_phi_ebj, MyWeight*scalFac_b);
      b_deltaR_wenu->Fill (DR_ebj, MyWeight*scalFac_b);
      b_mass_wenu_blepton->Fill(belectron.mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
      c_pt_W_wenu->Fill (wenu_pt, MyWeight*scalFac_b);
      c_eta_W_wenu->Fill (wenu_eta, MyWeight*scalFac_b);
      c_delta_wenu->Fill (delta_phi_ebj, MyWeight*scalFac_b);
      c_deltaR_wenu->Fill (DR_ebj, MyWeight*scalFac_b);
      c_mass_wenu_blepton->Fill(belectron.mass(), MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      h_mt_wenu_b->Fill (mt_wenu);
      w_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b); 
      w_pt_W_wenu_b->Fill (wenu_pt, MyWeight*scalFac_b);
      w_eta_W_wenu_b->Fill (wenu_eta, MyWeight*scalFac_b);
      w_delta_wenu_b->Fill (delta_phi_ebj, MyWeight*scalFac_b);
      w_deltaR_wenu_b->Fill (DR_ebj, MyWeight*scalFac_b);
      w_mass_wenu_blepton_b->Fill(belectron.mass(), MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b);
	t_pt_W_wenu_b->Fill (wenu_pt, MyWeight*scalFac_b);
	t_eta_W_wenu_b->Fill (wenu_eta, MyWeight*scalFac_b);
	t_delta_wenu_b->Fill (delta_phi_ebj, MyWeight*scalFac_b);
	t_deltaR_wenu_b->Fill (DR_ebj, MyWeight*scalFac_b);
	t_mass_wenu_blepton_b->Fill(belectron.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b);
	b_pt_W_wenu_b->Fill (wenu_pt, MyWeight*scalFac_b);
	b_eta_W_wenu_b->Fill (wenu_eta, MyWeight*scalFac_b);
	b_delta_wenu_b->Fill (delta_phi_ebj, MyWeight*scalFac_b);
	b_deltaR_wenu_b->Fill (DR_ebj, MyWeight*scalFac_b);
	b_mass_wenu_blepton_b->Fill(belectron.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b);
	c_pt_W_wenu_b->Fill (wenu_pt, MyWeight*scalFac_b);
	c_eta_W_wenu_b->Fill (wenu_eta, MyWeight*scalFac_b);
	c_delta_wenu_b->Fill (delta_phi_ebj, MyWeight*scalFac_b);
	c_deltaR_wenu_b->Fill (DR_ebj, MyWeight*scalFac_b);
	c_mass_wenu_blepton_b->Fill(belectron.mass(), MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      h_mt_wenu_bb->Fill (mt_wenu);
      w_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
      w_pt_W_wenu_bb->Fill (wenu_pt, MyWeight*scalFac_b);
      w_eta_W_wenu_bb->Fill (wenu_eta, MyWeight*scalFac_b);
      w_delta_wenu_bb->Fill (delta_phi_ebj, MyWeight*scalFac_b);
      w_deltaR_wenu_bb->Fill (DR_ebj, MyWeight*scalFac_b);
      delta_phi_ebjbj = fabs(vect_bjets[0].phi() - vect_bjets[1].phi());
      delta_eta_ebjbj = fabs(vect_bjets[0].eta() - vect_bjets[1].eta());
      if (delta_phi_ebjbj > acos (-1)) delta_phi_ebjbj = 2 * acos (-1) - delta_phi_ebjbj;
      DR_ebjbj = TMath::Sqrt(delta_phi_ebjbj*delta_phi_ebjbj + delta_eta_ebjbj*delta_eta_ebjbj);
      w_delta_wenu_2b->Fill (delta_phi_ebjbj, MyWeight*scalFac_b);
      w_deltaR_wenu_2b->Fill (DR_ebjbj, MyWeight*scalFac_b);
      w_mass_wenu_blepton_bb->Fill(belectron.mass(), MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
	t_pt_W_wenu_bb->Fill (wenu_pt, MyWeight*scalFac_b);
	t_eta_W_wenu_bb->Fill (wenu_eta, MyWeight*scalFac_b);
	t_delta_wenu_bb->Fill (delta_phi_ebj, MyWeight*scalFac_b);
	t_delta_wenu_2b->Fill (delta_phi_ebjbj, MyWeight*scalFac_b);
	t_deltaR_wenu_bb->Fill (DR_ebj, MyWeight*scalFac_b);
	t_deltaR_wenu_2b->Fill (DR_ebjbj, MyWeight*scalFac_b);
	t_mass_wenu_blepton_bb->Fill(belectron.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
	b_pt_W_wenu_bb->Fill (wenu_pt, MyWeight*scalFac_b);
	b_eta_W_wenu_bb->Fill (wenu_eta, MyWeight*scalFac_b);
	b_delta_wenu_bb->Fill (delta_phi_ebj, MyWeight*scalFac_b);
	b_delta_wenu_2b->Fill (delta_phi_ebjbj, MyWeight*scalFac_b);
	b_deltaR_wenu_bb->Fill (DR_ebj, MyWeight*scalFac_b);
	b_deltaR_wenu_2b->Fill (DR_ebjbj, MyWeight*scalFac_b);
	b_mass_wenu_blepton_bb->Fill(belectron.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
	c_pt_W_wenu_bb->Fill (wenu_pt, MyWeight*scalFac_b);
	c_eta_W_wenu_bb->Fill (wenu_eta, MyWeight*scalFac_b);
	c_delta_wenu_bb->Fill (delta_phi_ebj, MyWeight*scalFac_b);
	c_delta_wenu_2b->Fill (delta_phi_ebjbj, MyWeight*scalFac_b);
	c_deltaR_wenu_bb->Fill (DR_ebj, MyWeight*scalFac_b);
	c_deltaR_wenu_2b->Fill (DR_ebjbj, MyWeight*scalFac_b);
	c_mass_wenu_blepton_bb->Fill(belectron.mass(), MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ Zee PLOTS

  if ((lepton_=="electron" || lepton_ == "electronQCD") && iele1!=-1 && Nj > 1 && Nb > 0 && vtx_cut && !met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_mass_ee_wide->Fill (diele_mass, MyWeight*scalFac_b);
    if (ist) {
      t_mass_ee_wide->Fill (diele_mass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mass_ee_wide->Fill (diele_mass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mass_ee_wide->Fill (diele_mass, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      if (ist) {
        t_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
    }
  }

  if (ee_event && met_cut && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi());
    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
    h_mass_ee->Fill (diele_mass);
    w_mass_ee->Fill (diele_mass, MyWeight*scalFac_b);
    w_pt_Z_ee->Fill (diele_pt, MyWeight*scalFac_b);
    w_delta_ee->Fill (delta_phi_ee, MyWeight*scalFac_b);
    math::XYZTLorentzVector zj_ee_p = vect_jets[0].p4() + z_ee;
    double zj_ee_mass = zj_ee_p.mass();
    w_mass_Zj_ee->Fill (zj_ee_mass, MyWeight*scalFac_b);
    if (ist) {
      t_mass_ee->Fill (diele_mass, MyWeight*scalFac_b);
      t_pt_Z_ee->Fill (diele_pt, MyWeight*scalFac_b);
      t_mass_Zj_ee->Fill (zj_ee_mass, MyWeight*scalFac_b);
      t_delta_ee->Fill (delta_phi_ee, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mass_ee->Fill (diele_mass, MyWeight*scalFac_b);
      b_pt_Z_ee->Fill (diele_pt, MyWeight*scalFac_b);
      b_mass_Zj_ee->Fill (zj_ee_mass, MyWeight*scalFac_b);
      b_delta_ee->Fill (delta_phi_ee, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mass_ee->Fill (diele_mass, MyWeight*scalFac_b);
      c_pt_Z_ee->Fill (diele_pt, MyWeight*scalFac_b);
      c_mass_Zj_ee->Fill (zj_ee_mass, MyWeight*scalFac_b);
      c_delta_ee->Fill (delta_phi_ee, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
      w_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
      math::XYZTLorentzVector zb_ee_p = vect_jets[0].p4() + z_ee;
      double zb_ee_mass = zb_ee_p.mass();
      w_mass_Zj_ee_b->Fill (zb_ee_mass, MyWeight*scalFac_b);
      double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi());
      if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
      w_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
      if (Nj == 1) {
        w_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        w_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
      }
      if (ist) {
        t_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
        t_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        t_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        t_mass_Zj_ee_b->Fill (zb_ee_mass, MyWeight*scalFac_b);
	if (Nj == 1) {
          t_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
	  t_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
        b_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        b_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        b_mass_Zj_ee_b->Fill (zb_ee_mass, MyWeight*scalFac_b);
	if (Nj == 1) {
          b_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
	  b_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
        c_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        c_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        c_mass_Zj_ee_b->Fill (zb_ee_mass, MyWeight*scalFac_b);
        if (Nj == 1) {
          c_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
	  c_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        }
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      w_mass_ee_bb->Fill (diele_mass, MyWeight*scalFac_b);
      w_pt_Z_ee_bb->Fill (diele_pt, MyWeight*scalFac_b);
      math::XYZTLorentzVector zb_ee_p = vect_jets[0].p4() + z_ee;
      double zb_ee_mass = zb_ee_p.mass();
      w_mass_Zj_ee_bb->Fill (zb_ee_mass, MyWeight*scalFac_b);
      double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi());
      if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
      w_delta_ee_bb->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
      if (ist) {
        t_mass_ee_bb->Fill (diele_mass, MyWeight*scalFac_b);
        t_pt_Z_ee_bb->Fill (diele_pt, MyWeight*scalFac_b);
        t_delta_ee_bb->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        t_mass_Zj_ee_bb->Fill (zb_ee_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mass_ee_bb->Fill (diele_mass, MyWeight*scalFac_b);
        b_pt_Z_ee_bb->Fill (diele_pt, MyWeight*scalFac_b);
        b_delta_ee_bb->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        b_mass_Zj_ee_bb->Fill (zb_ee_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mass_ee_bb->Fill (diele_mass, MyWeight*scalFac_b);
        c_pt_Z_ee_bb->Fill (diele_pt, MyWeight*scalFac_b);
        c_delta_ee_bb->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        c_mass_Zj_ee_bb->Fill (zb_ee_mass, MyWeight*scalFac_b);
      }
    }
  }


  // ++++++++ WmNu PLOTS

  if (wmnu_event && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_mt_wmnu_wide->Fill (mt_wmnu);
    w_mt_wmnu_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
    if (ist) {
      t_mt_wmnu_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mt_wmnu_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mt_wmnu_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      h_mt_wmnu_b_wide->Fill (mt_wmnu);
      w_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      h_mt_wmnu_bb_wide->Fill (mt_wmnu);
      w_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
    }
  }


  double delta_phi_mj=0;
  double delta_eta_mj=0;
  double DR_mj=0;
  double delta_phi_mbj=0;
  double delta_eta_mbj=0;
  double DR_mbj=0;
  double delta_phi_mbjbj=0;
  double delta_eta_mbjbj=0;
  double DR_mbjbj=0;

  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_mt_wmnu->Fill (mt_wmnu);
    w_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
    w_pt_W_wmnu->Fill (wmnu_pt, MyWeight*scalFac_b);
    w_eta_W_wmnu->Fill (wmnu_eta, MyWeight*scalFac_b);
    delta_phi_mj = fabs(vect_muon[0].phi() - vect_jets[0].phi());
    delta_eta_mj = fabs(vect_muon[0].eta() - vect_jets[0].eta());
    if (delta_phi_mj > acos (-1)) delta_phi_mj = 2 * acos (-1) - delta_phi_mj;
    DR_mj = TMath::Sqrt(delta_phi_mj*delta_phi_mj + delta_eta_mj*delta_eta_mj);
    delta_phi_mbj = fabs(vect_muon[0].phi() - vect_bjets[0].phi());
    delta_eta_mbj = fabs(vect_muon[0].eta() - vect_bjets[0].eta());
    if (delta_phi_mbj > acos (-1)) delta_phi_mbj = 2 * acos (-1) - delta_phi_mbj;
    DR_mbj = TMath::Sqrt(delta_phi_mbj*delta_phi_mbj + delta_eta_mbj*delta_eta_mbj);
    w_delta_wmnu->Fill (delta_phi_mbj, MyWeight*scalFac_b);
    w_deltaR_wmnu->Fill (DR_mbj, MyWeight*scalFac_b);
    math::XYZTLorentzVector bmuon;
    bmuon = vect_muon[0].p4() + vect_bjets[0].p4();
    w_mass_wmnu_blepton->Fill(bmuon.mass(), MyWeight*scalFac_b);
    if (ist) {
      t_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
      t_pt_W_wmnu->Fill (wmnu_pt, MyWeight*scalFac_b);
      t_eta_W_wmnu->Fill (wmnu_eta, MyWeight*scalFac_b);
      t_delta_wmnu->Fill (delta_phi_mbj, MyWeight*scalFac_b);
      t_deltaR_wmnu->Fill (DR_mbj, MyWeight*scalFac_b);
      t_mass_wmnu_blepton->Fill(bmuon.mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
      b_pt_W_wmnu->Fill (wmnu_pt, MyWeight*scalFac_b);
      b_eta_W_wmnu->Fill (wmnu_eta, MyWeight*scalFac_b);
      b_delta_wmnu->Fill (delta_phi_mbj, MyWeight*scalFac_b);
      b_deltaR_wmnu->Fill (DR_mbj, MyWeight*scalFac_b);
      b_mass_wmnu_blepton->Fill(bmuon.mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
      c_pt_W_wmnu->Fill (wmnu_pt, MyWeight*scalFac_b);
      c_eta_W_wmnu->Fill (wmnu_eta, MyWeight*scalFac_b);
      c_delta_wmnu->Fill (delta_phi_mbj, MyWeight*scalFac_b);
      c_deltaR_wmnu->Fill (DR_mbj, MyWeight*scalFac_b);
      c_mass_wmnu_blepton->Fill(bmuon.mass(), MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      h_mt_wmnu_b->Fill (mt_wmnu);
      w_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
      w_pt_W_wmnu_b->Fill (wmnu_pt, MyWeight*scalFac_b);
      w_eta_W_wmnu_b->Fill (wmnu_eta, MyWeight*scalFac_b);
      w_delta_wmnu_b->Fill (delta_phi_mbj, MyWeight*scalFac_b);
      w_deltaR_wmnu_b->Fill (DR_mbj, MyWeight*scalFac_b);
      w_mass_wmnu_blepton_b->Fill(bmuon.mass(), MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
	t_pt_W_wmnu_b->Fill (wmnu_pt, MyWeight*scalFac_b);
	t_eta_W_wmnu_b->Fill (wmnu_eta, MyWeight*scalFac_b);
	t_delta_wmnu_b->Fill (delta_phi_mbj, MyWeight*scalFac_b);
	t_deltaR_wmnu_b->Fill (DR_mbj, MyWeight*scalFac_b);
	t_mass_wmnu_blepton_b->Fill(bmuon.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
	b_pt_W_wmnu_b->Fill (wmnu_pt, MyWeight*scalFac_b);
	b_eta_W_wmnu_b->Fill (wmnu_eta, MyWeight*scalFac_b);
	b_delta_wmnu_b->Fill (delta_phi_mbj, MyWeight*scalFac_b);
	b_deltaR_wmnu_b->Fill (DR_mbj, MyWeight*scalFac_b);
	b_mass_wmnu_blepton_b->Fill(bmuon.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
	c_pt_W_wmnu_b->Fill (wmnu_pt, MyWeight*scalFac_b);
	c_eta_W_wmnu_b->Fill (wmnu_eta, MyWeight*scalFac_b);
	c_delta_wmnu_b->Fill (delta_phi_mbj, MyWeight*scalFac_b);
	c_deltaR_wmnu_b->Fill (DR_mbj, MyWeight*scalFac_b);
	c_mass_wmnu_blepton_b->Fill(bmuon.mass(), MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      h_mt_wmnu_bb->Fill (mt_wmnu);
      w_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
      w_pt_W_wmnu_bb->Fill (wmnu_pt, MyWeight*scalFac_b);
      w_eta_W_wmnu_bb->Fill (wmnu_eta, MyWeight*scalFac_b);
      w_delta_wmnu_bb->Fill (delta_phi_mbj, MyWeight*scalFac_b);
      w_deltaR_wmnu_bb->Fill (DR_mbj, MyWeight*scalFac_b);
      delta_phi_mbjbj = fabs(vect_bjets[0].phi() - vect_bjets[1].phi());
      delta_eta_mbjbj = fabs(vect_bjets[0].eta() - vect_bjets[1].eta());
      if (delta_phi_mbjbj > acos (-1)) delta_phi_mbjbj = 2 * acos (-1) - delta_phi_mbjbj;
      DR_mbjbj = TMath::Sqrt(delta_phi_mbjbj*delta_phi_mbjbj + delta_eta_mbjbj*delta_eta_mbjbj);
      w_delta_wmnu_2b->Fill (delta_phi_mbjbj, MyWeight*scalFac_b);
      w_deltaR_wmnu_2b->Fill (DR_mbjbj, MyWeight*scalFac_b);
      w_mass_wmnu_blepton_bb->Fill(bmuon.mass(), MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
	t_pt_W_wmnu_bb->Fill (wmnu_pt, MyWeight*scalFac_b);
	t_eta_W_wmnu_bb->Fill (wmnu_eta, MyWeight*scalFac_b);
	t_delta_wmnu_bb->Fill (delta_phi_mbj, MyWeight*scalFac_b);
	t_delta_wmnu_2b->Fill (delta_phi_mbjbj, MyWeight*scalFac_b);
	t_deltaR_wmnu_bb->Fill (DR_mbj, MyWeight*scalFac_b);
	t_deltaR_wmnu_2b->Fill (DR_mbjbj, MyWeight*scalFac_b);
	t_mass_wmnu_blepton_bb->Fill(bmuon.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
	b_pt_W_wmnu_bb->Fill (wmnu_pt, MyWeight*scalFac_b);
	b_eta_W_wmnu_bb->Fill (wmnu_eta, MyWeight*scalFac_b);
	b_delta_wmnu_bb->Fill (delta_phi_mbj, MyWeight*scalFac_b);
	b_delta_wmnu_2b->Fill (delta_phi_mbjbj, MyWeight*scalFac_b);
	b_deltaR_wmnu_bb->Fill (DR_mbj, MyWeight*scalFac_b);
	b_deltaR_wmnu_2b->Fill (DR_mbjbj, MyWeight*scalFac_b);
	b_mass_wmnu_blepton_bb->Fill(bmuon.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
	c_pt_W_wmnu_bb->Fill (wmnu_pt, MyWeight*scalFac_b);
	c_eta_W_wmnu_bb->Fill (wmnu_eta, MyWeight*scalFac_b);
	c_delta_wmnu_bb->Fill (delta_phi_mbj, MyWeight*scalFac_b);
	c_delta_wmnu_2b->Fill (delta_phi_mbjbj, MyWeight*scalFac_b);
	c_deltaR_wmnu_bb->Fill (DR_mbj, MyWeight*scalFac_b);
	c_deltaR_wmnu_2b->Fill (DR_mbjbj, MyWeight*scalFac_b);
	c_mass_wmnu_blepton_bb->Fill(bmuon.mass(), MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ Zmm PLOTS

  if ((lepton_=="muon" || lepton_ == "muonQCD") && imuon1!=-1 && Nj > 1 && Nb>0 && vtx_cut && !met_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_mass_mm_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
    if (ist) {
      t_mass_mm_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mass_mm_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mass_mm_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
    }
    if (Nb > 0) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      if (ist) {
        t_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
    }
  }

  if (mm_event && met_cut && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi());
    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
    h_mass_mm->Fill (dimuon_mass);
    w_mass_mm->Fill (dimuon_mass, MyWeight);
    w_pt_Z_mm->Fill (dimuon_pt, MyWeight);
    w_delta_mm->Fill (delta_phi_mm, MyWeight);
    math::XYZTLorentzVector zj_mm_p = vect_jets[0].p4() + z_mm;
    double zj_mm_mass = zj_mm_p.mass();
    w_mass_Zj_mm->Fill (zj_mm_mass, MyWeight);
    if (ist) {
      t_mass_mm->Fill (dimuon_mass, MyWeight);
      t_pt_Z_mm->Fill (dimuon_pt, MyWeight);
      t_mass_Zj_mm->Fill (zj_mm_mass, MyWeight);
      t_delta_mm->Fill (delta_phi_mm, MyWeight);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mass_mm->Fill (dimuon_mass, MyWeight);
      b_pt_Z_mm->Fill (dimuon_pt, MyWeight);
      b_mass_Zj_mm->Fill (zj_mm_mass, MyWeight);
      b_delta_mm->Fill (delta_phi_mm, MyWeight);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mass_mm->Fill (dimuon_mass, MyWeight);
      c_pt_Z_mm->Fill (dimuon_pt, MyWeight);
      c_mass_Zj_mm->Fill (zj_mm_mass, MyWeight);
      c_delta_mm->Fill (delta_phi_mm, MyWeight);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
      w_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
      math::XYZTLorentzVector zb_mm_p = vect_jets[0].p4() + z_mm;
      double zb_mm_mass = zb_mm_p.mass();
      w_mass_Zj_mm_b->Fill (zb_mm_mass, MyWeight*scalFac_b);
      double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi());
      if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
      w_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
      if (Nj == 1) {
        w_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        w_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
      }
      if (ist) {
        t_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
        t_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        t_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        t_mass_Zj_mm_b->Fill (zb_mm_mass, MyWeight*scalFac_b);
        if (Nj == 1) {
          t_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
          t_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
        b_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        b_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        b_mass_Zj_mm_b->Fill (zb_mm_mass, MyWeight*scalFac_b);
        if (Nj == 1) {
          b_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
          b_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
        c_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        c_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        c_mass_Zj_mm_b->Fill (zb_mm_mass, MyWeight*scalFac_b);
        if (Nj == 1) {
          c_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
          c_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        }
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      w_mass_mm_bb->Fill (dimuon_mass, MyWeight*scalFac_b);
      w_pt_Z_mm_bb->Fill (dimuon_pt, MyWeight*scalFac_b);
      math::XYZTLorentzVector zb_mm_p = vect_jets[0].p4() + z_mm;
      double zb_mm_mass = zb_mm_p.mass();
      w_mass_Zj_mm_bb->Fill (zb_mm_mass, MyWeight*scalFac_b);
      double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi());
      if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
      w_delta_mm_bb->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
      if (ist) {
        t_mass_mm_bb->Fill (dimuon_mass, MyWeight*scalFac_b);
        t_pt_Z_mm_bb->Fill (dimuon_pt, MyWeight*scalFac_b);
        t_delta_mm_bb->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        t_mass_Zj_mm_bb->Fill (zb_mm_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_mass_mm_bb->Fill (dimuon_mass, MyWeight*scalFac_b);
        b_pt_Z_mm_bb->Fill (dimuon_pt, MyWeight*scalFac_b);
        b_delta_mm_bb->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        b_mass_Zj_mm_bb->Fill (zb_mm_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_mass_mm_bb->Fill (dimuon_mass, MyWeight*scalFac_b);
        c_pt_Z_mm_bb->Fill (dimuon_pt, MyWeight*scalFac_b);
        c_delta_mm_bb->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        c_mass_Zj_mm_bb->Fill (zb_mm_mass, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ MISC PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_pu_weights->Fill (MyWeight*scalFac_b);
    h_tracks->Fill (tracks->size());
    w_tracks->Fill (tracks->size(), MyWeight*scalFac_b);
    h_recoVTX->Fill (NVtx);
    w_recoVTX->Fill (NVtx, MyWeight*scalFac_b);
  }

  // ++++++++  ELECTRONS PLOTS

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_first_ele_pt->Fill (vect_ele[0].pt());
    w_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
    w_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
    w_first_ele_iso->Fill ((vect_ele[0].chargedHadronIso() + fmax(vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso(),0))/vect_ele[0].et(), MyWeight*scalFac_b);
    if (ee_event) { // temporary !!!
      w_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
      w_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
    }
    if (ist) {
      t_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
      t_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
      t_first_ele_iso->Fill ((vect_ele[0].chargedHadronIso() + fmax(vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso(),0))/vect_ele[0].et(), MyWeight*scalFac_b);
      if (ee_event) { // temporary !!!
        t_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
        t_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
      b_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
      b_first_ele_iso->Fill ((vect_ele[0].chargedHadronIso() + fmax(vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso(),0))/vect_ele[0].et(), MyWeight*scalFac_b);
      if (ee_event) { // temporary !!!
        b_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
        b_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
      c_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
      c_first_ele_iso->Fill ((vect_ele[0].chargedHadronIso() + fmax(vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso(),0))/vect_ele[0].et(), MyWeight*scalFac_b);
      if (ee_event) { // temporary !!!
        c_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
        c_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
      }
    }
  }

  if (wenu_event && mt_cut_wenu && vtx_cut && Nb == 1) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_first_ele_pt_b->Fill (vect_ele[0].pt());
    w_first_ele_pt_b->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  }

  if (wenu_event && mt_cut_wenu && vtx_cut && Nb > 1) {
    scalFac_b = btagSF(isMC, vect_bjets, 2);
    h_first_ele_pt_bb->Fill (vect_ele[0].pt());
    w_first_ele_pt_bb->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  }

  // ++++++++ MUONS PLOTS

  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_first_muon_pt->Fill (vect_muon[0].pt());
    w_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
    w_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
    w_first_muon_iso->Fill ((vect_muon[0].chargedHadronIso() + fmax(vect_muon[0].neutralHadronIso() + vect_muon[0].photonIso() - 0.5*vect_muon[0].puChargedHadronIso(),0))/vect_muon[0].pt(), MyWeight*scalFac_b);
    if (mm_event) { // temporary !!!
      w_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
      w_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
    }
    if (ist) {
      t_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
      t_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
      t_first_muon_iso->Fill ((vect_muon[0].chargedHadronIso() + fmax(vect_muon[0].neutralHadronIso() + vect_muon[0].photonIso() - 0.5*vect_muon[0].puChargedHadronIso(),0))/vect_muon[0].pt(), MyWeight*scalFac_b);
      if (mm_event) { // temporary !!!
        t_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
        t_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
      b_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
      b_first_muon_iso->Fill ((vect_muon[0].chargedHadronIso() + fmax(vect_muon[0].neutralHadronIso() + vect_muon[0].photonIso() - 0.5*vect_muon[0].puChargedHadronIso(),0))/vect_muon[0].pt(), MyWeight*scalFac_b);
      if (mm_event) { // temporary !!!
        b_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
        b_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
      c_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
      c_first_muon_iso->Fill ((vect_muon[0].chargedHadronIso() + fmax(vect_muon[0].neutralHadronIso() + vect_muon[0].photonIso() - 0.5*vect_muon[0].puChargedHadronIso(),0))/vect_muon[0].pt(), MyWeight*scalFac_b);
      if (mm_event) { // temporary !!!
        c_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
        c_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
      }
    }
  }

  if (wmnu_event && mt_cut_wmnu && vtx_cut && Nb == 1) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_first_muon_pt_b ->Fill (vect_muon[0].pt());
    w_first_muon_pt_b ->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
  }

  if (wmnu_event && mt_cut_wmnu && vtx_cut && Nb > 1) {
    scalFac_b = btagSF(isMC, vect_bjets, 2);
    h_first_muon_pt_bb ->Fill (vect_muon[0].pt());
    w_first_muon_pt_bb ->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
  }

  // ++++++++ SVTX MASS PLOTS

  double sumVertexMassJet = 0.0;
  double sumVertexMassTrk = 0.0;
  double sumVertexMass = 0.0;

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {

    reco::SecondaryVertexTagInfo const * svTagInfos = vect_bjets[0].tagInfoSecondaryVertex("secondaryVertex");

    if (svTagInfos && svTagInfos->nVertices() > 0) {
      ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVecJet;
      for (reco::Vertex::trackRef_iterator track = svTagInfos->secondaryVertex(0).tracks_begin(); track != svTagInfos->secondaryVertex(0).tracks_end(); ++track) {
        const double kPionMass = 0.13957018;
        ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > vec;
        vec.SetPx( (*track)->px() );
        vec.SetPy( (*track)->py() );
        vec.SetPz( (*track)->pz() );
        vec.SetM (kPionMass);
        sumVecJet += vec;
      }
      sumVertexMassJet = sumVecJet.M();
      //cout << "VTX mass JET = " << sumVecJet.M() << endl;
    }

    ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVecTrk;
    for (size_t itrack=0; itrack < vect_bjets[0].associatedTracks().size(); ++itrack) {
      const double kPionMass = 0.13957018;
      ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > vec;
      vec.SetPx( vect_bjets[0].associatedTracks()[itrack]->px() );
      vec.SetPy( vect_bjets[0].associatedTracks()[itrack]->py() );
      vec.SetPz( vect_bjets[0].associatedTracks()[itrack]->pz() );
      vec.SetM (kPionMass);
      sumVecTrk += vec;
    }
    sumVertexMassTrk = sumVecTrk.M();
    //cout << "VTX mass TRK = " << sumVecTrk.M() << endl;

    if (svTagInfos && svTagInfos->nVertices() > 0) {
      const reco::Vertex &vertex = svTagInfos->secondaryVertex(0);
      reco::TrackKinematics vertexKinematics(vertex);
      math::XYZTLorentzVector vertexSum = vertexKinematics.weightedVectorSum();
      sumVertexMass = vertexSum.M();
      //cout << "VTX mass NEW = " << vertexSum.M() << endl;
    }

    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
    w_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
    w_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    if (ist) {
      t_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      t_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      t_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      b_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      b_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      c_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      c_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_SVTX_mass_jet_b->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      w_SVTX_mass_trk_b->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      w_SVTX_mass_b->Fill (sumVertexMass, MyWeight*scalFac_b);
      if (ist) {
	t_SVTX_mass_jet_b->Fill (sumVertexMassJet, MyWeight*scalFac_b);
	t_SVTX_mass_trk_b->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
	t_SVTX_mass_b->Fill (sumVertexMass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
	b_SVTX_mass_jet_b->Fill (sumVertexMassJet, MyWeight*scalFac_b);
	b_SVTX_mass_trk_b->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
	b_SVTX_mass_b->Fill (sumVertexMass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
	c_SVTX_mass_jet_b->Fill (sumVertexMassJet, MyWeight*scalFac_b);
	c_SVTX_mass_trk_b->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
	c_SVTX_mass_b->Fill (sumVertexMass, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      w_SVTX_mass_jet_bb->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      w_SVTX_mass_trk_bb->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      w_SVTX_mass_bb->Fill (sumVertexMass, MyWeight*scalFac_b);
      if (ist) {
	t_SVTX_mass_jet_bb->Fill (sumVertexMassJet, MyWeight*scalFac_b);
	t_SVTX_mass_trk_bb->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
	t_SVTX_mass_bb->Fill (sumVertexMass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
	b_SVTX_mass_jet_bb->Fill (sumVertexMassJet, MyWeight*scalFac_b);
	b_SVTX_mass_trk_bb->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
	b_SVTX_mass_bb->Fill (sumVertexMass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
	c_SVTX_mass_jet_bb->Fill (sumVertexMassJet, MyWeight*scalFac_b);
	c_SVTX_mass_trk_bb->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
	c_SVTX_mass_bb->Fill (sumVertexMass, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ CSV PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    double discrSVTX = vect_bjets[0].bDiscriminator("combinedSecondaryVertexBJetTags");
    w_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    if (ist) {
      t_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    }
    if (sumVertexMass > 0.0 ) {
      w_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      if (ist) {
        t_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
    }
    if (sumVertexMass == 0.0 ) {
      w_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      if (ist) {
        t_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
        c_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ BJP/JBP PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    double discrBJP = vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags");
    double discrJBP = vect_bjets[0].bDiscriminator("jetProbabilityBJetTags");

    w_BJP->Fill (discrBJP, MyWeight*scalFac_b);
    w_JBP->Fill (discrJBP, MyWeight*scalFac_b);
    if (sumVertexMass > 0.0) {
      w_BJP_mass->Fill (discrBJP, MyWeight*scalFac_b);
      w_JBP_mass->Fill (discrJBP, MyWeight*scalFac_b);
    }
    if (sumVertexMass == 0.0) {
      w_BJP_nomass->Fill (discrBJP, MyWeight*scalFac_b);
      w_JBP_nomass->Fill (discrJBP, MyWeight*scalFac_b);
    }
    if (ist) {
      t_BJP->Fill (discrBJP, MyWeight*scalFac_b);
      t_JBP->Fill (discrJBP, MyWeight*scalFac_b);
      if (sumVertexMass > 0.0) {
        t_BJP_mass->Fill (discrBJP, MyWeight*scalFac_b);
        t_JBP_mass->Fill (discrJBP, MyWeight*scalFac_b);
      }
      if (sumVertexMass == 0.0) {
        t_BJP_nomass->Fill (discrBJP, MyWeight*scalFac_b);
        t_JBP_nomass->Fill (discrJBP, MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_BJP->Fill (discrBJP, MyWeight*scalFac_b);
      b_JBP->Fill (discrJBP, MyWeight*scalFac_b);
      if (sumVertexMass > 0.0) {
        b_BJP_mass->Fill (discrBJP, MyWeight*scalFac_b);
        b_JBP_mass->Fill (discrJBP, MyWeight*scalFac_b);
      }
      if (sumVertexMass == 0.0) {
        b_BJP_nomass->Fill (discrBJP, MyWeight*scalFac_b);
        b_JBP_nomass->Fill (discrJBP, MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_BJP->Fill (discrBJP, MyWeight*scalFac_b);
      c_JBP->Fill (discrJBP, MyWeight*scalFac_b);
      if (sumVertexMass > 0.0) {
        c_BJP_mass->Fill (discrBJP, MyWeight*scalFac_b);
        c_JBP_mass->Fill (discrJBP, MyWeight*scalFac_b);
      }
      if (sumVertexMass == 0.0) {
        c_BJP_nomass->Fill (discrBJP, MyWeight*scalFac_b);
        c_JBP_nomass->Fill (discrJBP, MyWeight*scalFac_b);
      }
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_BJP_b->Fill (discrBJP, MyWeight*scalFac_b);
      w_JBP_b->Fill (discrJBP, MyWeight*scalFac_b);
      if (sumVertexMass > 0.0) {
	w_BJP_mass_b->Fill (discrBJP, MyWeight*scalFac_b);
	w_JBP_mass_b->Fill (discrJBP, MyWeight*scalFac_b);
      }
      if (sumVertexMass == 0.0) {
	w_BJP_nomass_b->Fill (discrBJP, MyWeight*scalFac_b);
	w_JBP_nomass_b->Fill (discrJBP, MyWeight*scalFac_b);
      }
      if (ist) {
	t_BJP_b->Fill (discrBJP, MyWeight*scalFac_b);
	t_JBP_b->Fill (discrJBP, MyWeight*scalFac_b);
	if (sumVertexMass > 0.0) {
	  t_BJP_mass_b->Fill (discrBJP, MyWeight*scalFac_b);
	  t_JBP_mass_b->Fill (discrJBP, MyWeight*scalFac_b);
	}
	if (sumVertexMass == 0.0) {
	  t_BJP_nomass_b->Fill (discrBJP, MyWeight*scalFac_b);
	  t_JBP_nomass_b->Fill (discrJBP, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
	b_BJP_b->Fill (discrBJP, MyWeight*scalFac_b);
	b_JBP_b->Fill (discrJBP, MyWeight*scalFac_b);
	if (sumVertexMass > 0.0) {
	  b_BJP_mass_b->Fill (discrBJP, MyWeight*scalFac_b);
	  b_JBP_mass_b->Fill (discrJBP, MyWeight*scalFac_b);
	}
	if (sumVertexMass == 0.0) {
	  b_BJP_nomass_b->Fill (discrBJP, MyWeight*scalFac_b);
	  b_JBP_nomass_b->Fill (discrJBP, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
	c_BJP_b->Fill (discrBJP, MyWeight*scalFac_b);
	c_JBP_b->Fill (discrJBP, MyWeight*scalFac_b);
	if (sumVertexMass > 0.0) {
	  c_BJP_mass_b->Fill (discrBJP, MyWeight*scalFac_b);
	  c_JBP_mass_b->Fill (discrJBP, MyWeight*scalFac_b);
	}
	if (sumVertexMass == 0.0) {
	  c_BJP_nomass_b->Fill (discrBJP, MyWeight*scalFac_b);
	  c_JBP_nomass_b->Fill (discrJBP, MyWeight*scalFac_b);
	}
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      w_BJP_bb->Fill (discrBJP, MyWeight*scalFac_b);
      w_JBP_bb->Fill (discrJBP, MyWeight*scalFac_b);
      if (sumVertexMass > 0.0) {
	w_BJP_mass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	w_JBP_mass_bb->Fill (discrJBP, MyWeight*scalFac_b);
      }
      if (sumVertexMass == 0.0) {
	w_BJP_nomass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	w_JBP_nomass_bb->Fill (discrJBP, MyWeight*scalFac_b);
      }
      if (ist) {
	t_BJP_bb->Fill (discrBJP, MyWeight*scalFac_b);
	t_JBP_bb->Fill (discrJBP, MyWeight*scalFac_b);
	if (sumVertexMass > 0.0) {
	  t_BJP_mass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	  t_JBP_mass_bb->Fill (discrJBP, MyWeight*scalFac_b);
	}
	if (sumVertexMass == 0.0) {
	  t_BJP_nomass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	  t_JBP_nomass_bb->Fill (discrJBP, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
	b_BJP_bb->Fill (discrBJP, MyWeight*scalFac_b);
	b_JBP_bb->Fill (discrJBP, MyWeight*scalFac_b);
	if (sumVertexMass > 0.0) {
	  b_BJP_mass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	  b_JBP_mass_bb->Fill (discrJBP, MyWeight*scalFac_b);
	}
	if (sumVertexMass == 0.0) {
	  b_BJP_nomass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	  b_JBP_nomass_bb->Fill (discrJBP, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
	c_BJP_bb->Fill (discrBJP, MyWeight*scalFac_b);
	c_JBP_bb->Fill (discrJBP, MyWeight*scalFac_b);
	if (sumVertexMass > 0.0) {
	  c_BJP_mass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	  c_JBP_mass_bb->Fill (discrJBP, MyWeight*scalFac_b);
	}
	if (sumVertexMass == 0.0) {
	  c_BJP_nomass_bb->Fill (discrBJP, MyWeight*scalFac_b);
	  c_JBP_nomass_bb->Fill (discrJBP, MyWeight*scalFac_b);
	}
      }
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb == 1 && vtx_cut) {
    double discrBJP = vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags");
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    if (ist) {
      b_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 1 && vtx_cut) {
    double discrBJP = vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags");
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    if (ist) {
      t_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 1 && vtx_cut) {
    double discrBJP2 = vect_bjets[1].bDiscriminator("jetBProbabilityBJetTags");
    scalFac_b = btagSF(isMC, vect_bjets, 2);
    w_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    if (ist) {
      t_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    }
  }

  // ++++++++ JETS PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_jetmultiplicity->Fill (Nj);
    w_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
    h_first_jet_pt->Fill (vect_jets[0].pt());
    w_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
    w_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    w_first_jet_mass->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
      t_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      t_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      t_first_jet_mass->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
      b_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      b_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      b_first_jet_mass->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
      c_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      c_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      c_first_jet_mass->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    h_second_jet_pt->Fill (vect_jets[1].pt());
    w_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
    w_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    w_second_jet_mass->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      t_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      t_second_jet_mass->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      b_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      b_second_jet_mass->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      c_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      c_second_jet_mass->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb==1) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    h_first_jet_pt_b->Fill (vect_jets[0].pt());
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    w_first_jet_mass_b->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      t_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      t_first_jet_mass_b->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      b_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      b_first_jet_mass_b->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      c_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      c_first_jet_mass_b->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb==1) {
    h_second_jet_pt_b->Fill (vect_jets[1].pt());
    w_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
    w_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    w_second_jet_mass_b->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      t_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      t_second_jet_mass_b->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      b_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      b_second_jet_mass_b->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      c_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      c_second_jet_mass_b->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb>1) {
    scalFac_b = btagSF(isMC, vect_bjets, 2);
    h_first_jet_pt_bb->Fill (vect_jets[0].pt());
    w_first_jet_pt_bb->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
    w_first_jet_eta_bb->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    w_first_jet_mass_bb->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_first_jet_pt_bb->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      t_first_jet_eta_bb->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      t_first_jet_mass_bb->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_first_jet_pt_bb->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      b_first_jet_eta_bb->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      b_first_jet_mass_bb->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_first_jet_pt_bb->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      c_first_jet_eta_bb->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      c_first_jet_mass_bb->Fill (vect_jets[0].mass(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb>1) {
    scalFac_b = btagSF(isMC, vect_bjets, 2);
    h_second_jet_pt_bb->Fill (vect_jets[1].pt());
    w_second_jet_pt_bb->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
    w_second_jet_eta_bb->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    w_second_jet_mass_bb->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_second_jet_pt_bb->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      t_second_jet_eta_bb->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      t_second_jet_mass_bb->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_second_jet_pt_bb->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      b_second_jet_eta_bb->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      b_second_jet_mass_bb->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_second_jet_pt_bb->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      c_second_jet_eta_bb->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      c_second_jet_mass_bb->Fill (vect_jets[1].mass(), MyWeight*scalFac_b);
    }
  }

  math::XYZTLorentzVector dijet;
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    dijet = vect_jets[0].p4() + vect_jets[1].p4();
    h_dijet_pt->Fill (dijet.pt());
    w_dijet_pt->Fill (dijet.pt(), MyWeight*scalFac_b);
    w_dijet_eta->Fill (dijet.eta(), MyWeight*scalFac_b);
    w_dijet_mass->Fill (dijet.mass(), MyWeight*scalFac_b);
    if (ist) {
      t_dijet_pt->Fill (dijet.pt(), MyWeight*scalFac_b);
      t_dijet_eta->Fill (dijet.eta(), MyWeight*scalFac_b);
      t_dijet_mass->Fill (dijet.mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_dijet_pt->Fill (dijet.pt(), MyWeight*scalFac_b);
      b_dijet_eta->Fill (dijet.eta(), MyWeight*scalFac_b);
      b_dijet_mass->Fill (dijet.mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_dijet_pt->Fill (dijet.pt(), MyWeight*scalFac_b);
      c_dijet_eta->Fill (dijet.eta(), MyWeight*scalFac_b);
      c_dijet_mass->Fill (dijet.mass(), MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      w_dijet_pt_b->Fill (dijet.pt(), MyWeight*scalFac_b);
      w_dijet_eta_b->Fill (dijet.eta(), MyWeight*scalFac_b);
      w_dijet_mass_b->Fill (dijet.mass(), MyWeight*scalFac_b);
      if (ist) {
	t_dijet_pt_b->Fill (dijet.pt(), MyWeight*scalFac_b);
	t_dijet_eta_b->Fill (dijet.eta(), MyWeight*scalFac_b);
	t_dijet_mass_b->Fill (dijet.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
	b_dijet_pt_b->Fill (dijet.pt(), MyWeight*scalFac_b);
	b_dijet_eta_b->Fill (dijet.eta(), MyWeight*scalFac_b);
	b_dijet_mass_b->Fill (dijet.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
	c_dijet_pt_b->Fill (dijet.pt(), MyWeight*scalFac_b);
	c_dijet_eta_b->Fill (dijet.eta(), MyWeight*scalFac_b);
	c_dijet_mass_b->Fill (dijet.mass(), MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      w_dijet_pt_bb->Fill (dijet.pt(), MyWeight*scalFac_b);
      w_dijet_eta_bb->Fill (dijet.eta(), MyWeight*scalFac_b);
      w_dijet_mass_bb->Fill (dijet.mass(), MyWeight*scalFac_b);
      if (ist) {
	t_dijet_pt_bb->Fill (dijet.pt(), MyWeight*scalFac_b);
	t_dijet_eta_bb->Fill (dijet.eta(), MyWeight*scalFac_b);
	t_dijet_mass_bb->Fill (dijet.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
	b_dijet_pt_bb->Fill (dijet.pt(), MyWeight*scalFac_b);
	b_dijet_eta_bb->Fill (dijet.eta(), MyWeight*scalFac_b);
	b_dijet_mass_bb->Fill (dijet.mass(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
	c_dijet_pt_bb->Fill (dijet.pt(), MyWeight*scalFac_b);
	c_dijet_eta_bb->Fill (dijet.eta(), MyWeight*scalFac_b);
	c_dijet_mass_bb->Fill (dijet.mass(), MyWeight*scalFac_b);
      }
    }
  }


  // ++++++++ B JETS PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
    h_first_bjet_pt->Fill (vect_bjets[0].pt());
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    w_first_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      t_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      t_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
      t_first_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      b_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      b_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
      b_first_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      c_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      c_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
      c_first_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 1 && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 2);
    h_second_bjet_pt->Fill (vect_bjets[1].pt());
    w_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
    w_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    w_second_bjet_mass->Fill (vect_bjets[1].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      t_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
      t_second_bjet_mass->Fill (vect_bjets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      b_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
      b_second_bjet_mass->Fill (vect_bjets[1].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      c_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
      c_second_bjet_mass->Fill (vect_bjets[1].mass(), MyWeight*scalFac_b);
    }
  }

  // ++++++++ SINGLE BJET

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb == 1 && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    w_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
    w_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    w_single_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    if (ist) {
      t_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      t_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
      t_single_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      b_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
      b_single_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      c_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
      c_single_bjet_mass->Fill (vect_bjets[0].mass(), MyWeight*scalFac_b);
    }
  }

  // ++++++++ EXTRA PLOTS

  int Nf = 0;
  int Nbk = 0;

  double Afb = 0;

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_bjets, 1);
    if (fabs (vect_bjets[0].eta()) > 0) Nf++;
    if (fabs (vect_bjets[0].eta()) < 0) Nbk++;
    if ((Nf+Nbk) != 0) Afb = (Nf - Nbk) / (Nf + Nbk);
    w_Afb->Fill (Afb, MyWeight*scalFac_b);
  }

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    h_scaleFactor_first_ele->Fill (scalFac_first_e);
    h_scaleFactor_second_ele->Fill (scalFac_second_e);
    if (Nb == 1) {
      if (isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_scaleFactor_first_ele->Fill (scalFac_first_e);
        b_scaleFactor_second_ele->Fill (scalFac_second_e);
      }
    }
  }
  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    h_scaleFactor_first_muon->Fill (scalFac_first_m);
    h_scaleFactor_second_muon->Fill (scalFac_second_m);
    if (Nb == 1) {
      if (isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
        b_scaleFactor_first_muon->Fill (scalFac_first_m);
        b_scaleFactor_second_muon->Fill (scalFac_second_m);
      }
    }
  }

  // ++++++++ OUTPUT COLLECTIONS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    if (Nb < 2) {
      scalFac_b = btagSF(isMC, vect_bjets, 1);
      myEventWeight->push_back(MyWeight*scalFac_b);
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_bjets, 2);
      myEventWeight->push_back(MyWeight*scalFac_b);
    }
  }

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[0].px(),vect_ele[0].py(),vect_ele[0].pz(),vect_ele[0].energy()));
    myDeltaPhiEJ->push_back(delta_phi_ej);
    myDeltaPhiEBJ->push_back(delta_phi_ebj);
    if (Nb > 1) myDeltaPhiEBJBJ->push_back(delta_phi_ebjbj);
    myDeltaREJ->push_back(DR_ej);
    myDeltaREBJ->push_back(DR_ebj);
    if (Nb > 1) myDeltaREBJBJ->push_back(DR_ebjbj);
    myWenuPt->push_back(wenu_pt);
    myWenuEta->push_back(wenu_eta);
  }

  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[0].px(),vect_muon[0].py(),vect_muon[0].pz(),vect_muon[0].energy()));
    myDeltaPhiMJ->push_back(delta_phi_mj);
    myDeltaPhiMBJ->push_back(delta_phi_mbj);
    if (Nb > 1) myDeltaPhiMBJBJ->push_back(delta_phi_mbjbj);
    myDeltaRMJ->push_back(DR_mj);
    myDeltaRMBJ->push_back(DR_mbj);
    if (Nb > 1) myDeltaRMBJBJ->push_back(DR_mbjbj);
    myWmnuPt->push_back(wmnu_pt);
    myWmnuEta->push_back(wmnu_eta);
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    for (unsigned int i=0; i<vect_jets.size(); ++i) {
      myJets->push_back(math::XYZTLorentzVector(vect_jets[i].px(),vect_jets[i].py(),vect_jets[i].pz(),vect_jets[i].energy()));
    }
    for (unsigned int i=0; i<vect_bjets.size(); ++i) {
      myBJets->push_back(math::XYZTLorentzVector(vect_bjets[i].px(),vect_bjets[i].py(),vect_bjets[i].pz(),vect_bjets[i].energy()));
      if (Nb <= 1) {
	scalFac_b = btagSF(isMC, vect_bjets, 1);
	myBJetsWeights->push_back(scalFac_b);
      }
      if (Nb >= 2) {
	scalFac_b = btagSF(isMC, vect_bjets, 2);
	myBJetsWeights->push_back(scalFac_b);
      }
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    myHt->push_back(Ht);
    myDijetPt->push_back(dijet.pt());
    myDijetEta->push_back(dijet.eta());
    myDijetMass->push_back(dijet.mass());
  }

  iEvent.put( myEventWeight, "myEventWeight" );

  iEvent.put( myElectrons, "myElectrons" );
  iEvent.put( myMuons, "myMuons" );

  iEvent.put( myJets, "myJets" );

  iEvent.put( myHt, "myHt" );

  iEvent.put( myBJetsWeights, "myBJetsWeights" );

  iEvent.put( myBJets, "myBJets" );

  iEvent.put( myWenuPt, "myWenuPt" );
  iEvent.put( myWmnuPt, "myWmnuPt" );
  iEvent.put( myWenuEta, "myWenuEta" );
  iEvent.put( myWmnuEta, "myWmnuEta" );

  iEvent.put( myDijetPt, "myDijetPt" );
  iEvent.put( myDijetEta, "myDijetEta" );
  iEvent.put( myDijetMass, "myDijetMass" );

  iEvent.put( myDeltaPhiEJ, "myDeltaPhiEJ" );
  iEvent.put( myDeltaPhiEBJ, "myDeltaPhiEBJ" );
  iEvent.put( myDeltaPhiEBJBJ, "myDeltaPhiEBJBJ" );
  iEvent.put( myDeltaPhiMJ, "myDeltaPhiMJ" );
  iEvent.put( myDeltaPhiMBJ, "myDeltaPhiMBJ" );
  iEvent.put( myDeltaPhiMBJBJ, "myDeltaPhiMBJBJ" );

  iEvent.put( myDeltaREJ, "myDeltaREJ" );
  iEvent.put( myDeltaREBJ, "myDeltaREBJ" );
  iEvent.put( myDeltaREBJBJ, "myDeltaREBJBJ" );
  iEvent.put( myDeltaRMJ, "myDeltaRMJ" );
  iEvent.put( myDeltaRMBJ, "myDeltaRMBJ" );
  iEvent.put( myDeltaRMBJBJ, "myDeltaRMBJBJ" );

}

// ------------ method called once each job just before starting event loop ------------
void WbAnalyzer::beginJob () {
  jetCorrectionUncertainty_ = new JetCorrectionUncertainty(path_ + "/" + "Summer13_V4_DATA_Uncertainty_AK5PFchs.txt");
  LumiWeights_ = edm::LumiReWeighting(path_ + "/" + "pileup_" + pileupMC_ + ".root", path_ + "/" + "pileup_2012_" + pileupDT_ + ".root", "pileup", "pileup");

  ElSF_  = new table(path_ + "/" + "ele_sc_id_iso.txt");
  ElSF2_ = new table(path_ + "/" + "ele_sc_hlt.txt");
  MuSF_  = new table(path_ + "/" + "muon_sc_hlt.txt");
  MuSF2_ = new table(path_ + "/" + "muon_sc_id.txt");
  MuSF3_ = new table(path_ + "/" + "muon_sc_iso.txt");
  BtSF_  = new table(path_ + "/" + "btag_eff.txt");   //btagging scale factors SFb = SFc
  LtSF_  = new table(path_ + "/" + "light_eff.txt");  //light flavour scale factors

}

// ------------ method called once each job just after ending the event loop ------------
void WbAnalyzer::endJob () {
  delete jetCorrectionUncertainty_;

  delete ElSF_;
  delete ElSF2_;
  delete MuSF_;
  delete MuSF2_;
  delete MuSF3_;
  delete BtSF_;
  delete LtSF_;

}

// ------------ method called when starting to processes a run ------------
void WbAnalyzer::beginRun (edm::Run & iRun, edm::EventSetup const & iSetup) {
  nprup = 0;

  edm::Handle<LHERunInfoProduct> run;
  if (iRun.getByLabel ("source", run)) {
    const lhef::HEPRUP heprup = run->heprup();
    nprup = heprup.NPRUP;
  }

}

// ------------ method called when ending the processing of a run ------------
void WbAnalyzer::endRun (edm::Run &, edm::EventSetup const &) {
}

// ------------ method called when starting to processes a luminosity block ------------
void WbAnalyzer::beginLuminosityBlock (edm::LuminosityBlock &, edm::EventSetup const &) {
}

// ------------ method called when ending the processing of a luminosity block ------------
void WbAnalyzer::endLuminosityBlock (edm::LuminosityBlock &, edm::EventSetup const &) {
}

// define this as a plug-in
DEFINE_FWK_MODULE (WbAnalyzer);
