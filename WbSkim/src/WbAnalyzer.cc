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
        if ((fabs(jets[j].partonFlavour()) == 5 || fabs(jets[j].partonFlavour()) == 4)) {
          if (i==j) {
            w = w * BtSF_->Val(jets[j].pt(), jets[j].eta());
          } else {
            w = w * (1.0 - BtSF_->Val(jets[j].pt(), jets[j].eta()));
          }
        } else {
          if (i==j) {
            w = w * LtSF_->Val(jets[j].pt(), jets[j].eta());
          } else {
            w = w * (1.0 - LtSF_->Val(jets[j].pt(), jets[j].eta()));
          }
        }
      }
      w1n = w1n + w;
    }

    if (k==0) return (1.0-w0n);     // >= 1 b tagged jet
    if (k==1) return (w1n);         // == 1 b tagged jet
    if (k==2) return (1.0-w0n-w1n); // >= 2 b tagged jet
    if (k==3) return (1.0-w0n-w1n); /// FIXME //
    return (0);

  };

  double jetResolutionCorrection(double jetEta, double jetPt, double jetPtGen, int syst) {

    if (jetPt <= 0 || jetPtGen <= 0) return jetPt;

    double correctionFactor[5]     = {1.052, 1.057, 1.096, 1.134, 1.288};
    double correctionFactorUp[5]   = {1.115, 1.114, 1.161, 1.228, 1.488};
    double correctionFactorDown[5] = {0.990, 1.001, 1.032, 1.042, 1.089};

    int index = 0;

    if (                      fabs(jetEta) <= 0.5) index = 0;
    if (fabs(jetEta) > 0.5 && fabs(jetEta) <= 1.1) index = 1;
    if (fabs(jetEta) > 1.1 && fabs(jetEta) <= 1.7) index = 2;
    if (fabs(jetEta) > 1.7 && fabs(jetEta) <= 2.3) index = 3;
    if (fabs(jetEta) > 2.3 && fabs(jetEta) <= 3.0) index = 4;

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
  TH1F* h_nmult0;
  TH1F* h_nmult1;

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

  TH1F*     h_jetmultiplicity;
  TH1F*     h_eventYields;
  TH1F*     w_eventYields;

  TH1F*     h_pu_weights;
  TH1F*     h_tracks;
  TH1F*     w_tracks;
  TH1F*     h_recoVTX;
  TH1F*     w_recoVTX;

  TH1F*     w_jetmultiplicity;
  TH1F*     b_jetmultiplicity;
  TH1F*     c_jetmultiplicity;
  TH1F*     t_jetmultiplicity;

  TH1F*     w_first_jet_pt;	// leading jet of any type
  TH1F*     b_first_jet_pt;
  TH1F*     c_first_jet_pt;
  TH1F*     t_first_jet_pt;
  TH1F*     w_first_jet_eta;
  TH1F*     b_first_jet_eta;
  TH1F*     c_first_jet_eta;
  TH1F*     t_first_jet_eta;
  TH1F*     w_second_jet_pt;
  TH1F*     b_second_jet_pt;
  TH1F*     c_second_jet_pt;
  TH1F*     t_second_jet_pt;
  TH1F*     w_second_jet_eta;
  TH1F*     b_second_jet_eta;
  TH1F*     c_second_jet_eta;
  TH1F*     t_second_jet_eta;
  TH1F*     w_third_jet_pt;
  TH1F*     b_third_jet_pt;
  TH1F*     c_third_jet_pt;
  TH1F*     t_third_jet_pt;
  TH1F*     w_third_jet_eta;
  TH1F*     b_third_jet_eta;
  TH1F*     c_third_jet_eta;
  TH1F*     t_third_jet_eta;

  TH1F*     w_first_jet_pt_b;	// leading jet with at least one b jet in the event
  TH1F*     b_first_jet_pt_b;
  TH1F*     c_first_jet_pt_b;
  TH1F*     t_first_jet_pt_b;
  TH1F*     w_first_jet_eta_b;
  TH1F*     b_first_jet_eta_b;
  TH1F*     c_first_jet_eta_b;
  TH1F*     t_first_jet_eta_b;
  TH1F*     w_second_jet_pt_b;
  TH1F*     b_second_jet_pt_b;
  TH1F*     c_second_jet_pt_b;
  TH1F*     t_second_jet_pt_b;
  TH1F*     w_second_jet_eta_b;
  TH1F*     b_second_jet_eta_b;
  TH1F*     c_second_jet_eta_b;
  TH1F*     t_second_jet_eta_b;
  TH1F*     w_third_jet_pt_b;
  TH1F*     b_third_jet_pt_b;
  TH1F*     c_third_jet_pt_b;
  TH1F*     t_third_jet_pt_b;
  TH1F*     w_third_jet_eta_b;
  TH1F*     b_third_jet_eta_b;
  TH1F*     c_third_jet_eta_b;
  TH1F*     t_third_jet_eta_b;

  TH1F*     w_bjetmultiplicity;
  TH1F*     b_bjetmultiplicity;
  TH1F*     c_bjetmultiplicity;
  TH1F*     t_bjetmultiplicity;

  TH1F*     w_first_bjet_pt;	// leading b jet
  TH1F*     b_first_bjet_pt;
  TH1F*     c_first_bjet_pt;
  TH1F*     t_first_bjet_pt;
  TH1F*     w_first_bjet_eta;
  TH1F*     b_first_bjet_eta;
  TH1F*     c_first_bjet_eta;
  TH1F*     t_first_bjet_eta;

  TH1F*     w_single_bjet_pt;	// only 1 b jet
  TH1F*     b_single_bjet_pt;
  TH1F*     c_single_bjet_pt;
  TH1F*     t_single_bjet_pt;
  TH1F*     w_single_bjet_eta;
  TH1F*     b_single_bjet_eta;
  TH1F*     c_single_bjet_eta;
  TH1F*     t_single_bjet_eta;

  TH1F*     w_second_bjet_pt;
  TH1F*     b_second_bjet_pt;
  TH1F*     c_second_bjet_pt;
  TH1F*     t_second_bjet_pt;
  TH1F*     w_second_bjet_eta;
  TH1F*     b_second_bjet_eta;
  TH1F*     c_second_bjet_eta;
  TH1F*     t_second_bjet_eta;
  TH1F*     w_third_bjet_pt;
  TH1F*     b_third_bjet_pt;
  TH1F*     c_third_bjet_pt;
  TH1F*     t_third_bjet_pt;
  TH1F*     w_third_bjet_eta;
  TH1F*     b_third_bjet_eta;
  TH1F*     c_third_bjet_eta;
  TH1F*     t_third_bjet_eta;

  TH1F*     w_first_ele_pt;
  TH1F*     w_first_ele_pt_b;
  TH1F*     w_first_ele_pt_bb;
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

  //  TH1F*     w_mass_em_wide;
  //  TH1F*     b_mass_em_wide;
  //  TH1F*     c_mass_em_wide;
  //  TH1F*     t_mass_em_wide;

  TH1F*     w_mt_wenu_wide;
  TH1F*     b_mt_wenu_wide;
  TH1F*     c_mt_wenu_wide;
  TH1F*     t_mt_wenu_wide;

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

  //  TH1F*     h_mass_em;
  //  TH1F*     w_mass_em;
  //  TH1F*     b_mass_em;
  //  TH1F*     c_mass_em;
  //  TH1F*     t_mass_em;

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

  TH1F*     w_pt_Z_ee;
  TH1F*     b_pt_Z_ee;
  TH1F*     c_pt_Z_ee;
  TH1F*     t_pt_Z_ee;
  
  TH1F*     w_pt_Z_mm;
  TH1F*     b_pt_Z_mm;
  TH1F*     c_pt_Z_mm;
  TH1F*     t_pt_Z_mm;

  TH1F*     w_pt_Z_em;
  TH1F*     b_pt_Z_em;
  TH1F*     c_pt_Z_em;
  TH1F*     t_pt_Z_em;

  TH1F*     w_single_pt_Z_ee_b;
  TH1F*     b_single_pt_Z_ee_b;
  TH1F*     c_single_pt_Z_ee_b;
  TH1F*     t_single_pt_Z_ee_b;
  
  TH1F*     w_single_pt_Z_mm_b;
  TH1F*     b_single_pt_Z_mm_b;
  TH1F*     c_single_pt_Z_mm_b;
  TH1F*     t_single_pt_Z_mm_b;

  TH1F*     w_single_pt_Z_em_b;
  TH1F*     b_single_pt_Z_em_b;
  TH1F*     c_single_pt_Z_em_b;
  TH1F*     t_single_pt_Z_em_b;

  TH1F*     w_mass_ee_b_wide;	// at least one b jet in the event
  TH1F*     b_mass_ee_b_wide;
  TH1F*     c_mass_ee_b_wide;
  TH1F*     t_mass_ee_b_wide;

  TH1F*     w_mass_mm_b_wide;
  TH1F*     b_mass_mm_b_wide;
  TH1F*     c_mass_mm_b_wide;
  TH1F*     t_mass_mm_b_wide;

  //  TH1F*     w_mass_em_b_wide;
  //  TH1F*     b_mass_em_b_wide;
  //  TH1F*     c_mass_em_b_wide;
  //  TH1F*     t_mass_em_b_wide;

  TH1F*     w_mt_wenu_b_wide;	// at least one b jet in the event
  TH1F*     b_mt_wenu_b_wide;
  TH1F*     c_mt_wenu_b_wide;
  TH1F*     t_mt_wenu_b_wide;

  TH1F*     w_mt_wmnu_b_wide;
  TH1F*     b_mt_wmnu_b_wide;
  TH1F*     c_mt_wmnu_b_wide;
  TH1F*     t_mt_wmnu_b_wide;

  TH1F*     w_mt_wenu_bb_wide;	// at least one b jet in the event
  TH1F*     b_mt_wenu_bb_wide;
  TH1F*     c_mt_wenu_bb_wide;
  TH1F*     t_mt_wenu_bb_wide;

  TH1F*     w_mt_wmnu_bb_wide;
  TH1F*     b_mt_wmnu_bb_wide;
  TH1F*     c_mt_wmnu_bb_wide;
  TH1F*     t_mt_wmnu_bb_wide;

  //  TH1F*     w_mass_em_bb_wide;
  //  TH1F*     b_mass_em_bb_wide;
  //  TH1F*     c_mass_em_bb_wide;
  //  TH1F*     t_mass_em_bb_wide;

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

  //  TH1F*     w_mass_em_b;
  //  TH1F*     b_mass_em_b;
  //  TH1F*     c_mass_em_b;
  //  TH1F*     t_mass_em_b;

  TH1F*     w_mt_wenu_b;	// at least one b jet in the event
  TH1F*     b_mt_wenu_b;
  TH1F*     c_mt_wenu_b;
  TH1F*     t_mt_wenu_b;

  TH1F*     w_mt_wmnu_b;
  TH1F*     b_mt_wmnu_b;
  TH1F*     c_mt_wmnu_b;
  TH1F*     t_mt_wmnu_b;

  TH1F*     w_mt_wenu_bb;	// at least one b jet in the event
  TH1F*     b_mt_wenu_bb;
  TH1F*     c_mt_wenu_bb;
  TH1F*     t_mt_wenu_bb;

  TH1F*     w_mt_wmnu_bb;
  TH1F*     b_mt_wmnu_bb;
  TH1F*     c_mt_wmnu_bb;
  TH1F*     t_mt_wmnu_bb;

  //  TH1F*     w_mass_em_bb;
  //  TH1F*     b_mass_em_bb;
  //  TH1F*     c_mass_em_bb;
  //  TH1F*     t_mass_em_bb;

  TH1F*     w_mass_Zj_ee;
  TH1F*     b_mass_Zj_ee;
  TH1F*     c_mass_Zj_ee;
  TH1F*     t_mass_Zj_ee;
  
  TH1F*     w_mass_Zj_mm;
  TH1F*     b_mass_Zj_mm;
  TH1F*     c_mass_Zj_mm;
  TH1F*     t_mass_Zj_mm;

  TH1F*     w_mass_Zj_em;
  TH1F*     b_mass_Zj_em;
  TH1F*     c_mass_Zj_em;
  TH1F*     t_mass_Zj_em;

  TH1F*     w_mass_Zj_ee_b;
  TH1F*     b_mass_Zj_ee_b;
  TH1F*     c_mass_Zj_ee_b;
  TH1F*     t_mass_Zj_ee_b;
  
  TH1F*     w_mass_Zj_mm_b;
  TH1F*     b_mass_Zj_mm_b;
  TH1F*     c_mass_Zj_mm_b;
  TH1F*     t_mass_Zj_mm_b;

  TH1F*     w_mass_Zj_em_b;
  TH1F*     b_mass_Zj_em_b;
  TH1F*     c_mass_Zj_em_b;
  TH1F*     t_mass_Zj_em_b;

  TH1F*     w_mass_Zj_ee_bb;
  TH1F*     b_mass_Zj_ee_bb;
  TH1F*     c_mass_Zj_ee_bb;
  TH1F*     t_mass_Zj_ee_bb;
  
  TH1F*     w_mass_Zj_mm_bb;
  TH1F*     b_mass_Zj_mm_bb;
  TH1F*     c_mass_Zj_mm_bb;
  TH1F*     t_mass_Zj_mm_bb;

//  TH1F*     w_mass_Zj_em_bb;
//  TH1F*     b_mass_Zj_em_bb;
//  TH1F*     c_mass_Zj_em_bb;
//  TH1F*     t_mass_Zj_em_bb;

  TH1F*     w_pt_Z_ee_b;
  TH1F*     b_pt_Z_ee_b;
  TH1F*     c_pt_Z_ee_b;
  TH1F*     t_pt_Z_ee_b;
  
  TH1F*     w_pt_Z_mm_b;
  TH1F*     b_pt_Z_mm_b;
  TH1F*     c_pt_Z_mm_b;
  TH1F*     t_pt_Z_mm_b;

//  TH1F*     w_pt_Z_em_b;
//  TH1F*     b_pt_Z_em_b;
//  TH1F*     c_pt_Z_em_b;
//  TH1F*     t_pt_Z_em_b;

  TH1F*     w_pt_Z_ee_bb;
  TH1F*     b_pt_Z_ee_bb;
  TH1F*     c_pt_Z_ee_bb;
  TH1F*     t_pt_Z_ee_bb;
  
  TH1F*     w_pt_Z_mm_bb;
  TH1F*     b_pt_Z_mm_bb;
  TH1F*     c_pt_Z_mm_bb;
  TH1F*     t_pt_Z_mm_bb;

//  TH1F*     w_pt_Z_em_bb;
//  TH1F*     b_pt_Z_em_bb;
//  TH1F*     c_pt_Z_em_bb;
//  TH1F*     t_pt_Z_em_bb;

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

  TH1F*     w_delta_wenuMet;
  TH1F*     b_delta_wenuMet;
  TH1F*     c_delta_wenuMet;
  TH1F*     t_delta_wenuMet;
  TH1F*     w_delta_wenuMet_b;
  TH1F*     b_delta_wenuMet_b;
  TH1F*     c_delta_wenuMet_b;
  TH1F*     t_delta_wenuMet_b;
  TH1F*     w_delta_wenuMet_bb;
  TH1F*     b_delta_wenuMet_bb;
  TH1F*     c_delta_wenuMet_bb;
  TH1F*     t_delta_wenuMet_bb;

  TH1F*     w_delta_wmnuMet;
  TH1F*     b_delta_wmnuMet;
  TH1F*     c_delta_wmnuMet;
  TH1F*     t_delta_wmnuMet;
  TH1F*     w_delta_wmnuMet_b;
  TH1F*     b_delta_wmnuMet_b;
  TH1F*     c_delta_wmnuMet_b;
  TH1F*     t_delta_wmnuMet_b;
  TH1F*     w_delta_wmnuMet_bb;
  TH1F*     b_delta_wmnuMet_bb;
  TH1F*     c_delta_wmnuMet_bb;
  TH1F*     t_delta_wmnuMet_bb;

  TH1F*     w_single_delta_ee_b;
  TH1F*     b_single_delta_ee_b;
  TH1F*     c_single_delta_ee_b;
  TH1F*     t_single_delta_ee_b;

  TH1F*     w_single_delta_mm_b;
  TH1F*     b_single_delta_mm_b;
  TH1F*     c_single_delta_mm_b;
  TH1F*     t_single_delta_mm_b;

  //  TH1F*     w_single_delta_em_b;
  //  TH1F*     b_single_delta_em_b;
  //  TH1F*     c_single_delta_em_b;
  //  TH1F*     t_single_delta_em_b;

  TH1F*     w_single_delta_wenuMet_b;
  TH1F*     b_single_delta_wenuMet_b;
  TH1F*     c_single_delta_wenuMet_b;
  TH1F*     t_single_delta_wenuMet_b;

  TH1F*     w_single_delta_wmnuMet_b;
  TH1F*     b_single_delta_wmnuMet_b;
  TH1F*     c_single_delta_wmnuMet_b;
  TH1F*     t_single_delta_wmnuMet_b;

  TH1F*     w_single_delta_wenuMet_bb;
  TH1F*     b_single_delta_wenuMet_bb;
  TH1F*     c_single_delta_wenuMet_bb;
  TH1F*     t_single_delta_wenuMet_bb;

  TH1F*     w_single_delta_wmnuMet_bb;
  TH1F*     b_single_delta_wmnuMet_bb;
  TH1F*     c_single_delta_wmnuMet_bb;
  TH1F*     t_single_delta_wmnuMet_bb;

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

  TH1F*     w_SVTX_mass;
  TH1F*     b_SVTX_mass;
  TH1F*     c_SVTX_mass;
  TH1F*     t_SVTX_mass;

  TH1F*     w_BJP;
  TH1F*     b_BJP;
  TH1F*     c_BJP;
  TH1F*     t_BJP;

  TH1F*     w_JBP;
  TH1F*     b_JBP;
  TH1F*     c_JBP;
  TH1F*     t_JBP;

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

  TH1F*     w_MET;
  TH1F*     b_MET;
  TH1F*     c_MET;
  TH1F*     t_MET;
  TH1F*     w_MET_sign;
  TH1F*     b_MET_sign;
  TH1F*     c_MET_sign;
  TH1F*     t_MET_sign;

  TH1F*     w_MET_b;
  TH1F*     b_MET_b;
  TH1F*     c_MET_b;
  TH1F*     t_MET_b;
  TH1F*     w_MET_sign_b;
  TH1F*     b_MET_sign_b;
  TH1F*     c_MET_sign_b;
  TH1F*     t_MET_sign_b;

  TH1F*     w_MET_bb;
  TH1F*     b_MET_bb;
  TH1F*     c_MET_bb;
  TH1F*     t_MET_bb;
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

  h_pu_weights =        fs->make < TH1F > ("h_pu_weights",      "h_pu_weights;PU weight", 10, 0, 10);

  h_tracks =            fs->make < TH1F > ("h_tracks",          "h_tracks;N_tracks", 100, 0, 2500);
  w_tracks = 	        fs->make < TH1F > ("w_tracks",          "w_tracks;N_tracks", 100, 0, 2500);
  h_recoVTX =           fs->make < TH1F > ("h_recoVTX",         "h_recoVTX;N_vtx", 40, 0., 40.);
  w_recoVTX =           fs->make < TH1F > ("w_recoVTX",         "w_recoVTX;N_vtx", 40, 0., 40.);

  w_jetmultiplicity =   fs->make < TH1F > ("w_jetmultiplicity", "w_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  b_jetmultiplicity =   fs->make < TH1F > ("b_jetmultiplicity", "b_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  c_jetmultiplicity =   fs->make < TH1F > ("c_jetmultiplicity", "c_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  t_jetmultiplicity =   fs->make < TH1F > ("t_jetmultiplicity", "t_jetmultiplicity;N_jets", 8, 0.5, 8.5);
  w_first_jet_pt =      fs->make < TH1F > ("w_first_jet_pt",    "w_first_jet_pt;P_t [GeV]", 50, 30., 700.);
  b_first_jet_pt =      fs->make < TH1F > ("b_first_jet_pt",    "b_first_jet_pt;P_t [GeV]", 50, 30., 700.);
  c_first_jet_pt =      fs->make < TH1F > ("c_first_jet_pt",    "c_first_jet_pt;P_t [GeV]", 50, 30., 700.);
  t_first_jet_pt =      fs->make < TH1F > ("t_first_jet_pt",    "t_first_jet_pt;P_t [GeV]", 50, 30., 700.);
  w_first_jet_eta =     fs->make < TH1F > ("w_first_jet_eta",   "w_first_jet_eta;Eta", 16, -2.5, 2.5);
  b_first_jet_eta =     fs->make < TH1F > ("b_first_jet_eta",   "b_first_jet_eta;Eta", 16, -2.5, 2.5);
  c_first_jet_eta =     fs->make < TH1F > ("c_first_jet_eta",   "c_first_jet_eta;Eta", 16, -2.5, 2.5);
  t_first_jet_eta =     fs->make < TH1F > ("t_first_jet_eta",   "t_first_jet_eta;Eta", 16, -2.5, 2.5);
  w_second_jet_pt =     fs->make < TH1F > ("w_second_jet_pt",   "w_second_jet_pt;P_t [GeV]", 50, 30., 500.);
  b_second_jet_pt =     fs->make < TH1F > ("b_second_jet_pt",   "b_second_jet_pt;P_t [GeV]", 50, 30., 500.);
  c_second_jet_pt =     fs->make < TH1F > ("c_second_jet_pt",   "c_second_jet_pt;P_t [GeV]", 50, 30., 500.);
  t_second_jet_pt =     fs->make < TH1F > ("t_second_jet_pt",   "t_second_jet_pt;P_t [GeV]", 50, 30., 500.);
  w_second_jet_eta =    fs->make < TH1F > ("w_second_jet_eta",  "w_second_jet_eta;Eta", 16, -2.5, 2.5);
  b_second_jet_eta =    fs->make < TH1F > ("b_second_jet_eta",  "b_second_jet_eta;Eta", 16, -2.5, 2.5);
  c_second_jet_eta =    fs->make < TH1F > ("c_second_jet_eta",  "c_second_jet_eta;Eta", 16, -2.5, 2.5);
  t_second_jet_eta =    fs->make < TH1F > ("t_second_jet_eta",  "t_second_jet_eta;Eta", 16, -2.5, 2.5);
  w_third_jet_pt =      fs->make < TH1F > ("w_third_jet_pt",    "w_third_jet_pt;P_t [GeV]", 50, 30., 200.);
  b_third_jet_pt =      fs->make < TH1F > ("b_third_jet_pt",    "b_third_jet_pt;P_t [GeV]", 50, 30., 200.);
  c_third_jet_pt =      fs->make < TH1F > ("c_third_jet_pt",    "c_third_jet_pt;P_t [GeV]", 50, 30., 200.);
  t_third_jet_pt =      fs->make < TH1F > ("t_third_jet_pt",    "t_third_jet_pt;P_t [GeV]", 50, 30., 200.);
  w_third_jet_eta =     fs->make < TH1F > ("w_third_jet_eta",   "w_third_jet_eta;Eta", 16, -2.5, 2.5);
  b_third_jet_eta =     fs->make < TH1F > ("b_third_jet_eta",   "b_third_jet_eta;Eta", 16, -2.5, 2.5);
  c_third_jet_eta =     fs->make < TH1F > ("c_third_jet_eta",   "c_third_jet_eta;Eta", 16, -2.5, 2.5);
  t_third_jet_eta =     fs->make < TH1F > ("t_third_jet_eta",   "t_third_jet_eta;Eta", 16, -2.5, 2.5);

  w_first_jet_pt_b =    fs->make < TH1F > ("w_first_jet_pt_b",   "w_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  b_first_jet_pt_b =    fs->make < TH1F > ("b_first_jet_pt_b",   "b_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  c_first_jet_pt_b =    fs->make < TH1F > ("c_first_jet_pt_b",   "c_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  t_first_jet_pt_b =    fs->make < TH1F > ("t_first_jet_pt_b",   "t_first_jet_pt_b;P_t [GeV]", 50, 30., 700.);
  w_first_jet_eta_b =   fs->make < TH1F > ("w_first_jet_eta_b",  "w_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  b_first_jet_eta_b =   fs->make < TH1F > ("b_first_jet_eta_b",  "b_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  c_first_jet_eta_b =   fs->make < TH1F > ("c_first_jet_eta_b",  "c_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  t_first_jet_eta_b =   fs->make < TH1F > ("t_first_jet_eta_b",  "t_first_jet_eta_b;Eta", 16, -2.5, 2.5);
  w_second_jet_pt_b =   fs->make < TH1F > ("w_second_jet_pt_b",  "w_second_jet_pt_b;P_t [GeV]", 50, 30., 500.);
  b_second_jet_pt_b =   fs->make < TH1F > ("b_second_jet_pt_b",  "b_second_jet_pt_b;P_t [GeV]", 50, 30., 500.);
  c_second_jet_pt_b =   fs->make < TH1F > ("c_second_jet_pt_b",  "c_second_jet_pt_b;P_t [GeV]", 50, 30., 500.);
  t_second_jet_pt_b =   fs->make < TH1F > ("t_second_jet_pt_b",  "t_second_jet_pt_b;P_t [GeV]", 50, 30., 500.);
  w_second_jet_eta_b =  fs->make < TH1F > ("w_second_jet_eta_b", "w_second_jet_eta_b;Eta", 16, -2.5, 2.5);
  b_second_jet_eta_b =  fs->make < TH1F > ("b_second_jet_eta_b", "b_second_jet_eta_b;Eta", 16, -2.5, 2.5);
  c_second_jet_eta_b =  fs->make < TH1F > ("c_second_jet_eta_b", "c_second_jet_eta_b;Eta", 16, -2.5, 2.5);
  t_second_jet_eta_b =  fs->make < TH1F > ("t_second_jet_eta_b", "t_second_jet_eta_b;Eta", 16, -2.5, 2.5);
  w_third_jet_pt_b =    fs->make < TH1F > ("w_third_jet_pt_b",   "w_third_jet_pt_b;P_t [GeV]", 50, 30., 200.);
  b_third_jet_pt_b =    fs->make < TH1F > ("b_third_jet_pt_b",   "b_third_jet_pt_b;P_t [GeV]", 50, 30., 200.);
  c_third_jet_pt_b =    fs->make < TH1F > ("c_third_jet_pt_b",   "c_third_jet_pt_b;P_t [GeV]", 50, 30., 200.);
  t_third_jet_pt_b =    fs->make < TH1F > ("t_third_jet_pt_b",   "t_third_jet_pt_b;P_t [GeV]", 50, 30., 200.);
  w_third_jet_eta_b =   fs->make < TH1F > ("w_third_jet_eta_b",  "w_third_jet_eta_b;Eta", 16, -2.5, 2.5);
  b_third_jet_eta_b =   fs->make < TH1F > ("b_third_jet_eta_b",  "b_third_jet_eta_b;Eta", 16, -2.5, 2.5);
  c_third_jet_eta_b =   fs->make < TH1F > ("c_third_jet_eta_b",  "c_third_jet_eta_b;Eta", 16, -2.5, 2.5);
  t_third_jet_eta_b =   fs->make < TH1F > ("t_third_jet_eta_b",  "t_third_jet_eta_b;Eta", 16, -2.5, 2.5);

  w_bjetmultiplicity =  fs->make < TH1F > ("w_bjetmultiplicity", "w_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  b_bjetmultiplicity =  fs->make < TH1F > ("b_bjetmultiplicity", "b_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  c_bjetmultiplicity =  fs->make < TH1F > ("c_bjetmultiplicity", "c_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  t_bjetmultiplicity =  fs->make < TH1F > ("t_bjetmultiplicity", "t_bjetmultiplicity;N_bjets", 5, 0.5, 5.5);
  w_first_bjet_pt =     fs->make < TH1F > ("w_first_bjet_pt",    "w_first_bjet_pt;P_t [GeV]", 50, 30., 700.);
  b_first_bjet_pt =     fs->make < TH1F > ("b_first_bjet_pt",    "b_first_bjet_pt;P_t [GeV]", 50, 30., 700.);
  c_first_bjet_pt =     fs->make < TH1F > ("c_first_bjet_pt",    "c_first_bjet_pt;P_t [GeV]", 50, 30., 700.);
  t_first_bjet_pt =     fs->make < TH1F > ("t_first_bjet_pt",    "t_first_bjet_pt;P_t [GeV]", 50, 30., 700.);
  w_first_bjet_eta =    fs->make < TH1F > ("w_first_bjet_eta",   "w_first_bjet_eta;Eta", 16, -2.5, 2.5);
  b_first_bjet_eta =    fs->make < TH1F > ("b_first_bjet_eta",   "b_first_bjet_eta;Eta", 16, -2.5, 2.5);
  c_first_bjet_eta =    fs->make < TH1F > ("c_first_bjet_eta",   "c_first_bjet_eta;Eta", 16, -2.5, 2.5);
  t_first_bjet_eta =    fs->make < TH1F > ("t_first_bjet_eta",   "t_first_bjet_eta;Eta", 16, -2.5, 2.5);
  w_single_bjet_pt =    fs->make < TH1F > ("w_single_bjet_pt",    "w_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  b_single_bjet_pt =    fs->make < TH1F > ("b_single_bjet_pt",    "b_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  c_single_bjet_pt =    fs->make < TH1F > ("c_single_bjet_pt",    "c_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  t_single_bjet_pt =    fs->make < TH1F > ("t_single_bjet_pt",    "t_single_bjet_pt;P_t [GeV]", 50, 30., 700.);
  w_single_bjet_eta =   fs->make < TH1F > ("w_single_bjet_eta",   "w_single_bjet_eta;Eta", 16, -2.5, 2.5);
  b_single_bjet_eta =   fs->make < TH1F > ("b_single_bjet_eta",   "b_single_bjet_eta;Eta", 16, -2.5, 2.5);
  c_single_bjet_eta =   fs->make < TH1F > ("c_single_bjet_eta",   "c_single_bjet_eta;Eta", 16, -2.5, 2.5);
  t_single_bjet_eta =   fs->make < TH1F > ("t_single_bjet_eta",   "t_single_bjet_eta;Eta", 16, -2.5, 2.5);
  w_second_bjet_pt =    fs->make < TH1F > ("w_second_bjet_pt",   "w_second_bjet_pt;P_t [GeV]", 50, 30., 500.);
  b_second_bjet_pt =    fs->make < TH1F > ("b_second_bjet_pt",   "b_second_bjet_pt;P_t [GeV]", 50, 30., 500.);
  c_second_bjet_pt =    fs->make < TH1F > ("c_second_bjet_pt",   "c_second_bjet_pt;P_t [GeV]", 50, 30., 500.);
  t_second_bjet_pt =    fs->make < TH1F > ("t_second_bjet_pt",   "t_second_bjet_pt;P_t [GeV]", 50, 30., 500.);
  w_second_bjet_eta =   fs->make < TH1F > ("w_second_bjet_eta",  "w_second_bjet_eta;Eta", 16, -2.5, 2.5);
  b_second_bjet_eta =   fs->make < TH1F > ("b_second_bjet_eta",  "b_second_bjet_eta;Eta", 16, -2.5, 2.5);
  c_second_bjet_eta =   fs->make < TH1F > ("c_second_bjet_eta",  "c_second_bjet_eta;Eta", 16, -2.5, 2.5);
  t_second_bjet_eta =   fs->make < TH1F > ("t_second_bjet_eta",  "t_second_bjet_eta;Eta", 16, -2.5, 2.5);
  w_third_bjet_pt =     fs->make < TH1F > ("w_third_bjet_pt",    "w_third_bjet_pt;P_t [GeV]", 50, 30., 200.);
  b_third_bjet_pt =     fs->make < TH1F > ("b_third_bjet_pt",    "b_third_bjet_pt;P_t [GeV]", 50, 30., 200.);
  c_third_bjet_pt =     fs->make < TH1F > ("c_third_bjet_pt",    "c_third_bjet_pt;P_t [GeV]", 50, 30., 200.);
  t_third_bjet_pt =     fs->make < TH1F > ("t_third_bjet_pt",    "t_third_bjet_pt;P_t [GeV]", 50, 30., 200.);
  w_third_bjet_eta =    fs->make < TH1F > ("w_third_bjet_eta",   "w_third_bjet_eta;Eta", 16, -2.5, 2.5);
  b_third_bjet_eta =    fs->make < TH1F > ("b_third_bjet_eta",   "b_third_bjet_eta;Eta", 16, -2.5, 2.5);
  c_third_bjet_eta =    fs->make < TH1F > ("c_third_bjet_eta",   "c_third_bjet_eta;Eta", 16, -2.5, 2.5);
  t_third_bjet_eta =    fs->make < TH1F > ("t_third_bjet_eta",   "t_third_bjet_eta;Eta", 16, -2.5, 2.5);

  w_first_ele_pt =      fs->make < TH1F > ("w_first_ele_pt",     "w_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_first_ele_pt_b =    fs->make < TH1F > ("w_first_ele_pt_b",   "w_first_ele_pt_b;P_t [GeV]", 50, 0., 450.);
  w_first_ele_pt_bb =   fs->make < TH1F > ("w_first_ele_pt_bb",  "w_first_ele_pt_bb;P_t [GeV]", 50, 0., 450.);
  b_first_ele_pt =      fs->make < TH1F > ("b_first_ele_pt",     "b_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  c_first_ele_pt =      fs->make < TH1F > ("c_first_ele_pt",     "c_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  t_first_ele_pt =      fs->make < TH1F > ("t_first_ele_pt",     "t_first_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_second_ele_pt =     fs->make < TH1F > ("w_second_ele_pt",    "w_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  b_second_ele_pt =     fs->make < TH1F > ("b_second_ele_pt",    "b_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  c_second_ele_pt =     fs->make < TH1F > ("c_second_ele_pt",    "c_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  t_second_ele_pt =     fs->make < TH1F > ("t_second_ele_pt",    "t_second_ele_pt;P_t [GeV]", 50, 0., 450.);
  w_first_muon_pt =     fs->make < TH1F > ("w_first_muon_pt",    "w_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  w_first_muon_pt_b =   fs->make < TH1F > ("w_first_muon_pt_b",  "w_first_muon_pt_b [GeV]", 50, 0., 450.);
  w_first_muon_pt_bb =  fs->make < TH1F > ("w_first_muon_pt_bb",  "w_first_muon_pt_bb [GeV]", 50, 0., 450.);
  b_first_muon_pt =     fs->make < TH1F > ("b_first_muon_pt",    "b_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  c_first_muon_pt =     fs->make < TH1F > ("c_first_muon_pt",    "c_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  t_first_muon_pt =     fs->make < TH1F > ("t_first_muon_pt",    "t_first_muon_pt;P_t [GeV]", 50, 0., 450.);
  w_second_muon_pt =    fs->make < TH1F > ("w_second_muon_pt",   "w_second_muon_pt;P_t [GeV]", 50, 0., 450.);
  b_second_muon_pt =    fs->make < TH1F > ("b_second_muon_pt",   "b_second_muon_pt;P_t [GeV]", 50, 0., 450.);
  c_second_muon_pt =    fs->make < TH1F > ("c_second_muon_pt",   "c_second_muon_pt;P_t [GeV]", 50, 0., 450.);
  t_second_muon_pt =    fs->make < TH1F > ("t_second_muon_pt",   "t_second_muon_pt;P_t [GeV]", 50, 0., 450.);
  w_first_ele_eta =     fs->make < TH1F > ("w_first_ele_eta",    "w_first_ele_eta;Eta", 16, -2.5, 2.5);
  b_first_ele_eta =     fs->make < TH1F > ("b_first_ele_eta",    "b_first_ele_eta;Eta", 16, -2.5, 2.5);
  c_first_ele_eta =     fs->make < TH1F > ("c_first_ele_eta",    "c_first_ele_eta;Eta", 16, -2.5, 2.5);
  t_first_ele_eta =     fs->make < TH1F > ("t_first_ele_eta",    "t_first_ele_eta;Eta", 16, -2.5, 2.5);
  w_second_ele_eta =    fs->make < TH1F > ("w_second_ele_eta",   "w_second_ele_eta;Eta", 16, -2.5, 2.5);
  b_second_ele_eta =    fs->make < TH1F > ("b_second_ele_eta",   "b_second_ele_eta;Eta", 16, -2.5, 2.5);
  c_second_ele_eta =    fs->make < TH1F > ("c_second_ele_eta",   "c_second_ele_eta;Eta", 16, -2.5, 2.5);
  t_second_ele_eta =    fs->make < TH1F > ("t_second_ele_eta",   "t_second_ele_eta;Eta", 16, -2.5, 2.5);
  w_first_ele_iso =     fs->make < TH1F > ("w_first_ele_iso",    "w_first_ele_iso;Iso", 50, 0., 1.0);
  b_first_ele_iso =     fs->make < TH1F > ("b_first_ele_iso",    "b_first_ele_iso;Iso", 50, 0., 1.0);
  c_first_ele_iso =     fs->make < TH1F > ("c_first_ele_iso",    "c_first_ele_iso;Iso", 50, 0., 1.0);
  t_first_ele_iso =     fs->make < TH1F > ("t_first_ele_iso",    "t_first_ele_iso;Iso", 50, 0., 1.0);
  w_first_muon_eta =    fs->make < TH1F > ("w_first_muon_eta",   "w_first_muon_eta;Eta", 16, -2.5, 2.5);
  b_first_muon_eta =    fs->make < TH1F > ("b_first_muon_eta",   "b_first_muon_eta;Eta", 16, -2.5, 2.5);
  c_first_muon_eta =    fs->make < TH1F > ("c_first_muon_eta",   "c_first_muon_eta;Eta", 16, -2.5, 2.5);
  t_first_muon_eta =    fs->make < TH1F > ("t_first_muon_eta",   "t_first_muon_eta;Eta", 16, -2.5, 2.5);
  w_second_muon_eta =   fs->make < TH1F > ("w_second_muon_eta",  "w_second_muon_eta;Eta", 16, -2.5, 2.5);
  b_second_muon_eta =   fs->make < TH1F > ("b_second_muon_eta",  "b_second_muon_eta;Eta", 16, -2.5, 2.5);
  c_second_muon_eta =   fs->make < TH1F > ("c_second_muon_eta",  "c_second_muon_eta;Eta", 16, -2.5, 2.5);
  t_second_muon_eta =   fs->make < TH1F > ("t_second_muon_eta",  "t_second_muon_eta;Eta", 16, -2.5, 2.5);
  w_first_muon_iso =    fs->make < TH1F > ("w_first_muon_iso",   "w_first_muon_iso;Iso", 50, 0., 1.0);
  b_first_muon_iso =    fs->make < TH1F > ("b_first_muon_iso",   "b_first_muon_iso;Iso", 50, 0., 1.0);
  c_first_muon_iso =    fs->make < TH1F > ("c_first_muon_iso",   "c_first_muon_iso;Iso", 50, 0., 1.0);
  t_first_muon_iso =    fs->make < TH1F > ("t_first_muon_iso",   "t_first_muon_iso;Iso", 50, 0., 1.0);

  w_mass_ee_wide =      fs->make < TH1F > ("w_mass_ee_wide",    "w_mass_ee_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_ee_wide =      fs->make < TH1F > ("b_mass_ee_wide",    "b_mass_ee_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_ee_wide =      fs->make < TH1F > ("c_mass_ee_wide",    "c_mass_ee_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_ee_wide =      fs->make < TH1F > ("t_mass_ee_wide",    "t_mass_ee_wide;Mass [GeV]", 40, 50., 200.);
  w_mass_mm_wide =      fs->make < TH1F > ("w_mass_mm_wide",    "w_mass_mm_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_mm_wide =      fs->make < TH1F > ("b_mass_mm_wide",    "b_mass_mm_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_mm_wide =      fs->make < TH1F > ("c_mass_mm_wide",    "c_mass_mm_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_mm_wide =      fs->make < TH1F > ("t_mass_mm_wide",    "t_mass_mm_wide;Mass [GeV]", 40, 50., 200.);

  //  w_mass_em_wide =      fs->make < TH1F > ("w_mass_em_wide",    "w_mass_em_wide;Mass [GeV]", 40, 0., 160.);
  //  b_mass_em_wide =      fs->make < TH1F > ("b_mass_em_wide",    "b_mass_em_wide;Mass [GeV]", 40, 0., 160.);
  //  c_mass_em_wide =      fs->make < TH1F > ("c_mass_em_wide",    "c_mass_em_wide;Mass [GeV]", 40, 0., 160.);
  //  t_mass_em_wide =      fs->make < TH1F > ("t_mass_em_wide",    "t_mass_em_wide;Mass [GeV]", 40, 0., 160.);

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

  //  h_mass_em =           fs->make < TH1F > ("h_mass_em",         "h_mass_em;Mass [GeV]", 40, 0., 160.);
  //  w_mass_em =           fs->make < TH1F > ("w_mass_em",         "w_mass_em;Mass [GeV]", 40, 0., 160.);
  //  b_mass_em =           fs->make < TH1F > ("b_mass_em",         "b_mass_em;Mass [GeV]", 40, 0., 160.);
  //  c_mass_em =           fs->make < TH1F > ("c_mass_em",         "c_mass_em;Mass [GeV]", 40, 0., 160.);
  //  t_mass_em =           fs->make < TH1F > ("t_mass_em",         "t_mass_em;Mass [GeV]", 40, 0., 160.);

  w_mt_wenu_wide =      fs->make < TH1F > ("w_mt_wenu_wide",    "w_mt_wenu_wide;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wenu_wide =      fs->make < TH1F > ("b_mt_wenu_wide",    "b_mt_wenu_wide;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wenu_wide =      fs->make < TH1F > ("c_mt_wenu_wide",    "c_mt_wenu_wide;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wenu_wide =      fs->make < TH1F > ("t_mt_wenu_wide",    "t_mt_wenu_wide;M_{t} [GeV]", 40, 0., 160.);
  w_mt_wmnu_wide =      fs->make < TH1F > ("w_mt_wmnu_wide",    "w_mt_wmnu_wide;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wmnu_wide =      fs->make < TH1F > ("b_mt_wmnu_wide",    "b_mt_wmnu_wide;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wmnu_wide =      fs->make < TH1F > ("c_mt_wmnu_wide",    "c_mt_wmnu_wide;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wmnu_wide =      fs->make < TH1F > ("t_mt_wmnu_wide",    "t_mt_wmnu_wide;M_{t} [GeV]", 40, 0., 160.);

  h_mt_wenu =           fs->make < TH1F > ("h_mt_wenu",         "h_mt_wenu;M_{t} [GeV]", 40, 0., 160.);
  w_mt_wenu =           fs->make < TH1F > ("w_mt_wenu",         "w_mt_wenu;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wenu =           fs->make < TH1F > ("b_mt_wenu",         "b_mt_wenu;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wenu =           fs->make < TH1F > ("c_mt_wenu",         "c_mt_wenu;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wenu =           fs->make < TH1F > ("t_mt_wenu",         "t_mt_wenu;M_{t} [GeV]", 40, 0., 160.);
  h_mt_wmnu =           fs->make < TH1F > ("h_mt_wmnu",         "h_mt_wmnu;M_{t} [GeV]", 40, 0., 160.);
  w_mt_wmnu =           fs->make < TH1F > ("w_mt_wmnu",         "w_mt_wmnu;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wmnu =           fs->make < TH1F > ("b_mt_wmnu",         "b_mt_wmnu;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wmnu =           fs->make < TH1F > ("c_mt_wmnu",         "c_mt_wmnu;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wmnu =           fs->make < TH1F > ("t_mt_wmnu",         "t_mt_wmnu;M_{t} [GeV]", 40, 0., 160.);

  w_mass_Zj_ee =        fs->make < TH1F > ("w_mass_Zj_ee",      "w_mass_Zj_ee", 60, 100., 330.);
  b_mass_Zj_ee =        fs->make < TH1F > ("b_mass_Zj_ee",      "b_mass_Zj_ee", 60, 100., 330.);
  c_mass_Zj_ee =        fs->make < TH1F > ("c_mass_Zj_ee",      "c_mass_Zj_ee", 60, 100., 330.);
  t_mass_Zj_ee =        fs->make < TH1F > ("t_mass_Zj_ee",      "t_mass_Zj_ee", 60, 100., 330.);
  
  w_mass_Zj_ee_b =      fs->make < TH1F > ("w_mass_Zj_ee_b",    "w_mass_Zj_ee_b", 60, 100., 330.);
  b_mass_Zj_ee_b =      fs->make < TH1F > ("b_mass_Zj_ee_b",    "b_mass_Zj_ee_b", 60, 100., 330.);
  c_mass_Zj_ee_b =      fs->make < TH1F > ("c_mass_Zj_ee_b",    "c_mass_Zj_ee_b", 60, 100., 330.);
  t_mass_Zj_ee_b =      fs->make < TH1F > ("t_mass_Zj_ee_b",    "t_mass_Zj_ee_b", 60, 100., 330.);

  w_mass_Zj_ee_bb =      fs->make < TH1F > ("w_mass_Zj_ee_bb",    "w_mass_Zj_ee_bb", 60, 100., 330.);
  b_mass_Zj_ee_bb =      fs->make < TH1F > ("b_mass_Zj_ee_bb",    "b_mass_Zj_ee_bb", 60, 100., 330.);
  c_mass_Zj_ee_bb =      fs->make < TH1F > ("c_mass_Zj_ee_bb",    "c_mass_Zj_ee_bb", 60, 100., 330.);
  t_mass_Zj_ee_bb =      fs->make < TH1F > ("t_mass_Zj_ee_bb",    "t_mass_Zj_ee_bb", 60, 100., 330.);
  
  w_mass_Zj_mm =        fs->make < TH1F > ("w_mass_Zj_mm",      "w_mass_Zj_mm", 60, 100., 330.);
  b_mass_Zj_mm =        fs->make < TH1F > ("b_mass_Zj_mm",      "b_mass_Zj_mm", 60, 100., 330.);
  c_mass_Zj_mm =        fs->make < TH1F > ("c_mass_Zj_mm",      "c_mass_Zj_mm", 60, 100., 330.);
  t_mass_Zj_mm =        fs->make < TH1F > ("t_mass_Zj_mm",      "t_mass_Zj_mm", 60, 100., 330.);
  
  w_mass_Zj_mm_b =      fs->make < TH1F > ("w_mass_Zj_mm_b",    "w_mass_Zj_mm_b", 60, 100., 330.);
  b_mass_Zj_mm_b =      fs->make < TH1F > ("b_mass_Zj_mm_b",    "b_mass_Zj_mm_b", 60, 100., 330.);
  c_mass_Zj_mm_b =      fs->make < TH1F > ("c_mass_Zj_mm_b",    "c_mass_Zj_mm_b", 60, 100., 330.);
  t_mass_Zj_mm_b =      fs->make < TH1F > ("t_mass_Zj_mm_b",    "t_mass_Zj_mm_b", 60, 100., 330.);

  w_mass_Zj_mm_bb =      fs->make < TH1F > ("w_mass_Zj_mm_bb",    "w_mass_Zj_mm_bb", 60, 100., 330.);
  b_mass_Zj_mm_bb =      fs->make < TH1F > ("b_mass_Zj_mm_bb",    "b_mass_Zj_mm_bb", 60, 100., 330.);
  c_mass_Zj_mm_bb =      fs->make < TH1F > ("c_mass_Zj_mm_bb",    "c_mass_Zj_mm_bb", 60, 100., 330.);
  t_mass_Zj_mm_bb =      fs->make < TH1F > ("t_mass_Zj_mm_bb",    "t_mass_Zj_mm_bb", 60, 100., 330.);

  //  w_mass_Zj_em =        fs->make < TH1F > ("w_mass_Zj_em",      "w_mass_Zj_em", 60, 100., 330.);
  //  b_mass_Zj_em =        fs->make < TH1F > ("b_mass_Zj_em",      "b_mass_Zj_em", 60, 100., 330.);
  //  c_mass_Zj_em =        fs->make < TH1F > ("c_mass_Zj_em",      "c_mass_Zj_em", 60, 100., 330.);
  //  t_mass_Zj_em =        fs->make < TH1F > ("t_mass_Zj_em",      "t_mass_Zj_em", 60, 100., 330.);

  //  w_mass_Zj_em_b =      fs->make < TH1F > ("w_mass_Zj_em_b",    "w_mass_Zj_em_b", 60, 100., 330.);
  //  b_mass_Zj_em_b =      fs->make < TH1F > ("b_mass_Zj_em_b",    "b_mass_Zj_em_b", 60, 100., 330.);
  //  c_mass_Zj_em_b =      fs->make < TH1F > ("c_mass_Zj_em_b",    "c_mass_Zj_em_b", 60, 100., 330.);
  //  t_mass_Zj_em_b =      fs->make < TH1F > ("t_mass_Zj_em_b",    "t_mass_Zj_em_b", 60, 100., 330.);

  w_pt_Z_ee =           fs->make < TH1F > ("w_pt_Z_ee",         "w_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_ee =           fs->make < TH1F > ("b_pt_Z_ee",         "b_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_ee =           fs->make < TH1F > ("c_pt_Z_ee",         "c_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_ee =           fs->make < TH1F > ("t_pt_Z_ee",         "t_pt_Z_ee;P_t [GeV]", 40, 0., 400.);
  w_pt_Z_mm =           fs->make < TH1F > ("w_pt_Z_mm",         "w_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  b_pt_Z_mm =           fs->make < TH1F > ("b_pt_Z_mm",         "b_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  c_pt_Z_mm =           fs->make < TH1F > ("c_pt_Z_mm",         "c_pt_Z_mm;P_t [GeV]", 40, 0., 400.);
  t_pt_Z_mm =           fs->make < TH1F > ("t_pt_Z_mm",         "t_pt_Z_mm;P_t [GeV]", 40, 0., 400.);

  //  w_pt_Z_em =           fs->make < TH1F > ("w_pt_Z_em",         "w_pt_Z_em;P_t [GeV]", 40, 0., 400.);
  //  b_pt_Z_em =           fs->make < TH1F > ("b_pt_Z_em",         "b_pt_Z_em;P_t [GeV]", 40, 0., 400.);
  //  c_pt_Z_em =           fs->make < TH1F > ("c_pt_Z_em",         "c_pt_Z_em;P_t [GeV]", 40, 0., 400.);
  //  t_pt_Z_em =           fs->make < TH1F > ("t_pt_Z_em",         "t_pt_Z_em;P_t [GeV]", 40, 0., 400.);

  w_mass_ee_b_wide =    fs->make < TH1F > ("w_mass_ee_b_wide",  "w_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_ee_b_wide =    fs->make < TH1F > ("b_mass_ee_b_wide",  "b_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_ee_b_wide =    fs->make < TH1F > ("c_mass_ee_b_wide",  "c_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_ee_b_wide =    fs->make < TH1F > ("t_mass_ee_b_wide",  "t_mass_ee_b_wide;Mass [GeV]", 40, 50., 200.);
  w_mass_mm_b_wide =    fs->make < TH1F > ("w_mass_mm_b_wide",  "w_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);
  b_mass_mm_b_wide =    fs->make < TH1F > ("b_mass_mm_b_wide",  "b_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);
  c_mass_mm_b_wide =    fs->make < TH1F > ("c_mass_mm_b_wide",  "c_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);
  t_mass_mm_b_wide =    fs->make < TH1F > ("t_mass_mm_b_wide",  "t_mass_mm_b_wide;Mass [GeV]", 40, 50., 200.);

  //  w_mass_em_b_wide =    fs->make < TH1F > ("w_mass_em_b_wide",  "w_mass_em_b_wide;Mass [GeV]", 40, 0., 160.);
  //  b_mass_em_b_wide =    fs->make < TH1F > ("b_mass_em_b_wide",  "b_mass_em_b_wide;Mass [GeV]", 40, 0., 160.);
  //  c_mass_em_b_wide =    fs->make < TH1F > ("c_mass_em_b_wide",  "c_mass_em_b_wide;Mass [GeV]", 40, 0., 160.);
  //  t_mass_em_b_wide =    fs->make < TH1F > ("t_mass_em_b_wide",  "t_mass_em_b_wide;Mass [GeV]", 40, 0., 160.);

  w_mt_wenu_b_wide =    fs->make < TH1F > ("w_mt_wenu_b_wide",  "w_mt_wenu_b_wide;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wenu_b_wide =    fs->make < TH1F > ("b_mt_wenu_b_wide",  "b_mt_wenu_b_wide;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wenu_b_wide =    fs->make < TH1F > ("c_mt_wenu_b_wide",  "c_mt_wenu_b_wide;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wenu_b_wide =    fs->make < TH1F > ("t_mt_wenu_b_wide",  "t_mt_wenu_b_wide;M_{t} [GeV]", 40, 0., 160.);
  w_mt_wmnu_b_wide =    fs->make < TH1F > ("w_mt_wmnu_b_wide",  "w_mt_wmnu_b_wide;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wmnu_b_wide =    fs->make < TH1F > ("b_mt_wmnu_b_wide",  "b_mt_wmnu_b_wide;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wmnu_b_wide =    fs->make < TH1F > ("c_mt_wmnu_b_wide",  "c_mt_wmnu_b_wide;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wmnu_b_wide =    fs->make < TH1F > ("t_mt_wmnu_b_wide",  "t_mt_wmnu_b_wide;M_{t} [GeV]", 40, 0., 160.);

  w_mt_wenu_bb_wide =   fs->make < TH1F > ("w_mt_wenu_bb_wide", "w_mt_wenu_bb_wide;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wenu_bb_wide =   fs->make < TH1F > ("b_mt_wenu_bb_wide", "b_mt_wenu_bb_wide;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wenu_bb_wide =   fs->make < TH1F > ("c_mt_wenu_bb_wide", "c_mt_wenu_bb_wide;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wenu_bb_wide =   fs->make < TH1F > ("t_mt_wenu_bb_wide", "t_mt_wenu_bb_wide;M_{t} [GeV]", 40, 0., 160.);
  w_mt_wmnu_bb_wide =   fs->make < TH1F > ("w_mt_wmnu_bb_wide", "w_mt_wmnu_bb_wide;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wmnu_bb_wide =   fs->make < TH1F > ("b_mt_wmnu_bb_wide", "b_mt_wmnu_bb_wide;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wmnu_bb_wide =   fs->make < TH1F > ("c_mt_wmnu_bb_wide", "c_mt_wmnu_bb_wide;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wmnu_bb_wide =   fs->make < TH1F > ("t_mt_wmnu_bb_wide", "t_mt_wmnu_bb_wide;M_{t} [GeV]", 40, 0., 160.);

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

  //  w_mass_em_b =         fs->make < TH1F > ("w_mass_em_b",       "w_mass_em_b;Mass [GeV]", 40, 0., 160.);
  //  b_mass_em_b =         fs->make < TH1F > ("b_mass_em_b",       "b_mass_em_b;Mass [GeV]", 40, 0., 160.);
  //  c_mass_em_b =         fs->make < TH1F > ("c_mass_em_b",       "c_mass_em_b;Mass [GeV]", 40, 0., 160.);
  //  t_mass_em_b =         fs->make < TH1F > ("t_mass_em_b",       "t_mass_em_b;Mass [GeV]", 40, 0., 160.);

  w_mt_wenu_b =         fs->make < TH1F > ("w_mt_wenu_b",       "w_mt_wenu_b;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wenu_b =         fs->make < TH1F > ("b_mt_wenu_b",       "b_mt_wenu_b;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wenu_b =         fs->make < TH1F > ("c_mt_wenu_b",       "c_mt_wenu_b;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wenu_b =         fs->make < TH1F > ("t_mt_wenu_b",       "t_mt_wenu_b;M_{t} [GeV]", 40, 0., 160.);
  w_mt_wmnu_b =         fs->make < TH1F > ("w_mt_wmnu_b",       "w_mt_wmnu_b;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wmnu_b =         fs->make < TH1F > ("b_mt_wmnu_b",       "b_mt_wmnu_b;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wmnu_b =         fs->make < TH1F > ("c_mt_wmnu_b",       "c_mt_wmnu_b;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wmnu_b =         fs->make < TH1F > ("t_mt_wmnu_b",       "t_mt_wmnu_b;M_{t} [GeV]", 40, 0., 160.);

  w_mt_wenu_bb =        fs->make < TH1F > ("w_mt_wenu_bb",      "w_mt_wenu_bb;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wenu_bb =        fs->make < TH1F > ("b_mt_wenu_bb",      "b_mt_wenu_bb;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wenu_bb =        fs->make < TH1F > ("c_mt_wenu_bb",      "c_mt_wenu_bb;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wenu_bb =        fs->make < TH1F > ("t_mt_wenu_bb",      "t_mt_wenu_bb;M_{t} [GeV]", 40, 0., 160.);
  w_mt_wmnu_bb =        fs->make < TH1F > ("w_mt_wmnu_bb",      "w_mt_wmnu_bb;M_{t} [GeV]", 40, 0., 160.);
  b_mt_wmnu_bb =        fs->make < TH1F > ("b_mt_wmnu_bb",      "b_mt_wmnu_bb;M_{t} [GeV]", 40, 0., 160.);
  c_mt_wmnu_bb =        fs->make < TH1F > ("c_mt_wmnu_bb",      "c_mt_wmnu_bb;M_{t} [GeV]", 40, 0., 160.);
  t_mt_wmnu_bb =        fs->make < TH1F > ("t_mt_wmnu_bb",      "t_mt_wmnu_bb;M_{t} [GeV]", 40, 0., 160.);

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

  //  w_pt_Z_em_b =         fs->make < TH1F > ("w_pt_Z_em_b",       "w_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);
  //  b_pt_Z_em_b =         fs->make < TH1F > ("b_pt_Z_em_b",       "b_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);
  //  c_pt_Z_em_b =         fs->make < TH1F > ("c_pt_Z_em_b",       "c_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);
  //  t_pt_Z_em_b =         fs->make < TH1F > ("t_pt_Z_em_b",       "t_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);

  w_single_pt_Z_ee_b =  fs->make < TH1F > ("w_single_pt_Z_ee_b",       "w_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  b_single_pt_Z_ee_b =  fs->make < TH1F > ("b_single_pt_Z_ee_b",       "b_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  c_single_pt_Z_ee_b =  fs->make < TH1F > ("c_single_pt_Z_ee_b",       "c_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  t_single_pt_Z_ee_b =  fs->make < TH1F > ("t_single_pt_Z_ee_b",       "t_single_pt_Z_ee_b;P_t [GeV]", 40, 0., 400.);
  w_single_pt_Z_mm_b =  fs->make < TH1F > ("w_single_pt_Z_mm_b",       "w_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  b_single_pt_Z_mm_b =  fs->make < TH1F > ("b_single_pt_Z_mm_b",       "b_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  c_single_pt_Z_mm_b =  fs->make < TH1F > ("c_single_pt_Z_mm_b",       "c_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);
  t_single_pt_Z_mm_b =  fs->make < TH1F > ("t_single_pt_Z_mm_b",       "t_single_pt_Z_mm_b;P_t [GeV]", 40, 0., 400.);

  //  w_single_pt_Z_em_b =  fs->make < TH1F > ("w_single_pt_Z_em_b",       "w_single_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);
  //  b_single_pt_Z_em_b =  fs->make < TH1F > ("b_single_pt_Z_em_b",       "b_single_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);
  //  c_single_pt_Z_em_b =  fs->make < TH1F > ("c_single_pt_Z_em_b",       "c_single_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);
  //  t_single_pt_Z_em_b =  fs->make < TH1F > ("t_single_pt_Z_em_b",       "t_single_pt_Z_em_b;P_t [GeV]", 40, 0., 400.);

  w_delta_ee =          fs->make < TH1F > ("w_delta_phi_ee",    "w_delta_phi_ee",   12, 0, TMath::Pi ());
  b_delta_ee =          fs->make < TH1F > ("b_delta_phi_ee",    "b_delta_phi_ee",   12, 0, TMath::Pi ());
  c_delta_ee =          fs->make < TH1F > ("c_delta_phi_ee",    "c_delta_phi_ee",   12, 0, TMath::Pi ());
  t_delta_ee =          fs->make < TH1F > ("t_delta_phi_ee",    "t_delta_phi_ee",   12, 0, TMath::Pi ());
  w_delta_ee_b =        fs->make < TH1F > ("w_delta_phi_ee_b",  "w_delta_phi_ee_b", 12, 0, TMath::Pi ());  
  b_delta_ee_b =        fs->make < TH1F > ("b_delta_phi_ee_b",  "b_delta_phi_ee_b", 12, 0, TMath::Pi ());
  c_delta_ee_b =        fs->make < TH1F > ("c_delta_phi_ee_b",  "c_delta_phi_ee_b", 12, 0, TMath::Pi ());
  t_delta_ee_b =        fs->make < TH1F > ("t_delta_phi_ee_b",  "t_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_delta_ee_bb =        fs->make < TH1F > ("w_delta_phi_ee_bb",  "w_delta_phi_ee_bb", 12, 0, TMath::Pi ());  
  b_delta_ee_bb =        fs->make < TH1F > ("b_delta_phi_ee_bb",  "b_delta_phi_ee_bb", 12, 0, TMath::Pi ());
  c_delta_ee_bb =        fs->make < TH1F > ("c_delta_phi_ee_bb",  "c_delta_phi_ee_bb", 12, 0, TMath::Pi ());
  t_delta_ee_bb =        fs->make < TH1F > ("t_delta_phi_ee_bb",  "t_delta_phi_ee_bb", 12, 0, TMath::Pi ());
  w_delta_mm =          fs->make < TH1F > ("w_delta_phi_mm",    "w_delta_phi_mm",   12, 0, TMath::Pi ());
  b_delta_mm =          fs->make < TH1F > ("b_delta_phi_mm",    "b_delta_phi_mm",   12, 0, TMath::Pi ());
  c_delta_mm =          fs->make < TH1F > ("c_delta_phi_mm",    "c_delta_phi_mm",   12, 0, TMath::Pi ());
  t_delta_mm =          fs->make < TH1F > ("t_delta_phi_mm",    "t_delta_phi_mm",   12, 0, TMath::Pi ());
  w_delta_mm_b =        fs->make < TH1F > ("w_delta_phi_mm_b",  "w_delta_phi_mm_b", 12, 0, TMath::Pi ());
  b_delta_mm_b =        fs->make < TH1F > ("b_delta_phi_mm_b",  "b_delta_phi_mm_b", 12, 0, TMath::Pi ());
  c_delta_mm_b =        fs->make < TH1F > ("c_delta_phi_mm_b",  "c_delta_phi_mm_b", 12, 0, TMath::Pi ());
  t_delta_mm_b =        fs->make < TH1F > ("t_delta_phi_mm_b",  "t_delta_phi_mm_b", 12, 0, TMath::Pi ());
  w_delta_mm_bb =        fs->make < TH1F > ("w_delta_phi_mm_bb",  "w_delta_phi_mm_bb", 12, 0, TMath::Pi ());
  b_delta_mm_bb =        fs->make < TH1F > ("b_delta_phi_mm_bb",  "b_delta_phi_mm_bb", 12, 0, TMath::Pi ());
  c_delta_mm_bb =        fs->make < TH1F > ("c_delta_phi_mm_bb",  "c_delta_phi_mm_bb", 12, 0, TMath::Pi ());
  t_delta_mm_bb =        fs->make < TH1F > ("t_delta_phi_mm_bb",  "t_delta_phi_mm_bb", 12, 0, TMath::Pi ());

  //  w_delta_em =          fs->make < TH1F > ("w_delta_phi_em",    "w_delta_phi_em", 12, 0, TMath::Pi ());
  //  w_delta_em_b =        fs->make < TH1F > ("w_delta_phi_em_b",  "w_delta_phi_em_b", 12, 0, TMath::Pi ());
  //  b_delta_em_b =        fs->make < TH1F > ("b_delta_phi_em_b",  "b_delta_phi_em_b", 12, 0, TMath::Pi ());
  //  c_delta_em_b =        fs->make < TH1F > ("c_delta_phi_em_b",  "c_delta_phi_em_b", 12, 0, TMath::Pi ());
  //  t_delta_em_b =        fs->make < TH1F > ("t_delta_phi_em_b",  "t_delta_phi_em_b", 12, 0, TMath::Pi ());

  w_delta_wenuMet =          fs->make < TH1F > ("w_delta_phi_wenuMet",     "w_delta_phi_wenuMet",    12, 0, TMath::Pi ());
  b_delta_wenuMet =          fs->make < TH1F > ("b_delta_phi_wenuMet",     "b_delta_phi_wenuMet",    12, 0, TMath::Pi ());
  c_delta_wenuMet =          fs->make < TH1F > ("c_delta_phi_wenuMet",     "c_delta_phi_wenuMet",    12, 0, TMath::Pi ());
  t_delta_wenuMet =          fs->make < TH1F > ("t_delta_phi_wenuMet",     "t_delta_phi_wenuMet",    12, 0, TMath::Pi ());
  w_delta_wenuMet_b =        fs->make < TH1F > ("w_delta_phi_wenuMet_b",   "w_delta_phi_wenuMet_b",  12, 0, TMath::Pi ());  
  b_delta_wenuMet_b =        fs->make < TH1F > ("b_delta_phi_wenuMet_b",   "b_delta_phi_wenuMet_b",  12, 0, TMath::Pi ());
  c_delta_wenuMet_b =        fs->make < TH1F > ("c_delta_phi_wenuMet_b",   "c_delta_phi_wenuMet_b",  12, 0, TMath::Pi ());
  t_delta_wenuMet_b =        fs->make < TH1F > ("t_delta_phi_wenuMet_b",   "t_delta_phi_wenuMet_b",  12, 0, TMath::Pi ());
  w_delta_wenuMet_bb =       fs->make < TH1F > ("w_delta_phi_wenuMet_bb",  "w_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());  
  b_delta_wenuMet_bb =       fs->make < TH1F > ("b_delta_phi_wenuMet_bb",  "b_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());
  c_delta_wenuMet_bb =       fs->make < TH1F > ("c_delta_phi_wenuMet_bb",  "c_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());
  t_delta_wenuMet_bb =       fs->make < TH1F > ("t_delta_phi_wenuMet_bb",  "t_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());
  w_delta_wmnuMet =          fs->make < TH1F > ("w_delta_phi_wmnuMet",     "w_delta_phi_wmnuMet",    12, 0, TMath::Pi ());
  b_delta_wmnuMet =          fs->make < TH1F > ("b_delta_phi_wmnuMet",     "b_delta_phi_wmnuMet",    12, 0, TMath::Pi ());
  c_delta_wmnuMet =          fs->make < TH1F > ("c_delta_phi_wmnuMet",     "c_delta_phi_wmnuMet",    12, 0, TMath::Pi ());
  t_delta_wmnuMet =          fs->make < TH1F > ("t_delta_phi_wmnuMet",     "t_delta_phi_wmnuMet",    12, 0, TMath::Pi ());
  w_delta_wmnuMet_b =        fs->make < TH1F > ("w_delta_phi_wmnuMet_b",   "w_delta_phi_wmnuMet_b",  12, 0, TMath::Pi ());
  b_delta_wmnuMet_b =        fs->make < TH1F > ("b_delta_phi_wmnuMet_b",   "b_delta_phi_wmnuMet_b",  12, 0, TMath::Pi ());
  c_delta_wmnuMet_b =        fs->make < TH1F > ("c_delta_phi_wmnuMet_b",   "c_delta_phi_wmnuMet_b",  12, 0, TMath::Pi ());
  t_delta_wmnuMet_b =        fs->make < TH1F > ("t_delta_phi_wmnuMet_b",   "t_delta_phi_wmnuMet_b",  12, 0, TMath::Pi ());
  w_delta_wmnuMet_bb =       fs->make < TH1F > ("w_delta_phi_wmnuMet_bb",  "w_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());
  b_delta_wmnuMet_bb =       fs->make < TH1F > ("b_delta_phi_wmnuMet_bb",  "b_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());
  c_delta_wmnuMet_bb =       fs->make < TH1F > ("c_delta_phi_wmnuMet_bb",  "c_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());
  t_delta_wmnuMet_bb =       fs->make < TH1F > ("t_delta_phi_wmnuMet_bb",  "t_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());

  w_single_delta_ee_b =      fs->make < TH1F > ("w_single_delta_phi_ee_b", "w_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  b_single_delta_ee_b =      fs->make < TH1F > ("b_single_delta_phi_ee_b", "b_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  c_single_delta_ee_b =      fs->make < TH1F > ("c_single_delta_phi_ee_b", "c_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  t_single_delta_ee_b =      fs->make < TH1F > ("t_single_delta_phi_ee_b", "t_single_delta_phi_ee_b", 12, 0, TMath::Pi ());
  w_single_delta_mm_b =      fs->make < TH1F > ("w_single_delta_phi_mm_b", "w_single_delta_phi_mm_b", 12, 0, TMath::Pi ());
  b_single_delta_mm_b =      fs->make < TH1F > ("b_single_delta_phi_mm_b", "b_single_delta_phi_mm_b", 12, 0, TMath::Pi ());
  c_single_delta_mm_b =      fs->make < TH1F > ("c_single_delta_phi_mm_b", "c_single_delta_phi_mm_b", 12, 0, TMath::Pi ());
  t_single_delta_mm_b =      fs->make < TH1F > ("t_single_delta_phi_mm_b", "t_single_delta_phi_mm_b", 12, 0, TMath::Pi ());

  //  w_single_delta_em_b =        fs->make < TH1F > ("w_single_delta_phi_em_b",  "w_single_delta_phi_em_b", 12, 0, TMath::Pi ());
  //  b_single_delta_em_b =        fs->make < TH1F > ("b_single_delta_phi_em_b",  "b_single_delta_phi_em_b", 12, 0, TMath::Pi ());
  //  c_single_delta_em_b =        fs->make < TH1F > ("c_single_delta_phi_em_b",  "c_single_delta_phi_em_b", 12, 0, TMath::Pi ());
  //  t_single_delta_em_b =        fs->make < TH1F > ("t_single_delta_phi_em_b",  "t_single_delta_phi_em_b", 12, 0, TMath::Pi ());

  w_single_delta_wenuMet_b =        fs->make < TH1F > ("w_single_delta_phi_wenuMet_b",  "w_single_delta_phi_wenuMet_b", 12, 0, TMath::Pi ());
  b_single_delta_wenuMet_b =        fs->make < TH1F > ("b_single_delta_phi_wenuMet_b",  "b_single_delta_phi_wenuMet_b", 12, 0, TMath::Pi ());
  c_single_delta_wenuMet_b =        fs->make < TH1F > ("c_single_delta_phi_wenuMet_b",  "c_single_delta_phi_wenuMet_b", 12, 0, TMath::Pi ());
  t_single_delta_wenuMet_b =        fs->make < TH1F > ("t_single_delta_phi_wenuMet_b",  "t_single_delta_phi_wenuMet_b", 12, 0, TMath::Pi ());
  w_single_delta_wmnuMet_b =        fs->make < TH1F > ("w_single_delta_phi_wmnuMet_b",  "w_single_delta_phi_wmnuMet_b", 12, 0, TMath::Pi ());
  b_single_delta_wmnuMet_b =        fs->make < TH1F > ("b_single_delta_phi_wmnuMet_b",  "b_single_delta_phi_wmnuMet_b", 12, 0, TMath::Pi ());
  c_single_delta_wmnuMet_b =        fs->make < TH1F > ("c_single_delta_phi_wmnuMet_b",  "c_single_delta_phi_wmnuMet_b", 12, 0, TMath::Pi ());
  t_single_delta_wmnuMet_b =        fs->make < TH1F > ("t_single_delta_phi_wmnuMet_b",  "t_single_delta_phi_wmnuMet_b", 12, 0, TMath::Pi ());

  w_single_delta_wenuMet_bb =       fs->make < TH1F > ("w_single_delta_phi_wenuMet_bb", "w_single_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());
  b_single_delta_wenuMet_bb =       fs->make < TH1F > ("b_single_delta_phi_wenuMet_bb", "b_single_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());
  c_single_delta_wenuMet_bb =       fs->make < TH1F > ("c_single_delta_phi_wenuMet_bb", "c_single_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());
  t_single_delta_wenuMet_bb =       fs->make < TH1F > ("t_single_delta_phi_wenuMet_bb", "t_single_delta_phi_wenuMet_bb", 12, 0, TMath::Pi ());
  w_single_delta_wmnuMet_bb =       fs->make < TH1F > ("w_single_delta_phi_wmnuMet_bb", "w_single_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());
  b_single_delta_wmnuMet_bb =       fs->make < TH1F > ("b_single_delta_phi_wmnuMet_bb", "b_single_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());
  c_single_delta_wmnuMet_bb =       fs->make < TH1F > ("c_single_delta_phi_wmnuMet_bb", "c_single_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());
  t_single_delta_wmnuMet_bb =       fs->make < TH1F > ("t_single_delta_phi_wmnuMet_bb", "t_single_delta_phi_wmnuMet_bb", 12, 0, TMath::Pi ());

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

  w_SVTX_mass     =     fs->make < TH1F > ("w_SVTX_mass",       "w_SVTX_mass;Mass [GeV]", 50, 0, 6);
  b_SVTX_mass     =     fs->make < TH1F > ("b_SVTX_mass",       "b_SVTX_mass;Mass [GeV]", 50, 0, 6);
  c_SVTX_mass     =     fs->make < TH1F > ("c_SVTX_mass",       "c_SVTX_mass;Mass [GeV]", 50, 0, 6);
  t_SVTX_mass     =     fs->make < TH1F > ("t_SVTX_mass",       "t_SVTX_mass;Mass [GeV]", 50, 0, 6);

  w_BJP       =     fs->make < TH1F > ("w_BJP",   "w_BJP", 50, 0, 10);
  b_BJP       =     fs->make < TH1F > ("b_BJP",   "b_BJP", 50, 0, 10);
  c_BJP       =     fs->make < TH1F > ("c_BJP",   "c_BJP", 50, 0, 10);
  t_BJP       =     fs->make < TH1F > ("t_BJP",   "t_BJP", 50, 0, 10);

  w_JBP       =     fs->make < TH1F > ("w_JBP",   "w_JBP", 50, 0, 3);
  b_JBP       =     fs->make < TH1F > ("b_JBP",   "b_JBP", 50, 0, 3);
  c_JBP       =     fs->make < TH1F > ("c_JBP",   "c_JBP", 50, 0, 3);
  t_JBP       =     fs->make < TH1F > ("t_JBP",   "t_JBP", 50, 0, 3);

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

  w_Ht =                fs->make < TH1F > ("w_Ht",              "w_Ht [GeV]", 50, 30., 1000.);
  b_Ht =                fs->make < TH1F > ("b_Ht",              "b_Ht [GeV]", 50, 30., 1000.);
  c_Ht =                fs->make < TH1F > ("c_Ht",              "c_Ht [GeV]", 50, 30., 1000.);
  t_Ht =                fs->make < TH1F > ("t_Ht",              "t_Ht [GeV]", 50, 30., 1000.);
  w_Ht_b =              fs->make < TH1F > ("w_Ht_b",            "w_Ht_b [GeV]", 50, 30., 1000.);
  b_Ht_b =              fs->make < TH1F > ("b_Ht_b",            "b_Ht_b [GeV]", 50, 30., 1000.);
  c_Ht_b =              fs->make < TH1F > ("c_Ht_b",            "c_Ht_b [GeV]", 50, 30., 1000.);
  t_Ht_b =              fs->make < TH1F > ("t_Ht_b",            "t_Ht_b [GeV]", 50, 30., 1000.);
  w_Ht_bb =             fs->make < TH1F > ("w_Ht_bb",           "w_Ht_bb [GeV]", 50, 30., 1000.);
  b_Ht_bb =             fs->make < TH1F > ("b_Ht_bb",           "b_Ht_bb [GeV]", 50, 30., 1000.);
  c_Ht_bb =             fs->make < TH1F > ("c_Ht_bb",           "c_Ht_bb [GeV]", 50, 30., 1000.);
  t_Ht_bb =             fs->make < TH1F > ("t_Ht_bb",           "t_Ht_bb [GeV]", 50, 30., 1000.);

  w_single_Ht_b =       fs->make < TH1F > ("w_single_Ht_b",            "w_single_Ht_b [GeV]", 50, 30., 1000.);
  b_single_Ht_b =       fs->make < TH1F > ("b_single_Ht_b",            "b_single_Ht_b [GeV]", 50, 30., 1000.);
  c_single_Ht_b =       fs->make < TH1F > ("c_single_Ht_b",            "c_single_Ht_b [GeV]", 50, 30., 1000.);
  t_single_Ht_b =       fs->make < TH1F > ("t_single_Ht_b",            "t_single_Ht_b [GeV]", 50, 30., 1000.);

  w_MET =               fs->make < TH1F > ("w_MET",             "w_MET;MET [GeV]", 50, 0., 250.);
  b_MET =               fs->make < TH1F > ("b_MET",             "b_MET;MET [GeV]", 50, 0., 250.);
  c_MET =               fs->make < TH1F > ("c_MET",             "c_MET;MET [GeV]", 50, 0., 250.);
  t_MET =               fs->make < TH1F > ("t_MET",             "t_MET;MET [GeV]", 50, 0., 250.);
  w_MET_sign = 	        fs->make < TH1F > ("w_MET_sign",        "w_MET_sign;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign = 	        fs->make < TH1F > ("b_MET_sign",        "b_MET_sign;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign = 	        fs->make < TH1F > ("c_MET_sign",        "c_MET_sign;MET significance [GeV]", 50, 0., 100.);
  t_MET_sign = 	        fs->make < TH1F > ("t_MET_sign",        "t_MET_sign;MET significance [GeV]", 50, 0., 100.);

  w_MET_b =             fs->make < TH1F > ("w_MET_b",         "w_MET_b;MET [GeV]", 50, 0., 250.);
  b_MET_b =             fs->make < TH1F > ("b_MET_b",         "b_MET_b;MET [GeV]", 50, 0., 250.);
  c_MET_b =             fs->make < TH1F > ("c_MET_b",         "c_MET_b;MET [GeV]", 50, 0., 250.);
  t_MET_b =             fs->make < TH1F > ("t_MET_b",         "t_MET_b;MET [GeV]", 50, 0., 250.);
  w_MET_sign_b = 	fs->make < TH1F > ("w_MET_sign_b",    "w_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign_b = 	fs->make < TH1F > ("b_MET_sign_b",    "b_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign_b = 	fs->make < TH1F > ("c_MET_sign_b",    "c_MET_sign_b;MET significance [GeV]", 50, 0., 100.);
  t_MET_sign_b = 	fs->make < TH1F > ("t_MET_sign_b",    "t_MET_sign_b;MET significance [GeV]", 50, 0., 100.);

  w_MET_bb =            fs->make < TH1F > ("w_MET_bb",         "w_MET_bb;MET [GeV]", 50, 0., 250.);
  b_MET_bb =            fs->make < TH1F > ("b_MET_bb",         "b_MET_bb;MET [GeV]", 50, 0., 250.);
  c_MET_bb =            fs->make < TH1F > ("c_MET_bb",         "c_MET_bb;MET [GeV]", 50, 0., 250.);
  t_MET_bb =            fs->make < TH1F > ("t_MET_bb",         "t_MET_bb;MET [GeV]", 50, 0., 250.);
  w_MET_sign_bb = 	fs->make < TH1F > ("w_MET_sign_bb",    "w_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);
  b_MET_sign_bb = 	fs->make < TH1F > ("b_MET_sign_bb",    "b_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);
  c_MET_sign_bb = 	fs->make < TH1F > ("c_MET_sign_bb",    "c_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);
  t_MET_sign_bb = 	fs->make < TH1F > ("t_MET_sign_bb",    "t_MET_sign_bb;MET significance [GeV]", 50, 0., 100.);

  w_Afb =               fs->make < TH1F > ("b_asymmetry",       "b_asymmetry", 10, -1, 1);

  h_JEC_uncert =        fs->make < TH1F > ("JEC uncert", "JEC uncert", 10, -0.5, 0.5);

  h_scaleFactor_first_ele =   fs->make < TH1F > ("h_scaleFactor_first_ele",   "h_scaleFactor_first_ele", 50, 0.95, 1.05);
  b_scaleFactor_first_ele =   fs->make < TH1F > ("b_scaleFactor_first_ele",   "b_scaleFactor_first_ele", 50, 0.95, 1.05);
  h_scaleFactor_first_muon =  fs->make < TH1F > ("h_scaleFactor_first_muon",  "h_scaleFactor_first_muon", 50, 0.95, 1.05);
  b_scaleFactor_first_muon =  fs->make < TH1F > ("b_scaleFactor_first_muon",  "b_scaleFactor_first_muon", 50, 0.95, 1.05);
  h_scaleFactor_second_ele =  fs->make < TH1F > ("h_scaleFactor_second_ele",  "h_scaleFactor_second_ele", 50, 0.95, 1.05);
  b_scaleFactor_second_ele =  fs->make < TH1F > ("b_scaleFactor_second_ele",  "b_scaleFactor_second_ele", 50, 0.95, 1.05);
  h_scaleFactor_second_muon = fs->make < TH1F > ("h_scaleFactor_second_muon", "h_scaleFactor_second_muon", 50, 0.95, 1.05);
  b_scaleFactor_second_muon = fs->make < TH1F > ("b_scaleFactor_second_muon", "b_scaleFactor_second_muon", 50, 0.95, 1.05);

  produces<std::vector<double>>("myEventWeight");

  produces<std::vector<math::XYZTLorentzVector>>("myElectrons");
  produces<std::vector<math::XYZTLorentzVector>>("myMuons");

  produces<std::vector<double>>("myPtZ");
  produces<std::vector<double>>("myPtZb");

  produces<std::vector<double>>("myMassZj");
  produces<std::vector<double>>("myMassZb");

  produces<std::vector<math::XYZTLorentzVector>>("myJets");
  produces<std::vector<double>>("myDeltaPhi");

  produces<std::vector<double>>("myHt");
  produces<std::vector<double>>("myHtb");

  produces<std::vector<double>>("myBJetsWeights");

  produces<std::vector<math::XYZTLorentzVector>>("myBJets");
  produces<std::vector<double>>("myBDeltaPhi");

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

  // Get muon collection
  edm::Handle < pat::MuonCollection > muons;
  iEvent.getByLabel ("matchedMuons", muons);

  // Get jet collection
  edm::Handle < vector < pat::Jet > > jets;
  iEvent.getByLabel ("goodJets", jets);

  // Get tracks
  edm::Handle < vector < reco::Track > > tracks;
  iEvent.getByLabel ("generalTracks", tracks);

  // Get METs
  edm::Handle < vector < pat::MET > > mets;
  iEvent.getByLabel (edm::InputTag ("patMETsPFlow"), mets);

  if (debug && mets->empty()) cout << "Warning: empty MET collection." << endl; 
  if (debug && mets->size()>1) cout << "Warning: MET collection size > 1." << endl; 

  edm::Handle<vector<reco::GenParticle> > genPart;
  iEvent.getByLabel ("genParticles", genPart);

  edm::Handle<vector<reco::GenJet> > gJets;
  iEvent.getByLabel(edm::InputTag("goodJets","genJets"), gJets);

  std::auto_ptr<std::vector<double>> myEventWeight( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myElectrons( new std::vector<math::XYZTLorentzVector> );
  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myMuons( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myPtZ( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myPtZb( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myMassZj( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myMassZb( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myJets( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myDeltaPhi( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myHt( new std::vector<double> );
  std::auto_ptr<std::vector<double>> myHtb( new std::vector<double> );

  std::auto_ptr<std::vector<double>> myBJetsWeights( new std::vector<double> );

  std::auto_ptr<std::vector<math::XYZTLorentzVector>> myBJets( new std::vector<math::XYZTLorentzVector> );

  std::auto_ptr<std::vector<double>> myBDeltaPhi( new std::vector<double> );

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

  for (pat::ElectronCollection::const_iterator ele = electrons->begin (); ele != electrons->end (); ++ele) {

    if (ele->pt()>30 && fabs(ele->eta())<2.1) {
	vect_ele.push_back (*ele);
    }

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


  // +++++++++ MUONS

  vector < pat::Muon > vect_muon;

  for (pat::MuonCollection::const_iterator muon = muons->begin (); muon != muons->end (); ++muon) {

    if (muon->pt()>25 && fabs(muon->eta())<2.1) {
      vect_muon.push_back (*muon);
    }

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


  // +++++++++ Decisions:

  if (vect_ele.size()==0 && vect_muon.size()==0) {if (debug) cout << "No isolated leptons!!" << endl;}
  if (vect_ele.size()==1 && vect_muon.size()==0) wenu_event = true;
  if (vect_muon.size()==1 && vect_ele.size()==0) wmnu_event = true;

  wenu_event = wenu_event && (lepton_ == "electron");
  wmnu_event = wmnu_event && (lepton_ == "muon");
  ee_event = ee_event && (lepton_ == "electron");
  mm_event = mm_event && (lepton_ == "muon");


  // +++++++++ SCALE FACTORS:

  if (isMC) {
    if (wenu_event) {
      scalFac_first_e  =  ElSF_->Val (vect_ele[0].pt(), vect_ele[0].eta());  // HLT SC: to be fixed!
      MyWeight = MyWeight * scalFac_first_e;
    }
    if (wmnu_event) {
      scalFac_first_m  = MuSF_->Val (vect_muon[0].pt(), vect_muon[0].eta()) * MuSF2_->Val (vect_muon[0].pt(), vect_muon[0].eta()) * MuSF3_->Val (vect_muon[0].pt(), vect_muon[0].eta());
      MyWeight = MyWeight * scalFac_first_m;
    }
  }

  if (isMC) {
    if (ee_event) {
      scalFac_first_e  =  ElSF_->Val (vect_ele[iele0].pt(), vect_ele[iele0].eta());
      scalFac_second_e =  ElSF_->Val (vect_ele[iele1].pt(), vect_ele[iele1].eta());
      MyWeight = MyWeight * scalFac_first_e * scalFac_second_e;
    }
    if (mm_event) {
      scalFac_first_m  = MuSF_->Val (vect_muon[imuon0].pt(), vect_muon[imuon0].eta()) * MuSF2_->Val (vect_muon[imuon0].pt(), vect_muon[imuon0].eta()) * MuSF3_->Val (vect_muon[imuon0].pt(), vect_muon[imuon0].eta());
      scalFac_second_m = MuSF2_->Val (vect_muon[imuon1].pt(), vect_muon[imuon1].eta()) * MuSF3_->Val (vect_muon[imuon1].pt(), vect_muon[imuon1].eta());
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

  vector < pat::Jet > vect_jets;
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

    if (fabs(jetNew.eta()) < 2.4 && jetNew.pt() > 25) {

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
  double DR_ej = 0;
  double DR_mj = 0;
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
        DR_ej   = sqrt(DEta_ej*DEta_ej + DPhi_ej*DPhi_ej);
        if (DR_ej < R) iflag_ee = false;
      }
      for (unsigned int k=0; k<vect_muon.size(); ++k) {
        DEta_mj = fabs(vect_muon[k].eta() - vect_jets[i].eta()); 
        DPhi_mj = fabs(vect_muon[k].phi() - vect_jets[i].phi());
        if (DPhi_mj > TMath::ACos(-1)) DPhi_mj= 2*TMath::ACos(-1) - DPhi_mj;
        DR_mj   = sqrt(DEta_mj*DEta_mj + DPhi_mj*DPhi_mj);
        if (DR_mj < R) iflag_mm = false;
      }
    }
  }
  //cout << "DR(e;j) ="<<DR_ej<<"   "<<"DR(m;j) ="<<DR_mj<< endl;

  if (Nb > 0 && pcut_) {
    double discrBJP = vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags");
    if (discrBJP <= 5.) {
      Nb = 0;
      vect_bjets.clear();
    }
  }

  ee_event = ee_event && Nb>0 && Nj>1 && iflag_ee;
  mm_event = mm_event && Nb>0 && Nj>1 && iflag_mm;

  //  wenu_event = wenu_event && iflag_ee;
  //  wmnu_event = wmnu_event && iflag_mm;

  wenu_event = wenu_event && Nb>0 && Nj>1;
  wmnu_event = wmnu_event && Nb>0 && Nj>1;

  if (debug && Nj<1) cout << "Warning: 0 Jets in the event!" << endl;
  if (debug && Nb<1) cout << "Warning: 0 b-Jets in the event!" << endl;

  // ++++++++ EVENT YIELDS:

  h_eventYields->Fill(1);
  if (wenu_event || wmnu_event) h_eventYields->Fill(2);
  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) h_eventYields->Fill(3);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) h_eventYields->Fill(4);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) h_eventYields->Fill(5);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb>1) h_eventYields->Fill(6);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb==1) h_eventYields->Fill(7);

  scalFac_b = btagSF(isMC, vect_jets, 0);
  w_eventYields->Fill(1, MyWeight*scalFac_b);
  if (wenu_event || wmnu_event) w_eventYields->Fill(2, MyWeight*scalFac_b);
  if ((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) w_eventYields->Fill(3, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) w_eventYields->Fill(4, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) w_eventYields->Fill(5, MyWeight*scalFac_b);
  w_eventYields->Fill(2, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb>1) w_eventYields->Fill(6, MyWeight*scalFac_b);
  w_eventYields->Fill(1, MyWeight*scalFac_b);
  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut && Nb==1) w_eventYields->Fill(7, MyWeight*scalFac_b);

  // ++++++++ MET PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    w_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
    w_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    if (ist) {
      t_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      t_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      b_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_MET->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      c_MET_sign->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      w_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      if (ist) {
	t_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
	t_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        b_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_MET_b->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        c_MET_sign_b->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_jets, 2);
      w_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
      w_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      if (ist) {
	t_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
	t_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        b_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_MET_bb->Fill (mets->empty() ? 0 : (*mets)[0].et(), MyWeight*scalFac_b);
        c_MET_sign_bb->Fill (mets->empty() ? 0 : (*mets)[0].significance(), MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ HT PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (Nj == 1) w_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      if (ist) {
        t_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) t_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) b_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_Ht_b->Fill (Ht, MyWeight*scalFac_b);
        if (Nj == 1) c_single_Ht_b->Fill (Ht, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_jets, 2);
      w_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      if (ist) {
        t_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_Ht_bb->Fill (Ht, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ WeNu PLOTS

  if (wenu_event && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wenu_b_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_jets, 2);
      w_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wenu_bb_wide->Fill (mt_wenu, MyWeight*scalFac_b);
      }
    }
  }

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    h_mt_wenu->Fill (mt_wenu);
    w_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
    w_delta_wenuMet->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
    if (ist) {
      t_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
      t_delta_wenuMet->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
      b_delta_wenuMet->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mt_wenu->Fill (mt_wenu, MyWeight*scalFac_b);
      c_delta_wenuMet->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b);
      w_delta_wenuMet_b->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b);
	t_delta_wenuMet_b->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b);
	b_delta_wenuMet_b->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wenu_b->Fill (mt_wenu, MyWeight*scalFac_b);
	c_delta_wenuMet_b->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_jets, 2);
      w_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
      w_delta_wenuMet_bb->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
	t_delta_wenuMet_bb->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
	b_delta_wenuMet_bb->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wenu_bb->Fill (mt_wenu, MyWeight*scalFac_b);
	c_delta_wenuMet_bb->Fill (deltaPhiMetEle, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ Zee PLOTS

  if (lepton_=="electron" && iele1!=-1 && Nj > 1 && Nb > 0 && vtx_cut && !met_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      if (ist) {
        t_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mass_ee_b_wide->Fill (diele_mass, MyWeight*scalFac_b);
      }
    }
  }

  if (ee_event && met_cut && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
      scalFac_b = btagSF(isMC, vect_jets, 1);
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
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mass_ee_b->Fill (diele_mass, MyWeight*scalFac_b);
        b_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
        b_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        b_mass_Zj_ee_b->Fill (zb_ee_mass, MyWeight*scalFac_b);
	if (Nj == 1) {
          b_single_pt_Z_ee_b->Fill (diele_pt, MyWeight*scalFac_b);
	  b_single_delta_ee_b->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
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
      scalFac_b = btagSF(isMC, vect_jets, 2);
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
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mass_ee_bb->Fill (diele_mass, MyWeight*scalFac_b);
        b_pt_Z_ee_bb->Fill (diele_pt, MyWeight*scalFac_b);
        b_delta_ee_bb->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        b_mass_Zj_ee_bb->Fill (zb_ee_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mass_ee_bb->Fill (diele_mass, MyWeight*scalFac_b);
        c_pt_Z_ee_bb->Fill (diele_pt, MyWeight*scalFac_b);
        c_delta_ee_bb->Fill (delta_phi_ee_b, MyWeight*scalFac_b);
        c_mass_Zj_ee_bb->Fill (zb_ee_mass, MyWeight*scalFac_b);
      }
    }
  }


  // ++++++++ WmNu PLOTS

  if (wmnu_event && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wmnu_b_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_jets, 2);
      w_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wmnu_bb_wide->Fill (mt_wmnu, MyWeight*scalFac_b);
      }
    }
  }

  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    h_mt_wmnu->Fill (mt_wmnu);
    w_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
    w_delta_wmnuMet->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
    if (ist) {
      t_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
      t_delta_wmnuMet->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
      b_delta_wmnuMet->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_mt_wmnu->Fill (mt_wmnu, MyWeight*scalFac_b);
      c_delta_wmnuMet->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
    }
    if (Nb == 1) {
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
      w_delta_wmnuMet_b->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
	t_delta_wmnuMet_b->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
	b_delta_wmnuMet_b->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wmnu_b->Fill (mt_wmnu, MyWeight*scalFac_b);
	c_delta_wmnuMet_b->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      }
    }
    if (Nb > 1) {
      scalFac_b = btagSF(isMC, vect_jets, 2);
      w_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
      w_delta_wmnuMet_bb->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      if (ist) {
        t_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
	t_delta_wmnuMet_bb->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
	b_delta_wmnuMet_bb->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mt_wmnu_bb->Fill (mt_wmnu, MyWeight*scalFac_b);
	c_delta_wmnuMet_bb->Fill (deltaPhiMetMuo, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ Zmm PLOTS

  if (lepton_=="muon" && imuon1!=-1 && Nj > 1 && Nb>0 && vtx_cut && !met_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
      scalFac_b = btagSF(isMC, vect_jets, 1);
      w_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      if (ist) {
        t_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mass_mm_b_wide->Fill (dimuon_mass, MyWeight*scalFac_b);
      }
    }
  }

  if (mm_event && met_cut && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
      scalFac_b = btagSF(isMC, vect_jets, 1);
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
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mass_mm_b->Fill (dimuon_mass, MyWeight*scalFac_b);
        b_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
        b_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        b_mass_Zj_mm_b->Fill (zb_mm_mass, MyWeight*scalFac_b);
        if (Nj == 1) {
          b_single_pt_Z_mm_b->Fill (dimuon_pt, MyWeight*scalFac_b);
          b_single_delta_mm_b->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
	}
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
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
      scalFac_b = btagSF(isMC, vect_jets, 2);
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
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_mass_mm_bb->Fill (dimuon_mass, MyWeight*scalFac_b);
        b_pt_Z_mm_bb->Fill (dimuon_pt, MyWeight*scalFac_b);
        b_delta_mm_bb->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        b_mass_Zj_mm_bb->Fill (zb_mm_mass, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_mass_mm_bb->Fill (dimuon_mass, MyWeight*scalFac_b);
        c_pt_Z_mm_bb->Fill (dimuon_pt, MyWeight*scalFac_b);
        c_delta_mm_bb->Fill (delta_phi_mm_b, MyWeight*scalFac_b);
        c_mass_Zj_mm_bb->Fill (zb_mm_mass, MyWeight*scalFac_b);
      }
    }
  }

  //  // ++++++++ ELECTRON+MUON Z PLOTS
  //
  //  if (lepton_=="electron+muon" && (iele1!=-1||imuon1!=-1) && Nj > 0 && vtx_cut && !met_cut) {
  //    w_mass_em_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //    if (ist) {
  //      t_mass_em_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //    }
  //    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
  //      b_mass_em_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //    }
  //    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
  //      c_mass_em_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //    }
  //    if (Nb > 0) {
  //      scalFac_b = btagSF(isMC, vect_jets, 0);
  //      w_mass_em_b_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      if (ist) {
  //        t_mass_em_b_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      }
  //      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
  //        b_mass_em_b_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      }
  //      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
  //        c_mass_em_b_wide->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      }
  //    }
  //  }
  //
  //  if (em_event && Nj > 0 && vtx_cut) {
  //    w_numberOfZ->Fill (zem->size(), MyWeight*scalFac_b);
  //    if (isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
  //      b_numberOfZ->Fill (zem->size(), MyWeight*scalFac_b);
  //    }
  //    if (isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
  //      c_numberOfZ->Fill (zem->size(), MyWeight*scalFac_b);
  //    }
  //  }
  //
  //  if (em_event && Nj > 0 && vtx_cut) {
  //    double delta_phi_em = fabs(dielemuon_phi - vect_jets[0].phi());
  //    if (delta_phi_em > acos (-1)) delta_phi_em = 2 * acos (-1) - delta_phi_em;
  //    h_mass_em->Fill (dielemuon_mass);
  //    w_mass_em->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //    w_pt_Z_em->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //    w_delta_em->Fill (delta_phi_em, MyWeight*scalFac_b);
  //    math::XYZTLorentzVector zj_em_p = vect_jets[0].p4() + z_em;
  //    double zj_em_mass = zj_em_p.mass();
  //    w_mass_Zj_em->Fill (zj_em_mass, MyWeight*scalFac_b);
  //    if (ist) {
  //      t_mass_em->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      t_pt_Z_em->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //      t_mass_Zj_em->Fill (zj_em_mass, MyWeight*scalFac_b);
  //    }
  //    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
  //      b_mass_em->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      b_pt_Z_em->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //      b_mass_Zj_em->Fill (zj_em_mass, MyWeight*scalFac_b);
  //    }
  //    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
  //      c_mass_em->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      c_pt_Z_em->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //      c_mass_Zj_em->Fill (zj_em_mass, MyWeight*scalFac_b);
  //    }
  //    if (Nb > 0 && met_cut) {
  //      scalFac_b = btagSF(isMC, vect_jets, 0);
  //      w_mass_em_b->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //      w_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //      math::XYZTLorentzVector zb_em_p = vect_jets[0].p4() + z_em;
  //      double zb_em_mass = zb_em_p.mass();
  //      w_mass_Zj_em_b->Fill (zb_em_mass, MyWeight*scalFac_b);
  //      double delta_phi_em_b = fabs(dielemuon_phi - vect_bjets[0].phi());
  //      if (delta_phi_em_b > acos (-1)) delta_phi_em_b = 2 * acos (-1) - delta_phi_em_b;
  //      w_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //      if (Nj == 1) {
  //        w_single_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //        w_single_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //      }
  //      if (ist) {
  //        t_mass_em_b->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //        t_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //        t_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //        t_mass_Zj_em_b->Fill (zb_em_mass, MyWeight*scalFac_b);
  //	if (Nj == 1) {
  //          t_single_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //	  t_single_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //	}
  //      }
  //      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
  //        b_mass_em_b->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //        b_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //        b_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //        b_mass_Zj_em_b->Fill (zb_em_mass, MyWeight*scalFac_b);
  //	if (Nj == 1) {
  //          b_single_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //	  b_single_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //	}
  //      }
  //      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
  //        c_mass_em_b->Fill (dielemuon_mass, MyWeight*scalFac_b);
  //        c_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //        c_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //        c_mass_Zj_em_b->Fill (zb_em_mass, MyWeight*scalFac_b);
  //        if (Nj == 1) {
  //          c_single_pt_Z_em_b->Fill (dielemuon_pt, MyWeight*scalFac_b);
  //	  c_single_delta_em_b->Fill (delta_phi_em_b, MyWeight*scalFac_b);
  //        }
  //      }
  //    }
  //  }

  // ++++++++ MISC PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    h_pu_weights->Fill (MyWeight*scalFac_b);
    h_tracks->Fill (tracks->size());
    w_tracks->Fill (tracks->size(), MyWeight*scalFac_b);
    h_recoVTX->Fill (NVtx);
    w_recoVTX->Fill (NVtx, MyWeight*scalFac_b);
  }

  // ++++++++  ELECTRONS PLOTS

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    w_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
    w_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
    w_first_ele_iso->Fill (vect_ele[0].chargedHadronIso + max((vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso()),0.0))/et, MyWeight*scalFac_b);
    if (ee_event) { // temporary !!!
      w_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
      w_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
    }
    if (ist) {
      t_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
      t_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
      t_first_ele_iso->Fill (vect_ele[0].chargedHadronIso + max((vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso()),0.0))/et, MyWeight*scalFac_b);
      if (ee_event) { // temporary !!!
        t_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
        t_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
      b_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
      b_first_ele_iso->Fill (vect_ele[0].chargedHadronIso + max((vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso()),0.0))/et, MyWeight*scalFac_b);
      if (ee_event) { // temporary !!!
        b_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
        b_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
      c_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
      c_first_ele_iso->Fill (vect_ele[0].chargedHadronIso + max((vect_ele[0].neutralHadronIso() + vect_ele[0].photonIso() - 0.5*vect_ele[0].puChargedHadronIso()),0.0))/vect_ele[0].et(), MyWeight*scalFac_b);
      if (ee_event) { // temporary !!!
        c_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
        c_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
      }
    }
  }

  if (wenu_event && mt_cut_wenu && vtx_cut && Nb == 1) {
    scalFac_b = btagSF(isMC, vect_jets, 1);
    w_first_ele_pt_b->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  }

  if (wenu_event && mt_cut_wenu && vtx_cut && Nb > 1) {
    scalFac_b = btagSF(isMC, vect_jets, 2);
    w_first_ele_pt_bb->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  }

  // ++++++++ MUONS PLOTS

  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    w_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
    w_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
    w_first_muon_iso->Fill (vect_muon[0].sumChargedHadronPt() + max(vect_muon[0].sumNeutralHadronEt() + vect_muon[0].sumPhotonEt() - 0.5*vect_muon[0].sumPUPt(),0.0))/vect_muon[0].pt(), MyWeight*scalFac_b);
    if (mm_event) { // temporary !!!
      w_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
      w_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
    }
    if (ist) {
      t_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
      t_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
      t_first_muon_iso->Fill (vect_muon[0].sumChargedHadronPt() + max(vect_muon[0].sumNeutralHadronEt() + vect_muon[0].sumPhotonEt() - 0.5*vect_muon[0].sumPUPt(),0.0))/vect_muon[0].pt(), MyWeight*scalFac_b);
      if (mm_event) { // temporary !!!
        t_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
        t_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
      b_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
      b_first_muon_iso->Fill (vect_muon[0].sumChargedHadronPt() + max(vect_muon[0].sumNeutralHadronEt() + vect_muon[0].sumPhotonEt() - 0.5*vect_muon[0].sumPUPt(),0.0))/vect_muon[0].pt(), MyWeight*scalFac_b);
      if (mm_event) { // temporary !!!
        b_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
        b_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
      }
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_first_muon_pt->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
      c_first_muon_eta->Fill (vect_muon[0].eta(), MyWeight*scalFac_b);
      c_first_muon_iso->Fill (vect_muon[0].sumChargedHadronPt() + max(vect_muon[0].sumNeutralHadronEt() + vect_muon[0].sumPhotonEt() - 0.5*vect_muon[0].sumPUPt(),0.0))/vect_muon[0].pt(), MyWeight*scalFac_b);
      if (mm_event) { // temporary !!!
        c_second_muon_pt->Fill (vect_muon[1].pt(), MyWeight*scalFac_b);
        c_second_muon_eta->Fill (vect_muon[1].eta(), MyWeight*scalFac_b);
      }
    }
  }

  if (wmnu_event && mt_cut_wmnu && vtx_cut && Nb == 1) {
    scalFac_b = btagSF(isMC, vect_jets, 1);
    w_first_muon_pt_b ->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
  }

  if (wmnu_event && mt_cut_wmnu && vtx_cut && Nb > 1) {
    scalFac_b = btagSF(isMC, vect_jets, 2);
    w_first_muon_pt_bb ->Fill (vect_muon[0].pt(), MyWeight*scalFac_b);
  }

  //  // ++++++++  ELECTRONS+MUONS PLOTS
  //
  //  if (em_event && Nj > 0 && vtx_cut) {
  //    if (iele1!=-1) {
  //      w_first_muon_pt->Fill (vect_muon[imuon0].pt(), MyWeight*scalFac_b);
  //      w_first_muon_eta->Fill (vect_muon[imuon0].eta(), MyWeight*scalFac_b);
  //      w_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
  //      w_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
  //    }
  //    if (imuon1!=-1) {
  //      w_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  //      w_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
  //      w_second_muon_pt->Fill (vect_muon[imuon1].pt(), MyWeight*scalFac_b);
  //      w_second_muon_eta->Fill (vect_muon[imuon1].eta(), MyWeight*scalFac_b);
  //    }
  //    if (ist) {
  //      if (iele1!=-1) {
  //        t_first_muon_pt->Fill (vect_muon[imuon0].pt(), MyWeight*scalFac_b);
  //        t_first_muon_eta->Fill (vect_muon[imuon0].eta(), MyWeight*scalFac_b);
  //        t_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
  //        t_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
  //      }
  //      if(imuon1!=-1) {
  //        t_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  //        t_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
  //        t_second_muon_pt->Fill (vect_muon[imuon1].pt(), MyWeight*scalFac_b);
  //        t_second_muon_eta->Fill (vect_muon[imuon1].eta(), MyWeight*scalFac_b);
  //      }
  //    }
  //    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
  //      if (iele1!=-1) {
  //        b_first_muon_pt->Fill (vect_muon[imuon0].pt(), MyWeight*scalFac_b);
  //        b_first_muon_eta->Fill (vect_muon[imuon0].eta(), MyWeight*scalFac_b);
  //        b_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
  //        b_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
  //      }
  //      if (imuon1!=-1) {
  //        b_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  //        b_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
  //        b_second_muon_pt->Fill (vect_muon[imuon1].pt(), MyWeight*scalFac_b);
  //        b_second_muon_eta->Fill (vect_muon[imuon1].eta(), MyWeight*scalFac_b);
  //      }
  //    }
  //    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
  //      if (iele1!=-1) {
  //        c_first_muon_pt->Fill (vect_muon[imuon0].pt(), MyWeight*scalFac_b);
  //        c_first_muon_eta->Fill (vect_muon[imuon0].eta(), MyWeight*scalFac_b);
  //        c_second_ele_pt->Fill (vect_ele[1].pt(), MyWeight*scalFac_b);
  //        c_second_ele_eta->Fill (vect_ele[1].eta(), MyWeight*scalFac_b);
  //      }
  //      if (imuon1!=-1) {
  //        c_first_ele_pt->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  //        c_first_ele_eta->Fill (vect_ele[0].eta(), MyWeight*scalFac_b);
  //        c_second_muon_pt->Fill (vect_muon[imuon1].pt(), MyWeight*scalFac_b);
  //        c_second_muon_eta->Fill (vect_muon[imuon1].eta(), MyWeight*scalFac_b);
  //      }
  //    }
  //  }
  //
  //  if (em_event && Nj > 0 && Nb > 0 && vtx_cut && met_cut) {
  //    scalFac_b = btagSF(isMC, vect_jets, 0);
  //    if (iele1!=-1) {
  //      w_first_muon_pt_b ->Fill (vect_muon[imuon0].pt(), MyWeight*scalFac_b);
  //    }
  //    if (imuon1!=-1) {
  //      w_first_ele_pt_b->Fill (vect_ele[0].pt(), MyWeight*scalFac_b);
  //    }
  //  }

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

    scalFac_b = btagSF(isMC, vect_jets, 0);
    w_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
    w_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
    w_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    if (ist) {
      t_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      t_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      t_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      b_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      b_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_SVTX_mass_jet->Fill (sumVertexMassJet, MyWeight*scalFac_b);
      c_SVTX_mass_trk->Fill (sumVertexMassTrk, MyWeight*scalFac_b);
      c_SVTX_mass->Fill (sumVertexMass, MyWeight*scalFac_b);
    }
  }

  // ++++++++ CSV PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    double discrSVTX = vect_bjets[0].bDiscriminator("combinedSecondaryVertexBJetTags");
    w_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    if (ist) {
      t_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_secondvtx_N_zoom->Fill (discrSVTX, MyWeight*scalFac_b);
    }
    if (sumVertexMass > 0.0 ) {
      w_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      if (ist) {
        t_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_secondvtx_N_mass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
    }
    if (sumVertexMass == 0.0 ) {
      w_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      if (ist) {
        t_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
      if (isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
        c_secondvtx_N_nomass->Fill (discrSVTX, MyWeight*scalFac_b);
      }
    }
  }

  // ++++++++ BJP/JBP PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
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
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
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
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
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
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb == 1 && vtx_cut) {
    double discrBJP = vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags");
    scalFac_b = btagSF(isMC, vect_jets, 0);
    w_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    if (ist) {
      b_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_BJP0->Fill (discrBJP, MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 1 && vtx_cut) {
    double discrBJP = vect_bjets[0].bDiscriminator("jetBProbabilityBJetTags");
    scalFac_b = btagSF(isMC, vect_jets, 0);
    w_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    if (ist) {
      t_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_BJP1->Fill (discrBJP, MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 1 && vtx_cut) {
    double discrBJP2 = vect_bjets[1].bDiscriminator("jetBProbabilityBJetTags");
    scalFac_b = btagSF(isMC, vect_jets, 2);
    w_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    if (ist) {
      t_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_BJP2->Fill (discrBJP2, MyWeight*scalFac_b);
    }
  }

  // ++++++++ JETS PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    h_jetmultiplicity->Fill (Nj);
    w_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
    w_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
    w_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    if (ist) {
      t_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
      t_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      t_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 5) {
      b_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
      b_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      b_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[0].partonFlavour()) == 4) {
      c_jetmultiplicity->Fill (Nj, MyWeight*scalFac_b);
      c_first_jet_pt->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      c_first_jet_eta->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    w_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
    w_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    if (ist) {
      t_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      t_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[1].partonFlavour()) == 5) {
      b_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      b_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[1].partonFlavour()) == 4) {
      c_second_jet_pt->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      c_second_jet_eta->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nj > 2 && vtx_cut) {
    w_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
    w_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    if (ist) {
      t_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      t_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[2].partonFlavour()) == 5) {
      b_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      b_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_jets[2].partonFlavour()) == 4) {
      c_third_jet_pt->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      c_third_jet_eta->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    }
  }

  // ++++++++ B JETS PLOTS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    w_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
    w_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
    w_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
    w_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
    w_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    if (ist) {
      t_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      t_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      t_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      t_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      t_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      b_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      b_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      b_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      b_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_bjetmultiplicity->Fill (Nb, MyWeight*scalFac_b);
      c_first_jet_pt_b->Fill (vect_jets[0].pt(), MyWeight*scalFac_b);
      c_first_jet_eta_b->Fill (vect_jets[0].eta(), MyWeight*scalFac_b);
      c_first_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      c_first_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 1 && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 2);
    w_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
    w_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
    w_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
    w_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    if (ist) {
      t_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      t_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      t_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      t_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[1].partonFlavour()) == 5) {
      b_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      b_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      b_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      b_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[1].partonFlavour()) == 4) {
      c_second_jet_pt_b->Fill (vect_jets[1].pt(), MyWeight*scalFac_b);
      c_second_jet_eta_b->Fill (vect_jets[1].eta(), MyWeight*scalFac_b);
      c_second_bjet_pt->Fill (vect_bjets[1].pt(), MyWeight*scalFac_b);
      c_second_bjet_eta->Fill (vect_bjets[1].eta(), MyWeight*scalFac_b);
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb > 2 && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 3);
    w_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
    w_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
    w_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight*scalFac_b);
    w_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight*scalFac_b);
    if (ist) {
      t_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      t_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
      t_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight*scalFac_b);
      t_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[2].partonFlavour()) == 5) {
      b_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      b_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
      b_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight*scalFac_b);
      b_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[2].partonFlavour()) == 4) {
      c_third_jet_pt_b->Fill (vect_jets[2].pt(), MyWeight*scalFac_b);
      c_third_jet_eta_b->Fill (vect_jets[2].eta(), MyWeight*scalFac_b);
      c_third_bjet_pt->Fill (vect_bjets[2].pt(), MyWeight*scalFac_b);
      c_third_bjet_eta->Fill (vect_bjets[2].eta(), MyWeight*scalFac_b);
    }
  }

  // ++++++++ SINGLE BJET

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && Nb == 1 && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 1);
    w_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
    w_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    if (ist) {
      t_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      t_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
      b_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      b_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
    if (!ist && isMC && fabs(vect_bjets[0].partonFlavour()) == 4) {
      c_single_bjet_pt->Fill (vect_bjets[0].pt(), MyWeight*scalFac_b);
      c_single_bjet_eta->Fill (vect_bjets[0].eta(), MyWeight*scalFac_b);
    }
  }

  // ++++++++ EXTRA PLOTS

  int Nf = 0;
  int Nbk = 0;

  double Afb = 0;

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    if (fabs (vect_bjets[0].eta()) > 0) Nf++;
    if (fabs (vect_bjets[0].eta()) < 0) Nbk++;
    if ((Nf+Nbk) != 0) Afb = (Nf - Nbk) / (Nf + Nbk);
    w_Afb->Fill (Afb, MyWeight*scalFac_b);
  }

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    h_scaleFactor_first_ele->Fill (scalFac_first_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    h_scaleFactor_second_ele->Fill (scalFac_second_e, MyWeight / (scalFac_first_e * scalFac_second_e));
    if (Nb == 1) {
      if (isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_scaleFactor_first_ele->Fill (scalFac_first_e, MyWeight / (scalFac_first_e * scalFac_second_e));
        b_scaleFactor_second_ele->Fill (scalFac_second_e, MyWeight / (scalFac_first_e * scalFac_second_e));
      }
    }
  }
  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    h_scaleFactor_first_muon->Fill (scalFac_first_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    h_scaleFactor_second_muon->Fill (scalFac_second_m, MyWeight / (scalFac_first_m * scalFac_second_m));
    if (Nb == 1) {
      if (isMC && fabs(vect_bjets[0].partonFlavour()) == 5) {
        b_scaleFactor_first_muon->Fill (scalFac_first_m, MyWeight / (scalFac_first_m * scalFac_second_m));
        b_scaleFactor_second_muon->Fill (scalFac_second_m, MyWeight / (scalFac_first_m * scalFac_second_m));
      }
    }
  }

  // ++++++++ OUTPUT COLLECTIONS

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    scalFac_b = btagSF(isMC, vect_jets, 0);
    myEventWeight->push_back(MyWeight*scalFac_b);
  }

  if (wenu_event && mt_cut_wenu && vtx_cut) {
    myElectrons->push_back(math::XYZTLorentzVector(vect_ele[0].px(),vect_ele[0].py(),vect_ele[0].pz(),vect_ele[0].energy()));
    //    if (iele1!=-1) {
    //      myElectrons->push_back(math::XYZTLorentzVector(vect_ele[1].px(),vect_ele[1].py(),vect_ele[1].pz(),vect_ele[1].energy()));
    //    }
    //    myPtZ->push_back(diele_pt);
    //    math::XYZTLorentzVector zj_ee_p = vect_jets[0].p4() + z_ee;
    //    double zj_ee_mass = zj_ee_p.mass();
    //    myMassZj->push_back(zj_ee_mass);
    //    if (Nb > 0) {
    //      myPtZb->push_back(diele_pt);
    //      math::XYZTLorentzVector zb_ee_p = vect_jets[0].p4() + z_ee;
    //      double zb_ee_mass = zb_ee_p.mass();
    //      myMassZb->push_back(zb_ee_mass);
    //    }
  }

  if (wmnu_event && mt_cut_wmnu && vtx_cut) {
    myMuons->push_back(math::XYZTLorentzVector(vect_muon[0].px(),vect_muon[0].py(),vect_muon[0].pz(),vect_muon[0].energy()));
    //    if (imuon1!=-1) {
    //      myMuons->push_back(math::XYZTLorentzVector(vect_muon[imuon1].px(),vect_muon[imuon1].py(),vect_muon[imuon1].pz(),vect_muon[imuon1].energy()));
    //    }
    //    myPtZ->push_back(dimuon_pt);
    //    math::XYZTLorentzVector zj_mm_p = vect_jets[0].p4() + z_mm;
    //    double zj_mm_mass = zj_mm_p.mass();
    //    myMassZj->push_back(zj_mm_mass);
    //    if (Nb > 0) {
    //      myPtZb->push_back(dimuon_pt);
    //      math::XYZTLorentzVector zb_mm_p = vect_jets[0].p4() + z_mm;
    //      double zb_mm_mass = zb_mm_p.mass();
    //      myMassZb->push_back(zb_mm_mass);
    //    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    for (unsigned int i=0; i<vect_jets.size(); ++i) {
      myJets->push_back(math::XYZTLorentzVector(vect_jets[i].px(),vect_jets[i].py(),vect_jets[i].pz(),vect_jets[i].energy()));
    }
    if (Nb == 1) {
      for (unsigned int i=0; i<vect_bjets.size(); ++i) {
        scalFac_b = btagSF(isMC, vect_jets, 1);
        myBJetsWeights->push_back(scalFac_b);
        myBJets->push_back(math::XYZTLorentzVector(vect_bjets[i].px(),vect_bjets[i].py(),vect_bjets[i].pz(),vect_bjets[i].energy()));
      }
    }
  }

  if (((wenu_event && mt_cut_wenu) || (wmnu_event && mt_cut_wmnu)) && vtx_cut) {
    myHt->push_back(Ht);
    if (Nb == 1) {
      myHtb->push_back(Ht);
    }
  }

  //  if (wenu_event && mt_cut_wenu && Nj > 0 && vtx_cut) {
  //    double delta_phi_ee = fabs(diele_phi - vect_jets[0].phi());
  //    if (delta_phi_ee > acos (-1)) delta_phi_ee = 2 * acos (-1) - delta_phi_ee;
  //    myDeltaPhi->push_back(delta_phi_ee);
  //    if (Nb > 0) {
  //      double delta_phi_ee_b = fabs(diele_phi - vect_bjets[0].phi());
  //      if (delta_phi_ee_b > acos (-1)) delta_phi_ee_b = 2 * acos (-1) - delta_phi_ee_b;
  //      myBDeltaPhi->push_back(delta_phi_ee_b);
  //    }
  //  }
  //
  //  if (wmnu_event && mt_cut_wmnu && Nj > 0 && vtx_cut) {
  //    double delta_phi_mm = fabs(dimuon_phi - vect_jets[0].phi());
  //    if (delta_phi_mm > acos (-1)) delta_phi_mm = 2 * acos (-1) - delta_phi_mm;
  //    myDeltaPhi->push_back(delta_phi_mm);
  //    if (Nb > 0) {
  //      double delta_phi_mm_b = fabs(dimuon_phi - vect_bjets[0].phi());
  //      if (delta_phi_mm_b > acos (-1)) delta_phi_mm_b = 2 * acos (-1) - delta_phi_mm_b;
  //      myBDeltaPhi->push_back(delta_phi_mm_b);
  //    }
  //  }

  iEvent.put( myEventWeight, "myEventWeight" );

  iEvent.put( myElectrons, "myElectrons" );
  iEvent.put( myMuons, "myMuons" );

  iEvent.put( myPtZ, "myPtZ" );
  iEvent.put( myPtZb, "myPtZb" );

  iEvent.put( myMassZj, "myMassZj" );
  iEvent.put( myMassZb, "myMassZb" );

  iEvent.put( myJets, "myJets" );
  iEvent.put( myDeltaPhi, "myDeltaPhi" );

  iEvent.put( myHt, "myHt" );
  iEvent.put( myHtb, "myHtb" );

  iEvent.put( myBJetsWeights, "myBJetsWeights" );

  iEvent.put( myBJets, "myBJets" );
  iEvent.put( myBDeltaPhi, "myBDeltaPhi" );

}

// ------------ method called once each job just before starting event loop ------------
void WbAnalyzer::beginJob () {
  jetCorrectionUncertainty_ = new JetCorrectionUncertainty(path_ + "/" + "Summer13_V4_DATA_Uncertainty_AK5PFchs.txt");
  LumiWeights_ = edm::LumiReWeighting(path_ + "/" + "pileup_" + pileupMC_ + ".root", path_ + "/" + "pileup_2012_" + pileupDT_ + ".root", "pileup", "pileup");

  ElSF_  = new table(path_ + "/" + "ele_eff.txt");
  ElSF2_ = new table(path_ + "/" + "ele_eff2.txt");
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
