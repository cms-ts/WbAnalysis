// -*- C++ -*-
//
// Package:    WbFilter
// Class:      WbFilter
// 
/**\class WbFilter WbFilter.cc Zbanalysis/WbFilter/src/WbFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Schizzi
//         Created:  Thu Nov  1 11:32:14 CET 2012
// $Id: WbFilter.cc,v 1.2 2013/04/17 07:23:07 dellaric Exp $
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
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TLorentzVector.h"

//
// class declaration
//

class WbFilter : public edm::EDFilter {

public:

  explicit WbFilter(const edm::ParameterSet&);
  ~WbFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  int count_;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

  edm::Service<TFileService> fs; 

//
// constructors and destructor
//
};

WbFilter::WbFilter(const edm::ParameterSet& iConfig) {

  count_ = 0;

}


WbFilter::~WbFilter() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool WbFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   using namespace edm;

   // get electron collection
   edm::Handle < pat::ElectronCollection > electrons;
   iEvent.getByLabel ("matchedElectrons", electrons);

   // get muon collection
   edm::Handle < pat::MuonCollection > muons;
   iEvent.getByLabel ("matchedMuons", muons);

   // get electron collection with inverted ISO cuts for QCD studies
   edm::Handle < pat::ElectronCollection > electronsQCD;
   iEvent.getByLabel ("matchedElectronsQCD", electronsQCD);

   // get muon collection with inverted ISO cuts for QCD studies
   edm::Handle < pat::MuonCollection > muonsQCD;
   iEvent.getByLabel ("matchedMuonsQCD", muonsQCD);

   //get jet collection
   edm::Handle<std::vector<pat::Jet>  > jets;
   iEvent.getByLabel("goodJets",jets);

   //std::cout<<"numero j="<<jets->size()<<std::endl;

   if (electrons->size()==0 && muons->size()==0 && electronsQCD->size()==0 && muonsQCD->size()==0) return false;

   if (jets->size()==0) return false;

   bool hasEle=false;

   for (pat::ElectronCollection::const_iterator ele = electrons->begin (); ele != electrons->end (); ++ele) {
     if (ele->triggerObjectMatches().size()>0) hasEle=true;
   }
   
   bool hasMuo=false;
   
   for (pat::MuonCollection::const_iterator muon = muons->begin (); muon != muons->end (); ++muon) {
     if (muon->triggerObjectMatches().size()>0) hasMuo=true;
   }

   if (!hasEle && !hasMuo && electronsQCD->size()==0 && muonsQCD->size()==0) return false;

   ++count_;

   int prescaleFactorEle=1;
   if (!hasEle && electronsQCD->size()!=0 && count_%prescaleFactorEle!=0) return false;

   int prescaleFactorMuo=20;
   if (!hasMuo && muonsQCD->size()!=0 && count_%prescaleFactorMuo!=0) return false;

   for (std::vector < pat::Jet >::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {

     if (fabs(jet->eta())<2.5) {
       double discrCSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
       //if (discrCSV > 0.244) return true; // CSVL
       if (discrCSV > 0.679) return true; // CSVM
       //if (discrCSV > 0.898) return true; // CSVT
     }

   }

   return false;

}
// ------------ method called once each job just before starting event loop  ------------
void WbFilter::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void WbFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool WbFilter::beginRun(edm::Run&, edm::EventSetup const&) { 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool WbFilter::endRun(edm::Run&, edm::EventSetup const&) {
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool WbFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool WbFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
  return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(WbFilter);
