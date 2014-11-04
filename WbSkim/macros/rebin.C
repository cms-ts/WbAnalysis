#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

TH1F* rebin(TH1F* old) {

  int ngroup=0;
  string name = old->GetName();

//  if (name.find("first_jet_eta")!=string::npos) ngroup = 2;
//  if (name.find("first_jet_mass")!=string::npos) ngroup = 2;
//  if (name.find("first_jet_pt")!=string::npos) ngroup = 2;
//  if (name.find("second_jet_eta")!=string::npos) ngroup = 2;
//  if (name.find("second_jet_mass")!=string::npos) ngroup = 2;
//  if (name.find("second_jet_pt")!=string::npos) ngroup = 2;
  if (name.find("w_mt_")!=string::npos) ngroup = 2;
  if (name.find("b_mt_")!=string::npos) ngroup = 2;
  if (name.find("c_mt_")!=string::npos) ngroup = 2;
  if (name.find("t_mt_")!=string::npos) ngroup = 2;
//  if (name.find("w_Ht_")!=string::npos) ngroup = 2;
//  if (name.find("b_Ht_")!=string::npos) ngroup = 2;
//  if (name.find("c_Ht_")!=string::npos) ngroup = 2;
//  if (name.find("t_Ht_")!=string::npos) ngroup = 2;
//  if (name.find("delta_wenu")!=string::npos) ngroup = 2;
//  if (name.find("delta_wmnu")!=string::npos) ngroup = 2;
//  if (name.find("deltaR_wenu")!=string::npos) ngroup = 2;
//  if (name.find("deltaR_wmnu")!=string::npos) ngroup = 2;

  if (ngroup>1) old->Rebin(ngroup);

  return old;
}

TH2F* rebin(TH2F* old) {

  int ngroup=0;
  string name = old->GetName();

  if (name.find("first_jet_eta")!=string::npos) ngroup = 2;
  if (name.find("first_jet_mass")!=string::npos) ngroup = 2;
  if (name.find("first_jet_pt")!=string::npos) ngroup = 2;
  if (name.find("second_jet_eta")!=string::npos) ngroup = 2;
  if (name.find("second_jet_mass")!=string::npos) ngroup = 2;
  if (name.find("second_jet_pt")!=string::npos) ngroup = 2;
  if (name.find("w_mt_")!=string::npos) ngroup = 2;
  if (name.find("b_mt_")!=string::npos) ngroup = 2;
  if (name.find("c_mt_")!=string::npos) ngroup = 2;
  if (name.find("t_mt_")!=string::npos) ngroup = 2;
  if (name.find("w_Ht_")!=string::npos) ngroup = 2;
  if (name.find("b_Ht_")!=string::npos) ngroup = 2;
  if (name.find("c_Ht_")!=string::npos) ngroup = 2;
  if (name.find("t_Ht_")!=string::npos) ngroup = 2;
  if (name.find("delta_wenu")!=string::npos) ngroup = 2;
  if (name.find("delta_wmnu")!=string::npos) ngroup = 2;
  if (name.find("deltaR_wenu")!=string::npos) ngroup = 2;
  if (name.find("deltaR_wmnu")!=string::npos) ngroup = 2;

  if (ngroup>1) old->Rebin2D(ngroup, ngroup);

  return old;
}
