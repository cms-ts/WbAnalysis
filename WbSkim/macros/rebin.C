#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

TH1F* rebin(TH1F* old) {

  int ngroup=0;
  string name = old->GetName();

  if (name.find("first_jet_eta")!=string::npos) ngroup = 1;

  if (ngroup>1) old->Rebin(ngroup);

  return old;
}

TH2F* rebin(TH2F* old) {

  int ngroup=0;
  string name = old->GetName();

  if (name.find("first_jet_eta")!=string::npos) ngroup = 1;

  if (ngroup>1) old->Rebin2D(ngroup);

  return old;
}
