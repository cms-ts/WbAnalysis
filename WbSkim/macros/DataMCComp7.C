#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v09.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

TH1F* read(string subdir, string title, int ilepton) {
  TH1F* hist;
  TFile* file=0;
// to be removed when using unfolded data
  string title_tmp = title;
  if (title_tmp.find("_bb")==string::npos) {
    title_tmp = title_tmp + "b";
  }
  if (ilepton==1) {
    //file = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.root").c_str());
    file = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + title_tmp + "_xsecs.root").c_str());
  }
  if (ilepton==2) {
    //file = TFile::Open((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.root").c_str());
    file = TFile::Open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + title_tmp + "_xsecs.root").c_str());
  }
  hist = (TH1F*)gDirectory->Get(title.c_str())->Clone();
  hist->SetDirectory(0);
  file->Close();
  return hist;
}

void DataMCComp7(string title="", int plot=0, int ilepton=1, int isratio=1) {

//int drawInclusive = 0; // do not plot the "inclusive" histogram
int drawInclusive = 1; // do plot the "inclusive" histogram

int useSysBfit2=0;
//int useSysBfit2=1; // include Bfit2 systematics

int useSysDR=0;
//int useSysDR=1; // include deltaR systematics

int useSysRMS=0;
//int useSysRMS=1; // include xsec RMS systematics

//int useSysUnfold=0;
int useSysUnfold=1; // include unfolding systematics

int useSysUnfoldSherpa=0;
//int useSysUnfoldSherpa=1; // use Sherpa for unfolding systematics

int useSysUnfoldMadGraph4FS=0;
//int useSysUnfoldMadGraph4FS=1; // use MadGraph 4FS for unfolding systematics

//int useSysUnfoldWeight=0;
int useSysUnfoldWeight=1; // use data weighted MC for unfolding systematics

int useMC=0;
//int useMC=1; // use MC prediction

double ele_eff_sys=0.015;
double muo_eff_sys=0.020;
double btag_sys=0.030;
double lumi_sys=0.026;

string subdir="0";

	if (gROOT->GetVersionInt() >= 53401) {
	  gROOT->GetColor(kRed)->SetAlpha(0.5);
	  if (!useMC) gROOT->GetColor(kRed)->SetAlpha(0.0);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	  gROOT->GetColor(kOrange+7)->SetAlpha(0.5);
	}

	/* purity */

	double c_b=1.0;
	double ec_b=0.0;
	double c_c=1.0;
	double ec_c=0.0;
	double c_uds=1.0;
	double ec_uds=0.0;

	ifstream in;
	if (ilepton==1) {
	  //in.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_BJP_doFit" + ".dat").c_str());
	}
	if (ilepton==2) {
	  //in.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_BJP_doFit" + ".dat").c_str());
	}
	//in >> c_uds >> ec_uds;
	//in >> c_b >> ec_b;
	//in >> c_c >> ec_c;
	//in.close();

	double Lumi2012=0;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

	double norm1 = ((Lumi2012 * Xsec_wj) / Ngen_wj);

	if (title.empty()) title = "w_jetmultiplicity";

	if (ilepton==1) {
	  if (title.find("muon")!=string::npos) return;
	  if (title.find("mm")!=string::npos) return;
	  if (title.find("wmnu")!=string::npos) return;
	}
	if (ilepton==2) {
	  if (title.find("ele")!=string::npos) return;
	  if (title.find("ee")!=string::npos) return;
	  if (title.find("wenu")!=string::npos) return;
	}

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "Wj_merge.root").c_str());
	TFile *mcg = TFile::Open((path + "/" + version + "/" + "Wj_gen_merge.root").c_str());

	string title_b = title;

	if (title.find("_bjet_")!=string::npos) {
	  title.erase(title.find("_bjet_")+1, 1);
	} else {
	  if (title_b.find("_b")==string::npos) {
	    title_b = title + "_b";
	  } else {
	    title_b = title + "b";
	  }
        }

	TH1F* h_data;
	TH1F* h_data_b;
        h_data = read(subdir, title, ilepton);
        h_data_b = read(subdir, title_b, ilepton);
	h_data->SetStats(0);
	h_data_b->SetStats(0);

	const int NMAX = 100;
	TH1F* h_data_scan[NMAX];
	TH1F* h_data_b_scan[NMAX];

	for (int i=0;i<NMAX;i++) {
	  h_data_scan[i] = 0;
	  h_data_b_scan[i] = 0;
	  if (i<=13) {
	    stringstream ss;
	    ss << i;
	    h_data_scan[i] = read(ss.str(), title, ilepton);
	    h_data_b_scan[i] = read(ss.str(), title_b, ilepton);
	  }
	}
	if (useSysUnfoldMadGraph4FS) {
	  h_data_scan[77] = read("77", title, ilepton);
	  h_data_b_scan[77] = read("77", title_b, ilepton);
	}
	if (useSysUnfoldWeight) {
	  h_data_scan[66] = read("66", title, ilepton);
	  h_data_b_scan[66] = read("66", title_b, ilepton);
	}
	if (useSysDR) {
	  h_data_scan[88] = read("88", title, ilepton);
	  h_data_b_scan[88] = read("88", title_b, ilepton);
	}
	if (useSysBfit2) {
	  h_data_scan[99] = read("99", title, ilepton);
	  h_data_b_scan[99] = read("99", title_b, ilepton);
	}

	if (ilepton==1) mc1->cd("demoEle");
	if (ilepton==2) mc1->cd("demoMuo");
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1b_b = (TH1F*)gDirectory->Get(("b"+title_b.substr(1)).c_str());

	if (ilepton==1) mcg->cd("demoEleGen");
	if (ilepton==2) mcg->cd("demoMuoGen");
	TH1F* h_mcg = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg_b = (TH1F*)gDirectory->Get(title_b.c_str());

	h_mc1->Sumw2();
	h_mcg->Sumw2();

	h_mc1b_b->Sumw2();
	h_mcg_b->Sumw2();

	h_mc1->Scale(norm1);
	h_mcg->Scale(norm1);

	h_mc1b_b->Scale(norm1);
	h_mcg_b->Scale(norm1);

	h_mc1b_b->Scale(c_b);
	for (int i=0; i<=h_mc1b_b->GetNbinsX()+1; i++) {
	  float e = TMath::Power(h_mc1b_b->GetBinError(i),2);
	  e = e + TMath::Power(h_mc1b_b->GetBinContent(i)*(ec_b/c_b),2);
	  h_mc1b_b->SetBinError(i, TMath::Sqrt(e));
	}

        if (ilepton==1) {
	  TFile f((path + "/electrons/" + version + "/" + subdir +"/efficiency/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	  TFile f_b((path + "/electrons/" + version + "/" + subdir +"/efficiency/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
// FIXME
if (f.IsOpen()&&f_b.IsOpen()) {
// FIXME
	  TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	  TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	  h->SetDirectory(0);
	  h_b->SetDirectory(0);
	  f.Close();
	  f_b.Close();
	  h_mc1->Divide(h);
	  h_mc1b_b->Divide(h_b);
// FIXME
}
// FIXME
        }
	if (ilepton==2) {
	  TFile f((path + "/muons/" + version + "/" + subdir +"/efficiency/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	  TFile f_b((path + "/muons/" + version + "/" + subdir +"/efficiency/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
// FIXME
if (f.IsOpen()&&f_b.IsOpen()) {
// FIXME
	  TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	  TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	  h->SetDirectory(0);
	  h_b->SetDirectory(0);
	  f.Close();
	  f_b.Close();
	  h_mc1->Divide(h);
	  h_mc1b_b->Divide(h_b);
// FIXME
}
// FIXME
        }

// to be restored when using unfolded data
//	h_data->Scale(1./Lumi2012, "width");
//	h_data_b->Scale(1./Lumi2012, "width");
//	for (int i=0;i<NMAX;i++) {
//	  if (h_data_scan[i]) h_data_scan[i]->Scale(1./Lumi2012, "width");
//	  if (h_data_b_scan[i]) h_data_b_scan[i]->Scale(1./Lumi2012, "width");
//	}
	h_mc1->Scale(1./Lumi2012, "width");
	h_mc1b_b->Scale(1./Lumi2012, "width");
	if (isratio==1) {
	  h_data_b->Divide(h_data);
	  h_data_b->Scale(100.);
	  for (int i=0;i<NMAX;i++) {
	    if (h_data_scan[i]) h_data_b_scan[i]->Divide(h_data_scan[i]);
	    if (h_data_b_scan[i]) h_data_b_scan[i]->Scale(100.);
	  }
	  h_mc1b_b->Divide(h_mc1);
	  h_mc1b_b->Scale(100.);
	}

	TH1F* stat_bkg = (TH1F*)h_data->Clone();
	TH1F* stat_b_bkg = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_scan[13]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))/(TMath::Power(1.1,2)-1));
	  stat_bkg->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data->GetBinError(i),2)-TMath::Power(stat_bkg->GetBinError(i),2)));
	  h_data->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_b_scan[13]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))/(TMath::Power(1.1,2)-1));
	  stat_b_bkg->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_b->GetBinError(i),2)-TMath::Power(stat_b_bkg->GetBinError(i),2)));
	  h_data_b->SetBinError(i, val);
	}

	TH1F* syst_eff = (TH1F*)h_data->Clone();
	TH1F* syst_b_eff = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (ilepton==1) val = ele_eff_sys * h_data_scan[0]->GetBinContent(i);
	  if (ilepton==2) val = muo_eff_sys * h_data_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_eff->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (ilepton==1) val = ele_eff_sys * h_data_b_scan[0]->GetBinContent(i);
	  if (ilepton==2) val = muo_eff_sys * h_data_b_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_b_eff->SetBinError(i, val);
	}

	TH1F* syst_jec = (TH1F*)h_data->Clone();
	TH1F* syst_b_jec = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[2]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_scan[1]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i)));
	  syst_jec->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[2]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_b_scan[1]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i)));
	  syst_b_jec->SetBinError(i, val);
	}

	TH1F* syst_jer = (TH1F*)h_data->Clone();
	TH1F* syst_b_jer = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[12]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_scan[11]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i)));
	  syst_jer->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[12]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_b_scan[11]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i)));
	  syst_b_jer->SetBinError(i, val);
	}

	TH1F* syst_pu = (TH1F*)h_data->Clone();
	TH1F* syst_b_pu = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[3]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_scan[4]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i)));
	  syst_pu->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[3]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_b_scan[4]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i)));
	  syst_b_pu->SetBinError(i, val);
	}

	TH1F* syst_dr = (TH1F*)h_data->Clone();
	TH1F* syst_b_dr = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysDR) {
	    val = TMath::Abs(h_data_scan[88]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.0,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_scan[88]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))));
	  }
	  syst_dr->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysDR) {
	    val = TMath::Abs(h_data_b_scan[88]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.0,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_b_scan[88]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  syst_b_dr->SetBinError(i, val);
	}

	TH1F* syst_bkg = (TH1F*)h_data->Clone();
	TH1F* syst_b_bkg = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[10]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  syst_bkg->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[10]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  syst_b_bkg->SetBinError(i, val);
	}

	TH1F* stat_top = (TH1F*)h_data->Clone();
	TH1F* stat_b_top = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[5]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  stat_top->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[5]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  stat_b_top->SetBinError(i, val);
	}

	TH1F* stat_bfit = (TH1F*)h_data->Clone();
	TH1F* stat_b_bfit = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[6]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  stat_bfit->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[6]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  stat_b_bfit->SetBinError(i, val);
	}

	TH1F* syst_bfit2 = (TH1F*)h_data->Clone();
	TH1F* syst_b_bfit2 = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysBfit2) {
	    val = TMath::Abs(h_data_scan[99]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_scan[0]->GetBinError(i),2)-TMath::Power(h_data_scan[99]->GetBinError(i),2))));
	  }
	  syst_bfit2->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysBfit2) {
	    val = TMath::Abs(h_data_b_scan[99]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_b_scan[0]->GetBinError(i),2)-TMath::Power(h_data_b_scan[99]->GetBinError(i),2))));
	  }
	  syst_b_bfit2->SetBinError(i, val);
	}

	TH1F* syst_btag = (TH1F*)h_data->Clone();
	TH1F* syst_b_btag = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (title.find("_b")!=string::npos) {
	    val = btag_sys * h_data_scan[0]->GetBinContent(i);
	  }
	  syst_btag->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = btag_sys * h_data_b_scan[0]->GetBinContent(i);
	  syst_b_btag->SetBinError(i, val);
	}

	TH1F* stat_unfold = (TH1F*)h_data->Clone();
	TH1F* stat_b_unfold = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_scan[7]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2)));
	  stat_unfold->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2)));
	  stat_b_unfold->SetBinError(i, val);
	}

	TH1F* syst_unfold = (TH1F*)h_data->Clone();
	TH1F* syst_b_unfold = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysUnfold) {
	    val = TMath::Abs(h_data_scan[9]->GetBinContent(i)-h_data_scan[7]->GetBinContent(i));
	    val = TMath::Power(val,2);
	    //val = val - (TMath::Power(h_data_scan[9]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2));
	    //val = val - (TMath::Power(h_data_scan[7]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2));
	    val = TMath::Sqrt(TMath::Max(0.,val));
	    if (useSysUnfoldSherpa) {
	      val = TMath::Abs(h_data_scan[8]->GetBinContent(i)-h_data_scan[7]->GetBinContent(i));
	      val = TMath::Power(val,2);
	      //val = val - (TMath::Power(h_data_scan[8]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2));
	      //val = val - (TMath::Power(h_data_scan[7]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2));
	      val = TMath::Sqrt(TMath::Max(0.,val));
	    }
	    if (useSysUnfoldWeight) {
	      val = TMath::Abs(h_data_scan[66]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	    }
	  }
	  syst_unfold->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysUnfold) {
	    val = TMath::Abs(h_data_b_scan[9]->GetBinContent(i)-h_data_b_scan[7]->GetBinContent(i));
	    val = TMath::Power(val,2);
	    //val = val - (TMath::Power(h_data_b_scan[9]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2));
	    //val = val - (TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2));
	    val = TMath::Sqrt(TMath::Max(0.,val));
	    if (useSysUnfoldSherpa) {
	      val = TMath::Abs(h_data_b_scan[8]->GetBinContent(i)-h_data_b_scan[7]->GetBinContent(i));
	      val = TMath::Power(val,2);
	      //val = val - (TMath::Power(h_data_b_scan[8]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2));
	      //val = val - (TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2));
	      val = TMath::Sqrt(TMath::Max(0.,val));
	    }
	    if (useSysUnfoldMadGraph4FS) {
	      val = TMath::Abs(h_data_b_scan[77]->GetBinContent(i)-h_data_b_scan[7]->GetBinContent(i));
	      val = TMath::Power(val,2);
	      //val = val - (TMath::Power(h_data_b_scan[77]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2));
	      //val = val - (TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2));
	      val = TMath::Sqrt(TMath::Max(0.,val));
	    }
	    if (useSysUnfoldWeight) {
	      val = TMath::Abs(h_data_b_scan[66]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	    }
	  }
	  syst_b_unfold->SetBinError(i, val);
	}

	TH1F* syst_lumi = (TH1F*)h_data->Clone();
	TH1F* syst_b_lumi = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = lumi_sys * h_data_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_lumi->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = lumi_sys * h_data_b_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_b_lumi->SetBinError(i, val);
	}

	float sum1, sum2, sum3, sum4, sum5;
	float sum1_b, sum2_b, sum3_b, sum4_b, sum5_b;
	ifstream in1, in2, in3, in4, in5;
	if (ilepton==1) {
	  //in1.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_mt_wenu_b_wide" + "_xsecs_unfolding.dat").c_str());
	  //in2.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_first_jet_pt_b" + "_xsecs_unfolding.dat").c_str());
	  //in3.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_first_jet_eta_b" + "_xsecs_unfolding.dat").c_str());
	  //in4.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_second_jet_pt_b" + "_xsecs_unfolding.dat").c_str());
	  //in5.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_second_jet_eta_b" + "_xsecs_unfolding.dat").c_str());
	  in1.open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + "w_mt_wenu_bb_wide" + "_xsecs.dat").c_str());
	  in2.open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + "w_first_jet_pt_bb" + "_xsecs.dat").c_str());
	  in3.open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + "w_first_jet_eta_bb" + "_xsecs.dat").c_str());
	  in4.open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + "w_second_jet_pt_bb" + "_xsecs.dat").c_str());
	  in5.open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + "w_second_jet_eta_bb" + "_xsecs.dat").c_str());
	}
	if (ilepton==2) {
	  //in1.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_mt_wmnu_b_wide" + "_xsecs_unfolding.dat").c_str());
	  //in2.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_first_jet_pt_b" + "_xsecs_unfolding.dat").c_str());
	  //in3.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_first_jet_eta_b" + "_xsecs_unfolding.dat").c_str());
	  //in4.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_second_jet_pt_b" + "_xsecs_unfolding.dat").c_str());
	  //in5.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/" + "w_second_jet_eta_b" + "_xsecs_unfolding.dat").c_str());
	  in1.open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + "w_mt_wmnu_bb_wide" + "_xsecs.dat").c_str());
	  in2.open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + "w_first_jet_pt_bb" + "_xsecs.dat").c_str());
	  in3.open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + "w_first_jet_eta_bb" + "_xsecs.dat").c_str());
	  in4.open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + "w_second_jet_pt_bb" + "_xsecs.dat").c_str());
	  in5.open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + "w_second_jet_eta_bb" + "_xsecs.dat").c_str());
	}
	in1 >> sum1; in1 >> sum1_b;
	in2 >> sum2; in2 >> sum2_b;
	in3 >> sum3; in3 >> sum3_b;
	in4 >> sum4; in4 >> sum4_b;
	in5 >> sum5; in5 >> sum5_b;
	in1.close();
	in2.close();
	in3.close();
	in4.close();
	in5.close();

	float tot = (sum1+sum2+sum3+sum4+sum5)/5.;
	float tot_b = (sum1_b+sum2_b+sum3_b+sum4_b+sum5_b)/5.;
	float rms = TMath::Sqrt((TMath::Power(sum1-tot,2)+TMath::Power(sum2-tot,2)+TMath::Power(sum3-tot,2)+TMath::Power(sum4-tot,2)+TMath::Power(sum5-tot,2))/(5-1));
	float rms_b = TMath::Sqrt((TMath::Power(sum1_b-tot_b,2)+TMath::Power(sum2_b-tot_b,2)+TMath::Power(sum3_b-tot_b,2)+TMath::Power(sum4_b-tot_b,2)+TMath::Power(sum5_b-tot,2))/(5-1));

	if (isratio==1) {
	  float tmp1 = (tot_b/tot);
	  float tmp2 = (tot_b/tot)*TMath::Sqrt(TMath::Power(rms/tot,2)+TMath::Power(rms_b/tot_b,2));
	  tot_b = tmp1;
	  rms_b = tmp2;
	}

	TH1F* h_data_stat = (TH1F*)h_data->Clone();
	TH1F* h_data_b_stat = (TH1F*)h_data_b->Clone();
	TH1F* h_data_syst = (TH1F*)h_data->Clone();
	TH1F* h_data_b_syst = (TH1F*)h_data_b->Clone();
	TH1F* h_data_tot = (TH1F*)h_data->Clone();
	TH1F* h_data_b_tot = (TH1F*)h_data_b->Clone();

	for (int i=0;i<=h_data_stat->GetNbinsX()+1;i++) {
	  h_data_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(stat_top->GetBinError(i),2)));
	  h_data_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(stat_bfit->GetBinError(i),2)));
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_eff->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_jec->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_jer->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_pu->GetBinError(i),2));
	  if (useSysDR) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_dr->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_bkg->GetBinError(i),2));
	  if (useSysBfit2) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_bfit2->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_btag->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_lumi->GetBinError(i),2));
	  if (useSysRMS) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_stat->GetBinContent(i)*rms/tot,2));
	  h_data_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(h_data_syst->GetBinError(i),2));
	  h_data_tot->SetBinError(i, val);
	}

	for (int i=0;i<=h_data_b_stat->GetNbinsX()+1;i++) {
	  h_data_b_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(stat_b_top->GetBinError(i),2)));
	  h_data_b_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(stat_b_bfit->GetBinError(i),2)));
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_eff->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_jec->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_jer->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_pu->GetBinError(i),2));
	  if (useSysDR) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_dr->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_bkg->GetBinError(i),2));
	  if (useSysBfit2) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_bfit2->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_btag->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_lumi->GetBinError(i),2));
	  if (useSysRMS) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_b_stat->GetBinContent(i)*rms_b/tot_b,2));
	  h_data_b_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(h_data_b_syst->GetBinError(i),2));
	  h_data_b_tot->SetBinError(i, val);
	}

	h_mcg->Scale(1./Lumi2012, "width");
	h_mcg_b->Scale(1./Lumi2012, "width");
	if (isratio==1) {
	  h_mcg_b->Divide(h_mcg);
	  h_mcg_b->Scale(100.);
	}

/*
	h_data = fixrange(h_data);
	h_data_b = fixrange(h_data_b);
	for (int i=0;i<NMAX;i++) {
	  if (h_data_scan[i]) h_data_scan[i] = fixrange(h_data_scan[i]);
	  if (h_data_b_scan[i]) h_data_b_scan[i] = fixrange(h_data_b_scan[i]);
	}
	h_data_stat = fixrange(h_data_stat);
	h_data_b_stat = fixrange(h_data_b_stat);
	h_data_syst = fixrange(h_data_syst);
	h_data_b_syst = fixrange(h_data_b_syst);
	h_data_tot = fixrange(h_data_tot);
	h_data_b_tot = fixrange(h_data_b_tot);

	h_mc1 = fixrange(h_mc1);
	h_mc1b_b = fixrange(h_mc1b_b);
	h_mcg = fixrange(h_mcg);
	h_mcg_b = fixrange(h_mcg_b);
*/

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	h_mc1b_b->SetTitle("");
	if (isratio==1) {
	  h_mc1b_b->GetYaxis()->SetTitle("#sigma(W+2b) / #sigma(W+1b) [%]");
	} else {
	  h_mc1b_b->GetYaxis()->SetTitle("#sigma [pb]");
	}

	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	pad1->SetTopMargin(0.115);
	pad1->SetBottomMargin(0.0001);
	pad1->Draw();
	pad1->cd();

	h_mc1b_b->SetLineColor(kRed);
	h_mc1b_b->SetLineWidth(2);
	h_mc1b_b->SetMarkerColor(kRed);
	h_mc1b_b->SetFillColor(kRed);
	h_mc1b_b->SetStats(0);
	if (isratio==1) {
	  h_mc1b_b->Draw("E5");
	}

	h_mcg_b->SetLineColor(kGreen+2);
	h_mcg_b->SetLineWidth(2);
	h_mcg_b->SetFillColor(kGreen+2);
	h_mcg_b->SetMarkerColor(kGreen+2);
	if (isratio==1) {
	  h_mcg_b->Draw("E5SAME");
	}

	if (isratio==1) {
	  h_data_b_tot->GetYaxis()->SetTitle("#sigma_{W+b-jets}/#sigma_{W+jets} [%]");
	}
	h_data_b_tot->SetMarkerColor(kRed+1);
	h_data_b_tot->SetLineColor(kRed+1);
	h_data_b_tot->SetMarkerStyle(24);
	h_data_b_tot->SetMarkerSize(0.7);
	h_data_b_tot->SetStats(0);
	h_data_b_stat->GetXaxis()->SetTitleOffset(0.7);
	h_data_b_stat->SetMarkerColor(kBlack);
	h_data_b_stat->SetLineColor(kBlack);
	h_data_b_stat->SetMarkerStyle(24);
	h_data_b_stat->SetMarkerSize(0.7);
	h_data_b_stat->SetStats(0);
	if (isratio==1) {
	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");
	}

	TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (isratio==0) {
	  pad1->SetLogy();

	  h_mc1b_b->SetMaximum(4*h_data_tot->GetMaximum());
	  h_mc1b_b->SetMinimum(TMath::Max(0.002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));
	  if (title.find("_mt")!=string::npos) h_mc1b_b->SetMinimum(TMath::Max(0.00002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));
	  if (title.find("_pt")!=string::npos) h_mc1b_b->SetMinimum(TMath::Max(0.000002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));
	  if (title.find("_mass")!=string::npos) h_mc1b_b->SetMinimum(TMath::Max(0.00002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));

	  h_mc1b_b->Draw("E5");
	  TH1F* tmp1 = (TH1F*)h_mc1b_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp1->GetMinimum()==0) tmp1->GetXaxis()->SetRangeUser(0, tmp1->GetBinCenter(tmp1->GetMinimumBin()-1));
	  }
	  tmp1->SetFillColor(0);
	  tmp1->DrawClone("HISTLSAME");

	  h_mcg_b->Draw("E5SAME");
	  TH1F* tmp2 = (TH1F*)h_mcg_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2->GetMinimum()==0) tmp2->GetXaxis()->SetRangeUser(0, tmp2->GetBinCenter(tmp2->GetMinimumBin()-1));
	  }
	  tmp2->SetFillColor(0);
	  tmp2->DrawClone("HISTLSAME");

	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");

	  h_mc1->SetLineColor(kRed);
	  h_mc1->SetLineWidth(2);
	  h_mc1->SetMarkerColor(kRed);
	  h_mc1->SetFillColor(kRed);
	  if (drawInclusive) h_mc1->Draw("E5SAME");
	  TH1F* tmp3 = (TH1F*)h_mc1->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp3->GetMinimum()==0) tmp3->GetXaxis()->SetRangeUser(0, tmp3->GetBinCenter(tmp3->GetMinimumBin()-1));
	  }
	  tmp3->SetFillColor(0);
	  if (drawInclusive) tmp3->DrawClone("HISTLSAME");

	  h_mcg->SetLineColor(kGreen+2);
	  h_mcg->SetLineWidth(2);
	  h_mcg->SetMarkerColor(kGreen+2);
	  h_mcg->SetFillColor(kGreen+2);
	  if (drawInclusive) h_mcg->Draw("E5SAME");
	  TH1F* tmp4 = (TH1F*)h_mcg->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4->GetMinimum()==0) tmp4->GetXaxis()->SetRangeUser(0, tmp4->GetBinCenter(tmp4->GetMinimumBin()-1));
	  }
	  tmp4->SetFillColor(0);
	  if (drawInclusive) tmp4->DrawClone("HISTLSAME");

	  h_data_tot->SetMarkerColor(kRed+1);
	  h_data_tot->SetLineColor(kRed+1);
	  h_data_tot->SetMarkerStyle(20);
	  h_data_tot->SetMarkerSize (0.7);
	  h_data_stat->SetLineColor(kBlack);
	  h_data_stat->SetMarkerColor(kBlack);
	  h_data_stat->SetMarkerStyle(20);
	  h_data_stat->SetMarkerSize (0.7);
	  if (drawInclusive) h_data_tot->Draw("E1PX0SAME");
	  if (drawInclusive) h_data_stat->Draw("E1PX0SAME");

	  if (ilepton==1) {
	    if (drawInclusive) leg->AddEntry(h_data_stat,"W(#rightarrow e#nu)+1b DATA","p");
	    leg->AddEntry(h_data_b_stat,"W(#rightarrow e#nu)+2b DATA","p");
	    if (useMC) leg->AddEntry(h_mc1,"W(#rightarrow e#nu) MC","l");
	    leg->AddEntry(h_mcg,"W(#rightarrow e#nu) MadGraph","l");
	  }
	  if (ilepton==2){
	    if (drawInclusive) leg->AddEntry(h_data_stat,"W(#rightarrow #mu#nu)+1b DATA","p");
	    leg->AddEntry(h_data_b_stat,"W(#rightarrow #mu#nu)+2b DATA","p");
	    if (useMC) leg->AddEntry(h_mc1,"W(#rightarrow #mu#nu) MC","l");
	    leg->AddEntry(h_mcg,"W(#rightarrow #mu#nu) MadGraph","l");
	  }
	}

	if (isratio==1) {
	  if (ilepton==1) {
	    leg->AddEntry(h_data_b_stat,"W(#rightarrow e#nu) DATA","p");
	    if (useMC) leg->AddEntry(h_mc1b_b,"W(#rightarrow e#nu) MC","l");
	    leg->AddEntry(h_mcg_b,"W(#rightarrow e#nu) MadGraph","l");
	  }
	  if (ilepton==2){
	    leg->AddEntry(h_data_b_stat,"W(#rightarrow #mu#nu) DATA","p");
	    if (useMC) leg->AddEntry(h_mc1b_b,"W(#rightarrow #mu#nu) MC","l");
	    leg->AddEntry(h_mcg_b,"W(#rightarrow #mu#nu) MadGraph","l");
	  }
	}

	leg->Draw();

	c1->cd();

 	TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
	latexLabel->Draw("same");

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();

	TH1F *h_M_tot = (TH1F*)h_data_b_tot->Clone();
	TH1F *h_M_stat = (TH1F*)h_data_b_stat->Clone();

	h_M_tot->Divide(h_mcg_b);
	h_M_stat->Divide(h_mcg_b);

	h_M_tot->SetTitle("");
	h_M_tot->SetStats(0);
	h_M_tot->GetXaxis()->SetTitleOffset(0.9);
	h_M_tot->GetXaxis()->SetTitleSize(0.1);
	h_M_tot->GetXaxis()->SetLabelFont(42);
	h_M_tot->GetXaxis()->SetLabelSize(0.08);
	h_M_tot->GetXaxis()->SetTitleFont(42);
	h_M_tot->GetYaxis()->SetTitle("Data / Theory");
	h_M_tot->GetYaxis()->SetNdivisions(013);
	h_M_tot->GetYaxis()->SetTitleSize(0.09);
	h_M_tot->GetYaxis()->SetLabelSize(0.08);
	h_M_tot->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_M_tot->GetYaxis()->SetTitleOffset(0.4);

	h_M_tot->SetMarkerStyle(24);
	h_M_tot->Draw("E1PX0");
	h_M_stat->SetMarkerStyle(24);
	h_M_stat->Draw("E1PX0SAME");

	if (isratio==0) {
	  TH1F *h_M2_tot= (TH1F*)h_data_tot->Clone();
	  TH1F *h_M2_stat= (TH1F*)h_data_stat->Clone();

	  h_M2_tot->Divide(h_mcg);
	  h_M2_stat->Divide(h_mcg);

	  TGraphErrors *g_M2_tot = new TGraphErrors(h_M2_tot);
	  TGraphErrors *g_M2_stat = new TGraphErrors(h_M2_stat);

	  float dx = 0.1*(g_M2_tot->GetXaxis()->GetXmax()-g_M2_tot->GetXaxis()->GetXmin())/g_M2_tot->GetN();
	  for (int i=0; i<g_M2_tot->GetN(); i++) {
	    g_M2_stat->SetPoint(i, g_M2_stat->GetX()[i]-dx, g_M2_stat->GetY()[i]);
	    g_M2_stat->SetPointError(i, 0, g_M2_stat->GetEY()[i]);
	    g_M2_tot->SetPoint(i, g_M2_tot->GetX()[i]-dx, g_M2_tot->GetY()[i]);
	    g_M2_tot->SetPointError(i, 0, g_M2_tot->GetEY()[i]);
	  }

	  g_M2_tot->SetMarkerStyle(20);
	  g_M2_tot->Draw("EP0SAME");
	  g_M2_stat->SetMarkerStyle(20);
	  g_M2_stat->Draw("EP0SAME");
	}

	TLatex *t2 = new TLatex();
	t2->SetTextSize(0.09);
	t2->SetTextFont(42);
	t2->SetLineWidth(2);
	t2->SetNDC();
	//t2->DrawLatex(0.15,0.9,"MadGraph");

	TLine *OLine2 = new TLine(h_M_tot->GetXaxis()->GetXmin(),1.,h_M_tot->GetXaxis()->GetXmax(),1.);
	OLine2->SetLineColor(kGreen+2);
	OLine2->SetLineWidth(2);
	OLine2->Draw();

	c1->cd();

	if (isratio==1) {
	  h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 10);
	}

	if (title_b=="w_first_jet_pt_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(W+2b) / #sigma(W+1b)] / dp_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_first_jet_eta_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading jet #eta");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(W+2b) / #sigma(W+1b)] / d#eta [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_pt") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(W+2b) / #sigma(W+1b)] / dp^{b}_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_first_bjet_eta") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading b-jet #eta");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(W+2b) / #sigma(W+1b)] / d#eta^{b} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	}

	if (plot) {
	  ofstream out, out1;
	  if (isratio==0) {
	    if (ilepton==1) {
	      //gSystem->mkdir((path + "/electrons/" + version + "/xsecs_unfolding/").c_str(), kTRUE);
	      //c1->SaveAs((path + "/electrons/" + version + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      //out.open((path + "/electrons/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.dat").c_str());
	      //out1.open((path + "/electrons/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.txt").c_str());
	      gSystem->mkdir((path + "/electrons/" + version + "/xsecs/").c_str(), kTRUE);
	      c1->SaveAs((path + "/electrons/" + version + "/xsecs/" + title_b + "_xsecs.pdf").c_str());
	      out.open((path + "/electrons/" + version + "/" + "/xsecs/" + title_b + "_xsecs.dat").c_str());
	      out1.open((path + "/electrons/" + version + "/" + "/xsecs/" + title_b + "_xsecs.txt").c_str());
	    }
	    if (ilepton==2) {
	      //gSystem->mkdir((path + "/muons/" + version + "/xsecs_unfolding/").c_str(), kTRUE);
	      //c1->SaveAs((path + "/muons/" + version + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      //out.open((path + "/muons/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.dat").c_str());
	      //out1.open((path + "/muons/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.txt").c_str());
	      gSystem->mkdir((path + "/muons/" + version + "/xsecs/").c_str(), kTRUE);
	      c1->SaveAs((path + "/muons/" + version + "/xsecs/" + title_b + "_xsecs.pdf").c_str());
	      out.open((path + "/muons/" + version + "/" + "/xsecs/" + title_b + "_xsecs.dat").c_str());
	      out1.open((path + "/muons/" + version + "/" + "/xsecs/" + title_b + "_xsecs.txt").c_str());
	    }
	  }
	  if (isratio==1) {
	    if (ilepton==1) {
	      //gSystem->mkdir((path + "/electrons/" + version + "/ratios_unfolding/").c_str(), kTRUE);
	      //c1->SaveAs((path + "/electrons/" + version + "/ratios_unfolding/" + title_b + "_ratio_unfolding.pdf").c_str());
	      //out.open((path + "/electrons/" + version + "/" + "/ratios_unfolding/" + title_b + "_ratio_unfolding.dat").c_str());
	      //out1.open((path + "/electrons/" + version + "/" + "/ratios_unfolding/" + title_b + "_ratio_unfolding.txt").c_str());
	      gSystem->mkdir((path + "/electrons/" + version + "/ratios/").c_str(), kTRUE);
	      c1->SaveAs((path + "/electrons/" + version + "/ratios/" + title_b + "_ratio.pdf").c_str());
	      out.open((path + "/electrons/" + version + "/" + "/ratios/" + title_b + "_ratio.dat").c_str());
	      out1.open((path + "/electrons/" + version + "/" + "/ratios/" + title_b + "_ratio.txt").c_str());
	    }
	    if (ilepton==2) {
	      //gSystem->mkdir((path + "/muons/" + version + "/ratios_unfolding/").c_str(), kTRUE);
	      //c1->SaveAs((path + "/muons/" + version + "/ratios_unfolding/" + title_b + "_ratio_unfolding.pdf").c_str());
	      //out.open((path + "/muons/" + version + "/" + "/ratios_unfolding/" + title_b + "_ratio_unfolding.dat").c_str());
	      //out1.open((path + "/muons/" + version + "/" + "/ratios_unfolding/" + title_b + "_ratio_unfolding.txt").c_str());
	      gSystem->mkdir((path + "/muons/" + version + "/ratios/").c_str(), kTRUE);
	      c1->SaveAs((path + "/muons/" + version + "/ratios/" + title_b + "_ratio.pdf").c_str());
	      out.open((path + "/muons/" + version + "/" + "/ratios/" + title_b + "_ratio.dat").c_str());
	      out1.open((path + "/muons/" + version + "/" + "/ratios/" + title_b + "_ratio.txt").c_str());
	    }
	  }
	  if (isratio==0) {
	    out << h_data->GetName();
	    out << std::fixed << std::setprecision(4);
	    out << " : average unfolded total cross section = " << tot << " +- " << rms << " pb (" << 100*(rms/tot) << " %)";
	    out << endl;
	    out << std::setw(25) << "data";
	    out << std::setw(12) << "bkg";
	    out << std::setw(12) << "eff";
	    out << std::setw(12) << "jec";
	    out << std::setw(12) << "jer";
	    out << std::setw(12) << "pu";
	    if (useSysDR) out << std::setw(12) << "DR";
	    out << std::setw(12) << "bkg";
	    out << std::setw(12) << "ttbar";
	    out << std::setw(12) << "bfit";
	    if (useSysBfit2) out << std::setw(12) << "bfit2";
	    out << std::setw(12) << "btag";
	    out << std::setw(12) << "unfold";
	    out << std::setw(12) << "unfold";
	    out << std::setw(12) << "lumi";
	    if (useSysRMS) out << std::setw(12) << "unfold";
	    out << std::setw(12) << "total";
	    out << std::setw(12) << "total";
	    out << std::setw(12) << "total";
	    out << endl;
	    out << std::setw(25) << "stat";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    if (useSysDR) out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "stat";
	    if (useSysBfit2) out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    if (useSysRMS) out << std::setw(12) << "rms";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "error";
	    out << std::setw(8) << "%";
	    out << endl;
	    for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      out << std::fixed;
	      out << std::setw(2) << i;
	      out << " ";
	      out << std::setprecision(6);
	      out << std::setw(10) << h_data->GetBinContent(i);
	      out << " +- ";
	      out << std::setw(8) << h_data->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_bkg->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_eff->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_jec->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_jer->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_pu->GetBinError(i);
	      if (useSysDR) {
	        out << " +- ";
	        out << std::setw(8) << syst_dr->GetBinError(i);
	      }
	      out << " +- ";
	      out << std::setw(8) << syst_bkg->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_top->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_bfit->GetBinError(i);
	      if (useSysBfit2) {
	        out << " +- ";
	        out << std::setw(8) << syst_bfit2->GetBinError(i);
	      }
	      out << " +- ";
	      out << std::setw(8) << syst_btag->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_lumi->GetBinError(i);
	      if (useSysRMS) {
	        out << " +- ";
	        out << std::setw(8) << h_data->GetBinContent(i)*(rms/tot);
	      }
	      out << " => ";
	      out << std::setw(8) << h_data_stat->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << h_data_syst->GetBinError(i);
	      out << " => ";
	      out << std::setw(8) << h_data_tot->GetBinError(i);
	      out << " => ";
	      out << std::setprecision(1);
	      out << std::setw(4) << 100.*(h_data_stat->GetBinContent(i)==0 ? 0 : h_data_tot->GetBinError(i)/h_data_stat->GetBinContent(i));
	      out << endl;
	    }
	  }
	  out << h_data_b->GetName();
	  if (isratio==0) {
	    out << std::fixed << std::setprecision(4);
	    out << " : average unfolded total cross section = " << tot_b << " +- " << rms_b << " pb (" << 100*(rms_b/tot_b) << " %)";
	  }
	  out << endl;
	  out << std::setw(25) << "data";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "eff";
	  out << std::setw(12) << "jec";
	  out << std::setw(12) << "jer";
	  out << std::setw(12) << "pu";
	  if (useSysDR) out << std::setw(12) << "DR";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "ttbar";
	  out << std::setw(12) << "bfit";
	  if (useSysBfit2) out << std::setw(12) << "bfit2";
	  out << std::setw(12) << "btag";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "lumi";
	  if (useSysRMS) out << std::setw(12) << "unfold";
	  out << std::setw(12) << "total";
	  out << std::setw(12) << "total";
	  out << std::setw(12) << "total";
	  out << endl;
	  out << std::setw(25) << "stat";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  if (useSysDR) out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "stat";
	  if (useSysBfit2) out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  if (useSysRMS) out << std::setw(12) << "rms";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "error";
	  out << std::setw(8) << "%";
	  out << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    out << std::fixed;
	    out << std::setw(2) << i;
	    out << " ";
	    out << std::setprecision(6);
	    out << std::setw(10) << h_data_b->GetBinContent(i);
	    out << " +- ";
	    out << std::setw(8) << h_data_b->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_eff->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_jec->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_jer->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_pu->GetBinError(i);
	    if (useSysDR) {
	      out << " +- ";
	      out << std::setw(8) << syst_b_dr->GetBinError(i);
	    }
	    out << " +- ";
	    out << std::setw(8) << syst_b_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_top->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_bfit->GetBinError(i);
	    if (useSysBfit2) {
	      out << " +- ";
	      out << std::setw(8) << syst_b_bfit2->GetBinError(i);
	    }
	    out << " +- ";
	    out << std::setw(8) << syst_b_btag->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_lumi->GetBinError(i);
	    if (useSysRMS) {
	      out << " +- ";
	      out << std::setw(8) << h_data_b->GetBinContent(i)*(rms_b/tot_b);
	    }
	    out << " => ";
	    out << std::setw(8) << h_data_b_stat->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << h_data_b_syst->GetBinError(i);
	    out << " => ";
	    out << std::setw(8) << h_data_b_tot->GetBinError(i);
	    out << " => ";
	    out << std::setprecision(1);
	    out << std::setw(4) << 100.*(h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    out << endl;
	  }
	  out.close();
	  if (isratio==0) {
	    out1 << h_data->GetName() << " - RELATIVE ERRORS";
	    out1 << endl;
	    out1 << std::setw(7) << "data";
	    out1 << std::setw(8) << "bkg";
	    out1 << std::setw(8) << "eff";
	    out1 << std::setw(8) << "jec";
	    out1 << std::setw(8) << "jer";
	    out1 << std::setw(8) << "pu";
	    if (useSysDR) out1 << std::setw(8) << "DR";
	    out1 << std::setw(8) << "bkg";
	    out1 << std::setw(8) << "ttbar";
	    out1 << std::setw(8) << "bfit";
	    if (useSysBfit2) out1 << std::setw(8) << "bfit2";
	    out1 << std::setw(8) << "btag";
	    out1 << std::setw(8) << "unfold";
	    out1 << std::setw(8) << "unfold";
	    out1 << std::setw(8) << "lumi";
	    if (useSysRMS) out1 << std::setw(8) << "unfold";
	    out1 << std::setw(8) << "total";
	    out1 << std::setw(8) << "total";
	    out1 << std::setw(8) << "total";
	    out1 << endl;
	    out1 << std::setw(7) << "stat";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    if (useSysDR) out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "stat";
	    if (useSysBfit2) out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    if (useSysRMS) out1 << std::setw(8) << "rms";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "error";
	    out1 << endl;
	    for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      double val = 100.*(h_data->GetBinContent(i)==0 ? 0 : 1./h_data->GetBinContent(i));
	      out1 << std::fixed;
	      out1 << std::setw(2) << i;
	      out1 << " ";
	      out1 << std::setprecision(1);
	      out1 << std::setw(4) << h_data->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_bkg->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_eff->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_jec->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_jer->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_pu->GetBinError(i)*val;
	      if (useSysDR) {
	        out1 << " +- ";
	        out1 << std::setw(4) << syst_dr->GetBinError(i)*val;
	      }
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_bkg->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_top->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_bfit->GetBinError(i)*val;
	      if (useSysBfit2) {
	        out1 << " +- ";
	        out1 << std::setw(4) << syst_bfit2->GetBinError(i)*val;
	      }
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_btag->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_lumi->GetBinError(i)*val;
	      if (useSysRMS) {
	        out1 << " +- ";
	        out1 << std::setw(4) << h_data->GetBinContent(i)*(rms/tot)*val;
	      }
	      out1 << " => ";
	      out1 << std::setw(4) << h_data_stat->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << h_data_syst->GetBinError(i)*val;
	      out1 << " => ";
	      out1 << std::setw(4) << h_data_tot->GetBinError(i)*val;
	      out1 << endl;
	    }
	  }
	  out1 << h_data_b->GetName() << " - RELATIVE ERRORS";
	  out1 << endl;
	  out1 << std::setw(7) << "data";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "eff";
	  out1 << std::setw(8) << "jec";
	  out1 << std::setw(8) << "jer";
	  out1 << std::setw(8) << "pu";
	  if (useSysDR) out1 << std::setw(8) << "DR";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "ttbar";
	  out1 << std::setw(8) << "bfit";
	  if (useSysBfit2) out1 << std::setw(8) << "bfit2";
	  out1 << std::setw(8) << "btag";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "lumi";
	  if (useSysRMS) out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "total";
	  out1 << std::setw(8) << "total";
	  out1 << std::setw(8) << "total";
	  out1 << endl;
	  out1 << std::setw(7) << "stat";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  if (useSysDR) out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "stat";
	  if (useSysBfit2) out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  if (useSysRMS) out1 << std::setw(8) << "rms";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "error";
	  out1 << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    double val = 100.*(h_data_b->GetBinContent(i)==0 ? 0 : 1./h_data_b->GetBinContent(i));
	    out1 << std::fixed;
	    out1 << std::setw(2) << i;
	    out1 << " ";
	    out1 << std::setprecision(1);
	    out1 << std::setw(4) << h_data_b->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_eff->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_jec->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_jer->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_pu->GetBinError(i)*val;
	    if (useSysDR) {
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_b_dr->GetBinError(i)*val;
	    }
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_top->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_bfit->GetBinError(i)*val;
	    if (useSysBfit2) {
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_b_bfit2->GetBinError(i)*val;
	    }
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_btag->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_lumi->GetBinError(i)*val;
	    if (useSysRMS) {
	      out1 << " +- ";
	      out1 << std::setw(4) << h_data_b->GetBinContent(i)*(rms_b/tot_b)*val;
	    }
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_b_stat->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << h_data_b_syst->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_b_tot->GetBinError(i)*val;
	    out1 << endl;
	  }
	  out1.close();
	}
}

