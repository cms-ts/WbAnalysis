#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v09.h"

#include "fixrange.C"
#include "rebin.C"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

int unfold=0; // use pre-unfolding distributions
//int unfold=1; // use unfolded distributions

TH1F* read(string subdir, string title, int ilepton, TFile* infile=0) {
  TH1F* hist=0;
  TFile* file = infile;
  string title_tmp = title;
  if (ilepton==1) {
    if (title=="w_mt_b") title_tmp="w_mt_wenu_b";
    if (title=="w_mt_bb") title_tmp="w_mt_wenu_bb";
    string title_tmp2 = title_tmp;
    if (!unfold) {
      if (title_tmp2.find("_bb")==string::npos) {
        title_tmp2 = title_tmp2 + "b";
      }
    }
    if (file) {
      file->cd("demoEleGen");
    } else {
      if (unfold) {
        file = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title_tmp2 + "_unfolding.root").c_str());
      } else {
        file = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + title_tmp2 + "_xsecs.root").c_str());
      }
    }
  }
  if (ilepton==2) {
    if (title=="w_mt_b") title_tmp="w_mt_wmnu_b";
    if (title=="w_mt_bb") title_tmp="w_mt_wmnu_bb";
    string title_tmp2 = title_tmp;
    if (!unfold) {
      if (title_tmp2.find("_bb")==string::npos) {
        title_tmp2 = title_tmp2 + "b";
      }
    }
    if (file) {
      file->cd("demoMuoGen");
    } else {
      if (unfold) {
        file = TFile::Open((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title_tmp2 + "_unfolding.root").c_str());
      } else {
        file = TFile::Open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + title_tmp2 + "_xsecs.root").c_str());
      }
    }
  }
  hist = (TH1F*)gDirectory->Get(title_tmp.c_str())->Clone();
  hist->SetDirectory(0);
  if (!infile) file->Close();
  return hist;
}

double calc(int iflag, double cont1, double cont2, double stat1, double stat2, double stat_bkg1, double stat_bkg2, double syst_eff1, double syst_eff2, double syst_jec1, double syst_jec2, double syst_jer1, double syst_jer2, double syst_pu1, double syst_pu2, double syst_bkg1, double syst_bkg2, double stat_top1, double stat_top2, double stat_bfit1, double stat_bfit2, double syst_btag1, double syst_btag2, double stat_unfold1, double stat_unfold2, double syst_unfold1, double syst_unfold2, double syst_lumi1, double syst_lumi2) {
  double val = 0.0;

  if (iflag == 0) {
    if (cont1*cont2 != 0) {
      val = (cont1/(TMath::Power(stat1,2)+TMath::Power(stat_bkg1,2)+TMath::Power(syst_eff1,2)+TMath::Power(stat_top1,2)+TMath::Power(stat_bfit1,2)+TMath::Power(stat_unfold1,2))+cont2/(TMath::Power(stat2,2)+TMath::Power(stat_bkg2,2)+TMath::Power(syst_eff2,2)+TMath::Power(stat_top2,2)+TMath::Power(stat_bfit2,2)+TMath::Power(stat_unfold2,2)))/(1./(TMath::Power(stat1,2)+TMath::Power(stat_bkg1,2)+TMath::Power(syst_eff1,2)+TMath::Power(stat_top1,2)+TMath::Power(stat_bfit1,2)+TMath::Power(stat_unfold1,2))+1./(TMath::Power(stat2,2)+TMath::Power(stat_bkg2,2)+TMath::Power(syst_eff2,2)+TMath::Power(stat_top2,2)+TMath::Power(stat_bfit2,2)+TMath::Power(stat_unfold2,2)));
    }
  }

  if (iflag == 1) {
    if (cont1*cont2 != 0) {
      val = TMath::Sqrt(TMath::Power((syst_jec1+syst_jec2)/2.,2)+TMath::Power((syst_jer1+syst_jer2)/2.,2)+TMath::Power((syst_pu1+syst_pu2)/2.,2)+TMath::Power((syst_bkg1+syst_bkg2)/2.,2)+TMath::Power((syst_btag1+syst_btag2)/2.,2)+TMath::Power((syst_unfold1+syst_unfold2)/2.,2)+TMath::Power((syst_lumi1+syst_lumi2)/2.,2));
      val = TMath::Sqrt(TMath::Power(val,2)+1./(1./(TMath::Power(stat1,2)+TMath::Power(stat_bkg1,2)+TMath::Power(syst_eff1,2)+TMath::Power(stat_top1,2)+TMath::Power(stat_bfit1,2)+TMath::Power(stat_unfold1,2))+1./(TMath::Power(stat2,2)+TMath::Power(stat_bkg2,2)+TMath::Power(syst_eff2,2)+TMath::Power(stat_top2,2)+TMath::Power(stat_bfit2,2)+TMath::Power(stat_unfold2,2))));
    }
  }

  if (TMath::IsNaN(val)) val = 0.0;

  return val;
}

void DataMCComp8(string title="", int plot=0) {

//int drawInclusive = 0; // do not plot the "inclusive" histogram
int drawInclusive = 1; // do plot the "inclusive" histogram

string subdir="0";

	if (gROOT->GetVersionInt() >= 53401) {
	  gROOT->GetColor(kRed)->SetAlpha(0.5);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	  gROOT->GetColor(kOrange+7)->SetAlpha(0.5);
	}

	double Lumi2012=0;

	Lumi2012 = (Lumi2012_ele+Lumi2012_muon)/2.;

	double norm1 = ((Lumi2012 * Xsec_wj) / Ngen_wj);

	if (title.empty()) title = "w_jetmultiplicity";

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

	TH1F* w_data[2];
	TH1F* w_data_b[2];
	for (int i=0; i<2; i++) {
	  w_data[i] = read(subdir, title, i+1);
	  w_data_b[i] = read(subdir, title_b, i+1);
	  if (unfold) {
	    w_data[i]->Scale(1./Lumi2012, "width");
	    w_data_b[i]->Scale(1./Lumi2012, "width");
	  }
	}

	TH1F* w_mcg[2];
	TH1F* w_mcg_b[2];
	for (int i=0; i<2; i++) {
	  w_mcg[i] = read(subdir, title, i+1, mcg);
	  w_mcg_b[i] = read(subdir, title_b, i+1, mcg);
	}

	TH1F* h_mcg = (TH1F*)w_mcg[0]->Clone();
	TH1F* h_mcg_b = (TH1F*)w_mcg_b[0]->Clone();

	h_mcg->Sumw2();

	h_mcg_b->Sumw2();

	for (int i=0;i<=h_mcg->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (w_mcg[0]->GetBinContent(i)*w_mcg[1]->GetBinContent(i) != 0) {
	    val = (w_mcg[0]->GetBinContent(i)/TMath::Power(w_mcg[0]->GetBinError(i),2)+w_mcg[1]->GetBinContent(i)/TMath::Power(w_mcg[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg[0]->GetBinError(i),2)+1./TMath::Power(w_mcg[1]->GetBinError(i),2));
	    h_mcg->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg[0]->GetBinError(i),2)+1./TMath::Power(w_mcg[1]->GetBinError(i),2)));
	    h_mcg->SetBinError(i, val);
	  }
	  if (w_mcg_b[0]->GetBinContent(i)*w_mcg_b[1]->GetBinContent(i) != 0) {
	    val = (w_mcg_b[0]->GetBinContent(i)/TMath::Power(w_mcg_b[0]->GetBinError(i),2)+w_mcg_b[1]->GetBinContent(i)/TMath::Power(w_mcg_b[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg_b[1]->GetBinError(i),2));
	    h_mcg_b->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg_b[1]->GetBinError(i),2)));
	    h_mcg_b->SetBinError(i, val);
	  }
	}

	h_mcg->Scale(norm1);

	h_mcg_b->Scale(norm1);

	h_mcg->Scale(1./Lumi2012, "width");
	h_mcg_b->Scale(1./Lumi2012, "width");

	TH1F* w_stat_bkg[2];
	TH1F* w_stat_b_bkg[2];

	TH1F* w_syst_eff[2];
	TH1F* w_syst_b_eff[2];

	TH1F* w_syst_jec[2];
	TH1F* w_syst_b_jec[2];

	TH1F* w_syst_jer[2];
	TH1F* w_syst_b_jer[2];

	TH1F* w_syst_pu[2];
	TH1F* w_syst_b_pu[2];

	TH1F* w_syst_bkg[2];
	TH1F* w_syst_b_bkg[2];

	TH1F* w_stat_top[2];
	TH1F* w_stat_b_top[2];

	TH1F* w_stat_bfit[2];
	TH1F* w_stat_b_bfit[2];

	TH1F* w_syst_btag[2];
	TH1F* w_syst_b_btag[2];

	TH1F* w_stat_unfold[2];
	TH1F* w_stat_b_unfold[2];

	TH1F* w_syst_unfold[2];
	TH1F* w_syst_b_unfold[2];

	TH1F* w_syst_lumi[2];
	TH1F* w_syst_b_lumi[2];

	TH1F* w_stat_tot[2];
	TH1F* w_stat_b_tot[2];

	TH1F* w_syst_tot[2];
	TH1F* w_syst_b_tot[2];

	for (int i=0; i<2; i++) {

	  w_stat_bkg[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_bkg[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_eff[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_eff[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_jec[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_jec[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_jer[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_jer[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_pu[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_pu[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_bkg[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_bkg[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_top[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_top[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_bfit[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_bfit[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_btag[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_btag[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_unfold[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_unfold[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_unfold[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_unfold[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_lumi[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_lumi[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_tot[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_tot[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_tot[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_tot[i] = (TH1F*)w_data_b[0]->Clone();

	  ifstream in;
	  string title_b_tmp = title_b;
	  if (i==0) {
	    if (title_b=="w_mt_b") title_b_tmp="w_mt_wenu_b";
	    if (title_b=="w_mt_bb") title_b_tmp="w_mt_wenu_bb";
	    if (unfold) {
	      in.open((path + "/electrons/" + version + "/" + "/xsecs_unfolding/" + title_b_tmp + "_xsecs_unfolding.dat").c_str());
	    } else {
	      in.open((path + "/electrons/" + version + "/" + "/xsecs/" + title_b_tmp + "_xsecs.dat").c_str());
	    }
	  }
	  if (i==1) {
	    if (title_b=="w_mt_b") title_b_tmp="w_mt_wmnu_b";
	    if (title_b=="w_mt_bb") title_b_tmp="w_mt_wmnu_bb";
	    if (unfold) {
              in.open((path + "/muons/" + version + "/" + "/xsecs_unfolding/" + title_b_tmp + "_xsecs_unfolding.dat").c_str());
	    } else {
              in.open((path + "/muons/" + version + "/" + "/xsecs/" + title_b_tmp + "_xsecs.dat").c_str());
	    }
	  }

	  string tmp;

	  getline(in, tmp);
	  getline(in, tmp);
	  getline(in, tmp);
	  for (int j=0; j<w_data[0]->GetNbinsX()+2; j++) {
	    in >> tmp;
	    double val = 0.0;
	    in >> val; w_data[i]->SetBinContent(j, val); in >> tmp;
	    in >> val; w_data[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_bkg[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_eff[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_jec[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_jer[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_pu[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_bkg[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_top[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_bfit[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_btag[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_lumi[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; in >> tmp; in >> val;
	    in.ignore();
	  }

	  getline(in, tmp);
	  getline(in, tmp);
	  getline(in, tmp);
	  for (int j=0; j<w_data_b[0]->GetNbinsX()+2; j++) {
	    in >> tmp;
	    double val = 0.0;
	    in >> val; w_data_b[i]->SetBinContent(j, val); in >> tmp;
	    in >> val; w_data_b[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_bkg[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_eff[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_jec[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_jer[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_pu[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_bkg[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_top[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_bfit[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_btag[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_lumi[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; in >> tmp; in >> val;
	    in.ignore();
	  }

	  in.close();

	}

	TH1F* h_data = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b = (TH1F*)w_data_b[0]->Clone();
	TH1F* h_data_stat = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b_stat = (TH1F*)w_data_b[0]->Clone();
	TH1F* h_data_syst = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b_syst = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_tot = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b_tot = (TH1F*)w_data[0]->Clone();

	TH1F* stat_bkg = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_bkg = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_eff = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_eff = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_jec = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_jec = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_jer = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_jer = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_pu = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_pu = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_bkg = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_bkg = (TH1F*)w_data_b[0]->Clone();

	TH1F* stat_top = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_top = (TH1F*)w_data_b[0]->Clone();

	TH1F* stat_bfit = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_bfit = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_btag = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_btag = (TH1F*)w_data_b[0]->Clone();

	TH1F* stat_unfold = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_unfold = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_unfold = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_unfold = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_lumi = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_lumi = (TH1F*)w_data_b[0]->Clone();

	for (int i=0;i<=h_data_stat->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = calc(0, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  h_data->SetBinContent(i, val);
	  h_data_stat->SetBinContent(i, val);
	  h_data_syst->SetBinContent(i, val);
	  h_data_tot->SetBinContent(i, val);
	  double ref = 0.0;
	  ref = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			1.1*w_data[0]->GetBinError(i), 1.1*w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  h_data->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			1.1*w_stat_bkg[0]->GetBinError(i), 1.1*w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_bkg->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			1.1*w_syst_eff[0]->GetBinError(i), 1.1*w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_eff->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			1.1*w_syst_jer[0]->GetBinError(i), 1.1*w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_jer->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			1.1*w_syst_jec[0]->GetBinError(i), 1.1*w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_jec->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			1.1*w_syst_pu[0]->GetBinError(i), 1.1*w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_pu->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			1.1*w_syst_bkg[0]->GetBinError(i), 1.1*w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_bkg->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			1.1*w_stat_top[0]->GetBinError(i), 1.1*w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_top->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			1.1*w_stat_bfit[0]->GetBinError(i), 1.1*w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_bfit->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			1.1*w_syst_btag[0]->GetBinError(i), 1.1*w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_btag->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			1.1*w_stat_unfold[0]->GetBinError(i), 1.1*w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_unfold->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			1.1*w_syst_unfold[0]->GetBinError(i), 1.1*w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_unfold->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			1.1*w_syst_lumi[0]->GetBinError(i), 1.1*w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_lumi->SetBinError(i, val);

	  val = TMath::Sqrt(TMath::Power(h_data->GetBinError(i),2)+TMath::Power(stat_top->GetBinError(i),2)+TMath::Power(stat_bfit->GetBinError(i),2)+TMath::Power(stat_unfold->GetBinError(i),2));
	  h_data_stat->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(stat_bkg->GetBinError(i),2)+TMath::Power(syst_eff->GetBinError(i),2)+TMath::Power(syst_jec->GetBinError(i),2)+TMath::Power(syst_jer->GetBinError(i),2)+TMath::Power(syst_pu->GetBinError(i),2)+TMath::Power(syst_bkg->GetBinError(i),2)+TMath::Power(syst_btag->GetBinError(i),2)+TMath::Power(syst_unfold->GetBinError(i),2)+TMath::Power(syst_lumi->GetBinError(i),2));
	  h_data_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(h_data_syst->GetBinError(i),2));
	  h_data_tot->SetBinError(i, val);
	}

	for (int i=0;i<=h_data_b_stat->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = calc(0, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  h_data_b->SetBinContent(i, val);
	  h_data_b_stat->SetBinContent(i, val);
	  h_data_b_syst->SetBinContent(i, val);
	  h_data_b_tot->SetBinContent(i, val);
	  double ref = 0.0;
	  ref = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			1.1*w_data_b[0]->GetBinError(i), 1.1*w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  h_data_b->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			1.1*w_stat_b_bkg[0]->GetBinError(i), 1.1*w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_bkg->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			1.1*w_syst_b_eff[0]->GetBinError(i), 1.1*w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_eff->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			1.1*w_syst_b_jer[0]->GetBinError(i), 1.1*w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_jer->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			1.1*w_syst_b_jec[0]->GetBinError(i), 1.1*w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_jec->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			1.1*w_syst_b_pu[0]->GetBinError(i), 1.1*w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_pu->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			1.1*w_syst_b_bkg[0]->GetBinError(i), 1.1*w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_bkg->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			1.1*w_stat_b_top[0]->GetBinError(i), 1.1*w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_top->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			1.1*w_stat_b_bfit[0]->GetBinError(i), 1.1*w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_bfit->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			1.1*w_syst_b_btag[0]->GetBinError(i), 1.1*w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_btag->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			1.1*w_stat_b_unfold[0]->GetBinError(i), 1.1*w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_unfold->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			1.1*w_syst_b_unfold[0]->GetBinError(i), 1.1*w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_unfold->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			1.1*w_syst_b_lumi[0]->GetBinError(i), 1.1*w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_lumi->SetBinError(i, val);

	  val = TMath::Sqrt(TMath::Power(h_data_b->GetBinError(i),2)+TMath::Power(stat_b_top->GetBinError(i),2)+TMath::Power(stat_b_bfit->GetBinError(i),2)+TMath::Power(stat_b_unfold->GetBinError(i),2));
	  h_data_b_stat->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(stat_b_bkg->GetBinError(i),2)+TMath::Power(syst_b_eff->GetBinError(i),2)+TMath::Power(syst_b_jec->GetBinError(i),2)+TMath::Power(syst_b_jer->GetBinError(i),2)+TMath::Power(syst_b_pu->GetBinError(i),2)+TMath::Power(syst_b_bkg->GetBinError(i),2)+TMath::Power(syst_b_btag->GetBinError(i),2)+TMath::Power(syst_b_unfold->GetBinError(i),2)+TMath::Power(syst_b_lumi->GetBinError(i),2));
	  h_data_b_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(h_data_b_syst->GetBinError(i),2));
	  h_data_b_tot->SetBinError(i, val);
	}

/*
	h_data = fixrange(h_data);
	h_data_b = fixrange(h_data_b);
	h_data_stat = fixrange(h_data_stat);
	h_data_b_stat = fixrange(h_data_b_stat);
	h_data_syst = fixrange(h_data_syst);
	h_data_b_syst = fixrange(h_data_b_syst);
	h_data_tot = fixrange(h_data_tot);
	h_data_b_tot = fixrange(h_data_b_tot);

	h_mcg = fixrange(h_mcg);
	h_mcg_b = fixrange(h_mcg_b);
*/

	h_data = rebin(h_data);
	h_data_b = rebin(h_data_b);
	h_data_stat = rebin(h_data_stat);
	h_data_b_stat = rebin(h_data_b_stat);
	h_data_syst = rebin(h_data_syst);
	h_data_b_syst = rebin(h_data_b_syst);
	h_data_tot = rebin(h_data_tot);
	h_data_b_tot = rebin(h_data_b_tot);

	h_mcg = rebin(h_mcg);
	h_mcg_b = rebin(h_mcg_b);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	h_mcg_b->SetTitle("");
	h_mcg_b->GetYaxis()->SetTitle("#sigma [pb]");

	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	pad1->SetTopMargin(0.115);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();

	h_mcg_b->SetLineColor(kGreen+2);
	h_mcg_b->SetLineWidth(2);
	h_mcg_b->SetFillColor(kGreen+2);
	h_mcg_b->SetMarkerColor(kGreen+2);
	h_mcg_b->SetStats(0);

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

	TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	pad1->SetLogy();

	h_mcg_b->SetMaximum(4*h_data_tot->GetMaximum());
	h_mcg_b->SetMinimum(TMath::Max(0.002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));
	if (title.find("_mt")!=string::npos) h_mcg_b->SetMinimum(TMath::Max(0.00002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));
	if (title.find("_pt")!=string::npos) h_mcg_b->SetMinimum(TMath::Max(0.000002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));
	if (title.find("_mass")!=string::npos) h_mcg_b->SetMinimum(TMath::Max(0.00002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));

	h_mcg_b->Draw("E5");
	TH1F* tmp2 = (TH1F*)h_mcg_b->Clone();
	if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	  if (tmp2->GetMinimum()==0) tmp2->GetXaxis()->SetRangeUser(0, tmp2->GetBinCenter(tmp2->GetMinimumBin()-1));
	}
	tmp2->SetFillColor(0);
	tmp2->DrawClone("HISTLSAME");

	h_data_b_tot->Draw("E1PX0SAME");
	h_data_b_stat->Draw("E1PX0SAME");

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

	if (drawInclusive) leg->AddEntry(h_data_stat,"W(#rightarrow l#nu)+1b DATA","p");
	leg->AddEntry(h_data_b_stat,"W(#rightarrow l#nu)+2b DATA","p");
	leg->AddEntry(h_mcg,"W(#rightarrow l#nu) MadGraph","l");

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

	if (title_b=="w_first_jet_pt_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	} else if (title_b=="w_first_jet_eta_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading jet #eta");
	} else if (title_b=="w_first_bjet_pt") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	} else if (title_b=="w_first_bjet_eta") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_M_tot->GetXaxis()->SetTitle("leading b-jet #eta");
	}

	if (plot) {
	  ofstream out, out1, out2;
	  if (unfold) {
	    gSystem->mkdir((path + "/combined/" + version + "/xsecs_unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/combined/" + version + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	    out.open((path + "/combined/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.dat").c_str());
	    out1.open((path + "/combined/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.txt").c_str());
	    out2.open((path + "/combined/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.tex").c_str());
	  } else {
	    gSystem->mkdir((path + "/combined/" + version + "/xsecs/").c_str(), kTRUE);
	    c1->SaveAs((path + "/combined/" + version + "/xsecs/" + title_b + "_xsecs.pdf").c_str());
	    out.open((path + "/combined/" + version + "/" + "/xsecs/" + title_b + "_xsecs.dat").c_str());
	    out1.open((path + "/combined/" + version + "/" + "/xsecs/" + title_b + "_xsecs.txt").c_str());
	    out2.open((path + "/combined/" + version + "/" + "/xsecs/" + title_b + "_xsecs.tex").c_str());
	  }
	  out << h_data->GetName();
	  out << endl;
	  out << std::setw(25) << "data";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "eff";
	  out << std::setw(12) << "jec";
	  out << std::setw(12) << "jer";
	  out << std::setw(12) << "pu";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "ttbar";
	  out << std::setw(12) << "bfit";
	  out << std::setw(12) << "btag";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "lumi";
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
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
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
	    out << " +- ";
	    out << std::setw(8) << syst_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_top->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_bfit->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_btag->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_lumi->GetBinError(i);
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
	  out << h_data_b->GetName();
	  out << endl;
	  out << std::setw(25) << "data";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "eff";
	  out << std::setw(12) << "jec";
	  out << std::setw(12) << "jer";
	  out << std::setw(12) << "pu";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "ttbar";
	  out << std::setw(12) << "bfit";
	  out << std::setw(12) << "btag";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "lumi";
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
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
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
	    out << " +- ";
	    out << std::setw(8) << syst_b_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_top->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_bfit->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_btag->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_lumi->GetBinError(i);
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
	  out1 << h_data->GetName() << " - RELATIVE ERRORS";
	  out1 << endl;
	  out1 << std::setw(7) << "data";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "eff";
	  out1 << std::setw(8) << "jec";
	  out1 << std::setw(8) << "jer";
	  out1 << std::setw(8) << "pu";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "ttbar";
	  out1 << std::setw(8) << "bfit";
	  out1 << std::setw(8) << "btag";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "lumi";
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
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
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
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_top->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_bfit->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_btag->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_lumi->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_stat->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << h_data_syst->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_tot->GetBinError(i)*val;
	    out1 << endl;
	  }
	  out1 << h_data_b->GetName() << " - RELATIVE ERRORS";
	  out1 << endl;
	  out1 << std::setw(7) << "data";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "eff";
	  out1 << std::setw(8) << "jec";
	  out1 << std::setw(8) << "jer";
	  out1 << std::setw(8) << "pu";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "ttbar";
	  out1 << std::setw(8) << "bfit";
	  out1 << std::setw(8) << "btag";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "lumi";
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
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
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
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_top->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_bfit->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_btag->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_lumi->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_b_stat->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << h_data_b_syst->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_b_tot->GetBinError(i)*val;
	    out1 << endl;
	  }
	  out1.close();
	  //out2 << h_data->GetName() << " - RELATIVE ERRORS";
	  //out2 << endl;
	  out2 << std::setw(7) << "\\textbf{data} &"  ;
	  out2 << std::setw(8) << "\\textbf{bkg} &"   ;
	  out2 << std::setw(8) << "\\textbf{eff} &"   ;
	  out2 << std::setw(8) << "\\textbf{jec} &"   ;
	  out2 << std::setw(8) << "\\textbf{jer} &"   ;
	  out2 << std::setw(8) << "\\textbf{pu} &"    ;
	  out2 << std::setw(8) << "\\textbf{bkg} &"   ;
	  out2 << std::setw(8) << "\\textbf{ttbar} &" ;
	  out2 << std::setw(8) << "\\textbf{bfit} &"  ;
	  out2 << std::setw(8) << "\\textbf{btag} &"  ;
	  out2 << std::setw(8) << "\\textbf{unfold} &";
	  out2 << std::setw(8) << "\\textbf{unfold} &";
	  out2 << std::setw(8) << "\\textbf{lumi} &"  ;
	  out2 << std::setw(8) << "\\textbf{total} &" ;
	  out2 << std::setw(8) << "\\textbf{total} &" ;
	  out2 << std::setw(8) << "\\textbf{total} &" ;
	  out2 << endl;
	  out2 << std::setw(7) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{error & ";
	  out2 << endl;
	  for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	    double val = 100.*(h_data->GetBinContent(i)==0 ? 0 : 1./h_data->GetBinContent(i));
	    out2 << std::fixed;
	    out2 << std::setw(2) << i;
	    out2 << " ";
	    out2 << std::setprecision(1);
	    out2 << std::setw(4) << h_data->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_bkg->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_eff->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_jec->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_jer->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_pu->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_bkg->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_top->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_bfit->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_btag->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_unfold->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_unfold->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_lumi->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_stat->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_syst->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_tot->GetBinError(i)*val;
	    out2 << std::setw(4) << "\\tabularnewline" << "   " << "\\hline";
	    out2 << endl;
	  }
	  //out2 << h_data_b->GetName() << " - RELATIVE ERRORS";
	  //out2 << endl;
	  out2 << std::setw(7) << "\\textbf{data} &";
	  out2 << std::setw(8) << "\\textbf{bkg} &";
	  out2 << std::setw(8) << "\\textbf{eff} &";
	  out2 << std::setw(8) << "\\textbf{jec} &";
	  out2 << std::setw(8) << "\\textbf{jer} &";
	  out2 << std::setw(8) << "\\textbf{pu} &";
	  out2 << std::setw(8) << "\\textbf{bkg} &";
	  out2 << std::setw(8) << "\\textbf{ttbar} &";
	  out2 << std::setw(8) << "\\textbf{bfit} &";
	  out2 << std::setw(8) << "\\textbf{btag} &";
	  out2 << std::setw(8) << "\\textbf{unfold} &";
	  out2 << std::setw(8) << "\\textbf{unfold} &";
	  out2 << std::setw(8) << "\\textbf{lumi} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << endl;
	  out2 << std::setw(7) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{error} &";
	  out2 << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    double val = 100.*(h_data_b->GetBinContent(i)==0 ? 0 : 1./h_data_b->GetBinContent(i));
	    out2 << std::fixed;
	    out2 << std::setw(2) << i;
	    out2 << " ";
	    out2 << std::setprecision(1);
	    out2 << std::setw(4) << h_data_b->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_bkg->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_eff->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_jec->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_jer->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_pu->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_bkg->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_top->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_bfit->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_btag->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_unfold->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_unfold->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_lumi->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_b_stat->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_b_syst->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_b_tot->GetBinError(i)*val;
	    out2 << std::setw(4) << "\\tabularnewline" << "   " << "\\hline";
	    out2 << endl;
	  }
           
	}
}

