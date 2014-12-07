#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v12.h"

#include "fixrange.C"
#include "rebin.C"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

TH1F* h_data_fit = 0;
TH1F* h_mc_fit0 = 0;
TH1F* h_mc_fit1 = 0;
TH1F* h_mc_fit2 = 0;

float mc_fit0 = 1.0;
float mc_fit1 = 1.0;
float mc_fit2 = 1.0;
float mc_fit3 = 1.0;

float e_mc_fit0 = 0.0;
float e_mc_fit1 = 0.0;
float e_mc_fit2 = 0.0;
float e_mc_fit3 = 0.0;

void fcn(int& npar, double* gin, double& fun, double* par, int iflag) {
  double chisq = 0.0;
  if (iflag) {};
  if (gin) {};
  for (int i=1; i<=h_data_fit->GetNbinsX(); i++) {
    double xn = h_data_fit->GetBinContent(i);
    double xd = TMath::Power(h_data_fit->GetBinError(i),2);
    double w = 1.0;
    if (npar>3) {
      w = par[3];
    }
    if (npar>0) {
      xn = xn - w*par[0]*h_mc_fit0->GetBinContent(i);
      xd = xd + TMath::Power(w*par[0]*h_mc_fit0->GetBinError(i),2);
    }
    if (npar>1) {
      xn = xn - w*par[1]*h_mc_fit1->GetBinContent(i);
      xd = xd + TMath::Power(w*par[1]*h_mc_fit1->GetBinError(i),2);
    }
    if (npar>2) {
      xn = xn - par[2]*h_mc_fit2->GetBinContent(i);
      xd = xd + TMath::Power(par[2]*h_mc_fit2->GetBinError(i),2);
    }
    if (xd!=0) chisq = chisq + (xn*xn)/xd;
  }
  if (e_mc_fit0!=0) chisq = chisq + TMath::Power((par[0]-mc_fit0)/e_mc_fit0,2);
  if (e_mc_fit1!=0) chisq = chisq + TMath::Power((par[1]-mc_fit1)/e_mc_fit1,2);
  if (e_mc_fit2!=0) chisq = chisq + TMath::Power((par[2]-mc_fit2)/e_mc_fit2,2);
  if (e_mc_fit3!=0) chisq = chisq + TMath::Power((par[3]-mc_fit3)/e_mc_fit3,2);
  fun = chisq;
}

void DataMCComp(int irun=0, string title="", int plot=0, int ilepton=1, int doBkg=0, int doFit=0) {

//int useFitResults = 0; // do not use fit results for c_t
int useFitResults = 1; // use fit results for c_t

//int useFitResults2=0;  // do not use constrained fit results for c_b, c_t, c_qcd, c_s
int useFitResults2=1;  // use constrained fit results for c_b, c_t, c_qcd, c_s

//int useWbb=0; // do not use the special Wbb MC sample
int useWbb=1; // use the special Wbb MC sample

int printYield=0; // do not print pre-fit and post-fit yields
//int printYield=1; // print pre-fit and post-fit yields

string subdir="0";
string postfix="";
if (irun==1) {             // irun==1 => JEC Up
  subdir="1";
  postfix="Up";
}
if (irun==2) {             // irun==2 => JEC Down
  subdir="2";
  postfix="Down";
}
if (irun==3) {             // irun==3 => PU Up
  subdir="3";
  postfix="Pup";
}
if (irun==4) {             // irun==4 => PU Down
  subdir="4";
  postfix="Pum";
}
if (irun==5) {             // irun==5 => top bkg
  subdir="5";
  postfix="";
}
if (irun==6) {             // irun==6 => b purity
  subdir="6";
  postfix="";
}
if (irun==7) {             // irun==7 => unfolding
  subdir="7";
  postfix="";
}
if (irun==8) {             // irun==8 => unfolding with Sherpa
  subdir="8";
  postfix="";
}
if (irun==9) {             // irun==9 => unfolding with Powheg
  subdir="9";
  postfix="";
}
if (irun==10) {            // irun==10 => bkg systematics
  subdir="10";
  postfix="";
}
if (irun==11) {            // irun==11 => JER Up
  subdir="11";
  postfix="JerUp";
}
if (irun==12) {            // irun==12 => JER Down
  subdir="12";
  postfix="JerDown";
}
if (irun==13) {            // irun==13 => bkg statistics
  subdir="13";
  postfix="";
}
if (irun==14) {            // irun==14 => unfolding with data weight
  subdir="14";
  postfix="";
}
if (irun==15) {            // irun==15 => qcd bkg
  subdir="15";
  postfix="";
}
if (irun==77) {            // irun==77 => unfolding with MadGraph 4FS
  subdir="77";
  postfix="";
}
if (irun==88) {            // irun==88 => deltaR
  subdir="88";
  postfix="DR";
}
if (irun==99) {            // irun==99 => pur
  subdir="99";
  postfix="Pur";
}

// skip variations for FWD (5,6) and TOP (7,8) samples
if (ilepton>=5 && ilepton<=8) postfix="";

	/* QCD background */

	double c1_qcd=1.0;
	double ec1_qcd=0.0;
	double c2_qcd=1.0;
	double ec2_qcd=0.0;
	double c3_qcd=1.0;
	double ec3_qcd=0.0;

	/* top background */

	double c1_t=1.0;
	double ec1_t=0.0;
	double c2_t=1.0;
	double ec2_t=0.0;
	double c3_t=1.0;
	double ec3_t=0.0;

	/* scale factor 1 */
	double c1_s=1.0;
	double ec1_s=0.0;
	double c2_s=1.0;
	double ec2_s=0.0;
	double c3_s=1.0;
	double ec3_s=0.0;

	/* scale factor 2 */
	double c1_r=1.0;
	double ec1_r=0.0;
	double c2_r=1.0;
	double ec2_r=0.0;
	double c3_r=1.0;
	double ec3_r=0.0;

	/* purity */

	double c1_b=1.0;
	double ec1_b=0.0;
	double c2_b=1.0;
	double ec2_b=0.0;
	double c3_b=1.0;
	double ec3_b=0.0;

	if (doFit==1 || doFit==2) {
	  if (title=="w_mt_wenu_wide") useFitResults=0;
	  if (title=="w_mt_wenu_b_wide") useFitResults=0;
	  if (title=="w_mt_wenu_bb_wide") useFitResults=0;
	  if (title=="w_mt_wmnu_wide") useFitResults=0;
	  if (title=="w_mt_wmnu_b_wide") useFitResults=0;
	  if (title=="w_mt_wmnu_bb_wide") useFitResults=0;
	}

	if (title=="w_mt_wenu_wide") useFitResults2=0;
	if (title=="w_mt_wenu_b_wide") useFitResults2=0;
	if (title=="w_mt_wenu_bb_wide") useFitResults2=0;
	if (title=="w_mt_wmnu_wide") useFitResults2=0;
	if (title=="w_mt_wmnu_b_wide") useFitResults2=0;
	if (title=="w_mt_wmnu_bb_wide") useFitResults2=0;

	if (ilepton==3 || ilepton==4) useFitResults=0;

	if (ilepton==3 || ilepton==4) useFitResults2=0;
	if (ilepton==5 || ilepton==6) useFitResults2=0;
	if (ilepton==7 || ilepton==8) useFitResults2=0;

	ifstream in1, in2, in3, in4, in5, in6, in10, in11, in12;
	if (ilepton==1) {
	  in1.open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	  in2.open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	  in3.open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  if (useFitResults) {
	    in4.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	    in5.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	    in6.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  }
	  if (useFitResults2) {
	    in10.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	    in11.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	    in12.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  }
	}
	if (ilepton==5) {
	  in1.open((path + "/electronsFWD/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	  in2.open((path + "/electronsFWD/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	  in3.open((path + "/electronsFWD/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  if (useFitResults) {
	    in4.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	    in5.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	    in6.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  }
	}
	if (ilepton==7) {
	  in1.open((path + "/electronsTOP/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	  in2.open((path + "/electronsTOP/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	  in3.open((path + "/electronsTOP/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  if (useFitResults) {
	    in4.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	    in5.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	    in6.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  }
	}
	if (ilepton==2) {
	  in1.open((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	  in2.open((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	  in3.open((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
	  if (useFitResults) {
	    in4.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	    in5.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	    in6.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
	  }
	  if (useFitResults2) {
	    in10.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	    in11.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	    in12.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
	  }
	}
	if (ilepton==6) {
	  in1.open((path + "/muonsFWD/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	  in2.open((path + "/muonsFWD/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	  in3.open((path + "/muonsFWD/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
	  if (useFitResults) {
	    in4.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	    in5.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	    in6.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
	  }
	}
	if (ilepton==8) {
	  in1.open((path + "/muonsTOP/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	  in2.open((path + "/muonsTOP/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	  in3.open((path + "/muonsTOP/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
	  if (useFitResults) {
	    in4.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	    in5.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	    in6.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
	  }
	}

	in1 >> c1_qcd >> ec1_qcd;
	in2 >> c2_qcd >> ec2_qcd;
	in3 >> c3_qcd >> ec3_qcd;
	in1.close();
	in2.close();
	in3.close();
	if (useFitResults) {
	  in4 >> c1_t >> ec1_t;
	  in5 >> c2_t >> ec2_t;
	  in6 >> c3_t >> ec3_t;
	  in4.close();
	  in5.close();
	  in6.close();
	}
	if (useFitResults2) {
	  double c;
	  double ec;
	  in10 >> c1_b >> ec1_b;
	  in10 >> c1_t >> ec1_t;
	  c = 1.0; ec = 0.0;
	  in10 >> c >> ec;
	  c1_qcd = c1_qcd * c; ec1_qcd = c1_qcd * ec;
	  in10 >> c1_s >> ec1_s;
	  in10 >> c1_r >> ec1_r;
	  in11 >> c2_b >> ec2_b;
	  in11 >> c2_t >> ec2_t;
	  c = 1.0; ec = 0.0;
	  in11 >> c >> ec;
	  c2_qcd = c2_qcd * c; ec2_qcd = c2_qcd * ec;
	  in11 >> c2_s >> ec2_s;
	  in11 >> c2_r >> ec2_r;
	  in12 >> c3_b >> ec3_b;
	  in12 >> c3_t >> ec3_t;
	  c = 1.0; ec = 0.0;
	  in12 >> c >> ec;
	  c3_qcd = c3_qcd * c; ec3_qcd = c3_qcd * ec;
	  in12 >> c3_s >> ec3_s;
	  in12 >> c3_r >> ec3_r;
	  in10.close();
	  in11.close();
	  in12.close();
	}

	c1_b = c1_b * c1_s; ec1_b = TMath::Sqrt(TMath::Power(c1_b*ec1_s,2)+TMath::Power(ec1_b*c1_s,2));
	c1_t = c1_t * c1_s; ec1_t = TMath::Sqrt(TMath::Power(c1_t*ec1_s,2)+TMath::Power(ec1_t*c1_s,2));
	c2_b = c2_b * c2_s; ec2_b = TMath::Sqrt(TMath::Power(c2_b*ec2_s,2)+TMath::Power(ec2_b*c2_s,2));
	c2_t = c2_t * c2_s; ec2_t = TMath::Sqrt(TMath::Power(c2_t*ec2_s,2)+TMath::Power(ec2_t*c2_s,2));
	c3_b = c3_b * c3_s; ec3_b = TMath::Sqrt(TMath::Power(c3_b*ec3_s,2)+TMath::Power(ec3_b*c3_s,2));
	c3_t = c3_t * c3_s; ec3_t = TMath::Sqrt(TMath::Power(c3_t*ec3_s,2)+TMath::Power(ec3_t*c3_s,2));

	c1_b = c1_b * c1_r; ec1_b = TMath::Sqrt(TMath::Power(c1_b*ec1_r,2)+TMath::Power(ec1_b*c1_r,2));
	c1_t = c1_t * c1_r; ec1_t = TMath::Sqrt(TMath::Power(c1_t*ec1_r,2)+TMath::Power(ec1_t*c1_r,2));
	c2_b = c2_b * c2_r; ec2_b = TMath::Sqrt(TMath::Power(c2_b*ec2_r,2)+TMath::Power(ec2_b*c2_r,2));
	c2_t = c2_t * c2_r; ec2_t = TMath::Sqrt(TMath::Power(c2_t*ec2_r,2)+TMath::Power(ec2_t*c2_r,2));
	c3_b = c3_b * c3_r; ec3_b = TMath::Sqrt(TMath::Power(c3_b*ec3_r,2)+TMath::Power(ec3_b*c3_r,2));
	c3_t = c3_t * c3_r; ec3_t = TMath::Sqrt(TMath::Power(c3_t*ec3_r,2)+TMath::Power(ec3_t*c3_r,2));

	double Lumi2012=0;

	if (ilepton==1 || ilepton==3 || ilepton==5 || ilepton==7) Lumi2012 = Lumi2012_ele;
	if (ilepton==2 || ilepton==4 || ilepton==6 || ilepton==8) Lumi2012 = Lumi2012_muon;

	double norm1 = ((Lumi2012 * Xsec_wj) / Ngen_wj);
	double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
	double norm3 = ((Lumi2012 * Xsec_zz) / Ngen_zz);
	double norm4 = ((Lumi2012 * Xsec_wz) / Ngen_wz);
//	double norm5 = ((Lumi2012 * Xsec_qcd) / Ngen_qcd);
	double norm5 = 1.0;
	double norm6 = ((Lumi2012 * Xsec_ww) / Ngen_ww);
	double norm7 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm8 = ((Lumi2012 * Xsec_tbar_t) / Ngen_tbar_t);

	double enorm1 = ((Lumi2012 * eXsec_wj) / Ngen_wj);
	double enorm2 = ((Lumi2012 * eXsec_tt) / Ngen_tt);
	if (useFitResults) enorm2 = 0.0;
	double enorm3 = ((Lumi2012 * eXsec_zz) / Ngen_zz);
	double enorm4 = ((Lumi2012 * eXsec_wz) / Ngen_wz);
//	double enorm5 = ((Lumi2012 * eXsec_qcd) / Ngen_qcd);
	double enorm5 = 0.0;
	double enorm6 = ((Lumi2012 * eXsec_ww) / Ngen_ww);
	double enorm7 = ((Lumi2012 * eXsec_dy) / Ngen_dy);
	double enorm8 = ((Lumi2012 * eXsec_tbar_t) / Ngen_tbar_t);

	if (title.empty()) title = "w_jetmultiplicity";

	if (ilepton==1 || ilepton==3 || ilepton==5 || ilepton==7) {
	  if (title.find("muon")!=string::npos) return;
	  if (title.find("mm")!=string::npos) return;
	  if (title.find("wmnu")!=string::npos) return;
	}
	if (ilepton==2 || ilepton==4 || ilepton==6 || ilepton==8) {
	  if (title.find("ele")!=string::npos) return;
	  if (title.find("ee")!=string::npos) return;
	  if (title.find("wenu")!=string::npos) return;
	}

	TFile *data=0;
	if (ilepton==1 || ilepton==3 || ilepton==5 || ilepton==7) data = TFile::Open((path + "/" + version + "/" + "SingleElectron_2012_merge.root").c_str());
	if (ilepton==2 || ilepton==4 || ilepton==6 || ilepton==8) data = TFile::Open((path + "/" + version + "/" + "SingleMu_2012_merge.root").c_str());

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "Wj_merge.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar_merge.root").c_str());
	TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
	TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
	TFile *mc5 = 0;
//	mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
	TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
	TFile *mc7 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	TFile *mc8 = TFile::Open((path + "/" + version + "/" + "T_merge.root").c_str());

	if (ilepton==1) data->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) data->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) data->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) data->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) data->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) data->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) data->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) data->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_data = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc1->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) mc1->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) mc1->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) mc1->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) mc1->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) mc1->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1b = (TH1F*)gDirectory->Get(("b"+title.substr(1)).c_str());
	TH1F* h_mc1c = (TH1F*)gDirectory->Get(("c"+title.substr(1)).c_str());
	TH1F* h_mc1t = (TH1F*)gDirectory->Get(("t"+title.substr(1)).c_str());

	if (useWbb) {
	  if (title.find("_bb")!=string::npos || title.find("_2b")!=string::npos) {
	    if (h_mc1b) {
	      float xval = h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1);
	      for (int i=0; i<=h_mc1b->GetNbinsX()+1; i++) {
	        h_mc1->SetBinError(i, TMath::Sqrt(TMath::Power(h_mc1->GetBinError(i),2)-TMath::Power(h_mc1b->GetBinError(i),2)));
	        h_mc1->SetBinContent(i, h_mc1->GetBinContent(i)-h_mc1b->GetBinContent(i));
	      }
	      mc1 = TFile::Open((path + "/" + version + "/" + "Wbb.root").c_str());
	      if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
	      if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
	      if (ilepton==3) mc1->cd(("demoEleQCD"+postfix).c_str());
	      if (ilepton==4) mc1->cd(("demoMuoQCD"+postfix).c_str());
	      if (ilepton==5) mc1->cd(("demoEleFWD"+postfix).c_str());
	      if (ilepton==6) mc1->cd(("demoMuoFWD"+postfix).c_str());
	      if (ilepton==7) mc1->cd(("demoEleTOP"+postfix).c_str());
	      if (ilepton==8) mc1->cd(("demoMuoTOP"+postfix).c_str());
	      h_mc1b = (TH1F*)gDirectory->Get(("b"+title.substr(1)).c_str());
	      h_mc1b->Sumw2();
	      xval = xval / h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1);
	      h_mc1b->Scale(xval);
	      h_mc1->Add(h_mc1b, +1.);
	    }
	  }
	}

	if (ilepton==1) mc2->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc2->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc2->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) mc2->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) mc2->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) mc2->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) mc2->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) mc2->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc3->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc3->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc3->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) mc3->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) mc3->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) mc3->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) mc3->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) mc3->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc4->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc4->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc4->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) mc4->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) mc4->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) mc4->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) mc4->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) mc4->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());

//	if (ilepton==1) mc5->cd(("demoEle"+postfix).c_str());
//	if (ilepton==2) mc5->cd(("demoMuo"+postfix).c_str());
//	if (ilepton==3) mc5->cd(("demoEleQCD"+postfix).c_str());
//	if (ilepton==4) mc5->cd(("demoMuoQCD"+postfix).c_str());
//	if (ilepton==5) mc5->cd(("demoEleFWD"+postfix).c_str());
//	if (ilepton==6) mc5->cd(("demoMuoFWD"+postfix).c_str());
//	if (ilepton==7) mc5->cd(("demoEleTOP"+postfix).c_str());
//	if (ilepton==8) mc5->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc5 = 0;
//	h_mc5 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1||ilepton==5||ilepton==7) {
	  mc5 = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str());
	  mc5->cd();
	  h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
	}
	if (ilepton==2||ilepton==6||ilepton==8) {
	  mc5 = TFile::Open((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str());
	  mc5->cd();
	  h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
	}

	if (ilepton==1) mc6->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc6->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc6->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) mc6->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) mc6->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) mc6->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) mc6->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) mc6->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc7->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc7->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc7->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) mc7->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) mc7->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) mc7->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) mc7->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) mc7->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc8->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc8->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc8->cd(("demoEleQCD"+postfix).c_str());
	if (ilepton==4) mc8->cd(("demoMuoQCD"+postfix).c_str());
	if (ilepton==5) mc8->cd(("demoEleFWD"+postfix).c_str());
	if (ilepton==6) mc8->cd(("demoMuoFWD"+postfix).c_str());
	if (ilepton==7) mc8->cd(("demoEleTOP"+postfix).c_str());
	if (ilepton==8) mc8->cd(("demoMuoTOP"+postfix).c_str());
	TH1F* h_mc8 = (TH1F*)gDirectory->Get(title.c_str());

	h_data -> Sumw2();

	h_mc1 -> Sumw2();
	h_mc1 -> SetLineColor(kBlack);
	h_mc1 -> SetFillColor(kYellow-4);
	//h_mc1 -> SetFillStyle(3004);

	if (h_mc1b) {
	  h_mc1b -> Sumw2();
	  h_mc1b->SetLineColor(kBlack);
	  h_mc1b->SetFillColor(kYellow-4);
	  h_mc1b->SetFillStyle(3254);
	}
	if (h_mc1c) {
	  h_mc1c -> Sumw2();
	  h_mc1c->SetLineColor(kBlack);
	  h_mc1c->SetFillColor(kOrange);
	  h_mc1c->SetFillStyle(3245);
	}
	if (h_mc1t) {
	  h_mc1t -> Sumw2();
	  h_mc1t->SetLineColor(kBlack);
	  h_mc1t->SetFillColor(kOrange-4);
	  //h_mc1t->SetFillStyle(3004);
	}

	h_mc2 -> Sumw2();
	h_mc2 -> SetLineColor(kBlack);
	h_mc2 -> SetFillColor(kBlue);
	//h_mc2 -> SetFillStyle(3004);

	h_mc3 -> Sumw2();
	h_mc3 -> SetLineColor(kBlack);
	h_mc3 -> SetFillColor(kGray+2);
	//h_mc3 -> SetFillStyle(3004);

	h_mc4 -> Sumw2();
	h_mc4 -> SetLineColor(kBlack);
	h_mc4 -> SetFillColor(kGray+3);
	//h_mc4 -> SetFillStyle(3004);

	if (h_mc5) {
	  h_mc5 -> Sumw2();
	  h_mc5 -> SetLineColor(kBlack);
	  h_mc5 -> SetFillColor(kGreen-3);
	  //h_mc5 -> SetFillStyle(3004);
	}

	h_mc6 -> Sumw2();
	h_mc6 -> SetLineColor(kBlack);
	h_mc6 -> SetFillColor(kRed+2);
	//h_mc6 -> SetFillStyle(3004);

	h_mc7 -> Sumw2();
	h_mc7 -> SetLineColor(kBlack);
	h_mc7 -> SetFillColor(kGray);
	//h_mc7 -> SetFillStyle(3004);

	h_mc8 -> Sumw2();
	h_mc8 -> SetLineColor(kBlack);
	h_mc8 -> SetFillColor(kAzure-3);
	//h_mc8 -> SetFillStyle(3004);

	if (irun==10) {
	  norm1 = norm1 + 0.1*enorm1;
	  norm2 = norm2 + 0.1*enorm2;
	  norm3 = norm3 + 0.1*enorm3;
	  norm4 = norm4 + 0.1*enorm4;
	  norm5 = norm5 + 0.1*enorm5;
	  norm6 = norm6 + 0.1*enorm6;
	  norm7 = norm7 + 0.1*enorm7;
	  norm8 = norm8 + 0.1*enorm8;
	}

	h_mc1->Scale(norm1);
	if (h_mc1b) h_mc1b->Scale(norm1);
	if (h_mc1c) h_mc1c->Scale(norm1);
	if (h_mc1t) h_mc1t->Scale(norm1);
	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
	if (h_mc5) h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);
	h_mc8->Scale(norm8);

	if (printYield) {
	  cout << "******************" << endl;
	  cout << "PRE FIT YIELDS:" << endl;
	  cout << "data: " << h_data->Integral(0,h_data->GetNbinsX()+1) << endl;
	  cout << "Wjets: " << h_mc1->Integral(0,h_mc1->GetNbinsX()+1) << endl;
	  if (h_mc1b) cout << "Wb: " << h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1) << endl;
	  if (h_mc1c) cout << "Wc: " << h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1) << endl;
	  if (h_mc1t) cout << "Wtau: " << h_mc1t->Integral(0,h_mc1t->GetNbinsX()+1) << endl;
	  if (h_mc1b && h_mc1c && h_mc1t) cout << "Wlight: " << h_mc1->Integral(0,h_mc1->GetNbinsX()+1)
					    - h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1)
					    - h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1)
					    - h_mc1t->Integral(0,h_mc1t->GetNbinsX()+1) << endl;
	  cout << "TTbar: " << h_mc2->Integral(0,h_mc2->GetNbinsX()+1) << endl;
	  cout << "ZZ: " << h_mc3->Integral(0,h_mc3->GetNbinsX()+1) << endl;
	  cout << "WZ: " << h_mc4->Integral(0,h_mc4->GetNbinsX()+1) << endl;
	  if (h_mc5) cout << "QCD: " << h_mc5->Integral(0,h_mc5->GetNbinsX()+1) << endl;
	  cout << "WW: " << h_mc6->Integral(0,h_mc6->GetNbinsX()+1) << endl;
	  cout << "VV: " << h_mc3->Integral(0,h_mc3->GetNbinsX()+1)
	    + h_mc4->Integral(0,h_mc4->GetNbinsX()+1)
	    + h_mc6->Integral(0,h_mc6->GetNbinsX()+1) << endl;
	  cout << "DY: " << h_mc7->Integral(0,h_mc7->GetNbinsX()+1) << endl;
	  cout << "T: " << h_mc8->Integral(0,h_mc8->GetNbinsX()+1) << endl;
	  cout << "******************" << endl;
	  TH1F *dump_data;
	  TH1F *dump_wjets;
	  TH1F *dump_wb;
	  TH1F *dump_wc;
	  TH1F *dump_wl;
	  TH1F *dump_wtau;
	  TH1F *dump_TTbar;
	  TH1F *dump_zz;
	  TH1F *dump_wz;
	  TH1F *dump_qcd;
	  TH1F *dump_ww;
	  TH1F *dump_vv;
	  TH1F *dump_dy;
	  TH1F *dump_T;
	  dump_data = (TH1F*)h_data->Clone("Data");
	  dump_wjets = (TH1F*)h_mc1->Clone("W+jetsInclusive");
	  dump_wb = 0;
	  if (h_mc1b) dump_wb = (TH1F*)h_mc1b->Clone("W+b");
	  dump_wc = 0;
	  if (h_mc1c) dump_wc = (TH1F*)h_mc1c->Clone("W+c");
	  dump_wl = (TH1F*)h_mc1->Clone("W+l");
	  if (h_mc1b) dump_wl->Add(h_mc1b, -1.);
	  if (h_mc1c) dump_wl->Add(h_mc1c, -1.);
	  dump_wtau = 0;
	  if (h_mc1t) dump_wtau = (TH1F*)h_mc1t->Clone("W+tau");
	  dump_TTbar = (TH1F*)h_mc2->Clone("TTbar");
	  dump_zz = (TH1F*)h_mc3->Clone("ZZ");
	  dump_wz = (TH1F*)h_mc4->Clone("WZ");
	  dump_qcd = 0;
	  if (h_mc5) {
	    dump_qcd = (TH1F*)h_mc5->Clone("QCD");
	    dump_qcd->Scale(c3_qcd);
	  }
	  dump_ww = (TH1F*)h_mc6->Clone("WW");
	  dump_vv = (TH1F*)h_mc3->Clone("VV");
	  dump_vv->Add(h_mc4, +1.);
	  dump_vv->Add(h_mc6, +1.);
	  dump_dy = (TH1F*)h_mc7->Clone("DY");
	  dump_T = (TH1F*)h_mc8->Clone("T");
	  TFile dumphistos_file("dumped_histos.root","RECREATE");
	  dumphistos_file.cd();
	  dump_data->Write();
	  dump_wjets->Write();
	  if (dump_wb) dump_wb->Write();
	  if (dump_wc) dump_wc->Write();
	  dump_wl->Write();
	  if (dump_wtau) dump_wtau->Write();
	  dump_TTbar->Write();
	  dump_zz->Write();
	  dump_wz->Write();
	  if (dump_qcd) dump_qcd->Write();
	  dump_ww->Write();
	  dump_vv->Write();
	  dump_dy->Write();
	  dump_T->Write();
	  dumphistos_file.Close();
	}

	if (h_mc1b) h_mc1->Add(h_mc1b, -1.);
	if (h_mc1c) h_mc1->Add(h_mc1c, -1.);
	if (h_mc1t) h_mc1->Add(h_mc1t, -1.);
	for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	  float e = TMath::Power(h_mc1->GetBinError(i),2);
	  if (h_mc1b) e = e - TMath::Power(h_mc1b->GetBinError(i),2);
	  if (h_mc1c) e = e - TMath::Power(h_mc1c->GetBinError(i),2);
	  if (h_mc1t) e = e - TMath::Power(h_mc1t->GetBinError(i),2);
	  h_mc1->SetBinError(i, TMath::Sqrt(e));
	}

	if (useFitResults && !doFit) {
	  if (title.find("_bb")!=string::npos) {
	    if (irun==5) {
	      h_mc2->Scale(c3_t+0.1*ec3_t);
	      h_mc1->Scale(c3_r+0.1*ec3_r);
	      if (h_mc1t) h_mc1t->Scale(c3_r+0.1*ec3_r);
	      h_mc3->Scale(c3_r+0.1*ec3_r);
	      h_mc4->Scale(c3_r+0.1*ec3_r);
	      h_mc6->Scale(c3_r+0.1*ec3_r);
	      h_mc7->Scale(c3_r+0.1*ec3_r);
	      h_mc8->Scale(c3_r+0.1*ec3_r);
	    } else {
	      h_mc2->Scale(c3_t);
	      h_mc1->Scale(c3_r);
	      if (h_mc1t) h_mc1t->Scale(c3_r);
	      h_mc3->Scale(c3_r);
	      h_mc4->Scale(c3_r);
	      h_mc6->Scale(c3_r);
	      h_mc7->Scale(c3_r);
	      h_mc8->Scale(c3_r);
	    }
	  } else if (title.find("_b")!=string::npos) {
	    if (irun==5) {
	      h_mc2->Scale(c2_t+0.1*ec2_t);
	      h_mc1->Scale(c2_r+0.1*ec2_r);
	      if (h_mc1t) h_mc1t->Scale(c2_r+0.1*ec2_r);
	      h_mc3->Scale(c2_r+0.1*ec2_r);
	      h_mc4->Scale(c2_r+0.1*ec2_r);
	      h_mc6->Scale(c2_r+0.1*ec2_r);
	      h_mc7->Scale(c2_r+0.1*ec2_r);
	      h_mc8->Scale(c2_r+0.1*ec2_r);
	    } else {
	      h_mc2->Scale(c2_t);
	      h_mc1->Scale(c2_r);
	      if (h_mc1t) h_mc1t->Scale(c2_r);
	      h_mc3->Scale(c2_r);
	      h_mc4->Scale(c2_r);
	      h_mc6->Scale(c2_r);
	      h_mc7->Scale(c2_r);
	      h_mc8->Scale(c2_r);
	    }
	  } else {
	    if (irun==5) {
	      h_mc2->Scale(c1_t+0.1*ec1_t);
	      h_mc1->Scale(c1_r+0.1*ec1_r);
	      if (h_mc1t) h_mc1t->Scale(c1_r+0.1*ec1_r);
	      h_mc3->Scale(c1_r+0.1*ec1_r);
	      h_mc4->Scale(c1_r+0.1*ec1_r);
	      h_mc6->Scale(c1_r+0.1*ec1_r);
	      h_mc7->Scale(c1_r+0.1*ec1_r);
	      h_mc8->Scale(c1_r+0.1*ec1_r);
	    } else {
	      h_mc2->Scale(c1_t);
	      h_mc1->Scale(c1_r);
	      if (h_mc1t) h_mc1t->Scale(c1_r);
	      h_mc3->Scale(c1_r);
	      h_mc4->Scale(c1_r);
	      h_mc6->Scale(c1_r);
	      h_mc7->Scale(c1_r);
	      h_mc8->Scale(c1_r);
	    }
	  }
	}

	if (h_mc5) {
	  if (title.find("_bb")!=string::npos) {
	    if (irun==15) {
	      h_mc5->Scale(c3_qcd+0.1*ec3_qcd);
	    } else {
	      h_mc5->Scale(c3_qcd);
	    }
	  } else if (title.find("_b")!=string::npos) {
	    if (irun==15) {
	      h_mc5->Scale(c2_qcd+0.1*ec2_qcd);
	    } else {
	      h_mc5->Scale(c2_qcd);
	    }
	  } else {
	    if (irun==15) {
	      h_mc5->Scale(c1_qcd+0.1*ec1_qcd);
	    } else {
	      h_mc5->Scale(c1_qcd);
	    }
	  }
	}

	if (irun==13) {
	  for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	    h_mc1->SetBinError(i, 1.1*h_mc1->GetBinError(i));
	    if (h_mc1b) h_mc1b->SetBinError(i, 1.1*h_mc1b->GetBinError(i));
	    if (h_mc1c) h_mc1c->SetBinError(i, 1.1*h_mc1c->GetBinError(i));
	    if (h_mc1t) h_mc1t->SetBinError(i, 1.1*h_mc1t->GetBinError(i));
	    h_mc2->SetBinError(i, 1.1*h_mc2->GetBinError(i));
	    h_mc3->SetBinError(i, 1.1*h_mc3->GetBinError(i));
	    h_mc4->SetBinError(i, 1.1*h_mc4->GetBinError(i));
	    if (h_mc5) h_mc5->SetBinError(i, 1.1*h_mc5->GetBinError(i));
	    h_mc6->SetBinError(i, 1.1*h_mc6->GetBinError(i));
	    h_mc7->SetBinError(i, 1.1*h_mc7->GetBinError(i));
	    h_mc8->SetBinError(i, 1.1*h_mc8->GetBinError(i));
	  }
	}

	if (doBkg) {
	  h_data->Add(h_mc8, -1.);
	  h_data->Add(h_mc7, -1.);
	  h_data->Add(h_mc6, -1.);
	  if (h_mc5) h_data->Add(h_mc5, -1.);
	  h_data->Add(h_mc4, -1.);
	  h_data->Add(h_mc3, -1.);
	  h_data->Add(h_mc2, -1.);
	  h_data->Add(h_mc1t, -1.);
	}

	if (useFitResults && !doFit) {
	  if (title.find("_bb")!=string::npos) {
	    if (irun==6) {
	      if (h_mc1b) h_mc1b->Scale(c3_b+0.1*ec3_b);
	      if (h_mc1c) h_mc1c->Scale(c3_r+0.1*ec3_r);
	    } else {
	      if (h_mc1b) h_mc1b->Scale(c3_b);
	      if (h_mc1c) h_mc1c->Scale(c3_r);
	    }
	  } else if (title.find("_b")!=string::npos) {
	    if (irun==6) {
	      if (h_mc1b) h_mc1b->Scale(c2_b+0.1*ec2_b);
	      if (h_mc1c) h_mc1c->Scale(c2_r+0.1*ec2_r);
	    } else {
	      if (h_mc1b) h_mc1b->Scale(c2_b);
	      if (h_mc1c) h_mc1c->Scale(c2_r);
	    }
	  } else {
	    if (irun==6) {
	      if (h_mc1b) h_mc1b->Scale(c1_b+0.1*ec1_b);
	      if (h_mc1c) h_mc1c->Scale(c1_r+0.1*ec1_r);
	    } else {
	      if (h_mc1b) h_mc1b->Scale(c1_b);
	      if (h_mc1c) h_mc1c->Scale(c1_r);
	    }
	  }
	}

	h_data = rebin(h_data);
	h_mc1 = rebin(h_mc1);
	if (h_mc1b) h_mc1b = rebin(h_mc1b);
	if (h_mc1c) h_mc1c = rebin(h_mc1c);
	if (h_mc1t) h_mc1t = rebin(h_mc1t);
	h_mc2 = rebin(h_mc2);
	h_mc3 = rebin(h_mc3);
	h_mc4 = rebin(h_mc4);
	if (h_mc5) h_mc5 = rebin(h_mc5);
	h_mc6 = rebin(h_mc6);
	h_mc7 = rebin(h_mc7);
	h_mc8 = rebin(h_mc8);

	TVirtualFitter::SetDefaultFitter("Minuit2");
	TVirtualFitter* fitter=0;
	if (doFit==1) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc8, -1.);
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
	    if (h_mc5) h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc1t, -1.);
	  }
	  h_data_fit->Add(h_mc1, -1.);
	  if (h_mc1b) h_data_fit->Add(h_mc1b, -1.);
	  if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
	  h_mc_fit0 = h_mc2;
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    float e = TMath::Power(h_data_fit->GetBinError(i),2);
	    h_data_fit->SetBinError(i, TMath::Sqrt(e));
	  }
	  fitter = TVirtualFitter::Fitter(0, 1);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(ttbar)", 1.00, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	}
	if (doFit==2) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc8, -1.);
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc1t, -1.);
	  }
	  h_data_fit->Add(h_mc1, -1.);
	  if (h_mc1b) h_data_fit->Add(h_mc1b, -1.);
	  if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
	  h_mc_fit0 = h_mc2;
	  h_mc_fit1 = h_mc5;
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    float e = TMath::Power(h_data_fit->GetBinError(i),2);
	    h_data_fit->SetBinError(i, TMath::Sqrt(e));
	  }
	  fitter = TVirtualFitter::Fitter(0, 2);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(ttbar)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(qcd)", 1.00, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, -1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, -1.);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  if (h_mc1b) h_mc1b->Scale(fitter->GetParameter(0));
	  if (h_mc1c) h_mc1c->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
	}
	if (doFit==3) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc8, -1.);
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
	    if (h_mc5) h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc2, -1.);
	    h_data_fit->Add(h_mc1t, -1.);
	  }
	  h_mc_fit0 = h_mc1;
	  h_mc_fit1 = h_mc1b;
	  h_mc_fit2 = h_mc1c;
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    bool skip = false;
	    if (title=="w_SVTX_mass") {
	      if (h_data_fit->GetXaxis()->GetBinCenter(i) < 0.2) {
	        skip = true;
	      }
	    }
	    if (skip) {
	      h_data->SetBinContent(i, 0);
	      h_data->SetBinError(i, 0);
	      h_data_fit->SetBinContent(i, 0);
	      h_data_fit->SetBinError(i, 0);
	      h_mc1->SetBinContent(i, 0);
	      h_mc1->SetBinError(i, 0);
	      if (h_mc1b) h_mc1b->SetBinContent(i, 0);
	      if (h_mc1b) h_mc1b->SetBinError(i, 0);
	      if (h_mc1c) h_mc1c->SetBinContent(i, 0);
	      if (h_mc1c) h_mc1c->SetBinError(i, 0);
	      if (h_mc1t) h_mc1t->SetBinContent(i, 0);
	      if (h_mc1t) h_mc1t->SetBinError(i, 0);
	      h_mc2->SetBinContent(i, 0);
	      h_mc2->SetBinError(i, 0);
	      h_mc3->SetBinContent(i, 0);
	      h_mc3->SetBinError(i, 0);
	      h_mc4->SetBinContent(i, 0);
	      h_mc4->SetBinError(i, 0);
	      if (h_mc5) {
	        h_mc5->SetBinContent(i, 0);
	        h_mc5->SetBinError(i, 0);
	      }
	      h_mc6->SetBinContent(i, 0);
	      h_mc6->SetBinError(i, 0);
	      h_mc7->SetBinContent(i, 0);
	      h_mc7->SetBinError(i, 0);
	      h_mc8->SetBinContent(i, 0);
	      h_mc8->SetBinError(i, 0);
	    }
	  }
	  fitter = TVirtualFitter::Fitter(0, 3);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(uds)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(b)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(2, "c(c)", 1.00, 0.01, 0.00, 100.00);
	  if (irun==99) {
	    fitter->FixParameter(0);
	    fitter->FixParameter(2);
	  }
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
	  h_mc_fit2->Scale(fitter->GetParameter(2));
	}
	if (doFit==4) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc8, -1.);
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
	    if (h_mc5) h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc2, -1.);
	    h_data_fit->Add(h_mc1t, -1.);
	  }
	  h_mc_fit0 = h_mc1;
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, 1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, 1.);
	  fitter = TVirtualFitter::Fitter(0, 1);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(W+jets)", 1.00, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, -1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, -1.);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  if (h_mc1b) h_mc1b->Scale(fitter->GetParameter(0));
	  if (h_mc1c) h_mc1c->Scale(fitter->GetParameter(0));
	}
	if (doFit==5) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc8, -1.);
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc1t, -1.);
	  }
	  h_data_fit->Add(h_mc1, -1.);
	  if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
	  h_mc_fit0 = h_mc1b;
	  h_mc_fit1 = h_mc2;
	  h_mc_fit2 = h_mc5;
	  if (title.find("_bb")!=string::npos) {
	    mc_fit1 = c3_t;
	    e_mc_fit1 = ec3_t;
	  } else if (title.find("_b")!=string::npos) {
	    mc_fit1 = c2_t;
	    e_mc_fit1 = ec2_t;
	  } else {
	    mc_fit1 = c1_t;
	    e_mc_fit1 = ec1_t;
	  }
	  fitter = TVirtualFitter::Fitter(0, 3);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(W+b-jets)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(ttbar)", mc_fit1, 0.01, 0.00, 100.00);
	  fitter->SetParameter(2, "c(qcd)", 1.00, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
	  h_mc_fit2->Scale(fitter->GetParameter(2));
	}
	if (doFit==6) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc8, -1.);
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc1t, -1.);
	  }
	  h_data_fit->Add(h_mc1, -1.);
	  if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
	  h_mc_fit0 = h_mc1b;
	  h_mc_fit1 = h_mc2;
	  h_mc_fit2 = h_mc5;
	  mc_fit1 = 1.0;
	  e_mc_fit1 = eXsec_tt/Xsec_tt;
	  if (title.find("_bb")!=string::npos) {
	    mc_fit3 = c3_t;
	    e_mc_fit3 = ec3_t;
	  } else if (title.find("_b")!=string::npos) {
	    mc_fit3 = c2_t;
	    e_mc_fit3 = ec2_t;
	  } else {
	    mc_fit3 = c1_t;
	    e_mc_fit3 = ec1_t;
	  }
	  fitter = TVirtualFitter::Fitter(0, 4);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(W+b-jets)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(ttbar)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(2, "c(qcd)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(3, "c(scale)", mc_fit3, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0)*fitter->GetParameter(3));
	  h_mc_fit1->Scale(fitter->GetParameter(1)*fitter->GetParameter(3));
	  h_mc_fit2->Scale(fitter->GetParameter(2));
	}
	if (doFit==7) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  h_mc_fit0 = h_mc1b;
	  h_mc_fit1 = h_mc5;
	  h_mc_fit2 = (TH1F*)h_mc1->Clone("h_mc_fit");
	  h_mc_fit2->Add(h_mc1, 1.);
	  if (h_mc1c) h_mc_fit2->Add(h_mc1c, 1.);
	  if (!doBkg) {
	    h_mc_fit2->Add(h_mc8, 1.);
	    h_mc_fit2->Add(h_mc7, 1.);
	    h_mc_fit2->Add(h_mc6, 1.);
	    h_mc_fit2->Add(h_mc4, 1.);
	    h_mc_fit2->Add(h_mc3, 1.);
	    h_mc_fit2->Add(h_mc2, 1.);
	    h_mc_fit2->Add(h_mc1t, 1.);
	  }
	  if (title.find("_bb")!=string::npos) {
	    mc_fit2 = c3_t;
	    e_mc_fit2 = ec3_t;
	  } else if (title.find("_b")!=string::npos) {
	    mc_fit2 = c2_t;
	    e_mc_fit2 = ec2_t;
	  } else {
	    mc_fit2 = c1_t;
	    e_mc_fit2 = ec1_t;
	  }
	  fitter = TVirtualFitter::Fitter(0, 3);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(W+b-jets)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(qcd)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(2, "c(scale)", mc_fit2, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc1b->Scale(fitter->GetParameter(0)*fitter->GetParameter(2));
	  h_mc5->Scale(fitter->GetParameter(1));
	  h_mc1->Scale(fitter->GetParameter(2));
	  h_mc1t->Scale(fitter->GetParameter(2));
	  if (h_mc1c) h_mc1c->Scale(fitter->GetParameter(2));
	  h_mc2->Scale(fitter->GetParameter(2));
	  h_mc3->Scale(fitter->GetParameter(2));
	  h_mc4->Scale(fitter->GetParameter(2));
	  h_mc6->Scale(fitter->GetParameter(2));
	  h_mc7->Scale(fitter->GetParameter(2));
	  h_mc8->Scale(fitter->GetParameter(2));
	}

	if (printYield) {
	  cout << "******************" << endl;
	  cout << "POST FIT YIELDS:" << endl;
	  cout << "data: " << h_data->Integral(0,h_data->GetNbinsX()+1) << endl;
	  if (h_mc1b && h_mc1c && h_mc1t) cout << "Wjets: " << h_mc1->Integral(0,h_mc1->GetNbinsX()+1)
					    + h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1)
					    + h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1)
					    + h_mc1t->Integral(0,h_mc1t->GetNbinsX()+1) << endl;
	  if (h_mc1b) cout << "Wb: " << h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1) << endl;
	  if (h_mc1c) cout << "Wc: " << h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1) << endl;
	  if (h_mc1t) cout << "Wtau: " << h_mc1t->Integral(0,h_mc1t->GetNbinsX()+1) << endl;
	  if (h_mc1b && h_mc1c && h_mc1t) cout << "Wlight: " << h_mc1->Integral(0,h_mc1->GetNbinsX()+1) << endl;
	  cout << "TTbar: " << h_mc2->Integral(0,h_mc2->GetNbinsX()+1) << endl;
	  cout << "ZZ: " << h_mc3->Integral(0,h_mc3->GetNbinsX()+1) << endl;
	  cout << "WZ: " << h_mc4->Integral(0,h_mc4->GetNbinsX()+1) << endl;
	  if (h_mc5) cout << "QCD: " << h_mc5->Integral(0,h_mc5->GetNbinsX()+1) << endl;
	  cout << "WW: " << h_mc6->Integral(0,h_mc6->GetNbinsX()+1) << endl;
	  cout << "VV: " << h_mc3->Integral(0,h_mc3->GetNbinsX()+1)
	    + h_mc4->Integral(0,h_mc4->GetNbinsX()+1)
	    + h_mc6->Integral(0,h_mc6->GetNbinsX()+1) << endl;
	  cout << "DY: " << h_mc7->Integral(0,h_mc7->GetNbinsX()+1) << endl;
	  cout << "T: " << h_mc8->Integral(0,h_mc8->GetNbinsX()+1) << endl;
	  cout << "******************" << endl;
	}

	TH1F *ht = (TH1F*)h_mc1->Clone("ht");
	ht->Reset();
	if (h_mc1t) ht->Add(h_mc1t);
	if (!doBkg) {
	  ht->Add(h_mc8);
	  ht->Add(h_mc7);
	  ht->Add(h_mc6);
	  if (h_mc5) ht->Add(h_mc5);
	  ht->Add(h_mc4);
	  ht->Add(h_mc3);
	  ht->Add(h_mc2);
	}
	if (h_mc1b) ht->Add(h_mc1b);
	if (h_mc1c) ht->Add(h_mc1c);
	ht->Add(h_mc1);

	THStack *hs = new THStack("hs","");
	if (!doBkg) {
	  if (h_mc5) hs->Add(h_mc5);
	  hs->Add(h_mc7);
	  hs->Add(h_mc4);
	  hs->Add(h_mc3);
	  hs->Add(h_mc6);
	  hs->Add(h_mc8);
	  hs->Add(h_mc2);
	}
	if (h_mc1t) hs->Add(h_mc1t);
	hs->Add(h_mc1);
	if (h_mc1c) hs->Add(h_mc1c);
	if (h_mc1b) hs->Add(h_mc1b);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();
	//pad1->SetLogy();
	//if (title.find("MET")!=string::npos)   pad1->SetLogy(0);
	//if (title.find("mass_")!=string::npos) pad1->SetLogy(0);

	hs->Draw("HIST");
	hs->GetYaxis()->SetTitle("Events");
 	hs->GetXaxis()->SetLabelSize(0.08);
	hs->GetXaxis()->SetTitleOffset(0.7);
	hs->SetMinimum(8);

	h_data->Draw("EPX0SAMES");
	h_data->SetMarkerColor(kBlack);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize (1.0);
	//h_data->SetStats(0);

	hs->SetMaximum(1.2*TMath::Max(hs->GetMaximum(),h_data->GetMaximum()));

	TLegend *leg;
	if (doBkg) {
	  if (h_mc1c && h_mc1b) {
	    leg = new TLegend(0.62, 0.747, 0.88, 0.88);
	  } else if (h_mc1c || h_mc1b) {
	    leg = new TLegend(0.62, 0.780, 0.88, 0.88);
	  } else {
	    leg = new TLegend(0.62, 0.813, 0.88, 0.88);
	  }
	} else {
	  if (h_mc1c && h_mc1b) {
	    leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	  } else if (h_mc1c || h_mc1b) {
	    leg = new TLegend(0.62, 0.613, 0.88, 0.88);
	  } else {
	    leg = new TLegend(0.62, 0.647, 0.88, 0.88);
	  }
	}
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data,"W(#rightarrow e#nu)+jets","p");
	if (ilepton==2) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+jets","p");
	if (ilepton==3) leg->AddEntry(h_data,"W(#rightarrow e#nu)+jets [QCD]","p");
	if (ilepton==4) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+jets [QCD]","p");
	if (ilepton==5) leg->AddEntry(h_data,"W(#rightarrow e#nu)+jets [FWD]","p");
	if (ilepton==6) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+jets [FWD]","p");
	if (ilepton==7) leg->AddEntry(h_data,"W(#rightarrow e#nu)+jets [TOP]","p");
	if (ilepton==8) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+jets [TOP]","p");

	if (h_mc1b) leg->AddEntry(h_mc1b,"W+b-jets","f");
	if (h_mc1c) leg->AddEntry(h_mc1c,"W+c-jets","f");
	if (h_mc1c && h_mc1b) {
	  leg->AddEntry(h_mc1,"W+uds-jets","f");
	} else {
	  leg->AddEntry(h_mc1,"W+jets","f");
	}
	if (!doBkg) {
	  if (h_mc1t) leg->AddEntry(h_mc1t,"W(#rightarrow #tau#nu)+jets","f");
	  leg->AddEntry(h_mc2,"t#bar{t}","f");
	  leg->AddEntry(h_mc8,"t / #bar{t}", "f");
	  leg->AddEntry(h_mc6,"WW","f");
	  leg->AddEntry(h_mc3,"ZZ","f");
	  leg->AddEntry(h_mc4,"WZ","f");
	  leg->AddEntry(h_mc7,"Z+jets", "f");
	  if (h_mc5) leg->AddEntry(h_mc5,"QCD","f");
	}
	leg->Draw();

	pad1->Update();
	c1->Update();

	c1->cd();

	TH1F *h_ratio = (TH1F*)h_data->Clone("h_ratio");

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	h_ratio->SetTitle("");
	h_ratio->SetStats(0);

	if (title=="w_jetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("jet multiplicity");
	} else if (title=="w_mass_ee"||title=="w_mass_mm") {
	  h_ratio->GetXaxis ()->SetTitle("invariant mass [GeV/c^{2}]");
	} else if (title=="w_first_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="w_first_ele_pt") {
	  h_ratio->GetXaxis ()->SetTitle("electron p_{T} [GeV/c]");
	} else if (title=="w_first_jet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("jet p_{T} [GeV/c]");
	} else if (title=="w_secondvtx_N") {
	  h_ratio->GetXaxis ()->SetTitle("CSV discriminator");
	} else if (title=="h_recoVTX"||title=="w_recoVTX") {
	  h_ratio->GetXaxis ()->SetTitle("Number of offline vertices");
	} else if (title=="w_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="w_bjetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("b-jet multiplicity");
	} else if (title=="w_first_bjet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet p_{T}");
	} else if (title=="w_first_bjet_eta") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet #eta [rad]");
	} else if (title=="w_second_bjet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet p_{T}");
	} else if (title=="w_second_bjet_eta") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet #eta [rad]");
	} else if (title=="w_second_bjet_mass") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet mass [GeV/c^{2}");
	}
	if (title=="w_MET") {
	  h_ratio->GetXaxis ()->SetTitle("MET [GeV/c]");
	} else if (title=="w_MET_sign") {
	  h_ratio->GetXaxis ()->SetTitle("MET Significance [GeV/c]");
	} else if (title=="w_Ht") {
	  h_ratio->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	} else if (title=="w_mass_wenu_blepton"||title=="w_mass_wmnu_blepton") {
	  h_ratio->GetXaxis ()->SetTitle("lepton b-jet invariant mass [GeV/c^{2}]");
	} else if (title=="w_pt_W_wenu"||title=="w_pt_W_wmnu") {
	  h_ratio->GetXaxis ()->SetTitle("W boson p_{T} [GeV/c]");
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_mass_ee"||title=="w_mm_mass") {
	  h_ratio->GetXaxis ()->SetTitle("Z mass [GeV/c^{2}]");
	} else if (title=="w_first_jet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("leading jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta") {
	  h_ratio->GetXaxis ()->SetTitle("leading jet #eta");
	} else if (title=="w_first_jet_mass") {
	  h_ratio->GetXaxis ()->SetTitle("leading jet mass [GeV/c^{2}]");
	} else if (title=="w_second_jet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("subleading jet p_{T} [GeV/c]");
	} else if (title=="w_second_jet_eta") {
	  h_ratio->GetXaxis ()->SetTitle("subleading jet #eta");
	} else if (title=="w_second_jet_mass") {
	  h_ratio->GetXaxis ()->SetTitle("sub-leading jet mass [GeV/c^{2}]");
	} else if (title=="w_delta_wenu"||title=="w_delta_wenu") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(lj) [rad]");
	} else if (title=="w_deltaR_wenu"||title=="w_deltaR_wenu") {
	  h_ratio->GetXaxis ()->SetTitle("#DeltaR(lj) [rad]");
	} else if (title=="SVTX_mass_jet"||title=="SVTX_mass_trk"||title=="SVTX_mass") {
	  h_ratio->GetXaxis ()->SetTitle("SV mass [GeV/c^{2}]");
	} else if (title=="w_BJP"||title=="w_JBP") {
	  h_ratio->GetXaxis ()->SetTitle("JP Discriminator");
	}
	if (title=="w_MET_b") {
	  h_ratio->GetXaxis ()->SetTitle("MET [GeV/c]");
	} else if (title=="w_MET_sign_b") {
	  h_ratio->GetXaxis ()->SetTitle("MET Significance [GeV/c]");
	} else if (title=="w_Ht_b") {
	  h_ratio->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	} else if (title=="w_mass_wenu_blepton_b"||title=="w_mass_wmnu_blepton_b") {
	  h_ratio->GetXaxis ()->SetTitle("lepton b-jet invariant mass [GeV/c^{2}]");
	} else if (title=="w_pt_W_wenu_b"||title=="w_pt_W_wmnu_b") {
	  h_ratio->GetXaxis ()->SetTitle("W boson p_{T} [GeV/c]");
	} else if (title=="w_pt_Z_ee_b"||title=="w_pt_Z_mm_b") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_mass_ee_b"||title=="w_mm_mass_b") {
	  h_ratio->GetXaxis ()->SetTitle("Z mass [GeV/c^{2}]");
	} else if (title=="w_first_jet_pt_b") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta_b") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet #eta");
	} else if (title=="w_first_jet_mass_b") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet mass [GeV/c^{2}]");
	} else if (title=="w_second_jet_pt_b") {
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet p_{T} [GeV/c]");
	} else if (title=="w_second_jet_eta_b") {
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet #eta");
	} else if (title=="w_second_jet_mass_b") {
	  h_ratio->GetXaxis ()->SetTitle("sub-leading b-jet mass [GeV/c^{2}]");
	} else if (title=="w_delta_wenu_b"||title=="w_delta_wenu_b") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(lb) [rad]");
	} else if (title=="w_deltaR_wenu_b"||title=="w_deltaR_wenu_b") {
	  h_ratio->GetXaxis ()->SetTitle("#DeltaR(lb) [rad]");
	} else if (title=="SVTX_mass_jet_b"||title=="SVTX_mass_trk_b"||title=="SVTX_mass_b") {
	  h_ratio->GetXaxis ()->SetTitle("SV mass [GeV/c^{2}]");
	} else if (title=="w_BJP_b"||title=="w_JBP_b") {
	  h_ratio->GetXaxis ()->SetTitle("JP Discriminator");
	}
	if (title=="w_MET_bb") {
	  h_ratio->GetXaxis ()->SetTitle("MET [GeV/c]");
	} else if (title=="w_MET_sign_bb") {
	  h_ratio->GetXaxis ()->SetTitle("MET Significance [GeV/c]");
	} else if (title=="w_Ht_bb") {
	  h_ratio->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	} else if (title=="w_mass_wenu_blepton_bb"||title=="w_mass_wmnu_blepton_bb") {
	  h_ratio->GetXaxis ()->SetTitle("lepton b-jet invariant mass [GeV/c^{2}]");
	} else if (title=="w_pt_W_wenu_bb"||title=="w_pt_W_wmnu_bb") {
	  h_ratio->GetXaxis ()->SetTitle("W boson p_{T} [GeV/c]");
	} else if (title=="w_first_jet_pt_bb") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta_bb") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet #eta");
	} else if (title=="w_first_jet_mass_bb") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet mass [GeV/c^{2}]");
	} else if (title=="w_second_jet_pt_bb") {
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet p_{T} [GeV/c]");
	} else if (title=="w_second_jet_eta_bb") {
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet #eta");
	} else if (title=="w_second_jet_mass_bb") {
	  h_ratio->GetXaxis ()->SetTitle("sub-leading b-jet mass [GeV/c^{2}]");
	} else if (title=="w_delta_wenu_bb"||title=="w_delta_wenu_bb") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(lb) [rad]");
	} else if (title=="w_deltaR_wenu_bb"||title=="w_deltaR_wenu_bb") {
	  h_ratio->GetXaxis ()->SetTitle("#DeltaR(lb) [rad]");
	} else if (title=="SVTX_mass_jet_bb"||title=="SVTX_mass_trk_bb"||title=="SVTX_mass_bb") {
	  h_ratio->GetXaxis ()->SetTitle("SV mass [GeV/c^{2}]");
	} else if (title=="w_BJP_bb"||title=="w_JBP_bb") {
	  h_ratio->GetXaxis ()->SetTitle("JP Discriminator");
	}

	h_ratio->GetXaxis()->SetTitleOffset(0.9);
 	h_ratio->GetXaxis()->SetTitleSize(0.1);
	h_ratio->GetXaxis()->SetLabelFont(42);
 	h_ratio->GetXaxis()->SetLabelSize(0.08);
	h_ratio->GetXaxis()->SetTitleFont(42);
	h_ratio->GetYaxis()->SetTitle("Data/MC");
	h_ratio->GetYaxis()->SetNdivisions(505);
	h_ratio->GetYaxis()->SetTitleSize(0.09);
	h_ratio->GetYaxis()->SetLabelSize(0.08);
	h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
	h_ratio->GetYaxis()->SetTitleOffset(0.4);
	h_ratio->Divide(ht);
	h_ratio->SetMarkerStyle(20);
	h_ratio->Draw("EPX0");

	TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
	OLine->SetLineColor(kRed);
	OLine->SetLineWidth(2);
	OLine->Draw();

	c1->cd();

 	TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
	latexLabel->Draw("same");

	if (doFit) {
	  TLatex *fitLabel = new TLatex();
	  fitLabel->SetTextSize(0.0275);
	  fitLabel->SetTextFont(42);
	  fitLabel->SetLineWidth(2);
	  fitLabel->SetNDC();
	  char buff[100];
	  if (doFit==1) {
	    sprintf(buff, "c_{ttbar} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	  }
	  if (doFit==2) {
	    sprintf(buff, "c_{ttbar} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    sprintf(buff, "c_{qcd} = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	  }
	  if (doFit==3) {
	    sprintf(buff, "c_{uds} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.38, 0.48, buff);
	    sprintf(buff, "c_{b}   = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.38, 0.43, buff);
	    sprintf(buff, "c_{c}   = %5.3f #pm %5.3f", fitter->GetParameter(2), fitter->GetParError(2));
	    fitLabel->DrawLatex(0.38, 0.38, buff);
	    float f_uds = 100*h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)/(h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)+h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)+h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1));
	    float ef_uds = f_uds*(fitter->GetParError(0)/fitter->GetParameter(0));
	    sprintf(buff, "f_{uds} = %4.1f #pm %3.1f %%", f_uds, ef_uds);
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    float f_b = 100*h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)/(h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)+h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)+h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1));
	    float ef_b = f_b*(fitter->GetParError(1)/fitter->GetParameter(1));
	    sprintf(buff, "f_{b}   = %4.1f #pm %3.1f %%", f_b, ef_b);
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	    float f_c = 100*h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1)/(h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)+h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)+h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1));
	    float ef_c = f_c*(fitter->GetParError(2)/fitter->GetParameter(2));
	    sprintf(buff, "f_{c}   = %4.1f #pm %3.1f %%", f_c, ef_c);
	    fitLabel->DrawLatex(0.68, 0.38, buff);
	  }
	  if (doFit==4) {
	    sprintf(buff, "c_{W+jets} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	  }
	  if (doFit==5) {
	    sprintf(buff, "c_{W+b-jets} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    sprintf(buff, "c_{ttbar} = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	    sprintf(buff, "c_{qcd} = %5.3f #pm %5.3f", fitter->GetParameter(2), fitter->GetParError(2));
	    fitLabel->DrawLatex(0.68, 0.38, buff);
	  }
	  if (doFit==6) {
	    sprintf(buff, "c_{W+b-jets} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.53, buff);
	    sprintf(buff, "c_{ttbar} = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    sprintf(buff, "c_{qcd} = %5.3f #pm %5.3f", fitter->GetParameter(2), fitter->GetParError(2));
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	    sprintf(buff, "c_{scale} = %5.3f #pm %5.3f", fitter->GetParameter(3), fitter->GetParError(3));
	    fitLabel->DrawLatex(0.68, 0.38, buff);
	  }
	  if (doFit==7) {
	    sprintf(buff, "c_{W+b-jets} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    sprintf(buff, "c_{qcd} = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	    sprintf(buff, "c_{scale} = %5.3f #pm %5.3f", fitter->GetParameter(2), fitter->GetParError(2));
	    fitLabel->DrawLatex(0.68, 0.38, buff);
	  }
	}

	if (plot) {
	  if (doBkg) title = title + "_doBkg";
	  if (doFit) title = title + "_doFit";
	  ofstream out;
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==3) {
	    gSystem->mkdir((path + "/electronsQCD/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electronsQCD/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electronsQCD/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==4) {
	    gSystem->mkdir((path + "/muonsQCD/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muonsQCD/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/muonsQCD/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==5) {
	    gSystem->mkdir((path + "/electronsFWD/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electronsFWD/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electronsFWD/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==6) {
	    gSystem->mkdir((path + "/muonsFWD/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muonsFWD/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/muonsFWD/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==7) {
	    gSystem->mkdir((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==8) {
	    gSystem->mkdir((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/muonsTOP/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (doFit==1) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out.close();
	  }
	  if (doFit==2) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out << fitter->GetParameter(1) << " " << fitter->GetParError(1) << endl;
	    out.close();
	  }
	  if (doFit==3) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out << fitter->GetParameter(1) << " " << fitter->GetParError(1) << endl;
	    out << fitter->GetParameter(2) << " " << fitter->GetParError(2) << endl;
	    out.close();
	  }
	  if (doFit==4) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out.close();
	  }
	  if (doFit==5) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out << fitter->GetParameter(1) << " " << fitter->GetParError(1) << endl;
	    out << fitter->GetParameter(2) << " " << fitter->GetParError(2) << endl;
	    out << 1.0 << " " << 0.0 << endl;
	    out << 1.0 << " " << 0.0 << endl;
	    out.close();
	  }
	  if (doFit==6) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out << fitter->GetParameter(1) << " " << fitter->GetParError(1) << endl;
	    out << fitter->GetParameter(2) << " " << fitter->GetParError(2) << endl;
	    out << fitter->GetParameter(3) << " " << fitter->GetParError(3) << endl;
	    out << 1.0 << " " << 0.0 << endl;
	    out.close();
	  }
	  if (doFit==7) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out << 1.0 << " " << 0.0 << endl;
	    out << fitter->GetParameter(1) << " " << fitter->GetParError(1) << endl;
	    out << 1.0 << " " << 0.0 << endl;
	    out << fitter->GetParameter(2) << " " << fitter->GetParError(2) << endl;
	    out.close();
	  }
	}
}

