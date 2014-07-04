#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v09.h"

#include "fixrange.C"
#include "rebin.C"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

void DataMCComp2(int irun=0, string title="", int plot=0, int ilepton=1, int unfold=0) {

int useBinnedEfficiency=0; // use average efficiencies
//int useBinnedEfficiency=1; // use bin-by-bin efficiencies

//int useFitResults=0; // use MC predictions for c_b, c_c, c_uds, c_t
int useFitResults=1;  // use fit results for c_b, c_c, c_uds, c_t

//int useFitResults2=0;  // do not use constrained fit results for c_qcd, c_t
int useFitResults2=1;  // use constrained fit results for c_qcd, c_t

//int drawInclusive = 0; // do not plot the "inclusive" histogram
int drawInclusive = 1; // do plot the "inclusive" histogram

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

	if (gROOT->GetVersionInt() >= 53401) {
	  //gROOT->GetColor(kRed)->SetAlpha(0.5);
	  gROOT->GetColor(kRed)->SetAlpha(0.0);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	}

	/* efficiency */

	double e_W=1.0;
	double ee_W=0.0;
	double e_Wb=1.0;
	double ee_Wb=0.0;

	/* purity */

	double c1_b=1.0;
	double ec1_b=0.0;
	double c1_c=1.0;
	double ec1_c=0.0;
	double c1_uds=1.0;
	double ec1_uds=0.0;
	double c2_b=1.0;
	double ec2_b=0.0;
	double c2_c=1.0;
	double ec2_c=0.0;
	double c2_uds=1.0;
	double ec2_uds=0.0;
	double c3_b=1.0;
	double ec3_b=0.0;
	double c3_c=1.0;
	double ec3_c=0.0;
	double c3_uds=1.0;
	double ec3_uds=0.0;

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

	ifstream in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12;
	if (ilepton==1) {
	  in1.open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	  in2.open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	  in3.open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	  if (useFitResults) {
	    in4.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	    in5.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	    in6.open((path + "/electronsTOP/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
	    //in7.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_BJP_doFit" + ".dat").c_str());
	    //in8.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_BJP_b_doFit" + ".dat").c_str());
	    //in9.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_BJP_bb_doFit" + ".dat").c_str());
	  }
	  if (useFitResults2) {
	    in10.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_wide_doFit" + ".dat").c_str());
	    in11.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_b_wide_doFit" + ".dat").c_str());
	    in12.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_mt_wenu_bb_wide_doFit" + ".dat").c_str());
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
	    //in7.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_BJP_doFit" + ".dat").c_str());
	    //in8.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_BJP_b_doFit" + ".dat").c_str());
	    //in9.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_BJP_bb_doFit" + ".dat").c_str());
	  }
	  if (useFitResults2) {
	    in10.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_wide_doFit" + ".dat").c_str());
	    in11.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_b_wide_doFit" + ".dat").c_str());
	    in12.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_mt_wmnu_bb_wide_doFit" + ".dat").c_str());
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
	  //in7 >> c1_uds >> ec1_uds;
	  //in7 >> c1_b >> ec1_b;
	  //in7 >> c1_c >> ec1_c;
	  //in8 >> c2_uds >> ec2_uds;
	  //in8 >> c2_b >> ec2_b;
	  //in8 >> c2_c >> ec2_c;
	  //in9 >> c3_uds >> ec3_uds;
	  //in9 >> c3_b >> ec3_b;
	  //in9 >> c3_c >> ec3_c;
	  //in7.close();
	  //in8.close();
	  //in9.close();
	}
	if (useFitResults2) {
	  double c;
	  double ec;
	  c = 1.0;
	  ec = 0.0;
	  in10 >> c >> ec; // drop first line
	  in10 >> c >> ec;
	  c1_t = c1_t * c; ec1_t = ec;
	  in10 >> c >> ec;
	  c1_qcd = c1_qcd * c; ec1_qcd = ec;
	  c = 1.0;
	  ec = 0.0;
	  in11 >> c >> ec; // drop first line
	  in11 >> c >> ec;
	  c2_t = c2_t * c; ec2_t = ec;
	  in11 >> c >> ec;
	  c2_qcd = c2_qcd * c; ec2_qcd = ec;
	  c = 1.0;
	  ec = 0.0;
	  in12 >> c >> ec; // drop first line
	  in12 >> c >> ec;
	  c3_t = c3_t * c; ec3_t = ec;
	  in12 >> c >> ec;
	  c3_qcd = c3_qcd * c; ec3_qcd = ec;
	  in10.close();
	  in11.close();
	  in12.close();
	}

	double Lumi2012=0;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

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

	TFile *data=0;
	if (ilepton==1) data = TFile::Open((path + "/" + version + "/" + "SingleElectron_2012_merge.root").c_str());
	if (ilepton==2) data = TFile::Open((path + "/" + version + "/" + "SingleMu_2012_merge.root").c_str());

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "Wj_merge.root").c_str());
	TFile *mcg = TFile::Open((path + "/" + version + "/" + "Wj_gen_merge.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar_merge.root").c_str());
	TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
	TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
	TFile *mc5 = 0;
//	mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
	TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
	TFile *mc7 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	TFile *mc8 = TFile::Open((path + "/" + version + "/" + "T_merge.root").c_str());

	string title_b = title;

	if (title.find("_bjet_")!=string::npos) {
	  title.erase(title.find("_bjet_")+1, 1);
	} else {
	  title_b = title_b + "b";
        }

	if (ilepton==1) data->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) data->cd(("demoMuo"+postfix).c_str());
	TH1F* h_data=0;
	TH1F* h_data_b=0;
	if (unfold==0) {
	  h_data = (TH1F*)gDirectory->Get(title.c_str());
	  h_data_b = (TH1F*)gDirectory->Get(title_b.c_str());
	  if (!h_data_b) h_data_b = (TH1F*)h_data->Clone(); // protect against null pointers!
	}
	if (unfold==1) {
          if (ilepton==1) {
	    TFile f((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.root").c_str());
	    TFile f_b((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title_b + "_unfolding.root").c_str());
	    h_data = (TH1F*)f.Get(title.c_str())->Clone();
	    h_data_b = (TH1F*)f_b.Get(title_b.c_str())->Clone();
	    h_data->SetDirectory(0);
	    h_data_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
          }
          if (ilepton==2) {
	    TFile f((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.root").c_str());
	    TFile f_b((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title_b + "_unfolding.root").c_str());
	    h_data = (TH1F*)f.Get(title.c_str())->Clone();
	    h_data_b = (TH1F*)f_b.Get(title_b.c_str())->Clone();
	    h_data->SetDirectory(0);
	    h_data_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
          }
	  h_data->SetStats(0);
	  h_data_b->SetStats(0);
	}

	if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1b = (TH1F*)gDirectory->Get(("b"+title.substr(1)).c_str());
	TH1F* h_mc1c = (TH1F*)gDirectory->Get(("c"+title.substr(1)).c_str());
	TH1F* h_mc1t = (TH1F*)gDirectory->Get(("t"+title.substr(1)).c_str());
	TH1F* h_mc1_b = (TH1F*)gDirectory->Get(title_b.c_str());
	TH1F* h_mc1b_b = (TH1F*)gDirectory->Get(("b"+title_b.substr(1)).c_str());
	TH1F* h_mc1c_b = (TH1F*)gDirectory->Get(("c"+title_b.substr(1)).c_str());
	TH1F* h_mc1t_b = (TH1F*)gDirectory->Get(("t"+title_b.substr(1)).c_str());

	if (ilepton==1) mcg->cd("demoEleGen");
	if (ilepton==2) mcg->cd("demoMuoGen");
	TH1F* h_mcg = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc2->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc2->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc3->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc3->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc3_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc4->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc4->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc4_b = (TH1F*)gDirectory->Get(title_b.c_str());

//	if (ilepton==1) mc5->cd(("demoEle"+postfix).c_str());
//	if (ilepton==2) mc5->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc5 = 0;
//	h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc5_b = 0;
//	h_mc5_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) {
	  mc5 = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str());
	  mc5->cd();
	  h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
	  mc5 = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + title_b + ".root").c_str());
	  mc5->cd();
	  h_mc5_b = (TH1F*)gDirectory->Get(title_b.c_str());
	}
	if (ilepton==2) {
	  mc5 = TFile::Open((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str());
	  mc5->cd();
	  h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
	  mc5 = TFile::Open((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + title_b + ".root").c_str());
	  mc5->cd();
	  h_mc5_b = (TH1F*)gDirectory->Get(title_b.c_str());
	}

	if (ilepton==1) mc6->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc6->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc6_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc7->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc7->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc7_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc8->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc8->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc8 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc8_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (unfold==0) {
	  h_data->Sumw2();
	  h_data_b->Sumw2();
	}

	h_mc1->Sumw2();
	if (h_mc1t) h_mc1t->Sumw2();
	if (h_mc1b) h_mc1b->Sumw2();
	if (h_mc1c) h_mc1c->Sumw2();
	h_mcg->Sumw2();
	h_mc2->Sumw2();
	h_mc3->Sumw2();
	h_mc4->Sumw2();
	if (h_mc5) h_mc5->Sumw2();
	h_mc6->Sumw2();
	h_mc7->Sumw2();
	h_mc8->Sumw2();

	h_mc1_b->Sumw2();
	if (h_mc1b_b) h_mc1b_b->Sumw2();
	if (h_mc1c_b) h_mc1c_b->Sumw2();
	if (h_mc1t_b) h_mc1t_b->Sumw2();
	h_mcg_b->Sumw2();
	h_mc2_b->Sumw2();
	h_mc3_b->Sumw2();
	h_mc4_b->Sumw2();
	if (h_mc5_b) h_mc5_b->Sumw2();
	h_mc6_b->Sumw2();
	h_mc7_b->Sumw2();
	h_mc8_b->Sumw2();

	if (irun==10) {
	  norm1 = norm1 + enorm1;
	  norm2 = norm2 + enorm2;
	  norm3 = norm3 + enorm3;
	  norm4 = norm4 + enorm4;
	  norm5 = norm5 + enorm5;
	  norm6 = norm6 + enorm6;
	  norm7 = norm7 + enorm7;
	  norm8 = norm8 + enorm8;
	}

	h_mc1->Scale(norm1);
	if (h_mc1t) h_mc1t->Scale(norm1);
	if (h_mc1b) h_mc1b->Scale(norm1);
	if (h_mc1c) h_mc1c->Scale(norm1);
	h_mcg->Scale(norm1);
	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
	if (h_mc5) h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);
	h_mc8->Scale(norm8);

	h_mc1_b->Scale(norm1);
	if (h_mc1b_b) h_mc1b_b->Scale(norm1);
	if (h_mc1c_b) h_mc1c_b->Scale(norm1);
	if (h_mc1t_b) h_mc1t_b->Scale(norm1);
	h_mcg_b->Scale(norm1);
	h_mc2_b->Scale(norm2);
	h_mc3_b->Scale(norm3);
	h_mc4_b->Scale(norm4);
	if (h_mc5_b) h_mc5_b->Scale(norm5);
	h_mc6_b->Scale(norm6);
	h_mc7_b->Scale(norm7);
	h_mc8_b->Scale(norm8);

        if (useFitResults) {
          if (title.find("_b")!=string::npos) {
            h_mc2->Scale(c2_t);
            if (irun==5) h_mc2->Scale((c2_t+ec2_t)/c2_t);
          } else {
            h_mc2->Scale(c1_t);
            if (irun==5) h_mc2->Scale((c1_t+ec1_t)/c1_t);
          }
          if (title_b.find("_bb")!=string::npos) {
            h_mc2_b->Scale(c3_t);
            if (irun==5) h_mc2_b->Scale((c3_t+ec3_t)/c3_t);
          } else if (title_b.find("_b")!=string::npos) {
            h_mc2_b->Scale(c2_t);
            if (irun==5) h_mc2_b->Scale((c2_t+ec2_t)/c2_t);
          }
        }

	if (h_mc5) {
          if (title.find("_b")!=string::npos) {
            h_mc5->Scale(c2_qcd);
            if (irun==15) h_mc5->Scale((c2_qcd+ec2_qcd)/c2_qcd);
          } else {
            h_mc5->Scale(c1_qcd);
            if (irun==15) h_mc5->Scale((c1_qcd+ec1_qcd)/c1_qcd);
          }
          if (title_b.find("_bb")!=string::npos) {
            h_mc5_b->Scale(c3_qcd);
            if (irun==15) h_mc5_b->Scale((c3_qcd+ec3_qcd)/c3_qcd);
          } else if (title_b.find("_b")!=string::npos) {
            h_mc5_b->Scale(c2_qcd);
            if (irun==15) h_mc5_b->Scale((c2_qcd+ec2_qcd)/c2_qcd);
          }
	}

	if (irun==13) {
	  for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	    h_mc1->SetBinError(i, 1.1*h_mc1->GetBinError(i));
	    if (h_mc1t) h_mc1t->SetBinError(i, 1.1*h_mc1t->GetBinError(i));
	    if (h_mc1b) h_mc1b->SetBinError(i, 1.1*h_mc1b->GetBinError(i));
	    if (h_mc1c) h_mc1c->SetBinError(i, 1.1*h_mc1c->GetBinError(i));
	    h_mc1_b->SetBinError(i, 1.1*h_mc1_b->GetBinError(i));
	    if (h_mc1b_b) h_mc1b_b->SetBinError(i, 1.1*h_mc1b_b->GetBinError(i));
	    if (h_mc1c_b) h_mc1c_b->SetBinError(i, 1.1*h_mc1c_b->GetBinError(i));
	    if (h_mc1t_b) h_mc1t_b->SetBinError(i, 1.1*h_mc1t_b->GetBinError(i));
	    h_mc2->SetBinError(i, 1.1*h_mc2->GetBinError(i));
	    h_mc2_b->SetBinError(i, 1.1*h_mc2_b->GetBinError(i));
	    h_mc3->SetBinError(i, 1.1*h_mc3->GetBinError(i));
	    h_mc3_b->SetBinError(i, 1.1*h_mc3_b->GetBinError(i));
	    h_mc4->SetBinError(i, 1.1*h_mc4->GetBinError(i));
	    h_mc4_b->SetBinError(i, 1.1*h_mc4_b->GetBinError(i));
	    if (h_mc5) h_mc5->SetBinError(i, 1.1*h_mc5->GetBinError(i));
	    if (h_mc5_b) h_mc5_b->SetBinError(i, 1.1*h_mc5_b->GetBinError(i));
	    h_mc6->SetBinError(i, 1.1*h_mc6->GetBinError(i));
	    h_mc6_b->SetBinError(i, 1.1*h_mc6_b->GetBinError(i));
	    h_mc7->SetBinError(i, 1.1*h_mc7->GetBinError(i));
	    h_mc7_b->SetBinError(i, 1.1*h_mc7_b->GetBinError(i));
	    h_mc8->SetBinError(i, 1.1*h_mc8->GetBinError(i));
	    h_mc8_b->SetBinError(i, 1.1*h_mc8_b->GetBinError(i));
	  }
	}

	if (unfold==0) {
	  h_data->Add(h_mc8, -1.);
	  h_data->Add(h_mc7, -1.);
	  h_data->Add(h_mc6, -1.);
	  if (h_mc5) h_data->Add(h_mc5, -1.);
	  h_data->Add(h_mc4, -1.);
	  h_data->Add(h_mc3, -1.);
	  h_data->Add(h_mc2, -1.);
	  h_data->Add(h_mc1t, -1.);

	  h_data_b->Add(h_mc8_b, -1.);
	  h_data_b->Add(h_mc7_b, -1.);
	  h_data_b->Add(h_mc6_b, -1.);
	  if (h_mc5_b) h_data_b->Add(h_mc5_b, -1.);
	  h_data_b->Add(h_mc4_b, -1.);
	  h_data_b->Add(h_mc3_b, -1.);
	  h_data_b->Add(h_mc2_b, -1.);
	  h_data_b->Add(h_mc1t_b, -1.);
	}

	TH1F *h_mc1uds = (TH1F*)h_mc1->Clone("h_mc1uds");
	if (h_mc1b) h_mc1uds->Add(h_mc1b, -1);
	if (h_mc1c) h_mc1uds->Add(h_mc1c, -1);
	if (h_mc1t) h_mc1uds->Add(h_mc1t, -1);
	for (int i=0; i<=h_mc1uds->GetNbinsX()+1; i++) {
	  float e = TMath::Power(h_mc1uds->GetBinError(i),2);
	  if (h_mc1b) e = e - TMath::Power(h_mc1b->GetBinError(i),2);
	  if (h_mc1c) e = e - TMath::Power(h_mc1c->GetBinError(i),2);
	  if (h_mc1t) e = e - TMath::Power(h_mc1t->GetBinError(i),2);
	  h_mc1uds->SetBinError(i, TMath::Sqrt(e));
	}

	TH1F *h_mc1uds_b = (TH1F*)h_mc1_b->Clone("h_mc1uds_b");
	if (h_mc1b_b) h_mc1uds_b->Add(h_mc1b_b, -1);
	if (h_mc1c_b) h_mc1uds_b->Add(h_mc1c_b, -1);
	if (h_mc1t_b) h_mc1uds_b->Add(h_mc1t_b, -1);
	for (int i=0; i<=h_mc1uds_b->GetNbinsX()+1; i++) {
	  float e = TMath::Power(h_mc1uds_b->GetBinError(i),2);
	  if (h_mc1b_b) e = e - TMath::Power(h_mc1b_b->GetBinError(i),2);
	  if (h_mc1c_b) e = e - TMath::Power(h_mc1c_b->GetBinError(i),2);
	  if (h_mc1t_b) e = e - TMath::Power(h_mc1t_b->GetBinError(i),2);
	  h_mc1uds_b->SetBinError(i, TMath::Sqrt(e));
	}

        if (useFitResults) {
	  if (title.find("_b")!=string::npos) {
	    if (h_mc1uds) {
	      h_mc1uds->Scale(c2_uds);
	      if (irun==6) h_mc1uds->Scale((c2_uds+ec2_uds)/c2_uds);
	    }
	    if (h_mc1b) {
	      h_mc1b->Scale(c2_b);
	      if (irun==6) h_mc1b->Scale((c2_b+ec2_b)/c2_b);
	    }
	    if (h_mc1c) {
	      h_mc1c->Scale(c2_c);
	      if (irun==6) h_mc1c->Scale((c2_c+ec2_c)/c2_c);
            }
	  } else {
	    if (h_mc1uds) {
	      h_mc1uds->Scale(c1_uds);
	      if (irun==6) h_mc1uds->Scale((c1_uds+ec1_uds)/c1_uds);
	    }
	    if (h_mc1b) {
	      h_mc1b->Scale(c1_b);
	      if (irun==6) h_mc1b->Scale((c1_b+ec1_b)/c1_b);
	    }
	    if (h_mc1c) {
	      h_mc1c->Scale(c1_c);
	      if (irun==6) h_mc1c->Scale((c1_c+ec1_c)/c1_c);
            }
	  }
	}

        if (useFitResults) {
	  if (title_b.find("_bb")!=string::npos) {
	    if (h_mc1uds_b) {
	      h_mc1uds_b->Scale(c3_uds);
	      if (irun==6) h_mc1uds_b->Scale((c3_uds+ec3_uds)/c3_uds);
	    }
	    if (h_mc1b_b) {
	      h_mc1b_b->Scale(c3_b);
	      if (irun==6) h_mc1b_b->Scale((c3_b+ec3_b)/c3_b);
	    }
	    if (h_mc1c_b) {
	      h_mc1c_b->Scale(c3_c);
	      if (irun==6) h_mc1c_b->Scale((c3_c+ec3_c)/c3_c);
            }
          } else {
	    if (h_mc1uds_b) {
	      h_mc1uds_b->Scale(c2_uds);
	      if (irun==6) h_mc1uds_b->Scale((c2_uds+ec2_uds)/c2_uds);
	    }
	    if (h_mc1b_b) {
	      h_mc1b_b->Scale(c2_b);
	      if (irun==6) h_mc1b_b->Scale((c2_b+ec2_b)/c2_b);
	    }
	    if (h_mc1c_b) {
	      h_mc1c_b->Scale(c2_c);
	      if (irun==6) h_mc1c_b->Scale((c2_c+ec2_c)/c2_c);
            }
	  }
	}

	if (unfold==0) {
	  h_data->Add(h_mc1c, -1.);
	  h_data->Add(h_mc1uds, -1.);
	  h_data_b->Add(h_mc1c_b, -1.);
	  h_data_b->Add(h_mc1uds_b, -1.);
	}

	TH1F *h_data_raw=0;
	TH1F *h_data_b_raw=0;
	if (unfold==0) {
	  h_data_raw = (TH1F*)h_data->Clone();
	  h_data_b_raw = (TH1F*)h_data_b->Clone();
	}

	if (useBinnedEfficiency==0) {
	  ifstream in13, in14;
	  if (ilepton==1) {
	    in13.open((path + "/electrons/" + version + "/" + subdir + "/efficiency/" + string(h_data->GetName()) + "_efficiency.dat").c_str());
	    in14.open((path + "/electrons/" + version + "/" + subdir + "/efficiency/" + string(h_data_b->GetName()) + "_efficiency.dat").c_str());
	  }
	  if (ilepton==2) {
	    in13.open((path + "/muons/" + version + "/" + subdir + "/efficiency/" + string(h_data->GetName()) + "_efficiency.dat").c_str());
	    in14.open((path + "/muons/" + version + "/" + subdir + "/efficiency/" + string(h_data_b->GetName()) + "_efficiency.dat").c_str());
	  }
	  in13 >> e_W >> ee_W;
	  in14 >> e_Wb >> ee_Wb;
	  in13.close();
	  in14.close();
	  e_W = e_W>0 ? e_W : 1;
	  e_Wb = e_Wb>0 ? e_Wb : 1;
	  if (unfold==0) {
	    h_data->Scale(1./e_W);
	    h_data_b->Scale(1./e_Wb);
	  }
	  h_mc1b->Scale(1./e_W);
	  h_mc1b_b->Scale(1./e_Wb);
	}

	if (useBinnedEfficiency==1) {
          if (ilepton==1) {
	    TFile f((path + "/electrons/" + version + "/" + subdir + "/efficiency/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	    TFile f_b((path + "/electrons/" + version + "/" + subdir + "/efficiency/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
	    TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	    TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	    h->SetDirectory(0);
	    h_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
	    if (unfold==0) {
	      h_data->Divide(h);
	      h_data_b->Divide(h_b);
	    }
	    h_mc1b->Divide(h);
	    h_mc1b_b->Divide(h_b);
          }
	  if (ilepton==2) {
	    TFile f((path + "/muons/" + version + "/" + subdir + "/efficiency/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	    TFile f_b((path + "/muons/" + version + "/" + subdir + "/efficiency/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
	    TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	    TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	    h->SetDirectory(0);
	    h_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
	    if (unfold==0) {
	      h_data->Divide(h);
	      h_data_b->Divide(h_b);
	    }
	    h_mc1b->Divide(h);
	    h_mc1b_b->Divide(h_b);
          }
	}

	h_data->Scale(1./Lumi2012, "width");
	h_data_b->Scale(1./Lumi2012, "width");
	h_mc1b->Scale(1./Lumi2012, "width");
	h_mc1b_b->Scale(1./Lumi2012, "width");

	h_mcg->Scale(1./Lumi2012, "width");
	h_mcg_b->Scale(1./Lumi2012, "width");

	if (unfold==0) {
	  h_data = fixrange(h_data);
	  h_data_b = fixrange(h_data_b);
	}
	h_mc1 = fixrange(h_mc1);
	h_mc1b_b = fixrange(h_mc1b_b);
	h_mcg = fixrange(h_mcg);
	h_mcg_b = fixrange(h_mcg_b);

	if (unfold==0) {
	  h_data = rebin(h_data);
	  h_data_b = rebin(h_data_b);
	}
	h_mc1 = rebin(h_mc1);
	h_mc1b_b = rebin(h_mc1b_b);
	h_mcg = rebin(h_mcg);
	h_mcg_b = rebin(h_mcg_b);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	h_mc1b_b->SetTitle("");
	h_mc1b_b->GetYaxis()->SetTitle("#sigma [pb]");

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

	h_mcg_b->SetLineColor(kGreen+2);
	h_mcg_b->SetLineWidth(2);
	h_mcg_b->SetFillColor(kGreen+2);
	h_mcg_b->SetMarkerColor(kGreen+2);

	h_data_b->GetYaxis()->SetTitleOffset(1.2);
	h_data_b->GetXaxis()->SetTitleOffset(1.3);
	h_data_b->SetMarkerColor(kBlack);
	h_data_b->SetLineColor(kBlack);
	h_data_b->SetMarkerStyle(24);
	h_data_b->SetMarkerSize(0.7);
	h_data_b->SetStats(0);

	TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	pad1->SetLogy();

	if (drawInclusive) {
	  h_mc1b_b->SetMaximum(4*h_data->GetMaximum());
	} else {
	  h_mc1b_b->SetMaximum(4*h_data_b->GetMaximum());
	}
	h_mc1b_b->SetMinimum(TMath::Max(0.002,0.25*h_data_b->GetBinContent(h_data_b->GetMinimumBin())));
	if (title.find("_mt")!=string::npos) h_mc1b_b->SetMinimum(TMath::Max(0.00002,0.25*h_data_b->GetBinContent(h_data_b->GetMinimumBin())));
	if (title.find("_pt")!=string::npos) h_mc1b_b->SetMinimum(TMath::Max(0.000004,0.25*h_data_b->GetBinContent(h_data_b->GetMinimumBin())));
	if (title.find("_mass")!=string::npos) h_mc1b_b->SetMinimum(TMath::Max(0.0002,0.25*h_data_b->GetBinContent(h_data_b->GetMinimumBin())));

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

	h_data_b->Draw("E1PX0SAME");

	h_mc1b->SetLineColor(kRed);
	h_mc1b->SetLineWidth(2);
	h_mc1b->SetMarkerColor(kRed);
	h_mc1b->SetFillColor(kRed);
	if (drawInclusive) h_mc1b->Draw("E5SAME");
	TH1F* tmp3 = (TH1F*)h_mc1b->Clone();
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

	h_data->SetMarkerColor(kBlack);
	h_data->SetLineColor(kBlack);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize (0.7);
	if (drawInclusive) h_data->Draw("E1PX0SAME");

	if (ilepton==1) {
	  if (drawInclusive) leg->AddEntry(h_data,"W(#rightarrow e#nu)+1b DATA","p");
	  leg->AddEntry(h_data_b,"W(#rightarrow e#nu)+2b DATA","p");
	  //leg->AddEntry(h_mc1,"W(#rightarrow e#nu) MC","l");
	  leg->AddEntry(h_mcg,"W(#rightarrow e#nu) MadGraph","l");
	}
	if (ilepton==2){
	  if (drawInclusive) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+1b DATA","p");
	  leg->AddEntry(h_data_b,"W(#rightarrow #mu#nu)+2b DATA","p");
	  //leg->AddEntry(h_mc1,"W(#rightarrow #mu#nu) MC","l");
	  leg->AddEntry(h_mcg,"W(#rightarrow #mu#nu) MadGraph","l");
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

	TH1F *h_M = (TH1F*)h_data_b->Clone();
	h_M->Divide(h_mcg_b);

	h_M->SetTitle("");
	h_M->SetStats(0);
	h_M->GetXaxis()->SetTitleOffset(0.9);
	h_M->GetXaxis()->SetTitleSize(0.1);
	h_M->GetXaxis()->SetLabelFont(42);
	h_M->GetXaxis()->SetLabelSize(0.08);
	h_M->GetXaxis()->SetTitleFont(42);
	h_M->GetXaxis()->SetTickLength(0.1);
	h_M->GetYaxis()->SetTitle("Data / Theory");
	h_M->GetYaxis()->SetNdivisions(505);
	h_M->GetYaxis()->SetTitleSize(0.09);
	h_M->GetYaxis()->SetLabelSize(0.08);
	h_M->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_M->GetYaxis()->SetTitleOffset(0.4);
	h_M->GetYaxis()->SetTickLength(0.02);

	h_M->SetMarkerStyle(24);
	h_M->Draw("E1PX0");

	TH1F *h_M2= (TH1F*)h_data->Clone();
	h_M2->Divide(h_mcg);

	TGraphErrors *g_M2 = new TGraphErrors(h_M2);

	float dx = 0.1*(g_M2->GetXaxis()->GetXmax()-g_M2->GetXaxis()->GetXmin())/g_M2->GetN();
	for (int i=0; i<g_M2->GetN(); i++) {
	  g_M2->SetPoint(i, g_M2->GetX()[i]-dx, g_M2->GetY()[i]);
	  g_M2->SetPointError(i, 0, g_M2->GetEY()[i]);
	}

	g_M2->SetMarkerStyle(20);
	if (drawInclusive) g_M2->Draw("EP0SAME");

	TLine *OLine2 = new TLine(h_M->GetXaxis()->GetXmin(),1.,h_M->GetXaxis()->GetXmax(),1.);
	OLine2->SetLineColor(kGreen+2);
	OLine2->SetLineWidth(2);
	OLine2->Draw();

	c1->cd();

	if (title_b=="w_first_jet_pt_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb]");
	  h_M->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	} else if (title_b=="w_first_jet_eta_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta [pb]");
	  h_M->GetXaxis()->SetTitle("leading jet #eta");
	} else if (title_b=="w_first_bjet_pt") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_M->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	} else if (title_b=="w_first_bjet_eta") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_M->GetXaxis()->SetTitle("leading b-jet #eta");
	} else if (title_b=="w_Ht_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dH_{T} [pb]");
	  h_M->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	} else if (title_b=="w_delta_phi_ee_b" || title_b=="w_delta_phi_mm_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bZ} [pb]");
	  h_M->GetXaxis()->SetTitle("#Delta#phi(bZ) [rad]");
	}

	if (plot) {
	  if (unfold==0) {
	    ofstream out;
	    if (ilepton==1) {
	      gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/xsecs/").c_str(), kTRUE);
	      c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + title_b + "_xsecs.pdf").c_str());
	      out.open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + title_b + "_xsecs.dat").c_str());
	      TFile f((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + title_b + "_xsecs.root").c_str(),"RECREATE");
	      h_data_raw->Write((title+"_raw").c_str());
	      h_data_b_raw->Write((title_b+"_raw").c_str());
	      h_data->Write(title.c_str());
	      h_data_b->Write(title_b.c_str());
	      f.Close();
	    }
	    if (ilepton==2) {
	      gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/xsecs/").c_str(), kTRUE);
	      c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/xsecs/" + title_b + "_xsecs.pdf").c_str());
	      out.open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + title_b + "_xsecs.dat").c_str());
	      TFile f((path + "/muons/" + version + "/" + subdir + "/xsecs/" + title_b + "_xsecs.root").c_str(),"RECREATE");
	      h_data_raw->Write((title+"_raw").c_str());
	      h_data_b_raw->Write((title_b+"_raw").c_str());
	      h_data->Write(title.c_str());
	      h_data_b->Write(title_b.c_str());
	      f.Close();
	    }
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_data->Integral(0, h_data->GetNbinsX()+1, "width") << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_data_b->Integral(0, h_data_b->GetNbinsX()+1, "width") << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_mc1b->Integral(0, h_mc1b->GetNbinsX()+1, "width") << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_mc1b_b->Integral(0, h_mc1b_b->GetNbinsX()+1, "width") << endl;
	    out.close();
	  }
	  if (unfold==1) {
	    ofstream out;
	    if (ilepton==1) {
	      gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/").c_str(), kTRUE);
	      c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      out.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.dat").c_str());
	    }
	    if (ilepton==2) {
	      gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/").c_str(), kTRUE);
	      c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      out.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.dat").c_str());
	    }
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_data->Integral(0, h_data->GetNbinsX()+1, "width") << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_data_b->Integral(0, h_data_b->GetNbinsX()+1, "width") << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_mc1b->Integral(0, h_mc1b->GetNbinsX()+1, "width") << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	    out << h_mc1b_b->Integral(0, h_mc1b_b->GetNbinsX()+1, "width") << endl;
	    out.close();
	  }
	}
}

