#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

TH1F* h_data = 0;
TH1F* h_data_fit = 0;

void fcn(int& npar, double* gin, double& fun, double* par, int iflag) {
  double chisq = 0.0;
  if (iflag) {};
  if (gin) {};
  for (int i=1; i<=h_data->GetNbinsX(); i++) {
    double xn = h_data->GetBinContent(i);
    double xd = TMath::Power(h_data->GetBinError(i),2);
    if (npar>0) {
      xn = xn - par[0]*h_data_fit->GetBinContent(i);
      xd = xd + TMath::Power(par[0]*h_data_fit->GetBinError(i),2);
    }
    if (xd!=0) chisq = chisq + (xn*xn)/xd;
  }
  fun = chisq;
}

void DataMCComp5(int irun=0, string title="", int plot=0, int ilepton=1, int doFit=0) {

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

// skip variations for FWD (3,4) and TOP (5,6) samples
if (ilepton>=3 && ilepton<=6) postfix="";

      double Lumi2012=0;
      
      if (ilepton==1||ilepton==3||ilepton==5) Lumi2012 = Lumi2012_ele;
      if (ilepton==2||ilepton==4||ilepton==6) Lumi2012 = Lumi2012_muon;

      double norm1 = ((Lumi2012 * Xsec_wj) / Ngen_wj);
      double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
      double norm3 = ((Lumi2012 * Xsec_zz) / Ngen_zz);
      double norm4 = ((Lumi2012 * Xsec_wz) / Ngen_wz);
      double norm5 = ((Lumi2012 * Xsec_qcd) / Ngen_qcd);
      double norm6 = ((Lumi2012 * Xsec_ww) / Ngen_ww);
      double norm7 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
      double norm8 = ((Lumi2012 * Xsec_tbar_t) / Ngen_tbar_t);

      double enorm1 = ((Lumi2012 * eXsec_wj) / Ngen_wj);
      double enorm2 = ((Lumi2012 * eXsec_tt) / Ngen_tt);
      double enorm3 = ((Lumi2012 * eXsec_zz) / Ngen_zz);
      double enorm4 = ((Lumi2012 * eXsec_wz) / Ngen_wz);
      double enorm5 = ((Lumi2012 * eXsec_qcd) / Ngen_qcd);
      double enorm6 = ((Lumi2012 * eXsec_ww) / Ngen_ww);
      double enorm7 = ((Lumi2012 * eXsec_dy) / Ngen_dy);
      double enorm8 = ((Lumi2012 * eXsec_tbar_t) / Ngen_tbar_t);

      double norm1_fit = ((Lumi2012 * Xsec_wj) / Ngen_wj);
      double norm2_fit = ((Lumi2012 * Xsec_tt) / Ngen_tt);
      double norm3_fit = ((Lumi2012 * Xsec_zz) / Ngen_zz);
      double norm4_fit = ((Lumi2012 * Xsec_wz) / Ngen_wz);
      double norm5_fit = ((Lumi2012 * Xsec_qcd) / Ngen_qcd);
      double norm6_fit = ((Lumi2012 * Xsec_ww) / Ngen_ww);
      double norm7_fit = ((Lumi2012 * Xsec_dy) / Ngen_dy);
      double norm8_fit = ((Lumi2012 * Xsec_tbar_t) / Ngen_tbar_t);

      double enorm1_fit = ((Lumi2012 * eXsec_wj) / Ngen_wj);
      double enorm2_fit = ((Lumi2012 * eXsec_tt) / Ngen_tt);
      double enorm3_fit = ((Lumi2012 * eXsec_zz) / Ngen_zz);
      double enorm4_fit = ((Lumi2012 * eXsec_wz) / Ngen_wz);
      double enorm5_fit = ((Lumi2012 * eXsec_qcd) / Ngen_qcd);
      double enorm6_fit = ((Lumi2012 * eXsec_ww) / Ngen_ww);
      double enorm7_fit = ((Lumi2012 * eXsec_dy) / Ngen_dy);
      double enorm8_fit = ((Lumi2012 * eXsec_tbar_t) / Ngen_tbar_t);

      if (title.empty()) title = "w_jetmultiplicity";

      if (ilepton==1||ilepton==3||ilepton==5) {
        if (title.find("muon")!=string::npos) return;
	if (title.find("mm")!=string::npos) return;
	if (title.find("wmnu")!=string::npos) return;
      }
      if (ilepton==2||ilepton==4||ilepton==6) {
	if (title.find("ele")!=string::npos) return;
	if (title.find("ee")!=string::npos) return;
	if (title.find("wenu")!=string::npos) return;
      }

      TFile *data=0;
      if (ilepton==1||ilepton==3||ilepton==5) data = TFile::Open((path + "/" + version + "/" + "SingleElectron_2012_merge.root").c_str());
      if (ilepton==2||ilepton==4||ilepton==6) data = TFile::Open((path + "/" + version + "/" + "SingleMu_2012_merge.root").c_str());

      TFile *data_fit=data;

      TFile *mc1 = TFile::Open((path + "/" + version + "/" + "Wj_merge.root").c_str());
      TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar_merge.root").c_str());
      TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
      TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//    TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
      TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
      TFile *mc7 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
      TFile *mc8 = TFile::Open((path + "/" + version + "/" + "T_merge.root").c_str());

      string title_fit = title;

      if (ilepton==1) data->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) data->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) data->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) data->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) data->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) data->cd(("demoMuoTOP"+postfix).c_str());
      h_data = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) data_fit->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) data_fit->cd(("demoMuoQCD"+postfix).c_str());
      h_data_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) mc1->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) mc1->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) mc1->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) mc1->cd(("demoMuoTOP"+postfix).c_str());
      TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) mc1->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) mc1->cd(("demoMuoQCD"+postfix).c_str());
      TH1F* h_mc1_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc2->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc2->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) mc2->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) mc2->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) mc2->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) mc2->cd(("demoMuoTOP"+postfix).c_str());
      TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) mc2->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) mc2->cd(("demoMuoQCD"+postfix).c_str());
      TH1F* h_mc2_fit = (TH1F*)gDirectory->Get(title_fit.c_str());
      
      if (ilepton==1) mc3->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc3->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) mc3->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) mc3->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) mc3->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) mc3->cd(("demoMuoTOP"+postfix).c_str());
      TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) mc3->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) mc3->cd(("demoMuoQCD"+postfix).c_str());
      TH1F* h_mc3_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc4->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc4->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) mc4->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) mc4->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) mc4->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) mc4->cd(("demoMuoTOP"+postfix).c_str());
      TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) mc4->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) mc4->cd(("demoMuoQCD"+postfix).c_str());
      TH1F* h_mc4_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

//    if (ilepton==1) mc5->cd(("demoEle"+postfix).c_str());
//    if (ilepton==2) mc5->cd(("demoMuo"+postfix).c_str());
//    if (ilepton==3) mc5->cd(("demoEleFWD"+postfix).c_str());
//    if (ilepton==4) mc5->cd(("demoMuoFWD"+postfix).c_str());
//    if (ilepton==5) mc5->cd(("demoEleTOP"+postfix).c_str());
//    if (ilepton==6) mc5->cd(("demoMuoTOP"+postfix).c_str());
//    TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
//    if (ilepton==1) mc5->cd(("demoEleQCD"+postfix).c_str());
//    if (ilepton==2) mc5->cd(("demoMuoQCD"+postfix).c_str());
//    TH1F* h_mc5_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc6->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc6->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) mc6->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) mc6->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) mc6->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) mc6->cd(("demoMuoTOP"+postfix).c_str());
      TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) mc6->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) mc6->cd(("demoMuoQCD"+postfix).c_str());
      TH1F* h_mc6_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc7->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc7->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) mc7->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) mc7->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) mc7->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) mc7->cd(("demoMuoTOP"+postfix).c_str());
      TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) mc7->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) mc7->cd(("demoMuoQCD"+postfix).c_str());
      TH1F* h_mc7_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc8->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc8->cd(("demoMuo"+postfix).c_str());
      if (ilepton==3) mc8->cd(("demoEleFWD"+postfix).c_str());
      if (ilepton==4) mc8->cd(("demoMuoFWD"+postfix).c_str());
      if (ilepton==5) mc8->cd(("demoEleTOP"+postfix).c_str());
      if (ilepton==6) mc8->cd(("demoMuoTOP"+postfix).c_str());
      TH1F* h_mc8 = (TH1F*)gDirectory->Get(title.c_str());
      if (ilepton==1||ilepton==3||ilepton==5) mc8->cd(("demoEleQCD"+postfix).c_str());
      if (ilepton==2||ilepton==4||ilepton==6) mc8->cd(("demoMuoQCD"+postfix).c_str());
      TH1F* h_mc8_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      h_data->Sumw2();
      h_data_fit->Sumw2();

      h_mc1->Sumw2();
      h_mc2->Sumw2();
      h_mc3->Sumw2();
      h_mc4->Sumw2();
//    h_mc5->Sumw2();
      h_mc6->Sumw2();
      h_mc7->Sumw2();
      h_mc8->Sumw2();
      
      h_mc1_fit->Sumw2();
      h_mc2_fit->Sumw2();
      h_mc3_fit->Sumw2();
      h_mc4_fit->Sumw2();
//    h_mc5_fit->Sumw2();
      h_mc6_fit->Sumw2();
      h_mc7_fit->Sumw2();
      h_mc8_fit->Sumw2();

      if (irun==10) {
        norm1 = norm1 + enorm1;
        norm2 = norm2 + enorm2;
        norm3 = norm3 + enorm3;
        norm4 = norm4 + enorm4;
//      norm5 = norm5 + enorm5;
        norm6 = norm6 + enorm6;
        norm7 = norm7 + enorm7;
        norm8 = norm8 + enorm8;
      }

      h_mc1->Scale(norm1);
      h_mc2->Scale(norm2);
      h_mc3->Scale(norm3);
      h_mc4->Scale(norm4);
//    h_mc5->Scale(norm5);
      h_mc6->Scale(norm6);
      h_mc7->Scale(norm7);
      h_mc8->Scale(norm8);
      
      if (irun==10) {
        norm1_fit = norm1_fit + enorm1_fit;
        norm2_fit = norm2_fit + enorm2_fit;
        norm3_fit = norm3_fit + enorm3_fit;
        norm4_fit = norm4_fit + enorm4_fit;
//      norm5_fit = norm5_fit + enorm5_fit;
        norm6_fit = norm6_fit + enorm6_fit;
        norm7_fit = norm7_fit + enorm7_fit;
        norm8_fit = norm8_fit + enorm8_fit;
      }

      h_mc1_fit->Scale(norm1_fit);
      h_mc2_fit->Scale(norm2_fit);
      h_mc3_fit->Scale(norm3_fit);
      h_mc4_fit->Scale(norm4_fit);
//    h_mc5_fit->Scale(norm5_fit);
      h_mc6_fit->Scale(norm6_fit);
      h_mc7_fit->Scale(norm7_fit);
      h_mc8_fit->Scale(norm8_fit);

      if (irun==13) {
        for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
          h_mc1->SetBinError(i, 1.1*h_mc1->GetBinError(i));
          h_mc2->SetBinError(i, 1.1*h_mc2->GetBinError(i));
          h_mc3->SetBinError(i, 1.1*h_mc3->GetBinError(i));
          h_mc4->SetBinError(i, 1.1*h_mc4->GetBinError(i));
//        h_mc5->SetBinError(i, 1.1*h_mc5->GetBinError(i));
          h_mc6->SetBinError(i, 1.1*h_mc6->GetBinError(i));
          h_mc7->SetBinError(i, 1.1*h_mc7->GetBinError(i));
          h_mc8->SetBinError(i, 1.1*h_mc8->GetBinError(i));
        }
      }

      h_data->Add(h_mc8, -1.);
      h_data->Add(h_mc7, -1.);
      h_data->Add(h_mc6, -1.);
//    h_data->Add(h_mc5, -1.);
      h_data->Add(h_mc4, -1.);
      h_data->Add(h_mc3, -1.);
      h_data->Add(h_mc2, -1.);
      h_data->Add(h_mc1, -1.);

      h_data_fit->Add(h_mc8_fit, -1.);
      h_data_fit->Add(h_mc7_fit, -1.);
      h_data_fit->Add(h_mc6_fit, -1.);
//    h_data_fit->Add(h_mc5_fit, -1.);
      h_data_fit->Add(h_mc4_fit, -1.);
      h_data_fit->Add(h_mc3_fit, -1.);
      h_data_fit->Add(h_mc2_fit, -1.);
      h_data_fit->Add(h_mc1_fit, -1.);

      for (int i=0; i<=h_data->GetNbinsX()+1; i++) {
        if (h_data->GetBinContent(i) < 0) {
          h_data->SetBinContent(i, 0);
          h_data->SetBinError(i, 0);
        }
        if (h_data_fit->GetBinContent(i) < 0) {
          h_data_fit->SetBinContent(i, 0);
          h_data_fit->SetBinError(i, 0);
        }
      }

      TH1F* h_data_fit_raw = (TH1F*)h_data_fit->Clone();

      TVirtualFitter::SetDefaultFitter("Minuit2");
      TVirtualFitter* fitter=0;
      if (doFit==1) {
        for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	  bool skip = false;
          if (title.find("w_mt")!=string::npos) {
            if (h_data_fit->GetXaxis()->GetBinCenter(i)>20) {
	      skip = true;
            }
	  }
	  if (skip) {
  	    h_data->SetBinContent(i, 0);
	    h_data->SetBinError(i, 0);
	    h_data_fit->SetBinContent(i, 0);
	    h_data_fit->SetBinError(i, 0);
	  }
        }
	fitter = TVirtualFitter::Fitter(0, 1);
	fitter->SetFCN(fcn);
	double arglist[1] = {-1.0};
	fitter->ExecuteCommand("SET PRINT", arglist, 1);
	fitter->SetParameter(0, "c(QCD)", 1.00, 0.01, 0.00, 10000.00);
	fitter->ExecuteCommand("MIGRAD", arglist, 0);
	h_data_fit->Scale(fitter->GetParameter(0));
      }

      TCanvas* c1=0;
      if (doFit) {
        c1 = new TCanvas("c", "c", 800, 600);
        c1->cd();

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();

        h_data->SetMinimum(0);
        h_data->SetMaximum(-1111);

        h_data->Draw("EP");
        h_data_fit->Draw("HISTSAME");
        h_data->SetMarkerColor(kBlack);
        h_data_fit->SetMarkerColor(kRed);
        h_data_fit->SetLineColor(kRed);
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize (1.0);
        h_data->SetTitle("");
        //h_data->SetStats(0);

	pad1->Update();
	c1->Update();

	c1->cd();

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();

	TH1F* h_ratio = (TH1F*)h_data->Clone("h_ratio");
	h_ratio->Divide(h_data_fit);

	h_ratio->SetTitle("");
        h_ratio->SetStats(0);

	h_ratio->GetXaxis()->SetTitleOffset(0.9);
	h_ratio->GetXaxis()->SetTitleSize(0.1);
        h_ratio->GetXaxis()->SetLabelFont(42);
        h_ratio->GetXaxis()->SetLabelSize(0.08);
        h_ratio->GetXaxis()->SetTitleFont(42);
        h_ratio->GetYaxis()->SetTitle("ratio");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.09);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);

        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerColor(kBlack);

        h_ratio->Draw("EPX");

        TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kGreen);
        OLine->SetLineWidth(1);
        OLine->Draw();

        c1->cd();

	TLegend *leg;
	leg = new TLegend(0.5, 0.8, 0.65, 0.9);

	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data,"W(#rightarrow e#nu)+jets","p");
	if (ilepton==2) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+jets","p");
	if (ilepton==3) leg->AddEntry(h_data,"W(#rightarrow e#nu)+jets FWD","p");
	if (ilepton==4) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+jets FWD","p");
	if (ilepton==5) leg->AddEntry(h_data,"W(#rightarrow e#nu)+jets TOP","p");
	if (ilepton==6) leg->AddEntry(h_data,"W(#rightarrow #mu#nu)+jets TOP","p");

	if (ilepton==1||ilepton==3||ilepton==5) leg->AddEntry(h_data_fit,"W(#rightarrow e#nu)+jets [QCD]","l");
	if (ilepton==2||ilepton==4||ilepton==6) leg->AddEntry(h_data_fit,"W(#rightarrow #mu#nu)+jets [QCD]","l");

	leg->Draw();

        TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
        latexLabel->Draw("same");

        TLatex *fitLabel = new TLatex();
        fitLabel->SetTextSize(0.0275);
        fitLabel->SetTextFont(42);
        fitLabel->SetLineWidth(2);
        fitLabel->SetNDC();
        char buff[100];
        sprintf(buff, "c_{QCD} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
        fitLabel->DrawLatex(0.60, 0.68, buff);
        if (ilepton==1) sprintf(buff, "I_{e} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        if (ilepton==2) sprintf(buff, "I_{#mu} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        if (ilepton==3) sprintf(buff, "I_{e,FWD} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        if (ilepton==4) sprintf(buff, "I_{#mu,FWD} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        if (ilepton==5) sprintf(buff, "I_{e,TOP} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        if (ilepton==6) sprintf(buff, "I_{#mu,TOP} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        fitLabel->DrawLatex(0.60, 0.63, buff);
        if (ilepton==1||ilepton==3||ilepton==5) sprintf(buff, "I_{e,QCD} = %5.1f #pm %5.1f", h_data_fit->Integral(0,h_data_fit->GetNbinsX()+1), h_data_fit->Integral(0,h_data_fit->GetNbinsX()+1)*fitter->GetParError(0)/fitter->GetParameter(0));
        if (ilepton==2||ilepton==4||ilepton==6) sprintf(buff, "I_{#mu,QCD} = %5.1f #pm %5.1f", h_data_fit->Integral(0,h_data_fit->GetNbinsX()+1), h_data_fit->Integral(0,h_data_fit->GetNbinsX()+1)*fitter->GetParError(0)/fitter->GetParameter(0));
        fitLabel->DrawLatex(0.60, 0.58, buff);
      }

      if (plot) {
        if (doFit) title = title + "_doFit";
	ofstream out;
        if (ilepton==1) {
          gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/electrons/" + version + "/" + subdir + "/qcd_sub/" + title + ".dat").c_str());
        }
        if (ilepton==2) {
          gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/qcd_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/muons/" + version + "/" + subdir + "/qcd_sub/" + title + ".dat").c_str());
        }
        if (ilepton==3) {
          gSystem->mkdir((path + "/electronsFWD/" + version + "/" + subdir + "/qcd_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/electronsFWD/" + version + "/" + subdir + "/qcd_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/electronsFWD/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/electronsFWD/" + version + "/" + subdir + "/qcd_sub/" + title + ".dat").c_str());
        }
        if (ilepton==4) {
          gSystem->mkdir((path + "/muonsFWD/" + version + "/" + subdir + "/qcd_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/muonsFWD/" + version + "/" + subdir + "/qcd_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/muonsFWD/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/muonsFWD/" + version + "/" + subdir + "/qcd_sub/" + title + ".dat").c_str());
        }
        if (ilepton==5) {
          gSystem->mkdir((path + "/electronsTOP/" + version + "/" + subdir + "/qcd_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/electronsTOP/" + version + "/" + subdir + "/qcd_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/electronsTOP/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/electronsTOP/" + version + "/" + subdir + "/qcd_sub/" + title + ".dat").c_str());
        }
        if (ilepton==6) {
          gSystem->mkdir((path + "/muonsTOP/" + version + "/" + subdir + "/qcd_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/muonsTOP/" + version + "/" + subdir + "/qcd_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/muonsTOP/" + version + "/" + subdir + "/qcd_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/muonsTOP/" + version + "/" + subdir + "/qcd_sub/" + title + ".dat").c_str());
        }
	if (doFit==1) {
	  out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	  out.close();
	}
      }
}
