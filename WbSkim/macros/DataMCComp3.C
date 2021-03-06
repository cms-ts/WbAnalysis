#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v12.h"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

void DataMCComp3(int irun=0, string title="", int plot=0, int ilepton=1) {

//int useRecoFile=0; // use the GEN histograms as numerators when computing the efficiency
int useRecoFile=1; // use the REC histograms as numerators when computing the efficiency

//int useWbb=0; // do not use the special Wbb MC sample
int useWbb=1; // use the special Wbb MC sample

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

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "Wj_gen_merge.root").c_str());
	if (useRecoFile) mc1 = TFile::Open((path + "/" + version + "/" + "Wj_merge.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "Wj_gen_merge.root").c_str());

	if (useWbb) {
	  if (title.find("_bb")!=string::npos || title.find("_2b")!=string::npos) {
	    mc1 = TFile::Open((path + "/" + version + "/" + "Wbb_gen.root").c_str());
	    if (useRecoFile) mc1 = TFile::Open((path + "/" + version + "/" + "Wbb.root").c_str());
	    mc2 = TFile::Open((path + "/" + version + "/" + "Wbb_gen.root").c_str());
	  }
	}

/* efficiency:  e_W / e_Wb = e_W / e_W_1 * e_W_b */

int itype = 0; // e_W and e_Wb = e_W_1 * e_W_b
//int itype = 1; // e_W_1
//int itype = 2; // e_W_b

	string title_b = title;

	if (itype==0) title_b = "b"+title.substr(1);
	if (itype==2) title_b = "b"+title.substr(1);

	if (ilepton==1&&itype==0) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2&&itype==0) mc1->cd(("demoMuo"+postfix).c_str());
	if (ilepton==1&&itype==1) mc1->cd("demoEleBtag");
	if (ilepton==2&&itype==1) mc1->cd("demoMuoBtag");
	if (ilepton==1&&itype==2) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2&&itype==2) mc1->cd(("demoMuo"+postfix).c_str());
	TH1F* h_reco = (TH1F*)gDirectory->Get(title_b.c_str());
	if (ilepton==1&&itype==0) mc2->cd("demoEleGen");
	if (ilepton==2&&itype==0) mc2->cd("demoMuoGen");
	if (ilepton==1&&itype==1) mc2->cd("demoEleGen");
	if (ilepton==2&&itype==1) mc2->cd("demoMuoGen");
	if (ilepton==1&&itype==2) mc2->cd("demoEleBtag");
	if (ilepton==2&&itype==2) mc2->cd("demoMuoBtag");
	TH1F* h_gen = (TH1F*)gDirectory->Get(title.c_str());

	h_reco->Sumw2();
	h_gen->Sumw2();

	double N = 1.0;
	double errN = 0.0;
	N = h_reco->IntegralAndError(0,h_reco->GetNbinsX()+1,errN) / h_gen->Integral(0,h_gen->GetNbinsX()+1);
	errN = errN / h_gen->Integral(0,h_gen->GetNbinsX()+1);

//	TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
//	c1->cd();

//	h_gen->Draw("EPX0");
//	h_gen->SetMarkerStyle(20);
//	h_reco->Draw("histsame");

//    	c1->Update();

        TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
	c2->cd();

	h_reco->SetTitle("");
	h_reco->SetStats(0);

        if (title=="w_first_jet_pt_b") {
          h_reco->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
        } else if (title=="w_first_jet_eta_b") {
          h_reco->GetXaxis()->SetTitle("leading b-jet #eta [rad]");
        } else if (title=="w_first_jet_mass_b") {
          h_reco->GetXaxis()->SetTitle("leading b-jet mass [GeV/c^{2}]");
        } else if (title=="w_second_jet_pt_b") {
          h_reco->GetXaxis()->SetTitle("sub-leading b-jet p_{T} [GeV/c]");
        } else if (title=="w_second_jet_eta_b") {
          h_reco->GetXaxis()->SetTitle("sub-leading b-jet #eta [rad]");
        } else if (title=="w_second_jet_mass_b") {
          h_reco->GetXaxis()->SetTitle("sub-leading b-jet mass [GeV/c^{2}]");
        } else if (title=="w_Ht_b") {
          h_reco->GetXaxis()->SetTitle("H_{T} [GeV/c]");
        } else if (title=="w_delta_wenu_b"||title=="w_delta_wmnu_b") {
          h_reco->GetXaxis()->SetTitle("#Delta#phi (lepton b-jet) [rad]");
        } else if (title=="w_deltaR_wenu_b"||title=="w_deltaR_wmnu_b") {
          h_reco->GetXaxis()->SetTitle("#DeltaR (lepton b-jet) [rad]");
        }
        if (title=="w_first_jet_pt_bb") {
          h_reco->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
        } else if (title=="w_first_jet_eta_bb") {
          h_reco->GetXaxis()->SetTitle("leading b-jet #eta [rad]");
        } else if (title=="w_first_jet_mass_bb") {
          h_reco->GetXaxis()->SetTitle("leading b-jet mass [GeV/c^{2}]");
        } else if (title=="w_second_jet_pt_bb") {
          h_reco->GetXaxis()->SetTitle("sub-leading b-jet p_{T} [GeV/c]");
        } else if (title=="w_second_jet_eta_bb") {
          h_reco->GetXaxis()->SetTitle("sub-leading b-jet #eta [rad]");
        } else if (title=="w_second_jet_mass_bb") {
          h_reco->GetXaxis()->SetTitle("sub-leading b-jet mass [GeV/c^{2}]");
        } else if (title=="w_Ht_bb") {
          h_reco->GetXaxis()->SetTitle("H_{T} [GeV/c]");
        } else if (title=="w_delta_wenu_bb"||title=="w_delta_wmnu_bb") {
          h_reco->GetXaxis()->SetTitle("#Delta#phi (lepton b-jet) [rad]");
        } else if (title=="w_deltaR_wenu_bb"||title=="w_deltaR_wmnu_bb") {
          h_reco->GetXaxis()->SetTitle("#DeltaR (lepton b-jet) [rad]");
        }

	h_reco->GetXaxis()->SetTitleOffset(0.95);
	h_reco->GetXaxis()->SetTitleSize(0.04);
	h_reco->GetXaxis()->SetLabelFont(42);
	h_reco->GetXaxis()->SetLabelSize(0.04);
	h_reco->GetXaxis()->SetTitleFont(42);
	if (title_b == title) {
	  h_reco->GetYaxis()->SetTitle("#epsilon = N_{W}^{RECO} / N_{W}^{GEN}");
	} else {
	  h_reco->GetYaxis()->SetTitle("#epsilon = N_{W+b}^{RECO} / N_{W+b}^{GEN}");
	}
	h_reco->GetYaxis()->SetNdivisions(505);
	h_reco->GetYaxis()->SetTitleSize(0.04);
	h_reco->GetYaxis()->SetLabelSize(0.04);
	h_reco->GetYaxis()->SetRangeUser(0.0, 1.5);
	h_reco->GetYaxis()->SetTitleOffset(1.04);
	h_reco->SetMarkerStyle(20);

	h_reco->Divide(h_gen);
	h_reco->Draw("EPX0");

        TLatex *Label = new TLatex();
	Label->SetTextSize(0.0275);
	Label->SetTextFont(42);
	Label->SetLineWidth(2);
	Label->SetNDC();
	char buff[100];
	sprintf(buff, "< #epsilon > = #frac{#int RECO}{#int GEN} = %5.3f #pm %5.3f", N, errN);
	Label->DrawLatex(0.50, 0.68, buff);

	c2->Update();

	if (plot) {
	  ofstream out;
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/efficiency/").c_str(), kTRUE);
	    c2->SaveAs((path + "/electrons/" + version + "/" + subdir + "/efficiency/" + title + "_efficiency.pdf").c_str());
	    TFile f((path + "/electrons/" + version + "/" + subdir + "/efficiency/" + title + "_efficiency.root").c_str(),"RECREATE");
	    h_reco->Write(title.c_str());
	    f.Close();
	    out.open((path + "/electrons/" + version + "/" + subdir + "/efficiency/" + title + "_efficiency.dat").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/efficiency/").c_str(), kTRUE);
	    c2->SaveAs((path + "/muons/" + version + "/" + subdir + "/efficiency/" + title + "_efficiency.pdf").c_str());
	    TFile f((path + "/muons/" + version + "/" + subdir + "/efficiency/" + title + "_efficiency.root").c_str(),"RECREATE");
	    h_reco->Write(title.c_str());
	    f.Close();
	    out.open((path + "/muons/" + version + "/" + subdir + "/efficiency/" + title + "_efficiency.dat").c_str());
	  }
	  out << N << " " << errN << endl;
	  out.close();
	}
}

