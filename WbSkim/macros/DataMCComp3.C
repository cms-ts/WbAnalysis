#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v03.h"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

void DataMCComp3(int irun=0, string title="", int plot=0, int ilepton=1) {

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
if (irun==66) {            // irun==66 => unfolding with data weight
  subdir="66";
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
        }
        if (ilepton==2) {
          if (title.find("ele")!=string::npos) return;
          if (title.find("ee")!=string::npos) return;
        }

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "Wj_gen_merge.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "Wj_gen_merge.root").c_str());

/* efficiency:  e_Z / e_Zb = e_Z / e_Z_1 * e_Z_b */

int itype = 0; // e_Z and e_Zb = e_Z_1 * e_Z_b
//int itype = 1; // e_Z_1
//int itype = 2; // e_Z_b

	string title_b = title;

	if (title.find("_b")!=string::npos) {
	  if (itype==0) title_b = "b"+title.substr(1);
	  if (itype==2) title_b = "b"+title.substr(1);
	}

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
	double N = h_reco->GetEffectiveEntries() / h_gen->GetEffectiveEntries();
	double errN = TMath::Sqrt(h_reco->GetEffectiveEntries()) / h_gen->GetEffectiveEntries();

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

	if (title=="w_first_jet_pt") {
	  h_reco->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta") {
	  h_reco->GetXaxis()->SetTitle("leading jet #eta");
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm") {
	  h_reco->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee"||title=="w_delta_phi_mm") {
	  h_reco->GetXaxis()->SetTitle("#Delta #phi(Zj) [rad]");
	} else if (title=="w_Ht") {
          h_reco->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	}
	if (title=="w_first_bjet_pt") {
	  h_reco->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  h_reco->GetXaxis()->SetRangeUser(0, 200);
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	} else if (title=="w_first_bjet_eta") {
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	  h_reco->GetXaxis()->SetTitle("leading b-jet #eta");
	} else if (title=="w_pt_Z_ee_b"||title=="w_pt_Z_mm_b") {
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	  h_reco->GetXaxis()->SetRangeUser(0, 200);
	  h_reco->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee_b"||title=="w_delta_phi_mm_b") {
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	  h_reco->GetXaxis()->SetTitle("#Delta #phi(Zb) [rad]");
	} else if (title=="w_Ht_b") {
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	  h_reco->GetXaxis()->SetRangeUser(0, 200);
          h_reco->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	}

	h_reco->GetXaxis()->SetTitleOffset(0.95);
	h_reco->GetXaxis()->SetTitleSize(0.04);
	h_reco->GetXaxis()->SetLabelFont(42);
	h_reco->GetXaxis()->SetLabelSize(0.04);
	h_reco->GetXaxis()->SetTitleFont(42);
	if (title_b == title) {
	  h_reco->GetYaxis()->SetTitle("#epsilon = N_{Z}^{RECO} / N_{Z}^{GEN}");
	} else {
	  h_reco->GetYaxis()->SetTitle("#epsilon = N_{Z+b}^{RECO} / N_{Z+b}^{GEN}");
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
