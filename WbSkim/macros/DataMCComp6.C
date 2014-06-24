#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/schizzi/Wbb2012/test/data/";

void DataMCComp6(int irun=0, string title="", int plot=0) {

int unfold=0; // use pre-unfolding distributions
//int unfold=1; // use unfolded distributions

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

        string title_e = title;
        string title_e_b = title;
        string title_m = title;
        string title_m_b = title;

        if (title=="w_mt_b") {
	  title_e = "w_mt_wenu_b";
	  title_e_b = "w_mt_wenu_bb";
	  title_m = "w_mt_wmnu_b";
	  title_m_b = "w_mt_wmnu_bb";
	}

        if (title.find("_jet_")!=string::npos) {
	  title_e_b = title_e_b + "b";
	  title_m_b = title_m_b + "b";
	}

        if (title.find("_bjet_")!=string::npos) {
          title_e.erase(title_e.find("_bjet_")+1, 1);
          title_m.erase(title_m.find("_bjet_")+1, 1);
        }

        TFile* f_e;
        TFile* f_e_b;
        if (unfold) {
           f_e = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title_e + "_unfolding.root").c_str());
           f_e_b = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title_e_b + "_unfolding.root").c_str());
        } else {
          f_e = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + title_e_b + "_xsecs.root").c_str());
          f_e_b = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + title_e_b + "_xsecs.root").c_str());
        }
        TH1F* h_data_e = (TH1F*)f_e->Get(title_e.c_str())->Clone();
        TH1F* h_data_e_b = (TH1F*)f_e_b->Get(title_e_b.c_str())->Clone();
        h_data_e->SetDirectory(0);
        h_data_e_b->SetDirectory(0);
        f_e->Close();
        f_e_b->Close();

        TFile* f_m;
        TFile* f_m_b;
        if (unfold) {
          f_m = TFile::Open((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title_m + "_unfolding.root").c_str());
          f_m_b= TFile::Open((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title_m_b + "_unfolding.root").c_str());
        } else {
          f_m= TFile::Open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + title_m_b + "_xsecs.root").c_str());
          f_m_b= TFile::Open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + title_m_b + "_xsecs.root").c_str());
        }
        TH1F* h_data_m = (TH1F*)f_m->Get(title_m.c_str())->Clone();
        TH1F* h_data_m_b = (TH1F*)f_m_b->Get(title_m_b.c_str())->Clone();
        h_data_m->SetDirectory(0);
        h_data_m_b->SetDirectory(0);
        f_m->Close();
        f_m_b->Close();

        h_data_e->SetStats(0);
        h_data_e_b->SetStats(0);
        h_data_m->SetStats(0);
        h_data_m_b->SetStats(0);

	if (unfold) {
	  h_data_e->Scale(1./Lumi2012_ele, "width");
	  h_data_e_b->Scale(1./Lumi2012_ele, "width");
	  h_data_m->Scale(1./Lumi2012_muon, "width");
	  h_data_m_b->Scale(1./Lumi2012_muon, "width");
	}

        TCanvas* c1 = new TCanvas("c", "c", 800, 600);
        c1->cd();

        TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
        pad1->SetBottomMargin(0.001);
        pad1->Draw();
        pad1->cd();
        pad1->SetLogy();

        h_data_e->SetTitle("");
        h_data_e_b->SetTitle("");

        h_data_e->GetYaxis()->SetTitle("#sigma [pb]");

        h_data_e->SetMarkerStyle(20);
        h_data_e->SetMarkerSize(0.7);
        h_data_e->SetMarkerColor(kRed);
        h_data_e_b->SetMarkerStyle(24);
        h_data_e_b->SetMarkerSize(0.7);
        h_data_e_b->SetMarkerColor(kRed);
        h_data_m->SetMarkerStyle(20);
        h_data_m->SetMarkerSize(0.7);
        h_data_m->SetMarkerColor(kBlue);
        h_data_m_b->SetMarkerStyle(24);
        h_data_m_b->SetMarkerSize(0.7);
        h_data_m_b->SetMarkerColor(kBlue);

	if (drawInclusive) h_data_e_b->SetMaximum(4*h_data_e->GetMaximum());
        h_data_e_b->SetMinimum(TMath::Max(0.002,0.25*h_data_e_b->GetBinContent(h_data_e_b->GetMinimumBin())));
	if (title.find("_mt")!=string::npos) h_data_e_b->SetMinimum(TMath::Max(0.00002,0.25*h_data_e_b->GetBinContent(h_data_e_b->GetMinimumBin())));
	if (title.find("_pt")!=string::npos) h_data_e_b->SetMinimum(TMath::Max(0.000002,0.25*h_data_e_b->GetBinContent(h_data_e_b->GetMinimumBin())));
	if (title.find("_mass")!=string::npos) h_data_e_b->SetMinimum(TMath::Max(0.00002,0.25*h_data_e_b->GetBinContent(h_data_e_b->GetMinimumBin())));

        h_data_e_b->Draw("EPX");
        if (drawInclusive) h_data_e->Draw("EPXSAME");
        h_data_m_b->Draw("EPXSAME");
        if (drawInclusive) h_data_m->Draw("EPXSAME");

        TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetEntrySeparation(0.01);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);

        leg->AddEntry(h_data_e_b,"W(#rightarrow e#nu)+2b DATA","p");
        if (drawInclusive) leg->AddEntry(h_data_e,"W(#rightarrow e#nu)+1b DATA","p");
        leg->AddEntry(h_data_m_b,"W(#rightarrow #mu#nu)+2b DATA","p");
        if (drawInclusive) leg->AddEntry(h_data_m,"W(#rightarrow #mu#nu)+1b DATA","p");

        leg->Draw();

        pad1->Update();
        c1->Update();

        c1->cd();

        TLatex *latexLabel = CMSPrel((Lumi2012_ele+Lumi2012_muon)/2./1000.,"",0.15,0.94);
        latexLabel->Draw("same");

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();

        TH1F* h_ratio = (TH1F*)h_data_e->Clone("h_ratio");
        TH1F* h_ratio_b = (TH1F*)h_data_e_b->Clone("h_ratio_b");
        h_ratio->Divide(h_data_m);
        h_ratio_b->Divide(h_data_m_b);

        h_ratio_b->SetTitle("");

        h_ratio_b->GetXaxis()->SetTitleOffset(0.9);
        h_ratio_b->GetXaxis()->SetTitleSize(0.1);
        h_ratio_b->GetXaxis()->SetLabelFont(42);
        h_ratio_b->GetXaxis()->SetLabelSize(0.08);
        h_ratio_b->GetXaxis()->SetTitleFont(42);
        h_ratio_b->GetYaxis()->SetTitle("Electrons/Muons");
        h_ratio_b->GetYaxis()->SetNdivisions(505);
        h_ratio_b->GetYaxis()->SetTitleSize(0.09);
        h_ratio_b->GetYaxis()->SetLabelSize(0.08);
        h_ratio_b->GetYaxis()->SetRangeUser(0.5, 1.5);
        h_ratio_b->GetYaxis()->SetTitleOffset(0.4);

        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerColor(kBlack);
        h_ratio_b->SetMarkerStyle(24);
        h_ratio_b->SetMarkerColor(kBlack);

        h_ratio_b->Draw("EPX");
        if (drawInclusive) h_ratio->Draw("EPXSAME");

        TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kGreen);
        OLine->SetLineWidth(1);
        OLine->Draw();

        c1->cd();

	if (plot) {
          if (unfold) {
            gSystem->mkdir((path + "/combined/" + version + "/" + subdir + "/xsecs_unfolding/").c_str(), kTRUE);
            c1->SaveAs((path + "/combined/" + version + "/" + subdir + "/xsecs_unfolding/" + title + "b" + "_xsecs_unfolding.pdf").c_str());
          } else {
            gSystem->mkdir((path + "/combined/" + version + "/" + subdir + "/xsecs/").c_str(), kTRUE);
            c1->SaveAs((path + "/combined/" + version + "/" + subdir + "/xsecs/" + title + "b" + "_xsecs.pdf").c_str());
          }
	}

}
