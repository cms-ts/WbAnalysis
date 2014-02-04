
#include "LumiInfo_v02.h"

void hmerge_test() {

TFile f0(("data/" + version + "/Wj.root").c_str());
TFile f1(("data/" + version + "/W1j.root").c_str());
TFile f2(("data/" + version + "/W2j.root").c_str());
TFile f3(("data/" + version + "/W3j.root").c_str());
TFile f4(("data/" + version + "/W4j.root").c_str());
TFile f5(("data/" + version + "/Wj_merge.root").c_str());

f0.cd("demoEle");
TH1F* h0x = (TH1F*)gDirectory->Get("h_nmult0");
TH1F* h0 = (TH1F*)gDirectory->Get("h_nmult1");
f1.cd("demoEle");
TH1F* h1 = (TH1F*)gDirectory->Get("h_nmult1");
f2.cd("demoEle");
TH1F* h2 = (TH1F*)gDirectory->Get("h_nmult1");
f3.cd("demoEle");
TH1F* h3 = (TH1F*)gDirectory->Get("h_nmult1");
f4.cd("demoEle");
TH1F* h4 = (TH1F*)gDirectory->Get("h_nmult1");
f5.cd("demoEle");
TH1F* h5 = (TH1F*)gDirectory->Get("h_nmult1");

h0x->Sumw2();
h0->Sumw2();
h1->Sumw2();
h2->Sumw2();
h3->Sumw2();
h4->Sumw2();
h5->Sumw2();

h0x->Scale(Xsec_wj/Ngen_wj);
h0->Scale(Xsec_wj/Ngen_wj);

h1->Scale((Xsec_w1j/Ngen_w1j)*(Xsec_wj/30400.));
h2->Scale((Xsec_w2j/Ngen_w2j)*(Xsec_wj/30400.));
h3->Scale((Xsec_w3j/Ngen_w3j)*(Xsec_wj/30400.));
h4->Scale((Xsec_w4j/Ngen_w4j)*(Xsec_wj/30400.));

h5->Scale(Xsec_wj/Ngen_wj);

TCanvas* c = new TCanvas( "C", "MyCanvas C", 100, 0, 500, 500 );
c->cd();

h0x->SetFillColor(kCyan);
h0x->SetMarkerStyle(24);
h0x->GetXaxis()->SetRangeUser(0, 4);
h0x->SetMinimum(0);
h0x->SetMaximum(28);

h0x->DrawCopy("E2");

h1->DrawCopy("E1SAME");
h2->DrawCopy("E1SAME");
h3->DrawCopy("E1SAME");
h4->DrawCopy("E1SAME");

h5->SetLineColor(kRed);

h5->DrawCopy("E1SAME");

for (int i=1;i<=5;i++) {
cout << i << " ";
cout << h0x->GetBinContent(i) << "+-" << h0x->GetBinError(i) << " ";
cout << h0->GetBinContent(i) << "+-" << h0->GetBinError(i) << " ";
cout << h1->GetBinContent(i) << "+-" << h1->GetBinError(i) << " ";
cout << h2->GetBinContent(i) << "+-" << h2->GetBinError(i) << " ";
cout << h3->GetBinContent(i) << "+-" << h4->GetBinError(i) << " ";
cout << h4->GetBinContent(i) << "+-" << h4->GetBinError(i) << " ";
cout << h5->GetBinContent(i) << "+-" << h5->GetBinError(i) << " ";
cout << endl;
}

}
