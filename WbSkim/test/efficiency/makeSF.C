#include <iomanip>

void makeSF() {

  string pathFileData ="efficiency-data-WP70toHLT.root";
  string pathFileMC   ="efficiency-mc-WP70toHLT.root";
  
  TFile *histodata = TFile::Open (pathFileData.c_str ());
  histodata->cd ("WP70ToHLT/HLT_Ele27_WP80/fit_eff_plots");
  TCanvas *dataEff;
  gDirectory->GetObject ("probe_patEle_eta_probe_patEle_et_PLOT", dataEff);

  TFile *histomc = TFile::Open (pathFileMC.c_str ());
  histomc->cd ("WP70ToHLT/HLT_Ele27_WP80/fit_eff_plots");
  TCanvas *mcEff;
  gDirectory->GetObject ("probe_patEle_eta_probe_patEle_et_PLOT", mcEff);
  
  TH2D *hdata;
  hdata = (TH2D*) dataEff->GetPrimitive("probe_patEle_eta_probe_patEle_et_PLOT");
  hdata->Sumw2();

  TH2D *hmc;
  hmc = (TH2D*) mcEff->GetPrimitive("probe_patEle_eta_probe_patEle_et_PLOT");
  hmc->Sumw2();

  TH2D *hsf;
  hsf = (TH2D *) hdata->Clone();
  hsf->Sumw2();
  hsf->Divide(hmc);

  ofstream textfileData;
  textfileData.open("ele_HLT_data.txt");

  for (int i = 1; i<=hdata->GetNbinsX(); i++) {
    for (int j = 1; j<=hdata->GetNbinsY(); j++) {
      textfileData << std::fixed
	         << std::setprecision(1)
                 << std::setw(6) << hdata->GetYaxis()->GetBinLowEdge(j)
		 << std::setw(8) << hdata->GetYaxis()->GetBinUpEdge (j)
		 << std::setprecision(4)
		 << std::setw(7) << hdata->GetXaxis()->GetBinLowEdge(i)
		 << std::setw(9) << hdata->GetXaxis()->GetBinUpEdge (i)
                 << std::setprecision(4)
		 << std::setw(9) << hdata->GetBinContent(i,j)
                 << std::setw(9) << hdata->GetBinError(i,j)
		 << std::setw(9) << hdata->GetBinError(i,j) << endl;
    }
  }
  textfileData.close();


  ofstream textfileMC;
  textfileMC.open("ele_HLT_mc.txt");

  for (int i = 1; i<=hmc->GetNbinsX(); i++) {
    for (int j = 1; j<=hmc->GetNbinsY(); j++) {
      textfileMC << std::fixed
	         << std::setprecision(1)
                 << std::setw(6) << hmc->GetYaxis()->GetBinLowEdge(j)
		 << std::setw(8) << hmc->GetYaxis()->GetBinUpEdge (j)
		 << std::setprecision(4)
		 << std::setw(7) << hmc->GetXaxis()->GetBinLowEdge(i)
		 << std::setw(9) << hmc->GetXaxis()->GetBinUpEdge (i)
                 << std::setprecision(4)
		 << std::setw(9) << hmc->GetBinContent(i,j)
                 << std::setw(9) << hmc->GetBinError(i,j)
		 << std::setw(9) << hmc->GetBinError(i,j) << endl;
    }
  }
  textfileMC.close();


  ofstream textfileSF;
  textfileSF.open("ele_HLT_sf.txt");
  for (int i = 1; i<=hsf->GetNbinsX(); i++) {
    for (int j = 1; j<=hsf->GetNbinsY(); j++) {
      textfileSF << std::fixed
	         << std::setprecision(1)
                 << std::setw(6) << hsf->GetYaxis()->GetBinLowEdge(j)
		 << std::setw(8) << hsf->GetYaxis()->GetBinUpEdge (j)
		 << std::setprecision(4)
		 << std::setw(7) << hsf->GetXaxis()->GetBinLowEdge(i)
		 << std::setw(9) << hsf->GetXaxis()->GetBinUpEdge (i)
                 << std::setprecision(4)
		 << std::setw(9) << hsf->GetBinContent(i,j)
                 << std::setw(9) << hsf->GetBinError(i,j)
		 << std::setw(9) << hsf->GetBinError(i,j) << endl;
    }
  }
  textfileSF.close();
}
