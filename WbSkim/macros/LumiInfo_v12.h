string version = "v12";

//////////////////// Integrated luminosity in 1/picobarn

double Lumi2012_ele  = 19767.0; // full 2012 electrons luminosity
double Lumi2012_muon = 19783.0; // full 2012 muons luminosity

// single period lumi, same for electrons and muons dataset in 1/picobarn

double Lumi_2012A_22Jan =  889.362;
double Lumi_2012B_22Jan = 4422.000;
double Lumi_2012C_22Jan = 7137.000;
double Lumi_2012D_22Jan = 7318.000;

//////////////////////// DY MadGraph

double Ngen_dy = 30459503;
double Xsec_dy = 3*1177.3; // NNLO
double eXsec_dy = 3*39.1;

//////////////////////// TTbar (used also for the merged sample)

double Ngen_tt = 6923750;
double Xsec_tt = 240.6; // ATLAS+CMS https://cds.cern.ch/record/1950834
double eXsec_tt = 8.5;

//////////////////////// TTbar_FullLept

double Ngen_tt_fl = 12011428;
double Xsec_tt_fl = 25.8; // NNLO : 245.8*0.324*0.324
double eXsec_tt_fl = 0.039*Xsec_tt_fl;

//////////////////////// TTbar_SemiLept

double Ngen_tt_sl = 24963676;
double Xsec_tt_sl = 107.67; // NNLO : 245.8*0.324*0.676*2
double eXsec_tt_sl = 0.039*Xsec_tt_sl;

//////////////////////// TBar_s

double Ngen_tbar_s = 139974;
double Xsec_tbar_s = 1.76; // NNLO
double eXsec_tbar_s = 0.08;

//////////////////////// TBar_t

double Ngen_tbar_t = 1935072;
double Xsec_tbar_t = 29.74; // NNLO
double eXsec_tbar_t = 1.58;

//////////////////////// TBar_tW

double Ngen_tbar_tw = 493460;
double Xsec_tbar_tw = 11.1; // NNLO
double eXsec_tbar_tw = 0.8;

//////////////////////// T_s

double Ngen_t_s = 259961;
double Xsec_t_s = 3.79; // NNLO
double eXsec_t_s = 0.15;

//////////////////////// T_t

double Ngen_t_t = 3758227;
double Xsec_t_t = 54.87; // NNLO
double eXsec_t_t = 2.10;

//////////////////////// T_tW

double Ngen_t_tw = 497658;
double Xsec_t_tw = 11.1; // NNLO
double eXsec_t_tw = 0.8;

//////////////////////// T+Tbar (used for the merged sample)
double Ngen_t = Ngen_tbar_s+Ngen_tbar_t+Ngen_tbar_tw+Ngen_t_s+Ngen_t_t+Ngen_t_tw;
double Xsec_t = Xsec_tbar_s+Xsec_tbar_t+Xsec_tbar_tw+Xsec_t_s+Xsec_t_t+Xsec_t_tw;
double eXsec_t = TMath::Sqrt(TMath::Power(eXsec_tbar_t+eXsec_t_t,2)+TMath::Power(eXsec_tbar_s+eXsec_t_s,2)+TMath::Power(eXsec_tbar_tw+eXsec_t_tw,2));

//////////////////////// ZZ

double Ngen_zz = 9799908;
double Xsec_zz = 8.2; // NLO
double eXsec_zz = 0.4;

//////////////////////// WZ

double Ngen_wz = 10000283;
double Xsec_wz = 33.6; // NLO
double eXsec_wz = 1.3; 

//////////////////////// QCD

double Ngen_qcd = 7529312;
double Xsec_qcd = 3.64E8; // search the NLO !
double eXsec_qcd = 0.15*Xsec_qcd; 

//////////////////////// WW

double Ngen_ww = 10000431;
double Xsec_ww = 56.0; // NLO
double eXsec_ww = 3.0; 

//////////////////////// Wj

double Ngen_wj = 57709905;
double Xsec_wj = 3*12234.4; // NNLO
double eXsec_wj = 3*418.9; 

//////////////////////// W1j

double Ngen_w1j = 23141598;
double Xsec_w1j = 5400; // LO
double eXsec_w1j = 0; 

//////////////////////// W2j

double Ngen_w2j = 34044921;
double Xsec_w2j = 1750; // LO
double eXsec_w2j = 0; 

//////////////////////// W3j

double Ngen_w3j = 15539503;
double Xsec_w3j = 519; // LO
double eXsec_w3j = 0; 

//////////////////////// W4j

double Ngen_w4j = 13382803;
double Xsec_w4j = 214; // LO
double eXsec_w4j = 0; 

//////////////////////// Wbb

double Ngen_wbb = 20646001;
double Xsec_wbb = 3*46.3; // NLO
double eXsec_wbb = 0;

///////////////////////

