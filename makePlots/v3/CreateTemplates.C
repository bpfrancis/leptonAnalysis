#include "CreateTemplates.h"

using namespace std;

void CreateTemplates(TString input, TString variable, int controlRegion, TString channel, Int_t nBins, Float_t xlo, Float_t xhi, Float_t metCut = -1.0) {

  TFile * in = new TFile(input, "READ");

  TString sigName, qcdName;
  if(controlRegion == kSR1 || controlRegion == kSR2 || controlRegion == kCR0 || controlRegion == kAny) {
    sigName = channel+"_signalTree";
    qcdName = (channel.Contains("ele")) ? channel+"_eQCDTree" : channel+"_muQCDTree";
  }
  else {
    sigName = channel+"_fakeTree";
    qcdName = (channel.Contains("ele")) ? channel+"_eQCDfakeTree" : channel+"_muQCDfakeTree";
  }

  TTree * ggTree = (TTree*)in->Get(sigName);
  TTree * qcdTree = (TTree*)in->Get(qcdName);

  TemplateMaker * tMaker = new TemplateMaker(variable, channel, controlRegion, nBins, xlo, xhi, metCut);
  tMaker->LoadLeptonSFs("/eos/uscms/store/user/bfrancis/data/lepton_SF_8TeV_53x_baseline.root");
  tMaker->LoadPhotonSFs("/eos/uscms/store/user/bfrancis/data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root");

  bool loadSuccess = true;

  Double_t ttbar_hadronic_xsec = 245.8 * 0.457;
  Double_t ttbar_semiLep_xsec  = 245.8 * 0.438;
  Double_t ttbar_fullLep_xsec  = 245.8 * 0.105;

  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic", 
					  ttbar_hadronic_xsec, ttbar_hadronic_xsec * 0.025, ttbar_hadronic_xsec * 0.034, ttbar_hadronic_xsec * 0.026, ttbar_hadronic_xsec * 0.026,
					  false, true, true);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep", 
					  ttbar_semiLep_xsec, ttbar_semiLep_xsec * 0.025, ttbar_semiLep_xsec * 0.034, ttbar_semiLep_xsec * 0.026, ttbar_semiLep_xsec * 0.026,
					  false, true, true);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep", 
					  ttbar_fullLep_xsec, ttbar_fullLep_xsec * 0.025, ttbar_fullLep_xsec * 0.034, ttbar_fullLep_xsec * 0.026, ttbar_fullLep_xsec * 0.026,
					  false, true, true);

  double fix_wjets_xsec = 3. * 12234.4 / 37509.;

  double xsec_w1 = 6662. * fix_wjets_xsec;
  double xsec_w2 = 2159.9 * fix_wjets_xsec;
  double xsec_w3 = 640. * fix_wjets_xsec;
  double xsec_w4 = 264. * fix_wjets_xsec;

  double scaleUp_w1 = 79.0 * 6662. / 37509.;
  double scaleUp_w2 = 79.0 * 2159.9 / 37509.;
  double scaleUp_w3 = 79.0 * 640. / 37509.;
  double scaleUp_w4 = 79.0 * 264. / 37509.;

  double scaleDown_w1 = 39.7 * 6662. / 37509.;
  double scaleDown_w2 = 39.7 * 2159.9 / 37509.;
  double scaleDown_w3 = 39.7 * 640. / 37509.;
  double scaleDown_w4 = 39.7 * 264. / 37509.;

  double pdf_w1 = 414.7 * 6662. / 37509.;
  double pdf_w2 = 414.7 * 2159.9 / 37509.;
  double pdf_w3 = 414.7 * 640. / 37509.;
  double pdf_w4 = 414.7 * 264. / 37509.;

  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_W1JetsToLNu.root", "W1JetsToLNu", 
					  xsec_w1,
					  scaleUp_w1, scaleDown_w1,
					  pdf_w1, pdf_w1);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_W2JetsToLNu.root", "W2JetsToLNu", 
					  xsec_w2,
					  scaleUp_w2, scaleDown_w2,
					  pdf_w2, pdf_w2);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_W3JetsToLNu.root", "W3JetsToLNu", 
					  xsec_w3,
					  scaleUp_w3, scaleDown_w3,
					  pdf_w3, pdf_w3);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_W4JetsToLNu.root", "W4JetsToLNu", 
					  xsec_w4,
					  scaleUp_w4, scaleDown_w4,
					  pdf_w4, pdf_w4);

  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_dyJetsToLL.root", "dyJetsToLL", 
					  1177.3 * 3, 
					  5.9, 3.6, 
					  38.8, 38.8);
  
  double fix_zjets_xsec = 3. * 1177.3 / 3503.71;
  
  double xsec_dy1 = 666.7 * fix_zjets_xsec;
  double xsec_dy2 = 215.1 * fix_zjets_xsec;
  double xsec_dy3 = 66.07 * fix_zjets_xsec;
  double xsec_dy4 = 27.38 * fix_zjets_xsec;
  
  double scaleUp_dy1 = 5.9 * 666.7 / 3503.71;
  double scaleUp_dy2 = 5.9 * 215.1 / 3503.71;
  double scaleUp_dy3 = 5.9 * 66.07 / 3503.71;
  double scaleUp_dy4 = 5.9 * 27.38 / 3503.71;
  
  double scaleDown_dy1 = 3.6 * 666.7 / 3503.71;
  double scaleDown_dy2 = 3.6 * 215.1 / 3503.71;
  double scaleDown_dy3 = 3.6 * 66.07 / 3503.71;
  double scaleDown_dy4 = 3.6 * 27.38 / 3503.71;
  
  double pdf_dy1 = 38.8 * 666.7 / 3503.71;
  double pdf_dy2 = 38.8 * 215.1 / 3503.71;
  double pdf_dy3 = 38.8 * 66.07 / 3503.71;
  double pdf_dy4 = 38.8 * 27.38 / 3503.71;
  
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_dy1JetsToLL.root", "dy1JetsToLL", 
					  xsec_dy1,
					  scaleUp_dy1, scaleDown_dy1,
					  pdf_dy1, pdf_dy1);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_dy2JetsToLL.root", "dy2JetsToLL", 
					  xsec_dy2,
					  scaleUp_dy2, scaleDown_dy2,
					  pdf_dy2, pdf_dy2);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_dy3JetsToLL.root", "dy3JetsToLL", 
					  xsec_dy3,
					  scaleUp_dy3, scaleDown_dy3,
					  pdf_dy3, pdf_dy3);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_dy4JetsToLL.root", "dy4JetsToLL", 
					  xsec_dy4,
					  scaleUp_dy4, scaleDown_dy4,
					  pdf_dy4, pdf_dy4);

  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_TBar_s.root", "TBar_s", 
					  1.76, 0.01, 0.01, 0.08, 0.08);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_TBar_t.root", "TBar_t", 
					  30.7, 0.7, 0.7, 0.9, 1.1);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_TBar_tW.root", "TBar_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_T_s.root", "T_s", 
					  3.79, 0.07, 0.07, 0.13, 0.13);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_T_t.root", "T_t", 
					  56.4, 2.1, 0.3, 1.1, 1.1);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_T_tW.root", "T_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7);

  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_WW.root", "WW",
					  57.1097, 2.3, 2.3, 2.0, 2.0);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_WZ.root", "WZ",
					  32.3161, 1.3, 1.3, 1.3, 1.3);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_ZZ.root", "ZZ",
					  8.25561, 0.3, 0.3, 0.3, 0.3);
  
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_WGToLNuG.root", "WGToLNuG",
					  553.9, 0.5 * 553.9, 0.5 * 553.9, 0.5 * 553.9, 0.5 * 553.9);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_ZGToLLG.root", "ZGToLLG",
					  159.1, 0.5 * 159.1, 0.5 * 159.1, 0.5 * 159.1, 0.5 * 159.1);

  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_TTWJets.root", "TTWJets", 
					  0.232, 0.067, 0.067, 0.03, 0.03);
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_TTZJets.root", "TTZJets", 
					  0.2057, 0., 0., 0.019, 0.024);

  double ttgamma_xsec = 0.033 * 9 + 0.148 * 12 + 0.8; // 2l(NLO) + l+jets(NLO) + all_had (approx)
  loadSuccess |= tMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v8/signal_contamination_TTGamma.root", "TTGamma",
					  ttgamma_xsec, 0.5 * ttgamma_xsec, 0.5 * ttgamma_xsec, 0.076 * ttgamma_xsec, 0.099 * ttgamma_xsec);

  if(!loadSuccess) return;

  tMaker->SetTrees(ggTree, qcdTree);
  tMaker->SetAddresses();
  tMaker->BookTemplates();

  tMaker->FillTemplates();

  tMaker->SubtractMCFromQCD();

  tMaker->SaveOutput();

  in->Close();

}
