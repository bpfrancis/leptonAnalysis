#include "CreateHistograms.h"

using namespace std;

void CreateHistograms(TString input, int channel, double metCut, bool blinded, int controlRegion, int photonMode, TString metType, bool useNormalTopReweighting, bool useWhizard) {

  TFile * in = new TFile(input, "READ");

  TString sigName, qcdName;
  if(controlRegion == kSR1 || controlRegion == kSR2 || controlRegion == kCR0 || controlRegion == kAny) {
    sigName = channels[channel]+"_signalTree";
    qcdName = qcdChannels[channel];
  }
  else if(photonMode == kFake) {
    sigName = channels[channel]+"_fakeTree";
    qcdName = qcdChannels_fakePhotons[channel];
  }
  else {
    cout << "Invalid photonMode!" << endl;
    return;
  }

  TTree * ggTree = (TTree*)in->Get(sigName);
  TTree * qcdTree = (TTree*)in->Get(qcdName);

  TFile * fSigA = new TFile("/eos/uscms/store/user/bfrancis/inputs_v7/acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSigA->Get(sigName);

  TFile * fSigB = new TFile("/eos/uscms/store/user/bfrancis/inputs_v7/acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  TTree * sigbTree = (TTree*)fSigB->Get(sigName);

  HistogramMaker * hMaker = new HistogramMaker(channel, blinded, controlRegion, metCut, photonMode, metType, useNormalTopReweighting);
  hMaker->LoadLeptonSFs("/eos/uscms/store/user/bfrancis/data/lepton_SF_8TeV_53x_baseline.root");
  hMaker->LoadPhotonSFs("/eos/uscms/store/user/bfrancis/data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root");

  bool loadSuccess = true;

  Double_t ttbar_hadronic_xsec = 245.8 * 0.457;
  Double_t ttbar_semiLep_xsec  = 245.8 * 0.438;
  Double_t ttbar_fullLep_xsec  = 245.8 * 0.105;

  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic", 
					  ttbar_hadronic_xsec, ttbar_hadronic_xsec * 0.025, ttbar_hadronic_xsec * 0.034, ttbar_hadronic_xsec * 0.026, ttbar_hadronic_xsec * 0.026,
					  useWhizard, !useWhizard, true);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep", 
					  ttbar_semiLep_xsec, ttbar_semiLep_xsec * 0.025, ttbar_semiLep_xsec * 0.034, ttbar_semiLep_xsec * 0.026, ttbar_semiLep_xsec * 0.026,
					  useWhizard, !useWhizard, true);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep", 
					  ttbar_fullLep_xsec, ttbar_fullLep_xsec * 0.025, ttbar_fullLep_xsec * 0.034, ttbar_fullLep_xsec * 0.026, ttbar_fullLep_xsec * 0.026,
					  useWhizard, !useWhizard, true);

  /*
    loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_WJetsToLNu.root", "WJetsToLNu", 
    12234.4 * 3, 79.0, 39.7, 414.7, 414.7,
    false, false);
  */
  
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

  /* No longer using W1 and W2
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_W1JetsToLNu.root", "W1JetsToLNu", 
					  xsec_w1,
					  scaleUp_w1, scaleDown_w1,
					  pdf_w1, pdf_w1);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_W2JetsToLNu.root", "W2JetsToLNu", 
					  xsec_w2,
					  scaleUp_w2, scaleDown_w2,
					  pdf_w2, pdf_w2);
  */
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_W3JetsToLNu.root", "W3JetsToLNu", 
					  xsec_w3,
					  scaleUp_w3, scaleDown_w3,
					  pdf_w3, pdf_w3);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_W4JetsToLNu.root", "W4JetsToLNu", 
					  xsec_w4,
					  scaleUp_w4, scaleDown_w4,
					  pdf_w4, pdf_w4);

  /*
    loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_dyJetsToLL.root", "dyJetsToLL", 
    1177.3 * 3, 5.9, 3.6, 38.8, 38.8);
  */
  
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

  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_dy1JetsToLL.root", "dy1JetsToLL", 
					  xsec_dy1,
					  scaleUp_dy1, scaleDown_dy1,
					  pdf_dy1, pdf_dy1);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_dy2JetsToLL.root", "dy2JetsToLL", 
					  xsec_dy2,
					  scaleUp_dy2, scaleDown_dy2,
					  pdf_dy2, pdf_dy2);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_dy3JetsToLL.root", "dy3JetsToLL", 
					  xsec_dy3,
					  scaleUp_dy3, scaleDown_dy3,
					  pdf_dy3, pdf_dy3);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_dy4JetsToLL.root", "dy4JetsToLL", 
					  xsec_dy4,
					  scaleUp_dy4, scaleDown_dy4,
					  pdf_dy4, pdf_dy4);
  
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_TBar_s.root", "TBar_s", 
					  1.76, 0.01, 0.01, 0.08, 0.08);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_TBar_t.root", "TBar_t", 
					  30.7, 0.7, 0.7, 0.9, 1.1);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_TBar_tW.root", "TBar_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_T_s.root", "T_s", 
					  3.79, 0.07, 0.07, 0.13, 0.13);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_T_t.root", "T_t", 
					  56.4, 2.1, 0.3, 1.1, 1.1);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_T_tW.root", "T_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7);

  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_WW.root", "WW",
					  57.1097, 2.3, 2.3, 2.0, 2.0);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_WZ.root", "WZ",
					  32.3161, 1.3, 1.3, 1.3, 1.3);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_ZZ.root", "ZZ",
					  8.25561, 0.3, 0.3, 0.3, 0.3);
  
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_WGToLNuG.root", "WGToLNuG",
					  553.9, 0.5 * 553.9, 0.5 * 553.9, 0.5 * 553.9, 0.5 * 553.9);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_ZGToLLG.root", "ZGToLLG",
					  159.1, 0.5 * 159.1, 0.5 * 159.1, 0.5 * 159.1, 0.5 * 159.1);

  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_TTWJets.root", "TTWJets", 
					  0.232, 0.067, 0.067, 0.03, 0.03);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_TTZJets.root", "TTZJets", 
					  0.2057, 0., 0., 0.019, 0.024);

  // http://arxiv.org/abs/1102.1967
  //loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_ttgjets.root", "ttgjets", 
  //2.166, 2.166 * .25, 2.166 * .25, 2.166 * 0.076, 2.166 * 0.099);

  if(useWhizard) {
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/WhizardMCTeeTeeGamma#2_to_5_All_ttbar_decay_channels
    loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_ttA_2to5.root", "ttA_2to5", 
					    .9081 * 2, .9081 * .5, .9081 * .5, .9081 * 2 * 0.076, .9081 * 2 * 0.099);
  }
  else {
    // else use madgraph TTGamma 2to7
    double ttgamma_xsec = 0.033 * 9 + 0.148 * 12 + 0.8; // 2l(NLO) + l+jets(NLO) + all_had (approx)
    loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_TTGamma.root", "TTGamma",
					    ttgamma_xsec, 0.5 * ttgamma_xsec, 0.5 * ttgamma_xsec, 0.076 * ttgamma_xsec, 0.099 * ttgamma_xsec);
  }

  if(!loadSuccess) return;

  hMaker->SetTrees(ggTree, qcdTree, sigaTree, sigbTree);

  if(controlRegion == kCR0 || controlRegion == kAny) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("Ngamma", 4, 0., 4.);
    hMaker->BookHistogram("Nfake", 4, 0., 4.);
    hMaker->BookHistogram(metType, 200, 0., 2000.);
    hMaker->BookHistogram("HT", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("Njets", 20, 0., 20.);
    hMaker->BookHistogram("Nbtags", 20, 0., 20.);
    hMaker->BookHistogram("max_csv", 20, 0., 1.);
    hMaker->BookHistogram("submax_csv", 20, 0., 1.);
    hMaker->BookHistogram("HT_jets", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("hadronic_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet2_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet3_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("btag1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("w_mT", nKinematicBins_1g, xbins_kinematic_1g);    // 13
    hMaker->BookHistogram("m3", 200, 0., 2000.);
    hMaker->BookHistogram("ele_pt", nKinematicBins_1g, xbins_kinematic_1g);  // 15
    hMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_1g, xbins_kinematic_1g); // 17
    hMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);
    hMaker->BookHistogram("dPhi_met_l", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_met_ht", 63, -3.14159, 3.14159);

    hMaker->BookHistogram2D("Ngamma", "Nfake", 4, 0., 4., 4, 0., 4.);
  }

  if(controlRegion == kCR1 || controlRegion == kSR1 || controlRegion == kSigmaPlot) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("Ngamma", 4, 0., 4.);
    hMaker->BookHistogram("Nfake", 4, 0., 4.);
    hMaker->BookHistogram(metType, nMetBins_1g, xbins_met_1g);
    hMaker->BookHistogram("HT", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("Njets", 20, 0., 20.);
    hMaker->BookHistogram("Nbtags", 20, 0., 20.);
    hMaker->BookHistogram("max_csv", 20, 0., 1.);
    hMaker->BookHistogram("submax_csv", 20, 0., 1.);
    hMaker->BookHistogram("HT_jets", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("hadronic_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet2_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet3_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("btag1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("w_mT", nKinematicBins_1g, xbins_kinematic_1g);    // 13
    hMaker->BookHistogram("m3", 200, 0., 2000.);
    hMaker->BookHistogram("ele_pt", nKinematicBins_1g, xbins_kinematic_1g);  // 15
    hMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_1g, xbins_kinematic_1g); // 17
    hMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);
    hMaker->BookHistogram("dPhi_met_l", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_met_ht", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_leadPhoton_l", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_leadPhoton_b_min", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("leadPhotonEt", nKinematicBins_1g, xbins_kinematic_1g); // 19
    hMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);                  // 20
    hMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159); // 21
    hMaker->BookHistogram("leadSigmaIetaIeta", 80, 0., 0.04);      // 22
    hMaker->BookHistogram("leadChargedHadronIso", 110, -2.0, 20.0);    // 23
    hMaker->BookHistogram("leadNeutralHadronIso", 50, -20., 30.);
    hMaker->BookHistogram("leadPhotonIso", 100, 0., 100.);
    hMaker->BookHistogram("mLepGammaLead", 40, 0., 400.);

    hMaker->BookHistogram2D("leadSigmaIetaIeta", , 80, 0., 0.04, 20, 0., 350.);
    hMaker->BookHistogram2D("leadChargedHadronIso", metType, 110, -2.0, 20.0, 20, 0., 350.);
    hMaker->BookHistogram2D("Ngamma", "Nfake", 4, 0., 4., 4, 0., 4.);
  }
  
  if(controlRegion == kSR2 || controlRegion == kCR2 || controlRegion == kCR2a) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("Ngamma", 4, 0., 4.);
    hMaker->BookHistogram("Nfake", 4, 0., 4.);
    hMaker->BookHistogram(metType, nMetBins_2g, xbins_met_2g);
    hMaker->BookHistogram("HT", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("Njets", 20, 0., 20.);
    hMaker->BookHistogram("Nbtags", 20, 0., 20.);
    hMaker->BookHistogram("max_csv", 20, 0., 1.);
    hMaker->BookHistogram("submax_csv", 20, 0., 1.);
    hMaker->BookHistogram("HT_jets", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("hadronic_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("jet1_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("jet2_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("jet3_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("btag1_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("w_mT", nKinematicBins_2g, xbins_kinematic_2g);    // 13
    hMaker->BookHistogram("m3", 200, 0., 2000.);
    hMaker->BookHistogram("ele_pt", nKinematicBins_2g, xbins_kinematic_2g);  // 15
    hMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_2g, xbins_kinematic_2g); // 17
    hMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);
    hMaker->BookHistogram("dPhi_met_l", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_met_ht", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_leadPhoton_l", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_leadPhoton_b_min", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_trailPhoton_l", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("dPhi_trailPhoton_b_min", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("leadPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 19
    hMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);                  // 20
    hMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("leadSigmaIetaIeta", 80, 0., 0.04); // 22
    hMaker->BookHistogram("leadChargedHadronIso", 110, -2.0, 20.0);
    hMaker->BookHistogram("leadNeutralHadronIso", 50, -20., 30.);
    hMaker->BookHistogram("leadPhotonIso", 100, 0., 100.);
    hMaker->BookHistogram("mLepGammaLead", 40, 0., 400.);
    
    hMaker->BookHistogram2D("leadSigmaIetaIeta", metType, 80, 0., 0.04, 20, 0., 350.);
    hMaker->BookHistogram2D("leadChargedHadronIso", metType, 110, -2.0, 20.0, 20, 0., 350.);
    hMaker->BookHistogram2D("Ngamma", "Nfake", 4, 0., 4., 4, 0., 4.);

    hMaker->BookHistogram("trailPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 25
    hMaker->BookHistogram("trailPhotonPhi", 63, -3.14159, 3.14159);          // 26
    hMaker->BookHistogram("trailPhotonEta", 40, -1.5, 1.5);                  // 27
    hMaker->BookHistogram("trailSigmaIetaIeta", 80, 0., 0.04);                // 28
    hMaker->BookHistogram("trailChargedHadronIso", 110, -2.0, 20.0);
    hMaker->BookHistogram("trailNeutralHadronIso", 50, -20., 30.);
    hMaker->BookHistogram("trailPhotonIso", 100, 0., 100.);
    hMaker->BookHistogram("diEMpT", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("diJetPt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("photon_invmass", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("photon_dR", 50, 0., 5.);
    hMaker->BookHistogram("photon_dPhi", 35, 0., 3.14159);
    hMaker->BookHistogram("mLepGammaTrail", 40, 0., 400.);
    hMaker->BookHistogram("mLepGammaGamma", 40, 0., 400.);
  }

  hMaker->FillHistograms();
  hMaker->SubtractMCFromQCD();

  hMaker->SaveOutput();

  hMaker->CreateDatacards();

  in->Close();
  fSigA->Close();
  fSigB->Close();

}
