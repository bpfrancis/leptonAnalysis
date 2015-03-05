#include "CreateHistograms.h"

using namespace std;

void CreateHistograms(TString input, int channel, double metCut, bool blinded, int controlRegion) {

  TFile * in = new TFile(input, "READ");

  TString ggName = channels[channel]+"_superFakeTree";
  TString qcdName = qcdChannels_superFake[channel];
  if(controlRegion == kSR1 || controlRegion == kSR2) {
    ggName = channels[channel]+"_signalTree";
    qcdName = qcdChannels[channel];
  }

  TTree * ggTree = (TTree*)in->Get(ggName);
  TTree * qcdTree = (TTree*)in->Get(qcdName);

  TString sigName = channels[channel]+"_superFakeTree";
  if(controlRegion == kSR1 || controlRegion == kSR2) sigName = chanenls[channel]+"_signalTree";

  TFile * fSigA = new TFile("/eos/uscms/store/user/bfrancis/inputs_v5/acceptance_v2/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSigA->Get(sigName);

  TFile * fSigB = new TFile("/eos/uscms/store/user/bfrancis/inputs_v5/acceptance_v2/signal_contamination_mst_560_m1_325.root", "READ");
  TTree * sigbTree = (TTree*)fSigB->Get(sigName);

  HistogramMaker * hMaker = new HistogramMaker(channel, blinded, controlRegion, metCut);
  hMaker->LoadLeptonSFs("../../data/lepton_SF_8TeV_53x_baseline.root");
  hMaker->LoadPhotonSFs("../../data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root");

  bool loadSuccess = true;
  
  Double_t wjetsSF, wjetsSFerror;
  Double_t ttjetsSF, ttjetsSFerror;
  Double_t ttgammaSF, ttgammaSFerror;

  wjetsSF = -1.;
  wjetsSFerror = 0.;
  ttjetsSF = -1.;
  ttjetsSFerror = 0.;
  ttgammaSF = -1.;
  ttgammaSFerror = 0.;

  Double_t ttbar_hadronic_xsec = 245.8 * 0.457;
  Double_t ttbar_semiLep_xsec  = 245.8 * 0.438;
  Double_t ttbar_fullLep_xsec  = 245.8 * 0.105;

  bool reallyDoTopPt = true;

  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic", 
					  ttbar_hadronic_xsec, ttbar_hadronic_xsec * 0.025, ttbar_hadronic_xsec * 0.034, ttbar_hadronic_xsec * 0.026, ttbar_hadronic_xsec * 0.026,
					  true, reallyDoTopPt,
					  ttjetsSF, ttjetsSFerror);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep", 
					  ttbar_semiLep_xsec, ttbar_semiLep_xsec * 0.025, ttbar_semiLep_xsec * 0.034, ttbar_semiLep_xsec * 0.026, ttbar_semiLep_xsec * 0.026,
					  true, reallyDoTopPt,
					  ttjetsSF, ttjetsSFerror);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep", 
					  ttbar_fullLep_xsec, ttbar_fullLep_xsec * 0.025, ttbar_fullLep_xsec * 0.034, ttbar_fullLep_xsec * 0.026, ttbar_fullLep_xsec * 0.026,
					  true, reallyDoTopPt,
					  ttjetsSF, ttjetsSFerror);

  /*
    loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_WJetsToLNu.root", "WJetsToLNu", 
    12234.4 * 3, 79.0, 39.7, 414.7, 414.7,
    false, false,
    wjetsSF, wjetsSFerror);
  */
  
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_W1JetsToLNu.root", "W1JetsToLNu", 
					  12234.4 * 3 * 6662. / 37509., 79.0 * 3 * 6662. / 37509., 39.7 * 3 * 6662. / 37509., 414.7 * 3 * 6662. / 37509., 414.7 * 3 * 6662. / 37509.,
					  false, false,
					  wjetsSF, wjetsSFerror);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_W2JetsToLNu.root", "W2JetsToLNu", 
					  12234.4 * 3 * 2159. / 37509., 79.0 * 3 * 2159. / 37509., 39.7 * 3 * 2159. / 37509., 414.7 * 3 * 2159. / 37509., 414.7 * 3 * 2159. / 37509.,
					  false, false,
					  wjetsSF, wjetsSFerror);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_W3JetsToLNu.root", "W3JetsToLNu", 
					  12234.4 * 3 * 640. / 37509., 79.0 * 3 * 640. / 37509., 39.7 * 3 * 640. / 37509., 414.7 * 3 * 640. / 37509., 414.7 * 3 * 640. / 37509.,
					  false, false,
					  wjetsSF, wjetsSFerror);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_W4JetsToLNu.root", "W4JetsToLNu", 
					  12234.4 * 3 * 264. / 37509., 79.0 * 3 * 264. / 37509., 39.7 * 3 * 264. / 37509., 414.7 * 3 * 264. / 37509., 414.7 * 3 * 264. / 37509.,
					  false, false,
					  wjetsSF, wjetsSFerror);


  /*
    loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_dyJetsToLL.root", "dyJetsToLL", 
    1177.3 * 3, 5.9, 3.6, 38.8, 38.8,
    false, false);
  */
  
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_dy1JetsToLL.root", "dy1JetsToLL", 
					  666.7 * 1177.3 * 3 / 3503.71, 
					  5.9 * 3 * 1177.3 / 3503.71, 3.6 * 3, 
					  38.8 * 3, 38.8 * 3,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_dy2JetsToLL.root", "dy2JetsToLL", 
					  215.1 * 1177.3 * 3 / 3503.71, 5.9 * 3 * 215.1 / 666.7, 3.6 * 3 * 215.1 / 666.7, 38.8 * 3 * 215.1 / 666.7, 38.8 * 3 * 215.1 / 666.7,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_dy3JetsToLL.root", "dy3JetsToLL", 
					  66.07 * 1177.3 * 3 / 3503.71, 5.9 * 3 * 66.07 / 666.7, 3.6 * 3 * 66.07 / 666.7, 38.8 * 3 * 66.07 / 666.7, 38.8 * 3 * 66.07 / 666.7,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_dy4JetsToLL.root", "dy4JetsToLL", 
					  27.38 * 1177.3 * 3 / 3503.71, 5.9 * 3 * 27.38 / 666.7, 3.6 * 3 * 27.38 / 666.7, 38.8 * 3 * 27.38 / 666.7, 38.8 * 3 * 27.38 / 666.7,
					  false, false);
  
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_TBar_s.root", "TBar_s", 
					  1.76, 0.01, 0.01, 0.08, 0.08,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_TBar_t.root", "TBar_t", 
					  30.7, 0.7, 0.7, 0.9, 1.1,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_TBar_tW.root", "TBar_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_T_s.root", "T_s", 
					  3.79, 0.07, 0.07, 0.13, 0.13,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_T_t.root", "T_t", 
					  56.4, 2.1, 0.3, 1.1, 1.1,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_T_tW.root", "T_tW", 
					  11.1, 0.3, 0.3, 0.7, 0.7,
					  false, false);

  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_WW.root", "WW",
					  57.1097, 2.3, 2.3, 2.0, 2.0,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_WZ.root", "WZ",
					  32.3161, 1.3, 1.3, 1.3, 1.3,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_ZZ.root", "ZZ",
					  8.25561, 0.3, 0.3, 0.3, 0.3,
					  false, false);
  
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_TTWJets.root", "TTWJets", 
					  0.232, 0.067, 0.067, 0.03, 0.03,
					  false, false);
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_TTZJets.root", "TTZJets", 
					  0.2057, 0., 0., 0.019, 0.024,
					  false, false);

  // http://arxiv.org/abs/1102.1967
  //loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_ttgjets.root", "ttgjets", 
  //2.166, 2.166 * .25, 2.166 * .25, 2.166 * 0.076, 2.166 * 0.099,
  //false, true,
  //ttgammaSF, ttgammaSFerror);

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/WhizardMCTeeTeeGamma#2_to_5_All_ttbar_decay_channels
  loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_ttA_2to5.root", "ttA_2to5", 
					  .9081 * 2, .9081 * .5, .9081 * .5, .9081 * 2 * 0.076, .9081 * 2 * 0.099, 
					  false, reallyDoTopPt,
					  ttgammaSF, ttgammaSFerror);

  //loadSuccess |= hMaker->LoadMCBackground("/eos/uscms/store/user/bfrancis/inputs_v5/signal_contamination_ttGG.root", "ttGG", 0.146, channel, 6, kCyan+3, "t#bar{t} + #gamma#gamma");

  if(!loadSuccess) return;

  hMaker->SetTrees(ggTree, qcdTree, sigaTree, sigbTree);

  if(controlRegion == kCR0) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("Ngamma", 4, 0., 4.);
    hMaker->BookHistogram("Nfake", 4, 0., 4.);
    hMaker->BookHistogram("pfMET", 200, 0., 2000.);
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
    
    hMaker->BookHistogram2D("Ngamma", "Nfake", 4, 0., 4., 4, 0., 4.);
  }

  if(controlRegion == kCR1 || controlRegion == kSR1 || controlRegion == kSigmaPlot) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("Ngamma", 4, 0., 4.);
    hMaker->BookHistogram("Nfake", 4, 0., 4.);
    hMaker->BookHistogram("pfMET", nMetBins_2g, xbins_met_2g);
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
    hMaker->BookHistogram("m3", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("ele_pt", nKinematicBins_1g, xbins_kinematic_1g);  // 15
    hMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_1g, xbins_kinematic_1g); // 17
    hMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);
    hMaker->BookHistogram("leadPhotonEt", nKinematicBins_1g, xbins_kinematic_1g); // 19
    hMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);                  // 20
    hMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159); // 21
    hMaker->BookHistogram("leadSigmaIetaIeta", 40, 0., 0.02);      // 22
    hMaker->BookHistogram("leadChargedHadronIso", 35, 0, 15.0);    // 23
    hMaker->BookHistogram("leadNeutralHadronIso", 50, -20., 30.);
    hMaker->BookHistogram("leadPhotonIso", 100, 0., 100.);
    hMaker->BookHistogram("mLepGammaLead", nKinematicBins_1g, xbins_kinematic_1g); // 24

    hMaker->BookHistogram2D("leadSigmaIetaIeta", "pfMET", 80, 0., 0.04, 20, 0., 350.);
    hMaker->BookHistogram2D("leadChargedHadronIso", "pfMET", 70, 0., 15., 20, 0., 350.);
    hMaker->BookHistogram2D("Ngamma", "Nfake", 4, 0., 4., 4, 0., 4.);
  }
  
  if(controlRegion == kSR2 || controlRegion == kCR2 || controlRegion == kCR2a) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("Ngamma", 4, 0., 4.);
    hMaker->BookHistogram("Nfake", 4, 0., 4.);
    hMaker->BookHistogram("pfMET", nMetBins_2g, xbins_met_2g);
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
    hMaker->BookHistogram("m3", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("ele_pt", nKinematicBins_2g, xbins_kinematic_2g);  // 15
    hMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_2g, xbins_kinematic_2g); // 17
    hMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);
    hMaker->BookHistogram("leadPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 19
    hMaker->BookHistogram("leadPhotonEta", 40, -1.5, 1.5);                  // 20
    hMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("leadSigmaIetaIeta", 40, 0., 0.02); // 22
    hMaker->BookHistogram("leadChargedHadronIso", 35, 0, 15.0);
    hMaker->BookHistogram("leadNeutralHadronIso", 50, -20., 30.);
    hMaker->BookHistogram("leadPhotonIso", 100, 0., 100.);
    hMaker->BookHistogram("mLepGammaLead", nKinematicBins_2g, xbins_kinematic_2g);
    
    hMaker->BookHistogram2D("leadSigmaIetaIeta", "pfMET", 80, 0., 0.04, 20, 0., 350.);
    hMaker->BookHistogram2D("leadChargedHadronIso", "pfMET", 70, 0., 15., 20, 0., 350.);
    hMaker->BookHistogram2D("Ngamma", "Nfake", 4, 0., 4., 4, 0., 4.);

    hMaker->BookHistogram("trailPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 25
    hMaker->BookHistogram("trailPhotonPhi", 63, -3.14159, 3.14159);          // 26
    hMaker->BookHistogram("trailPhotonEta", 40, -1.5, 1.5);                  // 27
    hMaker->BookHistogram("trailSigmaIetaIeta", 40, 0, 0.02);                // 28
    hMaker->BookHistogram("trailChargedHadronIso", 35, 0, 15.0);
    hMaker->BookHistogram("trailNeutralHadronIso", 50, -20., 30.);
    hMaker->BookHistogram("trailPhotonIso", 100, 0., 100.);
    hMaker->BookHistogram("diEMpT", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("diJetPt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("photon_invmass", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("photon_dR", 50, 0., 5.);
    hMaker->BookHistogram("photon_dPhi", 35, 0., 3.14159);
    hMaker->BookHistogram("mLepGammaTrail", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("mLepGammaGamma", nKinematicBins_2g, xbins_kinematic_2g);
  }

  hMaker->FillHistograms();
  hMaker->SubtractMCFromQCD();

  hMaker->SaveOutput();

  hMaker->CreateDatacards();

  in->Close();
  fSigA->Close();
  fSigB->Close();

}
