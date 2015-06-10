enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};

void makeTemplates() {

  gROOT->LoadMacro("CreateTemplates.C+");

  TString input_ele = "/eos/uscms/store/user/bfrancis/inputs_v8/SingleElectron.root";
  TString input_muon = "/eos/uscms/store/user/bfrancis/inputs_v8/SingleMu.root";

  // zFit -- mLepGammaLead in SR1, no b-tag (jjj), z+vgamma matched to ele, all else
  CreateTemplates(input_ele, "mLepGammaLead", kSR1, "ele_jjj", 40, 0., 400.);
  CreateTemplates(input_muon, "mLepGammaLead", kSR1, "muon_jjj", 40, 0., 400.);

  // qcdFit -- pfMET_t01 in pre-selection, bjj, no matching needed
  CreateTemplates(input_ele, "pfMET_t01", kAny, "ele_bjj", 200, 0., 2000.);
  CreateTemplates(input_muon, "pfMET_t01", kAny, "muon_bjj", 200, 0., 2000.);

  // m3 fit -- M3 in pre-selection, bjj, no matching needed
  CreateTemplates(input_ele, "m3", kAny, "ele_bjj", 200, 0., 2000.);
  CreateTemplates(input_muon, "m3", kAny, "muon_bjj", 200, 0., 2000.);

  // sigma fit -- leadSigmaIetaIeta in 'sigmaPlot', bjj, matching, with and without MET cut of 50
  //Float_t metCut = -1.0, , bool cutOnSigma = false, bool cutOnChHadIso = false)
  CreateTemplates(input_ele, "leadSigmaIetaIeta", kSigmaPlot, "ele_bjj", 80, 0., 0.04, -1.0, false, true);
  CreateTemplates(input_muon, "leadSigmaIetaIeta", kSigmaPlot, "muon_bjj", 80, 0., 0.04, -1.0, false, true);

  CreateTemplates(input_ele, "leadSigmaIetaIeta", kSigmaPlot, "ele_bjj", 80, 0., 0.04, 50.0, false, true);
  CreateTemplates(input_muon, "leadSigmaIetaIeta", kSigmaPlot, "muon_bjj", 80, 0., 0.04, 50.0, false, true);

  // chHadIso fit -- leadChargedHadronIso in 'sigmaPlot', bjj, matching, with and without MET cut of 50
  //Float_t metCut = -1.0, , bool cutOnSigma = false, bool cutOnChHadIso = false) {
  CreateTemplates(input_ele, "leadChargedHadronIso", kSigmaPlot, "ele_bjj", 110, -2.0, 20.0, -1.0, true, false);
  CreateTemplates(input_muon, "leadChargedHadronIso", kSigmaPlot, "muon_bjj", 110, -2.0, 20.0, -1.0, true, false);

  CreateTemplates(input_ele, "leadChargedHadronIso", kSigmaPlot, "ele_bjj", 110, -2.0, 20.0, 50.0, true, false);
  CreateTemplates(input_muon, "leadChargedHadronIso", kSigmaPlot, "muon_bjj", 110, -2.0, 20.0, 50.0, true, false);

  // dis-entangling the purity scale factor
  // sigma ==> cut on chhadiso
  CreateTemplates(input_ele, "pfMET_t01", kSR1, "ele_bjj", 200, 0., 2000., -1.0, false, true);
  CreateTemplates(input_muon, "pfMET_t01", kSR1, "muon_bjj", 200, 0., 2000., -1.0, false, true);
  
  CreateTemplates(input_ele, "pfMET_t01", kSR2, "ele_bjj", 200, 0., 2000., -1.0, false, true);
  CreateTemplates(input_muon, "pfMET_t01", kSR2, "muon_bjj", 200, 0., 2000., -1.0, false, true);

  // chhadiso ==> cut on sigma
  CreateTemplates(input_ele, "pfMET_t01", kSR1, "ele_bjj", 200, 0., 2000., -1.0, true, false);
  CreateTemplates(input_muon, "pfMET_t01", kSR1, "muon_bjj", 200, 0., 2000., -1.0, true, false);
  
  CreateTemplates(input_ele, "pfMET_t01", kSR2, "ele_bjj", 200, 0., 2000., -1.0, true, false);
  CreateTemplates(input_muon, "pfMET_t01", kSR2, "muon_bjj", 200, 0., 2000., -1.0, true, false);


}
