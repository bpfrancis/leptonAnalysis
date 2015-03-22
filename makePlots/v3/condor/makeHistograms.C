enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kNumControlRegions};
enum photonModes {kSignal, kNoSigmaIetaIeta, kNoChargedHadronIso, kSuperFake, kNumPhotonModes};

void makeHistograms() {

  gROOT->Reset();
  gROOT->LoadMacro("CreateHistograms_C.so");

  TString input_ele = "/eos/uscms/store/user/bfrancis/inputs_v6/SingleElectron.root";
  TString input_muon = "/eos/uscms/store/user/bfrancis/inputs_v6/SingleMu.root";

  double metCut = -1.;

  int controlRegion = PHOTON_REGION;
  int photonMode = PHOTON_MODE;

  if(controlRegion == kSigmaPlot) metCut = 50.;

  bool blinded = false;

  const int nChannels = 4;
  TString channels[nChannels] = {"ele_jjj", "muon_jjj",
                                 "ele_bjj", "muon_bjj"};

  for(int i = 0; i < nChannels; i++) {
    if(controlRegion == kSigmaPlot && channels[i].Contains("jjj")) continue;

    if(channels[i].Contains("ele")) CreateHistograms(input_ele, i, metCut, blinded, controlRegion, photonMode);
    else CreateHistograms(input_muon, i, metCut, blinded, controlRegion, photonMode);
  }  

}
