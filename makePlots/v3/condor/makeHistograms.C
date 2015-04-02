enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};
enum photonModes {kSignal, kFake, kNumPhotonModes};

void makeHistograms() {

  gROOT->Reset();
  gROOT->LoadMacro("CreateHistograms_C.so");

  TString input_ele = "/eos/uscms/store/user/bfrancis/inputs_v7/SingleElectron.root";
  TString input_muon = "/eos/uscms/store/user/bfrancis/inputs_v7/SingleMu.root";

  double metCut = -1.;

  int controlRegion = PHOTON_REGION;
  int photonMode = PHOTON_MODE;

  if(controlRegion == kSigmaPlot) metCut = 50.;

  bool blinded = false;

  const int nChannels = 2;
  TString channels[nChannels] = {"ele_bjj", "muon_bjj"};

  for(int i = 0; i < nChannels; i++) {
    if(channels[i].Contains("ele")) CreateHistograms(input_ele, i, metCut, blinded, controlRegion, photonMode);
    else CreateHistograms(input_muon, i, metCut, blinded, controlRegion, photonMode);
  }  

}
