enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};

void makeHistograms() {

  gROOT->LoadMacro("CreateHistograms.C+");

  TString input_ele = "/eos/uscms/store/user/MYUSERNAME/stopBino_inputs/SingleElectron.root";
  TString input_muon = "/eos/uscms/store/user/MYUSERNAME/stopBino_inputs/SingleMu.root";

  TString metType = "pfMET_t01";
  bool useNormalTopReweighting = true;

  bool useWhizard = false; // false = use madgraph

  double metCut = -1.;

  int controlRegion = PHOTON_REGION;
  int photonMode = PHOTON_MODE;

  if(controlRegion == kSigmaPlot) metCut = 50.;

  bool blinded = false;

  const int nChannels = 4;
  TString channels[nChannels] = {"ele_bjj", "muon_bjj",
				 "ele_jjj", "muon_jjj"};

  for(int i = 0; i < nChannels; i++) {
    if(channels[i].Contains("ele")) CreateHistograms(input_ele, i, metCut, blinded, controlRegion, photonMode, metType, useNormalTopReweighting, useWhizard);
    else CreateHistograms(input_muon, i, metCut, blinded, controlRegion, photonMode, metType, useNormalTopReweighting, useWhizard);
  }  

}
