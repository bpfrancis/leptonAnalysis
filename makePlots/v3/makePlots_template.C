enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};

void makePlots() {

  gROOT->LoadMacro("CreatePlots.C+");

  int controlRegion = PHOTON_REGION;

  bool needsQCD = true;

  const int nChannels = 2;
  TString channels[nChannels] = {"ele_bjj", "muon_bjj"};

  for(int i = 0; i < nChannels; i++) {
    CreatePlots(i, controlRegion, needsQCD);
  }  

}

