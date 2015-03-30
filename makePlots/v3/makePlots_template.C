enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};

void makePlots() {

  gROOT->LoadMacro("CreatePlots.C+");

  TStopwatch ts;
  ts.Start();

  int controlRegion = PHOTON_REGION;

  bool needsQCD = true;

  const int nChannels = 4;
  TString channels[nChannels] = {"ele_jjj", "muon_jjj",
                                 "ele_bjj", "muon_bjj"};

  for(int i = 0; i < nChannels; i++) {
    if(controlRegion == kSigmaPlot && channels[i].Contains("jjj")) continue;

    CreatePlots(i, controlRegion, needsQCD);
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}

