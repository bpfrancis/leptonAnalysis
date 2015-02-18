enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kNumControlRegions};

void makePlots() {

  gROOT->LoadMacro("CreatePlots.C+");

  TStopwatch ts;
  ts.Start();

  int controlRegion = PHOTON_REGION;

  bool needsQCD = true;

  const int nChannels = 8;

  for(int i = 0; i < nChannels; i++) {
    CreatePlots(i, controlRegion, needsQCD);
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}

