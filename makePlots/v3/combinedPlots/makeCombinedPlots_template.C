enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};

void makePlots() {

  gROOT->LoadMacro("CreateCombinedPlots.C+");

  int controlRegion = kSR1;

  bool needsQCD = (controlRegion == kCR0 || controlRegion == kAny || controlRegion == kCR1);

  CreateCombinedPlots(controlRegion);

}

