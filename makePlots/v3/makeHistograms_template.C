enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kNumControlRegions};

void makeHistograms() {

  gROOT->LoadMacro("CreateHistograms.C+");

  TStopwatch ts;
  ts.Start();

  TString input_ele = "ELE_FILE_TO_RUN";
  TString input_muon = "MUON_FILE_TO_RUN";

  double metCut = -1.;

  int controlRegion = PHOTON_REGION;

  bool blinded = true;

  const int nChannels = 4;
  TString channels[nChannels] = {"ele_jjj", "muon_jjj",
                                 "ele_bjj", "muon_bjj"};

  for(int i = 0; i < nChannels; i++) {
    if(controlRegion == kSigmaPlot && channels[i].Contains("jjj")) continue;

    if(channels[i].Contains("ele")) CreateHistograms(input_ele, i, metCut, blinded, controlRegion);
    else CreateHistograms(input_muon, i, metCut, blinded, controlRegion);
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
