enum photonRegions {kCR1, kCR2, kSR1, kSR2};

void maxMET(TFile * input, TString treeName, int region, Float_t& max) {

  Float_t met, ngamma, nfake;

  TTree * tree = (TTree*)input->Get(treeName);
  tree->SetBranchAddress("pfMET_t01", &met);
  tree->SetBranchAddress("Ngamma", &ngamma);
  tree->SetBranchAddress("Nfake", &nfake);

  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(region == kCR1 && (ngamma != 0 || nfake != 1)) continue;
    if(region == kCR2 && (ngamma != 0 || nfake < 2)) continue;
    if(region == kSR1 && ngamma != 1) continue;
    if(region == kSR2 && ngamma < 2) continue;

    if(met > max) max = met;
  }

  tree->ResetBranchAddresses();

}

void findOverflow() {

  TString prefix = "/eos/uscms/store/user/MYUSERNAME/stopBino_inputs/";

  const int nBkgs = 23;

  TString backgrounds[nBkgs] = {"ttJetsSemiLep", "ttJetsFullLep", "ttJetsHadronic",
				"TBar_s", "TBar_tW", "TBar_t", "T_s", "T_t", "T_tW",
				"W3JetsToLNu", "W4JetsToLNu",
				"dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
				"ZGToLLG", "WGToLNuG",
				"WW", "WZ", "ZZ",
				"TTZJets", "TTWJets",
				"TTGamma"};

  TString channels[2] = {"ele_bjj", "muon_bjj"};

  // CR1, CR2, SR1, SR2
  TCut photonRegion[4] = {"Ngamma==0 && Nfake==1", "Ngamma==0 && Nfake>=2", "Ngamma==1", "Ngamma>=2"};

  Float_t met;

  TFile * fEle = new TFile(prefix+"SingleElectron.root");
  TFile * fMuon = new TFile(prefix+"SingleMu.root");

  TFile * fBkg;

  Float_t max = 0.;

  // CR1
  maxMET(fEle, "ele_bjj_fakeTree", kCR1, max);
  maxMET(fEle, "ele_bjj_eQCDfakeTree", kCR1, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "ele_bjj_fakeTree", kCR1, max);
  }
  cout << "ele CR1 max = " << max << endl;

  max = 0.;
  maxMET(fMuon, "muon_bjj_fakeTree", kCR1, max);
  maxMET(fMuon, "muon_bjj_muQCDfakeTree", kCR1, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "muon_bjj_fakeTree", kCR1, max);
  }
  cout << "muon CR1 max = " << max << endl;

  // CR2
  maxMET(fEle, "ele_bjj_fakeTree", kCR2, max);
  maxMET(fEle, "ele_bjj_eQCDfakeTree", kCR2, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "ele_bjj_fakeTree", kCR2, max);
  }
  cout << "ele CR2 max = " << max << endl;

  max = 0.;
  maxMET(fMuon, "muon_bjj_fakeTree", kCR2, max);
  maxMET(fMuon, "muon_bjj_muQCDfakeTree", kCR2, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "muon_bjj_fakeTree", kCR2, max);
  }
  cout << "muon CR2 max = " << max << endl;

  // SR1
  maxMET(fEle, "ele_bjj_fakeTree", kSR1, max);
  maxMET(fEle, "ele_bjj_eQCDfakeTree", kSR1, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "ele_bjj_fakeTree", kSR1, max);
  }
  cout << "ele SR1 max = " << max << endl;

  max = 0.;
  maxMET(fMuon, "muon_bjj_fakeTree", kSR1, max);
  maxMET(fMuon, "muon_bjj_muQCDfakeTree", kSR1, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "muon_bjj_fakeTree", kSR1, max);
  }
  cout << "muon SR1 max = " << max << endl;

  // SR2
  maxMET(fEle, "ele_bjj_fakeTree", kSR2, max);
  maxMET(fEle, "ele_bjj_eQCDfakeTree", kSR2, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "ele_bjj_fakeTree", kSR2, max);
  }
  cout << "ele SR2 max = " << max << endl;

  max = 0.;
  maxMET(fMuon, "muon_bjj_fakeTree", kSR2, max);
  maxMET(fMuon, "muon_bjj_muQCDfakeTree", kSR2, max);
  for(int i = 0; i < nBkgs; i++) {
    fBkg = new TFile(prefix+"signal_contamination_"+backgrounds[i]+".root");
    maxMET(fBkg, "muon_bjj_fakeTree", kSR2, max);
  }
  cout << "muon SR2 max = " << max << endl;

}
