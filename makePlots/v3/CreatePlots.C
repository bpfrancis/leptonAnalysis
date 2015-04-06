#include "CreatePlots.h"

using namespace std;

void CreatePlots(int channel, int controlRegion, bool needsQCD, TString metType, bool useWhizard) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  // derived from CR0
  //systematic    topSF           topSFerror              wjetsSF         wjetsSFerror    QCDSF           QCDSFerror              MCSF            MCSFerror
  //e             0.858155455758  0.00530696636038        1.72687920543   0.0272087469619 0.329063044627  0.00388312156924        1.06815211075   0.0018533136753
  //mu            0.901305990296  0.00512420403339        1.50657303747   0.0263282690732 0.0145787340776 0.000532393744709	  1.07209994177   0.00218894359086

  Float_t sf_wJets = (channels[channel].Contains("ele")) ? 1.72687920543 : 1.50657303747;
  Float_t sfError_wJets = (channels[channel].Contains("ele")) ? 0.0272087469619 : 0.0263282690732;

  Float_t sf_ttbar = (channels[channel].Contains("ele")) ? 0.858155455758 : 0.901305990296;
  Float_t sfError_ttbar = (channels[channel].Contains("ele")) ? 0.00530696636038 : 0.00512420403339;

  Float_t sf_mc = (channels[channel].Contains("ele")) ? 1.06815211075 : 1.07209994177;
  Float_t sfError_mc = (channels[channel].Contains("ele")) ? 0.0018533136753 : 0.00218894359086;

  //Float_t sf_qcd = (channels[channel].Contains("ele")) ? 0.329063044627 : 0.0145787340776;
  //Float_t sfError_qcd = (channels[channel].Contains("ele")) ? 0.00388312156924 : 0.000532393744709;

  // invert both sIetaIeta and chHadIso
  // fit sIetaIeta (0.006-0.02) --> central
  // also fit chHadIso (0-20) --> systematic
  // MET < 50

  Float_t sf_ttjets = (channels[channel].Contains("ele")) ? 0.948467461175 : 1.04147892559;
  Float_t sfError_ttjets = (channels[channel].Contains("ele")) ? 0.0213520487837 : 0.0219131353559;
  
  Float_t sf_ttgamma = (channels[channel].Contains("ele")) ? 0.802607703245 : 1.0182321871;
  Float_t sfError_ttgamma = (channels[channel].Contains("ele")) ? 0.156110926201 : 0.160640063318;

  sf_mc = -1.;
  sf_ttbar = -1.;

  PlotMaker * pMaker = new PlotMaker(channel, controlRegion, needsQCD, metType);

  vector<TString> ttJets;
  ttJets.push_back("ttJetsSemiLep");
  ttJets.push_back("ttJetsFullLep");
  ttJets.push_back("ttJetsHadronic");
  pMaker->BookMCLayer(ttJets, kGray, "ttjets", "t#bar{t} inclusive", kGG, kTTbar, sf_ttjets, sfError_ttjets);

  vector<TString> wJets;
  wJets.push_back("W1JetsToLNu");
  wJets.push_back("W2JetsToLNu");
  wJets.push_back("W3JetsToLNu");
  wJets.push_back("W4JetsToLNu");
  pMaker->BookMCLayer(wJets, kOrange-3, "wjets", "W + Jets", kQQ, kV, sf_wJets, sfError_wJets);

  vector<TString> zJets;
  zJets.push_back("dy1JetsToLL");
  zJets.push_back("dy2JetsToLL");
  zJets.push_back("dy3JetsToLL");
  zJets.push_back("dy4JetsToLL");
  pMaker->BookMCLayer(zJets, kYellow, "zjets", "Z/#gamma* + Jets", kQQ, kV, sf_mc, sfError_mc);

  vector<TString> singleTop;
  singleTop.push_back("TBar_s");
  singleTop.push_back("TBar_t");
  singleTop.push_back("TBar_tW");
  singleTop.push_back("T_s");
  singleTop.push_back("T_t");
  singleTop.push_back("T_tW");
  pMaker->BookMCLayer(singleTop, kRed, "singleTop", "Single top", kQG, kTTbar, sf_mc, sfError_mc);

  vector<TString> diboson;
  diboson.push_back("WW");
  diboson.push_back("WZ");
  diboson.push_back("ZZ");
  pMaker->BookMCLayer(diboson, kCyan, "diboson", "Diboson", kQQ, kVV, sf_mc, sfError_mc);

  vector<TString> vgamma;
  vgamma.push_back("WGToLNuG");
  vgamma.push_back("ZGToLLG");
  pMaker->BookMCLayer(vgamma, kRed, "vgamma", "V#gamma", kQQ, kVV, sf_mc, sfError_mc);

  vector<TString>  ttW;
  ttW.push_back("TTWJets");
  pMaker->BookMCLayer(ttW, kAzure-2, "ttW", "t#bar{t} + W", kQQ, kTTbar, sf_mc, sfError_mc);

  vector<TString>  ttZ;
  ttZ.push_back("TTZJets");
  pMaker->BookMCLayer(ttZ, kOrange-7, "ttZ", "t#bar{t} + Z", kGG, kTTbar, sf_mc, sfError_mc);

  vector<TString> ttgamma;
  if(useWhizard) ttgamma.push_back("ttA_2to5");
  else ttgamma.push_back("TTGamma");
  pMaker->BookMCLayer(ttgamma, 8, "ttgamma", "t#bar{t} + #gamma", kGG, kTTbar, sf_ttgamma, sfError_ttgamma);

  ///////////////////////////////////////////////////////

  pMaker->BookPlot(metType, true,
		   "#slash{E}_{T} (GeV)", "Number of Events / GeV",
		   0., 300., 7.e-3, 2.5e5,
		   0., 1.9,
		   true, true, true);
  
  //if(controlRegion == kCR1) pMaker->SetDoRebinMET(true);

  pMaker->BookPlot("m3", false,
		   "M3 (GeV/c^{2})", "Number of Events / GeV",
		   0., 300., 7.e-5, 2.e5,
		   0., 1.9,
		   true, true, true);
		   
  pMaker->BookPlot("Njets", false,
		   "Njets", "Number of Events",
		   2., 14., 2.e-3, 9.e5,
		   0., 1.9,
		   true, true, true);

  pMaker->BookPlot("Nbtags", false,
		   "Nbtags", "Number of Events",
		   0., 8., 2.e-3, 9.e5,
		   0., 1.9,
		   true, true, true);

  pMaker->BookPlot("HT_jets", true,
		   "HT_jets", "Number of Events / GeV",
		   0., 1200., 2.e-4, 4.e3,
		   0., 1.9,
		   true, true, true);

  pMaker->BookPlot("hadronic_pt", true, 
		   "MHT (GeV/c)", "Number of Events / GeV",
		   0, 1500, 2.e-5, 8.e3,
		   0., 2.1,
		   true, true, true);

  pMaker->BookPlot("HT", true,
		   "HT (GeV)", "Number of Events / GeV",
		   0, 2000, 2.e-4, 8.e3,
		   0., 2.1, 
		   true, true, true);

  pMaker->BookPlot("jet1_pt", true,
		   "P_{T} of leading jet", "Number of Events / GeV",
		   0, 1500, 2.e-4, 8.e3,
		   0., 2.1, 
		   true, true, true);

  pMaker->BookPlot("jet2_pt", true,
		   "P_{T} of sub-leading jet", "Number of Events / GeV",
		   0, 1200, 2.e-4, 8.e3,
		   0., 2.1, 
		   true, true, true);

  pMaker->BookPlot("jet3_pt", true,
		   "P_{T} of third-leading jet", "Number of Events / GeV",
		   0, 800, 2.e-4, 8.e3,
		   0., 2.1, 
		   true, true, true);

  pMaker->BookPlot("btag1_pt", true,
		   "P_{T} of leading btag", "Number of Events / GeV",
		   0, 1400, 2.e-4, 8.e3,
		   0., 2.1, 
		   true, true, true);

  pMaker->BookPlot("w_mT", true,
		   "W Transverse Mass", "Number of Events / GeV",
		   0, 1000, 2.e-4, 8.e3,
		   0., 2.1, 
		   true, true, true);

  pMaker->BookPlot("dPhi_met_l", false,
		   "#Delta#phi(#ell, #slash{E}_T)", "Number of Events",
		   -3.2, 3.2, 2.e-4, 8.e3,
		   0., 2.1,
		   true, false, false);

  pMaker->BookPlot("dPhi_met_ht", false,
		   "#Delta#phi(#ell, #vec{HT})", "Number of Events",
		   -3.2, 3.2, 2.e-4, 8.e3,
		   0., 2.1,
		   true, false, false);

  if(channels[channel].Contains("ele")) {
    pMaker->BookPlot("ele_pt", true,
		     "Electron P_{T}", "Number of Events / GeV",
		     0, 1500, 2.e-4, 8.e3,
		     0., 2.1, 
		     true, true, true);

    pMaker->BookPlot("ele_eta", false,
		     "Electron #Eta", "Number of Events / GeV",
		     -2.4, 2.4, 2.e-4, 8.e3,
		     0., 2.1, 
		     true, false, false);
  }
  else if(channels[channel].Contains("muon")) {
    pMaker->BookPlot("muon_pt", true,
		     "Muon P_{T}", "Number of Events / GeV",
		     0, 2000, 2.e-4, 8.e3,
		     0., 2.1, 
		     true, true, true);

    pMaker->BookPlot("muon_eta", false,
		     "Muon #Eta", "Number of Events / GeV",
		     -2.4, 2.4, 2.e-4, 8.e3,
		     0., 2.1, 
		     true, false, false);
  }

  if(controlRegion != kCR0 && controlRegion != kAny) {
    pMaker->BookPlot("leadSigmaIetaIeta", false,
		     "lead #sigma_{i#eta i#eta}", "Number of Events",
		     0, 0.035, 2.e-5, 8.e3,
		     0, 2.1,
		     true, true, true);

    pMaker->BookPlot("leadPhotonEta", false,
		     "#eta of leading #gamma", "Number of Events",
		     -1.5, 1.5, 2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false);

    pMaker->BookPlot("leadPhotonPhi", false,
		     "#phi of leading #gamma", "Number of Events",
		     -3.2, 3.2, 2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false);

    pMaker->BookPlot("leadPhotonEt", true,
		     "E_{T} of leading #gamma", "Number of Events / GeV",
		     0., 700., 2.e-4, 5.e2,
		     0., 5.1,
		     true, true, true);

    pMaker->BookPlot("mLepGammaLead", true,
		     "m_{#ell#gamma_{lead}}", "Number of Events / GeV",
		     0., 1200., 2.e-3, 5.e4,
		     0., 5.1, 
		     true, true, true);
  }

  if(controlRegion == kSR2 || controlRegion == kCR2 || controlRegion == kCR2a) {
    pMaker->BookPlot("trailSigmaIetaIeta", false,
		     "trail #sigma_{i#eta i#eta}", "Number of Events",
		     0, 0.035, 2.e-5, 8.e3,
		     0, 2.1,
		     true, true, true);

    pMaker->BookPlot("trailPhotonEta", false,
		     "#eta of trailing #gamma", "Number of Events",
		     -1.5, 1.5, 2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false);

    pMaker->BookPlot("trailPhotonPhi", false,
		     "#phi of trailing #gamma", "Number of Events",
		     -3.2, 3.2, 2.e-3, 3.e4,
		     0., 2.1,
		     false, false, false);

    pMaker->BookPlot("trailPhotonEt", true,
		     "E_{T} of trailing #gamma", "Number of Events / GeV",
		     0., 700., 2.e-4, 5.e2,
		     0., 5.1,
		     true, true, true);

    pMaker->BookPlot("mLepGammaTrail", true,
		     "m_{#ell#gamma_{trail}}", "Number of Events / GeV",
		     0., 1200., 2.e-3, 5.e4,
		     0., 5.1, 
		     true, true, true);

    pMaker->BookPlot("photon_dR", false,
		     "#DeltaR_{#gamma#gamma}", "Number of Events",
		     0.5, 5., 2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false);

    pMaker->BookPlot("photon_dPhi", false,
		     "#Delta#phi_{#gamma#gamma}", "Number of Events",
		     0., 3.14159, 2.e-2, 3.e5,
		     0., 2.1,
		     true, false, false);

    pMaker->BookPlot("photon_invmass", true,
		     "m_{#gamma#gamma} (GeV/c^{2})", "Number of Events",
		     0, 2000, 2.e-3, 3.e4,
		     0., 11.5,
		     true, true, true);

    pMaker->BookPlot("diEMpT", true,
		     "di-EM Pt", "Number of Events",
		     0, 1200, 2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true);
    
    pMaker->BookPlot("diJetPt", true,
		     "di-Jet Pt", "Number of Events",
		     0, 1400, 2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true);

    pMaker->BookPlot("mLepGammaGamma", true,
		     "m_{#ell#gamma#gamma}", "Number of Events",
		     0, 1200, 2.e-3, 5.e4,
		     0., 5.1,
		     true, true, true);
  }

  pMaker->CreatePlots();

  delete pMaker;

}
