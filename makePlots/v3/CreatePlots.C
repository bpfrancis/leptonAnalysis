#include "CreatePlots.h"

using namespace std;

void CreatePlots(int channel, int controlRegion, bool needsQCD) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  // derived from CR0
  //systematic    topSF           topSFerror              wjetsSF         wjetsSFerror    QCDSF           QCDSFerror              MCSF            MCSFerror
  //e             0.840742963543  0.00666019431831        1.8169554357    0.0365078500672 0.329381930897  0.00388688459937        1.06790795756   0.00185289005357
  //mu            0.885471810863  0.00652764702301        1.5925431947    0.0358831914025 0.0145930197407 0.000532915435935       1.07191340778   0.00218856273796

  //e             0.858155455758  0.00530696636038        1.72687920543   0.0272087469619 0.329063044627  0.00388312156924        1.06815211075   0.0018533136753
  //mu            0.901305990296  0.00512420403339        1.50657303747   0.0263282690732 0.0145787340776 0.000532393744709	  1.07209994177   0.00218894359086


  Float_t sf_wJets = (channel < 2) ? 1.72687920543 : 1.50657303747;
  Float_t sfError_wJets = (channel < 2) ? 0.0272087469619 : 0.0263282690732;

  Float_t sf_ttbar = (channel < 2) ? 0.858155455758 : 0.901305990296;
  Float_t sfError_ttbar = (channel < 2) ? 0.00530696636038 : 0.00512420403339;

  Float_t sf_mc = (channel < 2) ? 1.06815211075 : 1.07209994177;
  Float_t sfError_mc = (channel < 2) ? 0.0018533136753 : 0.00218894359086;

  //Float_t sf_qcd = (channel < 2) ? 0.329063044627 : 0.0145787340776;
  //Float_t sfError_qcd = (channel < 2) ? 0.00388312156924 : 0.000532393744709;

  //durp
  sf_mc = -1.;
  sf_ttbar = -1.;

  PlotMaker * pMaker = new PlotMaker(channel, controlRegion, needsQCD);

  vector<TString> ttJets;
  ttJets.push_back("ttJetsSemiLep");
  ttJets.push_back("ttJetsFullLep");
  ttJets.push_back("ttJetsHadronic");
  pMaker->BookMCLayer(ttJets, kGray, "ttjets", "t#bar{t} inclusive", kGG, kTTbar, sf_ttbar, sfError_ttbar);

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

  vector<TString>  ttW;
  ttW.push_back("TTWJets");
  pMaker->BookMCLayer(ttW, kAzure-2, "ttW", "t#bar{t} + W", kQQ, kTTbar, sf_mc, sfError_mc);

  vector<TString>  ttZ;
  ttZ.push_back("TTZJets");
  pMaker->BookMCLayer(ttZ, kOrange-7, "ttZ", "t#bar{t} + Z", kGG, kTTbar, sf_mc, sfError_mc);

  vector<TString> ttgamma;
  ttgamma.push_back("ttA_2to5");
  pMaker->BookMCLayer(ttgamma, 8, "ttgamma", "t#bar{t} + #gamma", kGG, kTTbar, sf_mc, sfError_mc);

  ///////////////////////////////////////////////////////

  pMaker->BookPlot("pfMET", false,
		   "#slash{E}_{T} (GeV)", "Number of Events / GeV",
		   0., 300., 7.e-3, 2.5e5,
		   0., 1.9,
		   true, true, true);
  
  //if(controlRegion == kCR1) pMaker->SetDoRebinMET(true);

  pMaker->BookPlot("m3", false,
		   "M3 (GeV/c^{2})", "Number of Events / GeV",
		   0., 500., 7.e-5, 2.e5,
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



  pMaker->CreatePlots();

  delete pMaker;

}
