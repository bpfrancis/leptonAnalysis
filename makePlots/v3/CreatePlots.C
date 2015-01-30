#include "CreatePlots.h"

using namespace std;

void CreatePlots(int channel, int controlRegion, bool needsQCD) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  PlotMaker * pMaker = new PlotMaker(channel, controlRegion, needsQCD);

  vector<TString> ttJets;
  ttJets.push_back("ttJetsSemiLep");
  ttJets.push_back("ttJetsFullLep");
  ttJets.push_back("ttJetsHadronic");
  pMaker->BookMCLayer(ttJets, kGray, "t#bar{t} inclusive");

  vector<TString> wJets;
  wJets.push_back("W1JetsToLNu");
  wJets.push_back("W2JetsToLNu");
  wJets.push_back("W3JetsToLNu");
  wJets.push_back("W4JetsToLNu");
  pMaker->BookMCLayer(wJets, kOrange-3, "W + Jets");

  vector<TString> zJets;
  zJets.push_back("dy1JetsToLL");
  zJets.push_back("dy2JetsToLL");
  zJets.push_back("dy3JetsToLL");
  zJets.push_back("dy4JetsToLL");
  pMaker->BookMCLayer(zJets, kYellow, "Z/#gamma* + Jets");

  vector<TString> singleTop;
  singleTop.push_back("TBar_s");
  singleTop.push_back("TBar_t");
  singleTop.push_back("TBar_tW");
  singleTop.push_back("T_s");
  singleTop.push_back("T_t");
  singleTop.push_back("T_tW");
  pMaker->BookMCLayer(singleTop, kRed, "Single top");

  vector<TString> diboson;
  diboson.push_back("WW");
  diboson.push_back("WZ");
  diboson.push_back("ZZ");
  pMaker->BookMCLayer(diboson, kCyan, "Diboson");

  vector<TString>  ttV;
  ttV.push_back("TTWJets");
  ttV.push_back("TTZJets");
  pMaker->BookMCLayer(ttV, kAzure-2, "t#bar{t} + W/Z");

  vector<TString> ttgamma;
  ttgamma.push_back("ttA_2to5");
  pMaker->BookMCLayer(ttgamma, 8, "t#bar{t} + #gamma");

  ///////////////////////////////////////////////////////

  pMaker->BookPlot("pfMET", true,
		   "#slash{E}_{T} (GeV)", "Number of Events / GeV",
		   0., 300., 7.e-3, 2.5e4,
		   0., 1.9,
		   true, true, true);

  pMaker->CreatePlots();

  delete pMaker;

}
