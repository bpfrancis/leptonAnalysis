#include "CreatePlots.h"

using namespace std;

void CreatePlots(int controlRegion) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  PlotMaker * pMaker = new PlotMaker(channel, controlRegion, needsQCD, metType, sf_qcd, sfError_qcd);

  pMaker->BookMCLayer("ttjets", kGray, "t#bar{t} + Jets");
  pMaker->BookMCLayer("wjets", kOrange-3, "W + Jets");
  pMaker->BookMCLayer("zjets", kYellow, "Z/#gamma* + Jets");
  pMaker->BookMCLayer("singleTop", kRed, "Single top");
  pMaker->BookMCLayer("diboson", kViolet-2, "VV");
  pMaker->BookMCLayer("vgamma", kAzure-2, "V#gamma");
  pMaker->BookMCLayer("ttW", kCyan, "t#bar{t} + V");
  pMaker->BookMCLayer("ttZ", kOrange-5, "t#bar{t} + V");
  pMaker->BookMCLayer("ttgamma", 8, "t#bar{t} + #gamma");

  ///////////////////////////////////////////////////////

  bool divideByWidth = true;
  TString xTitle = "#slash{E}_{T} (GeV)";
  TString yTitle = "Number of Events / GeV";

  Float_t xlo = 0.;
  Float_t xhi = 300.;

  Float_t ylo = 7.e-3;
  Float_t yhi = 2.5e6;

  Float_t ratiolo = 0.;
  Float_t ratiohi = 1.9;

  bool drawSignal = true;
  bool drawLegend = true;
  bool drawPrelim = true;
  bool usePasStyle = true;

  if(controlRegion == kCR1) {
    xhi = 800.;
    ylo = 2.e-4;
    yhi = 1.e2;
    ratiolo = 0.35;
    ratiohi = 1.65;
    drawSignal = false;
  }
  if(controlRegion == kCR2) {
    xhi = 799.9;
    ylo = 3.e-4;
    yhi = 40.0;
    ratiohi = 1.9;
    drawSignal = false;
  }
  if(controlRegion == kSR1) {
    xhi = 799.9;
    ylo = 7.e-4;
    yhi = 1.e2;
    ratiolo = 0.45;
    ratiohi = 1.55;
  }
  if(controlRegion == kSR2) {
    xhi = 799.9;
    ylo = 3.e-4;
    yhi = 2.0;
    ratiohi = 5.9;
  }

  pMaker->GoPlot(divideByWidth,
		 xTitle, yTitle,
		 xlo, xhi, ylo, yhi,
		 ratiolo, ratiohi,
		 drawSignal, drawLegend, drawPrelim,
		 usePasStyle);

  delete pMaker;

}
