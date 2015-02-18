#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <stdio.h>

#include "rootRoutines.h"
#include "Binning.h"

using namespace std;

const int nChannels = 8;
TString channels[nChannels] = {"ele_jj", "ele_jjj",
			       "ele_bj", "ele_bjj",
			       "muon_jj", "muon_jjj",
			       "muon_bj", "muon_bjj"};

TString channelLabels[nChannels] = {"XYZ e (jj)", "XYZ e (jjj)",
				    "XYZ e (bj)", "XYZ e",
				    "XYZ #mu (jj)", "XYZ #mu (jjj)",
				    "XYZ #mu (bj)", "XYZ #mu"};

enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kNumControlRegions};
TString crNames[kNumControlRegions] = {"SR1", "SR2", "CR1", "CR2", "CR2a", "CR0"};

class PlotMaker : public TObject {

  ClassDef(PlotMaker, 1);

 public:
  PlotMaker(int chanNo, int cr, bool useQCD);
  virtual ~PlotMaker();

  void BookMCLayer(vector<TString> newNames, int color, TString legendEntry, Float_t scale = -1., Float_t scaleErr = -1.) { 
    TH1D * h;
    mc.push_back(h);

    mc_btagWeightUp.push_back(h); 
    mc_btagWeightDown.push_back(h);
    mc_puWeightUp.push_back(h); 
    mc_puWeightDown.push_back(h);
    mc_scaleUp.push_back(h); 
    mc_scaleDown.push_back(h);
    mc_pdfUp.push_back(h);
    mc_pdfDown.push_back(h);
    mc_topPtUp.push_back(h);
    mc_topPtDown.push_back(h);
    mc_JECUp.push_back(h);
    mc_JECDown.push_back(h);
    mc_leptonSFUp.push_back(h);
    mc_leptonSFDown.push_back(h);
    mc_photonSFUp.push_back(h);
    mc_photonSFDown.push_back(h);

    layerNames.push_back(newNames);
    layerColors.push_back(color);
    layerLegends.push_back(legendEntry);

    fitScales.push_back(scale);
    fitScaleErrors.push_back(scaleErr);
  };

  void BookPlot(TString variable, bool divideByWidth,
		TString xaxisTitle, TString yaxisTitle,
		Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax,
		Float_t ratiomin, Float_t ratiomax,
		bool drawSignal, bool drawLegend, bool drawPrelim);

  void DivideWidth() {
    data = (TH1D*)DivideByBinWidth(data);
    qcd = (TH1D*)DivideByBinWidth(qcd);
    siga = (TH1D*)DivideByBinWidth(siga);
    sigb = (TH1D*)DivideByBinWidth(sigb);
    
    bkg = (TH1D*)DivideByBinWidth(bkg);
    bkg_btagWeightUp = (TH1D*)DivideByBinWidth(bkg_btagWeightUp);
    bkg_btagWeightDown = (TH1D*)DivideByBinWidth(bkg_btagWeightDown);
    bkg_puWeightUp = (TH1D*)DivideByBinWidth(bkg_puWeightUp);
    bkg_puWeightDown = (TH1D*)DivideByBinWidth(bkg_puWeightDown);
    bkg_scaleUp = (TH1D*)DivideByBinWidth(bkg_scaleUp);
    bkg_scaleDown = (TH1D*)DivideByBinWidth(bkg_scaleDown);
    bkg_pdfUp = (TH1D*)DivideByBinWidth(bkg_pdfUp);
    bkg_pdfDown = (TH1D*)DivideByBinWidth(bkg_pdfDown);
    bkg_topPtUp = (TH1D*)DivideByBinWidth(bkg_topPtUp);
    bkg_topPtDown = (TH1D*)DivideByBinWidth(bkg_topPtDown);
    bkg_JECUp = (TH1D*)DivideByBinWidth(bkg_JECUp);
    bkg_JECDown = (TH1D*)DivideByBinWidth(bkg_JECDown);
    bkg_leptonSFUp = (TH1D*)DivideByBinWidth(bkg_leptonSFUp);
    bkg_leptonSFDown = (TH1D*)DivideByBinWidth(bkg_leptonSFDown);
    bkg_photonSFUp = (TH1D*)DivideByBinWidth(bkg_photonSFUp);
    bkg_photonSFDown = (TH1D*)DivideByBinWidth(bkg_photonSFDown);

    for(unsigned int i = 0; i < mc.size(); i++) {
      mc[i] = (TH1D*)DivideByBinWidth(mc[i]);
      mc_btagWeightUp[i] = (TH1D*)DivideByBinWidth(mc_btagWeightUp[i]);
      mc_btagWeightDown[i] = (TH1D*)DivideByBinWidth(mc_btagWeightDown[i]);
      mc_puWeightUp[i] = (TH1D*)DivideByBinWidth(mc_puWeightUp[i]);
      mc_puWeightDown[i] = (TH1D*)DivideByBinWidth(mc_puWeightDown[i]);
      mc_scaleUp[i] = (TH1D*)DivideByBinWidth(mc_scaleUp[i]);
      mc_scaleDown[i] = (TH1D*)DivideByBinWidth(mc_scaleDown[i]);
      mc_pdfUp[i] = (TH1D*)DivideByBinWidth(mc_pdfUp[i]);
      mc_pdfDown[i] = (TH1D*)DivideByBinWidth(mc_pdfDown[i]);
      mc_topPtUp[i] = (TH1D*)DivideByBinWidth(mc_topPtUp[i]);
      mc_topPtDown[i] = (TH1D*)DivideByBinWidth(mc_topPtDown[i]);
      mc_JECUp[i] = (TH1D*)DivideByBinWidth(mc_JECUp[i]);
      mc_JECDown[i] = (TH1D*)DivideByBinWidth(mc_JECDown[i]);
      mc_leptonSFUp[i] = (TH1D*)DivideByBinWidth(mc_leptonSFUp[i]);
      mc_leptonSFDown[i] = (TH1D*)DivideByBinWidth(mc_leptonSFDown[i]);
      mc_photonSFUp[i] = (TH1D*)DivideByBinWidth(mc_photonSFUp[i]);
      mc_photonSFDown[i] = (TH1D*)DivideByBinWidth(mc_photonSFDown[i]);
    }

    errors_stat = (TH1D*)DivideByBinWidth(errors_stat);
    errors_sys = (TH1D*)DivideByBinWidth(errors_sys);

  };

  // done once for all variables

  void MakeCanvas() {
    can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
    padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
    padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
    padhi->SetLogy(true);
    padhi->SetTickx(true);
    padhi->SetTicky(true);
    //padhi->SetGridx(true);
    //padhi->SetGridy(true);
    padhi->SetBottomMargin(0);

    padlo->SetTopMargin(0);
    padlo->SetBottomMargin(0.2);

    padhi->Draw();
    padlo->Draw();
  };

  void MakeLegends();

  // done for each variable

  void GetHistograms(unsigned int n);
  void StackHistograms(unsigned int n);
  void CalculateRatio(unsigned int n);
  void SetStyles(unsigned int n);

  void ScaleByFit(unsigned int n, vector<TH1D*>& h);

  void CalculateQCDNormalization();
  void ScaleQCD();

  void CreatePlot(unsigned int n);
  void CreatePlots() {
    MakeCanvas();
    CalculateQCDNormalization();
    for(unsigned int i = 0; i < variables.size(); i++) CreatePlot(i);
  };
  
  void CreateSignalOutputs();
  void METDifference();
  void SaveLimitOutputs();

  bool inControlRegion(Float_t ngamma, Float_t nfake) {
    switch(controlRegion) {
    case kSR1:
      return ngamma == 1;
    case kSR2:
      return ngamma >= 2;
    case kCR1:
      return (ngamma == 0 && nfake == 1);
    case kCR2:
      return (ngamma == 0 && nfake >= 2);
    case kCR2a:
      return (ngamma == 1 && nfake == 1) || (ngamma == 0 && nfake >= 2);
    case kCR0:
      return ngamma == 0;
    default:
      return false;
    }
  };

 private:
  TFile * input;

  // only one copy of each below for all variables -- for each variable, copy over histograms
  
  Float_t qcdScale, qcdScale_defUp, qcdScale_defDown;
  Float_t qcdScaleError;

  TH1D * data;
  
  TH1D * siga;
  TH1D * sigb;

  TH1D * qcd;
  TH1D * qcd_defUp;
  TH1D * qcd_defDown;

  vector<TH1D*> mc;
  vector<TH1D*> mc_btagWeightUp, mc_btagWeightDown;
  vector<TH1D*> mc_puWeightUp, mc_puWeightDown;
  vector<TH1D*> mc_scaleUp, mc_scaleDown;
  vector<TH1D*> mc_pdfUp, mc_pdfDown;
  vector<TH1D*> mc_topPtUp, mc_topPtDown;
  vector<TH1D*> mc_JECUp, mc_JECDown;
  vector<TH1D*> mc_leptonSFUp, mc_leptonSFDown;
  vector<TH1D*> mc_photonSFUp, mc_photonSFDown;

  TH1D *bkg, 
    *bkg_btagWeightUp, *bkg_btagWeightDown, 
    *bkg_puWeightUp, *bkg_puWeightDown, 
    *bkg_scaleUp, *bkg_scaleDown, 
    *bkg_pdfUp, *bkg_pdfDown, 
    *bkg_topPtUp, *bkg_topPtDown, 
    *bkg_JECUp, *bkg_JECDown,
    *bkg_leptonSFUp, *bkg_leptonSFDown,
    *bkg_photonSFUp, *bkg_photonSFDown;

  TH1D * errors_stat;
  TH1D * errors_sys;
  
  TH1D * ratio;
  TH1D * ratio_stat;
  TH1D * ratio_sys;

  TLine * oneLine;

  TCanvas * can;
  TPad * padhi;
  TPad * padlo;

  TPaveText * reqText;
  TPaveText * lumiHeader;

  // these three stay the same for all MC and all variables

  vector< vector<TString> > layerNames;
  vector<int> layerColors;
  vector<TString> layerLegends;

  vector<Float_t> fitScales;
  vector<Float_t> fitScaleErrors;

  // below exist vectors of qualities for each variable -- these are what is looped over

  vector<TString> variables;
  vector<bool> divideByBinWidth;
  vector<TString> xaxisTitles, yaxisTitles;
  vector<Float_t> xMinimums, xMaximums, yMinimums, yMaximums;
  vector<Float_t> ratioMinimums, ratioMaximums;
  vector<bool> doDrawSignal, doDrawLegend, doDrawPrelim;

  // channel-wide qualities

  TString channel;
  TString channelLabel;
  int controlRegion;
  bool needsQCD;

  TLegend * leg;
  TLegend * legDrawSignal;
  TLegend * ratioLeg;

};

PlotMaker::PlotMaker(int chanNo, int cr, bool useQCD) {

  channel = channels[chanNo];
  channelLabel = channelLabels[chanNo];
  controlRegion = cr;
  needsQCD = useQCD;

  mc.clear();
  mc_btagWeightUp.clear();
  mc_btagWeightDown.clear();
  mc_puWeightUp.clear();
  mc_puWeightDown.clear();
  mc_scaleUp.clear();
  mc_scaleDown.clear();
  mc_pdfUp.clear();
  mc_pdfDown.clear();
  mc_topPtUp.clear();
  mc_topPtDown.clear();
  mc_JECUp.clear();
  mc_JECDown.clear();
  mc_leptonSFUp.clear();
  mc_leptonSFDown.clear();
  mc_photonSFUp.clear();
  mc_photonSFDown.clear();
  
  layerNames.clear();
  layerColors.clear();
  layerLegends.clear();

  fitScales.clear();
  fitScaleErrors.clear();

  variables.clear();
  divideByBinWidth.clear();
  xaxisTitles.clear();
  yaxisTitles.clear();
  xMinimums.clear();
  xMaximums.clear();
  yMinimums.clear();
  yMaximums.clear();
  ratioMinimums.clear();
  ratioMaximums.clear();
  doDrawSignal.clear();
  doDrawLegend.clear();
  doDrawPrelim.clear();

  input = new TFile("histograms_"+channel+"_"+crNames[controlRegion]+".root", "READ"); 

}

PlotMaker::~PlotMaker() {

  mc.clear();
  mc_btagWeightUp.clear();
  mc_btagWeightDown.clear();
  mc_puWeightUp.clear();
  mc_puWeightDown.clear();
  mc_scaleUp.clear();
  mc_scaleDown.clear();
  mc_pdfUp.clear();
  mc_pdfDown.clear();
  mc_topPtUp.clear();
  mc_topPtDown.clear();
  mc_JECUp.clear();
  mc_JECDown.clear();
  mc_leptonSFUp.clear();
  mc_leptonSFDown.clear();
  mc_photonSFUp.clear();
  mc_photonSFDown.clear();

  layerNames.clear();
  layerColors.clear();
  layerLegends.clear();

  fitScales.clear();
  fitScaleErrors.clear();

  variables.clear();
  divideByBinWidth.clear();
  xaxisTitles.clear();
  yaxisTitles.clear();
  xMinimums.clear();
  xMaximums.clear();
  yMinimums.clear();
  yMaximums.clear();
  ratioMinimums.clear();
  ratioMaximums.clear();
  doDrawSignal.clear();
  doDrawLegend.clear();
  doDrawPrelim.clear();

  delete leg;
  delete legDrawSignal;
  delete ratioLeg;

  delete can;

  input->Close();

}

void PlotMaker::GetHistograms(unsigned int n) {

  data = (TH1D*)input->Get(variables[n]+"_gg_"+channel);

  qcd = (TH1D*)input->Get(variables[n]+"_qcd_"+channel);
  qcd_defUp = (TH1D*)input->Get(variables[n]+"_qcd_10_"+channel);

  qcd_defDown = (TH1D*)qcd->Clone(variables[n]+"_qcd_down_"+channel);
  
  for(int ibin = 0; ibin < qcd_defUp->GetNbinsX(); ibin++) {
    // central - (up - central) = 2*central - up
    Float_t value = 2. * qcd->GetBinContent(ibin+1) - qcd_defUp->GetBinContent(ibin+1);
    qcd_defDown->SetBinContent(ibin+1, value);
  }

  siga = (TH1D*)input->Get(variables[n]+"_a_"+channel);
  sigb = (TH1D*)input->Get(variables[n]+"_b_"+channel);

  for(unsigned int i = 0; i < layerNames.size(); i++) {
    mc[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel);

    mc_btagWeightUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_btagWeightUp");
    mc_btagWeightDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_btagWeightDown");

    mc_puWeightUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_puWeightUp");
    mc_puWeightDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_puWeightDown");

    mc_scaleUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_scaleUp");
    mc_scaleDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_scaleDown");

    mc_pdfUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_pdfUp");
    mc_pdfDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_pdfDown");

    mc_topPtUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_topPtUp");
    mc_topPtDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_topPtDown");

    mc_JECUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_JECUp");
    mc_JECDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_JECDown");

    mc_leptonSFUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_leptonSFUp");
    mc_leptonSFDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_leptonSFDown");

    mc_photonSFUp[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_photonSFUp");
    mc_photonSFDown[i] = (TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_photonSFDown");

    for(unsigned int j = 1; j < layerNames[i].size(); j++) {
      mc[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel));

      mc_btagWeightUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_btagWeightUp"));
      mc_btagWeightDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_btagWeightDown"));
      
      mc_puWeightUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_puWeightUp"));
      mc_puWeightDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_puWeightDown"));
      
      mc_scaleUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_scaleUp"));
      mc_scaleDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_scaleDown"));
      
      mc_pdfUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_pdfUp"));
      mc_pdfDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_pdfDown"));
      
      mc_topPtUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_topPtUp"));
      mc_topPtDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_topPtDown"));
      
      mc_JECUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_JECUp"));
      mc_JECDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_JECDown"));
      
      mc_leptonSFUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_leptonSFUp"));
      mc_leptonSFDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_leptonSFDown"));
      
      mc_photonSFUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_photonSFUp"));
      mc_photonSFDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][0]+"_"+channel+"_photonSFDown"));

    }
  }

  for(unsigned int i = 0; i < mc.size(); i++) {
    if(fitScales[i] > 0) {
      ScaleByFit(i, mc);
      ScaleByFit(i, mc_btagWeightUp);
      ScaleByFit(i, mc_btagWeightDown);
      ScaleByFit(i, mc_puWeightUp);
      ScaleByFit(i, mc_puWeightDown);
      ScaleByFit(i, mc_scaleUp);
      ScaleByFit(i, mc_scaleDown);
      ScaleByFit(i, mc_pdfUp);
      ScaleByFit(i, mc_pdfDown);
      ScaleByFit(i, mc_topPtUp);
      ScaleByFit(i, mc_topPtDown);
      ScaleByFit(i, mc_JECUp);
      ScaleByFit(i, mc_JECDown);
      ScaleByFit(i, mc_leptonSFUp);
      ScaleByFit(i, mc_leptonSFDown);
      ScaleByFit(i, mc_photonSFUp);
      ScaleByFit(i, mc_photonSFDown);
    }
  }

}

void PlotMaker::BookPlot(TString variable, bool divideByWidth,
			 TString xaxisTitle, TString yaxisTitle,
			 Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax,
			 Float_t ratiomin, Float_t ratiomax,
			 bool drawSignal, bool drawLegend, bool drawPrelim) {

  variables.push_back(variable);
  divideByBinWidth.push_back(divideByWidth);
  xaxisTitles.push_back(xaxisTitle);
  yaxisTitles.push_back(yaxisTitle);
  xMinimums.push_back(xmin);
  xMaximums.push_back(xmax);
  yMinimums.push_back(ymin);
  yMaximums.push_back(ymax);
  ratioMinimums.push_back(ratiomin);
  ratioMaximums.push_back(ratiomax);
  doDrawSignal.push_back(drawSignal);
  doDrawLegend.push_back(drawLegend);
  doDrawPrelim.push_back(drawPrelim);

}

void PlotMaker::StackHistograms(unsigned int n) {

  if(needsQCD) {
    bkg = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]);
    bkg_btagWeightUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_btagWeightUp");
    bkg_btagWeightDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_btagWeightDown");
    bkg_puWeightUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_puWeightUp");
    bkg_puWeightDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_puWeightDown");
    bkg_scaleUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_scaleUp");
    bkg_scaleDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_scaleDown");
    bkg_pdfUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_pdfUp");
    bkg_pdfDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_pdfDown");
    bkg_topPtUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_topPtUp");
    bkg_topPtDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_topPtDown");
    bkg_JECUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_JECUp");
    bkg_JECDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_JECDown");
    bkg_leptonSFUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_leptonSFUp");
    bkg_leptonSFDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_leptonSFDown");
    bkg_photonSFUp = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_photonSFUp");
    bkg_photonSFDown = (TH1D*)qcd->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_photonSFDown");
  }
  
  else {
    bkg = (TH1D*)mc[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]);
    bkg_btagWeightUp = (TH1D*)mc_btagWeightUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_btagWeightUp");
    bkg_btagWeightDown = (TH1D*)mc_btagWeightDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_btagWeightDown");
    bkg_puWeightUp = (TH1D*)mc_puWeightUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_puWeightUp");
    bkg_puWeightDown = (TH1D*)mc_puWeightDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_puWeightDown");
    bkg_scaleUp = (TH1D*)mc_scaleUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_scaleUp");
    bkg_scaleDown = (TH1D*)mc_scaleDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_scaleDown");
    bkg_pdfUp = (TH1D*)mc_pdfUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_pdfUp");
    bkg_pdfDown = (TH1D*)mc_pdfDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_pdfDown");
    bkg_topPtUp = (TH1D*)mc_topPtUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_topPtUp");
    bkg_topPtDown = (TH1D*)mc_topPtDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_topPtDown");
    bkg_JECUp = (TH1D*)mc_JECUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_JECUp");
    bkg_JECDown = (TH1D*)mc_JECDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_JECDown");
    bkg_leptonSFUp = (TH1D*)mc_leptonSFUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_leptonSFUp");
    bkg_leptonSFDown = (TH1D*)mc_leptonSFDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_leptonSFDown");
    bkg_photonSFUp = (TH1D*)mc_photonSFUp[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_photonSFUp");
    bkg_photonSFDown = (TH1D*)mc_photonSFDown[0]->Clone(variables[n]+"_bkg_"+crNames[controlRegion]+"_photonSFDown");
  }

  for(unsigned int i = 0; i < mc.size(); i++) {

    if(!needsQCD && i == 0) continue;

    bkg->Add(mc[i]);
    bkg_btagWeightUp->Add(mc_btagWeightUp[i]);
    bkg_btagWeightDown->Add(mc_btagWeightDown[i]);
    bkg_puWeightUp->Add(mc_puWeightUp[i]);
    bkg_puWeightDown->Add(mc_puWeightDown[i]);
    bkg_scaleUp->Add(mc_scaleUp[i]);
    bkg_scaleDown->Add(mc_scaleDown[i]);
    bkg_pdfUp->Add(mc_pdfUp[i]);
    bkg_pdfDown->Add(mc_pdfDown[i]);
    bkg_topPtUp->Add(mc_topPtUp[i]);
    bkg_topPtDown->Add(mc_topPtDown[i]);
    bkg_JECUp->Add(mc_JECUp[i]);
    bkg_JECDown->Add(mc_JECDown[i]);
    bkg_leptonSFUp->Add(mc_leptonSFUp[i]);
    bkg_leptonSFDown->Add(mc_leptonSFDown[i]);
    bkg_photonSFUp->Add(mc_photonSFUp[i]);
    bkg_photonSFDown->Add(mc_photonSFDown[i]);
    
    for(unsigned int j = i + 1; j < mc.size(); j++) {
      mc[i]->Add(mc[j]);
    }
  }
  
}

void PlotMaker::CalculateRatio(unsigned int n) {

  errors_stat = (TH1D*)bkg->Clone("errors_stat_"+variables[n]);

  errors_sys = (TH1D*)bkg->Clone("errors_sys_"+variables[n]);

  ratio = (TH1D*)data->Clone("ratio_"+variables[n]);
  ratio->Reset();
  ratio->SetTitle("Data / Background");
  for(int i = 0; i < ratio->GetNbinsX(); i++) {
    if(bkg->GetBinContent(i+1) == 0.) continue;
    ratio->SetBinContent(i+1, data->GetBinContent(i+1) / bkg->GetBinContent(i+1));
    ratio->SetBinError(i+1, data->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  ratio_stat = (TH1D*)bkg->Clone("ratio_stat_"+variables[n]);
  for(int i = 0; i < ratio_stat->GetNbinsX(); i++) {
    ratio_stat->SetBinContent(i+1, 1.);
    if(bkg->GetBinContent(i+1) == 0.) ratio_stat->SetBinError(i+1, 0.);
    else ratio_stat->SetBinError(i+1, bkg->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  ratio_sys = (TH1D*)bkg->Clone("ratio_sys_"+variables[n]);

  for(int i = 0; i < errors_sys->GetNbinsX(); i++) {

    Double_t stat = bkg->GetBinError(i+1);

    Double_t qcdUp = qcd_defUp->GetBinContent(i+1);
    Double_t qcdDown = qcd_defDown->GetBinContent(i+1);
    Double_t qcd_sys = fabs(qcdUp - qcdDown) / 2.;

    Double_t btagUp = bkg_btagWeightUp->GetBinContent(i+1);
    Double_t btagDown = bkg_btagWeightDown->GetBinContent(i+1);
    Double_t btag_sys = fabs(btagUp - btagDown) / 2.;

    Double_t puUp = bkg_puWeightUp->GetBinContent(i+1);
    Double_t puDown = bkg_puWeightDown->GetBinContent(i+1);
    Double_t pu_sys = fabs(puUp - puDown) / 2.;

    Double_t scaleUp = bkg_scaleUp->GetBinContent(i+1);
    Double_t scaleDown = bkg_scaleDown->GetBinContent(i+1);
    Double_t scale_sys = fabs(scaleUp - scaleDown) / 2.;
    
    Double_t pdfUp = bkg_pdfUp->GetBinContent(i+1);
    Double_t pdfDown = bkg_pdfDown->GetBinContent(i+1);
    Double_t pdf_sys = fabs(pdfUp - pdfDown) / 2.;

    Double_t topPtUp = bkg_topPtUp->GetBinContent(i+1);
    Double_t topPtDown = bkg_topPtDown->GetBinContent(i+1);
    Double_t topPt_sys = fabs(topPtUp - topPtDown) / 2.;

    Double_t JECup = bkg_JECUp->GetBinContent(i+1);
    Double_t JECdown = bkg_JECDown->GetBinContent(i+1);
    Double_t JEC_sys = fabs(JECup - JECdown) / 2.;

    Double_t leptonSFup = bkg_leptonSFUp->GetBinContent(i+1);
    Double_t leptonSFdown = bkg_leptonSFDown->GetBinContent(i+1);
    Double_t leptonSF_sys = fabs(leptonSFup - leptonSFdown) / 2.;

    Double_t photonSFup = bkg_photonSFUp->GetBinContent(i+1);
    Double_t photonSFdown = bkg_photonSFDown->GetBinContent(i+1);
    Double_t photonSF_sys = fabs(photonSFup - photonSFdown) / 2.;

    Double_t totalError2 = stat*stat + 
      qcd_sys*qcd_sys +
      btag_sys*btag_sys +
      pu_sys*pu_sys +
      scale_sys*scale_sys + 
      pdf_sys*pdf_sys + 
      topPt_sys*topPt_sys + 
      JEC_sys*JEC_sys + 
      leptonSF_sys*leptonSF_sys + 
      photonSF_sys*photonSF_sys;

    if(bkg->GetBinContent(i+1) == 0.) {
      errors_sys->SetBinError(i+1, 0.);
      ratio_sys->SetBinError(i+1, 0.);
    }
    else {
      errors_sys->SetBinError(i+1, sqrt(totalError2));
      ratio_sys->SetBinError(i+1, sqrt(totalError2) / bkg->GetBinContent(i+1));
    }

    ratio_sys->SetBinContent(i+1, 1.);

  }

}

void PlotMaker::MakeLegends() {

  leg = new TLegend(0.45, 0.6, 0.85, 0.85, NULL, "brNDC");
  leg->SetNColumns(2);
  leg->AddEntry(data, "Data", "LP");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(errors_sys, "Stat. #oplus Syst. Errors", "F");
  leg->AddEntry((TObject*)0, "", "");
  if(needsQCD) leg->AddEntry(bkg, "QCD", "F");

  leg->AddEntry(mc[0], layerLegends[0], "F");
  for(unsigned int i = 1; i < mc.size(); i++) leg->AddEntry(mc[i], layerLegends[i], "F");
  if(!needsQCD) leg->AddEntry((TObject*)0, "", "");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  legDrawSignal = new TLegend(0.45, 0.6, 0.85, 0.85, NULL, "brNDC");
  legDrawSignal->SetNColumns(2);
  legDrawSignal->AddEntry(data, "Data", "LP");
  legDrawSignal->AddEntry((TObject*)0, "", "");
  legDrawSignal->AddEntry(errors_sys, "Stat. #oplus Syst. Errors", "F");
  legDrawSignal->AddEntry((TObject*)0, "", "");
  if(needsQCD) legDrawSignal->AddEntry(bkg, "QCD", "F");

  legDrawSignal->AddEntry(mc[0], layerLegends[0], "F");
  for(unsigned int i = 1; i < mc.size(); i++) legDrawSignal->AddEntry(mc[i], layerLegends[i], "F");
  if(!needsQCD) legDrawSignal->AddEntry((TObject*)0, "", "");
  legDrawSignal->SetFillColor(0);
  legDrawSignal->SetTextSize(0.028);
  legDrawSignal->AddEntry(siga, "GGM (460_175)", "L");
  legDrawSignal->AddEntry(sigb, "GGM (560_325)", "L");

  ratioLeg = new TLegend(0.78, 0.7, 0.88, 0.95, NULL, "brNDC");
  ratioLeg->AddEntry(ratio_stat, "Stat.", "F");
  ratioLeg->AddEntry(ratio_sys, "Stat. #oplus Syst.", "F");
  ratioLeg->SetFillColor(0);
  ratioLeg->SetTextSize(0.032);

  reqText = new TPaveText(0.45, 0.47, 0.85, 0.57, "NDC");
  reqText->SetFillColor(0);
  reqText->SetFillStyle(0);
  reqText->SetLineColor(0);
  reqText->AddText(channelLabel.ReplaceAll("XYZ", crNames[controlRegion]));

  lumiHeader = new TPaveText(0.1, 0.901, 0.9, 0.94, "NDC");
  lumiHeader->SetFillColor(0);
  lumiHeader->SetFillStyle(0);
  lumiHeader->SetLineColor(0);
  lumiHeader->AddText("CMS Preliminary 2015     #sqrt{s} = 8 TeV     #intL = 19.7 fb^{-1}");

}

void PlotMaker::SetStyles(unsigned int n) {

  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.5);
  
  errors_stat->SetFillColor(kOrange+10);
  errors_stat->SetFillStyle(3154);
  errors_stat->SetMarkerSize(0);
  
  errors_sys->SetFillColor(kOrange+10);
  errors_sys->SetFillStyle(3154);
  errors_sys->SetMarkerSize(0);

  ratio_stat->SetFillStyle(1001);
  ratio_stat->SetFillColor(kGray+1);
  ratio_stat->SetLineColor(kGray+1);
  ratio_stat->SetMarkerColor(kGray+1);

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);
  
  if(needsQCD) bkg->SetFillColor(kSpring-6);
  else bkg->SetFillColor(layerColors[0]);
  bkg->SetMarkerSize(0);
  bkg->SetLineColor(1);

  for(unsigned int i = 0; i < mc.size(); i++) {
    mc[i]->SetFillColor(layerColors[i]);
    mc[i]->SetMarkerSize(0);
    mc[i]->SetLineColor(1);
  }

  bkg->SetTitle(variables[n]);

  bkg->GetXaxis()->SetTitle(xaxisTitles[n]);
  ratio->GetXaxis()->SetTitle(xaxisTitles[n]);
  ratio->GetXaxis()->SetLabelFont(63);
  ratio->GetXaxis()->SetLabelSize(48);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleOffset(0.6);

  if(xMaximums[n] > xMinimums[n]) {
    bkg->GetXaxis()->SetRangeUser(xMinimums[n], xMaximums[n]);
    ratio->GetXaxis()->SetRangeUser(xMinimums[n], xMaximums[n]);
  }

  bkg->GetYaxis()->SetTitle(yaxisTitles[n]);
  ratio->GetYaxis()->SetTitle("Data / Background");
  ratio->GetYaxis()->SetLabelFont(63);
  ratio->GetYaxis()->SetLabelSize(48);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetYaxis()->SetNdivisions(508);

  bkg->GetYaxis()->SetRangeUser(yMinimums[n], yMaximums[n]);
  ratio->GetYaxis()->SetRangeUser(ratioMinimums[n], ratioMaximums[n]);

  oneLine = new TLine(xMinimums[n], 1, xMaximums[n], 1);
  oneLine->SetLineStyle(2);

  siga->SetLineColor(kMagenta);
  siga->SetLineWidth(3);
  
  sigb->SetLineColor(kBlue);
  sigb->SetLineWidth(3);

}

void PlotMaker::ScaleByFit(unsigned int n, vector<TH1D*>& h) {

  if(n >= fitScales.size() ||
     n >= fitScaleErrors.size() ||
     n >= h.size()) return;

  Float_t sf = fitScales[n];
  Float_t sfError = fitScaleErrors[n];

  Float_t olderror, newerror;
  Float_t oldcontent;

  for(Int_t i = 0; i < h[n]->GetNbinsX(); i++) {
    olderror = h[n]->GetBinError(i+1);
    oldcontent = h[n]->GetBinContent(i+1);

    if(olderror == 0) continue;

    newerror = sfError*sfError / sf / sf;
    newerror += olderror*olderror / oldcontent / oldcontent;
    newerror = sf * oldcontent * sqrt(newerror);
    
    h[n]->SetBinContent(i+1, sf * oldcontent);
    h[n]->SetBinError(i+1, newerror);
  }

}

void PlotMaker::ScaleQCD() {

  Float_t olderror, newerror;
  Float_t oldcontent;

  for(Int_t i = 0; i < qcd->GetNbinsX(); i++) {
    olderror = qcd->GetBinError(i+1);
    oldcontent = qcd->GetBinContent(i+1);

    if(olderror == 0) continue;
    
    newerror = qcdScaleError*qcdScaleError / qcdScale / qcdScale;
    newerror += olderror*olderror / oldcontent / oldcontent;
    newerror = qcdScale * oldcontent * sqrt(newerror);

    qcd->SetBinContent(i+1, qcdScale * oldcontent);
    qcd->SetBinError(i+1, newerror);
  }

  qcd_defUp->Scale(qcdScale_defUp);
  qcd_defDown->Scale(qcdScale_defDown);

}

void PlotMaker::CalculateQCDNormalization() {

  unsigned int met_index = 0;
  bool foundMET = false;

  for(unsigned int i = 0; i < variables.size(); i++) {
    if(variables[i] == "pfMET") {
      met_index = i;
      foundMET = true;
      break;
    }
  }

  if(!foundMET) {
    cout << endl << endl << "Can't normalize QCD in pfMET if you don't plot pfMET!" << endl << endl;
    return;
  }

  GetHistograms(met_index);

  const int endBin = data->GetXaxis()->FindBin(20) - 1;

  double n_data = data->Integral(1, endBin);
  double n_qcd = qcd->Integral(1, endBin);

  double n_qcd_defUp = qcd_defUp->Integral(1, endBin);
  double n_qcd_defDown = qcd_defDown->Integral(1, endBin);

  if(n_qcd < 1) {
    qcdScale = 0.0;
    qcdScale_defUp = 0.0;
    qcdScale_defDown = 0.0;
    qcdScaleError = 0.0;
    return;
  }

  double n_mc = 0;
  for(unsigned int i = 0; i < mc.size(); i++) n_mc += mc[i]->Integral(1, endBin);

  double sigma_data = 0;
  double sigma_qcd = 0;
  double sigma_mc = 0;
  
  for(int i = 0; i < endBin; i++) {
    sigma_data += data->GetBinError(i+1) * data->GetBinError(i+1);
    sigma_qcd += qcd->GetBinError(i+1) * qcd->GetBinError(i+1);

    for(unsigned int j = 0; j < mc.size(); j++) sigma_mc +=mc[j]->GetBinError(i+1) * mc[j]->GetBinError(i+1);
  }

  sigma_data = sqrt(sigma_data);
  sigma_qcd = sqrt(sigma_qcd);
  sigma_mc = sqrt(sigma_mc);

  qcdScale = (n_data - n_mc) / n_qcd;
  if(qcdScale < 0) {
    qcdScale = 0.0;
    qcdScale_defUp = 0.0;
    qcdScale_defDown = 0.0;
    qcdScaleError = 0.0;
    return;
  }

  qcdScale_defUp = (n_data - n_mc) / n_qcd_defUp;
  qcdScale_defDown = (n_data - n_mc) / n_qcd_defDown;

  qcdScaleError = sigma_data*sigma_data + sigma_mc*sigma_mc;
  qcdScaleError /= (n_data - n_mc) * (n_data - n_mc);
  qcdScaleError += sigma_qcd*sigma_qcd / n_qcd / n_qcd;
  qcdScaleError = qcdScale * sqrt(qcdScaleError);

  cout << endl << "CalculateQCDNormalization(): " << qcdScale << " +- " << qcdScaleError << endl << endl;

  return;

}

void PlotMaker::CreatePlot(unsigned int n) {

  if(n > 0) GetHistograms(n);
  ScaleQCD();
  StackHistograms(n);
  if(n == 0) {
    METDifference();
    CreateSignalOutputs();
  }
  CalculateRatio(n);
  if(n == 0)  MakeLegends();

  if(divideByBinWidth[n]) DivideWidth();

  SetStyles(n);

  padhi->cd();

  bkg->Draw("hist");
  for(unsigned int i = 0; i < mc.size(); i++) {
    if(!needsQCD && i == 0) continue;
    mc[i]->Draw("same hist");
  }
  //errors_stat->Draw("same e2");
  errors_sys->Draw("same e2");
  data->Draw("same e1");
  bkg->Draw("same axis");

  if(doDrawSignal[n]) {
    siga->Draw("same hist");
    sigb->Draw("same hist");
  }

  lumiHeader->Draw("same");
  if(doDrawLegend[n]) {
    if(doDrawSignal[n]) legDrawSignal->Draw("same");
    else leg->Draw("same");
  }
  if(doDrawPrelim[n] && doDrawLegend[n]) reqText->Draw("same");

  padlo->cd();

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio_stat->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  ratioLeg->Draw("same");

  oneLine->Draw();  

  can->SaveAs(variables[n]+"_"+channel+"_"+crNames[controlRegion]+".pdf");

}

void PlotMaker::METDifference() {

  TH1D * h = (TH1D*)data->Clone(channel+"_"+crNames[controlRegion]);
  h->Add(bkg, -1.);

  TFile * output = new TFile("met_differences.root", "UPDATE");

  h->Write();
  output->Close();
}

void PlotMaker::SaveLimitOutputs() {

  // durp
  
}

// durp
void PlotMaker::CreateSignalOutputs() {

  Double_t xbins[31];
  xbins[0] = 0;
  xbins[1] = 55;
  for(int i = 1; i < 29; i++) xbins[i+1] = (mst[i] + mst[i-1])/2.;
  xbins[30] = 6510;

  Double_t ybins[33];
  ybins[0] = 0;
  ybins[1] = 12.5;
  for(int i = 1; i < 31; i++) ybins[i+1] = (mBino[i] + mBino[i-1])/2.;
  ybins[32] = 2175;

  char code[100];
  int index1, index2;

  TFile * f_xsec = new TFile("../../data/stop-bino_xsecs.root", "READ");
  TH2D * h_xsec = (TH2D*)f_xsec->Get("real_xsec");
  TH2D * h_xsec_errors = (TH2D*)f_xsec->Get("real_errors");

  TH2D * h_acc = new TH2D("acc_"+req, "acc_"+req, 30, xbins, 32, ybins);
  TH2D * h_contamination = new TH2D("contamination_"+req, "contamination_"+req, 30, xbins, 32, ybins);

  TFile * fSignalOut = new TFile("limitInputs.root", "UPDATE");
  if(req.Contains("ele")) {
    fSignalOut->mkdir("ele");
    fSignalOut->cd("ele");
  }
  else {
    fSignalOut->mkdir("muon");
    fSignalOut->cd("muon");
  }

  for(int imass = 0; imass < 899; imass++) {

    index1 = mst[int(imass)/31];
    index2 = mBino[int(imass)%31];

    if(index1 < index2) continue;

    sprintf(code, "_mst_%d_m1_%d", index1, index2);
    TString code_t = code;

    TFile * f = new TFile("../../acceptance_v2/signal_contamination"+code_t+".root", "READ");
    if(f->IsZombie()) {
      f->Close();
      continue;
    }
    
    TTree * tree = (TTree*)f->Get(channel+"_noSigmaIetaIetaTree");
    TTree * tree_JECup = (TTree*)f->Get(channel+"_noSigmaIetaIetaTree_JECup");
    TTree * tree_JECdown = (TTree*)f->Get(channel+"_noSigmaIetaIetaTree_JECdown");

    TTree * tree_contam;
    if(channel.Contains("ele")) tree_contam = (TTree*)f->Get("ele_jjj_veto_eQCDnoSigmaIetaIetaTree");
    else if(channel.Contains("muon")) tree_contam = (TTree*)f->Get("muon_jjj_veto_muQCDnoSigmaIetaIetaTree");

    if(!tree || !tree_JECup || !tree_JECdown || !tree_contam) {
      f->Close();
      continue;
    }

    Float_t met, ngamma, nfake, nphotons;
    Float_t puWeight, btagWeight;
    Float_t puWeightErr, btagWeightErr;
    Float_t puWeightUp, puWeightDown, btagWeightUp, btagWeightDown;
    Float_t overlaps_ttA;
    Float_t topPtReweighting;

    Float_t lepton_pt, lepton_eta;
    Float_t lead_photon_et, lead_photon_eta;
    Float_t trail_photon_et, trail_photon_eta;

    tree->SetBranchAddress("pfMET", &met);
    tree_JECup->SetBranchAddress("pfMET", &met);
    tree_JECdown->SetBranchAddress("pfMET", &met);
    tree_contam->SetBranchAddress("pfMET", &met);

    tree->SetBranchAddress("Nphotons", &nphotons);
    tree_JECup->SetBranchAddress("Nphotons", &nphotons);
    tree_JECdown->SetBranchAddress("Nphotons", &nphotons);
    tree_contam->SetBranchAddress("Nphotons", &nphotons);

    tree->SetBranchAddress("Ngamma", &ngamma);
    tree_JECup->SetBranchAddress("Ngamma", &ngamma);
    tree_JECdown->SetBranchAddress("Ngamma", &ngamma);
    tree_contam->SetBranchAddress("Ngamma", &ngamma);

    tree->SetBranchAddress("Nfake", &nfake);
    tree_JECup->SetBranchAddress("Nfake", &nfake);
    tree_JECdown->SetBranchAddress("Nfake", &nfake);
    tree_contam->SetBranchAddress("Nfake", &nfake);

    tree->SetBranchAddress("leadPhotonEt", &lead_photon_et);
    tree_JECup->SetBranchAddress("leadPhotonEt", &lead_photon_et);
    tree_JECdown->SetBranchAddress("leadPhotonEt", &lead_photon_et);
    tree_contam->SetBranchAddress("leadPhotonEt", &lead_photon_et);

    tree->SetBranchAddress("leadPhotonEta", &lead_photon_eta);
    tree_JECup->SetBranchAddress("leadPhotonEta", &lead_photon_eta);
    tree_JECdown->SetBranchAddress("leadPhotonEta", &lead_photon_eta);
    tree_contam->SetBranchAddress("leadPhotonEta", &lead_photon_eta);

    tree->SetBranchAddress("trailPhotonEta", &trail_photon_eta);
    tree_JECup->SetBranchAddress("trailPhotonEta", &trail_photon_eta);
    tree_JECdown->SetBranchAddress("trailPhotonEta", &trail_photon_eta);
    tree_contam->SetBranchAddress("trailPhotonEta", &trail_photon_eta);

    tree->SetBranchAddress("trailPhotonEt", &trail_photon_et);
    tree_JECup->SetBranchAddress("trailPhotonEt", &trail_photon_et);
    tree_JECdown->SetBranchAddress("trailPhotonEt", &trail_photon_et);
    tree_contam->SetBranchAddress("trailPhotonEt", &trail_photon_et);

    if(channel.Contains("ele")) {
      tree->SetBranchAddress("ele_pt", &lepton_pt);
      tree_JECup->SetBranchAddress("ele_pt", &lepton_pt);
      tree_JECdown->SetBranchAddress("ele_pt", &lepton_pt);
      tree_contam->SetBranchAddress("ele_pt", &lepton_pt);

      tree->SetBranchAddress("ele_eta", &lepton_eta);
      tree_JECup->SetBranchAddress("ele_eta", &lepton_eta);
      tree_JECdown->SetBranchAddress("ele_eta", &lepton_eta);
      tree_contam->SetBranchAddress("ele_eta", &lepton_eta);
    }
    else if(channel.Contains("muon")) {
      tree->SetBranchAddress("muon_pt", &lepton_pt);
      tree_JECup->SetBranchAddress("muon_pt", &lepton_pt);
      tree_JECdown->SetBranchAddress("muon_pt", &lepton_pt);
      tree_contam->SetBranchAddress("muon_pt", &lepton_pt);

      tree->SetBranchAddress("muon_eta", &lepton_eta);
      tree_JECup->SetBranchAddress("muon_eta", &lepton_eta);
      tree_JECdown->SetBranchAddress("muon_eta", &lepton_eta);
      tree_contam->SetBranchAddress("muon_eta", &lepton_eta);
    }

    tree->SetBranchAddress("pileupWeight", &puWeight);
    tree->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree->SetBranchAddress("btagWeight", &btagWeight);
    tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    
    tree_JECup->SetBranchAddress("pileupWeight", &puWeight);
    tree_JECup->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree_JECup->SetBranchAddress("btagWeight", &btagWeight);
    tree_JECup->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree_JECup->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree_JECup->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree_JECup->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree_JECup->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree_JECup->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    
    tree_JECdown->SetBranchAddress("pileupWeight", &puWeight);
    tree_JECdown->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree_JECdown->SetBranchAddress("btagWeight", &btagWeight);
    tree_JECdown->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree_JECdown->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree_JECdown->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree_JECdown->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree_JECdown->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree_JECdown->SetBranchAddress("TopPtReweighting", &topPtReweighting);

    tree_contam->SetBranchAddress("pileupWeight", &puWeight);
    tree_contam->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree_contam->SetBranchAddress("btagWeight", &btagWeight);
    tree_contam->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree_contam->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree_contam->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree_contam->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree_contam->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree_contam->SetBranchAddress("TopPtReweighting", &topPtReweighting);

    TH1D * h = new TH1D("signal"+code_t, "signal"+code_t, nMetBins, xbins_met); h->Sumw2();

    TH1D * h_btagWeightUp = new TH1D("signal"+code_t+"_btagWeightUp", "signal"+code_t+"_btagWeightUp", nMetBins, xbins_met); h_btagWeightUp->Sumw2();
    TH1D * h_btagWeightDown = new TH1D("signal"+code_t+"_btagWeightDown", "signal"+code_t+"_btagWeightDown", nMetBins, xbins_met); h_btagWeightDown->Sumw2();

    TH1D * h_puWeightUp = new TH1D("signal"+code_t+"_puWeightUp", "signal"+code_t+"_puWeightUp", nMetBins, xbins_met); h_puWeightUp->Sumw2();
    TH1D * h_puWeightDown = new TH1D("signal"+code_t+"_puWeightDown", "signal"+code_t+"_puWeightDown", nMetBins, xbins_met); h_puWeightDown->Sumw2();

    TH1D * h_topPtUp = new TH1D("signal"+code_t+"_topPtUp", "signal"+code_t+"_topPtUp", nMetBins, xbins_met); h_topPtUp->Sumw2();
    TH1D * h_topPtDown = new TH1D("signal"+code_t+"_topPtDown", "signal"+code_t+"_topPtDown", nMetBins, xbins_met); h_topPtDown->Sumw2();

    TH1D * h_JECup = new TH1D("signal"+code_t+"_JECUp", "signal"+code_t+"_JECUp", nMetBins, xbins_met); h_JECup->Sumw2();
    TH1D * h_JECdown = new TH1D("signal"+code_t+"_JECDown", "signal"+code_t+"_JECDown", nMetBins, xbins_met); h_JECdown->Sumw2();

    TH1D * h_leptonSFup = new TH1D("signal"+code_t+"_leptonSFUp", "signal"+code_t+"_leptonSFUp", nMetBins, xbins_met); h_leptonSFup->Sumw2();
    TH1D * h_leptonSFdown = new TH1D("signal"+code_t+"_leptonSFDown", "signal"+code_t+"_leptonSFDown", nMetBins, xbins_met); h_leptonSFdown->Sumw2();

    TH1D * h_photonSFup = new TH1D("signal"+code_t+"_photonSFUp", "signal"+code_t+"_photonSFUp", nMetBins, xbins_met); h_photonSFup->Sumw2();
    TH1D * h_photonSFdown = new TH1D("signal"+code_t+"_photonSFDown", "signal"+code_t+"_photonSFDown", nMetBins, xbins_met); h_photonSFdown->Sumw2();

    for(int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);

      if(!inControlRegion(ngamma, nfake)) continue;

      if(!channel.Contains("b")) {
	btagWeight = 1.;
	btagWeightErr = 0.;
	btagWeightUp = 1.;
	btagWeightDown = 1.;
      }
      
      if(btagWeightUp < 0) btagWeightUp = 0.;
      if(btagWeightDown < 0) btagWeightDown = 0;
      
      if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
      
      if(btagWeight != btagWeight) continue;

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;

      if(totalWeight < 0) continue;

      Float_t olderror = h->GetBinError(h->FindBin(met));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h->Fill(met, totalWeight);
      h->SetBinError(h->FindBin(met), newerror);

      totalWeight = puWeight * btagWeightUp * leptonSF * photonSF * topPtReweighting;
      olderror = h_btagWeightUp->GetBinError(h_btagWeightUp->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_btagWeightUp->Fill(met, totalWeight);
      h_btagWeightUp->SetBinError(h_btagWeightUp->FindBin(met), newerror);

      totalWeight = puWeight * btagWeightDown * leptonSF * photonSF * topPtReweighting;
      olderror = h_btagWeightDown->GetBinError(h_btagWeightDown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_btagWeightDown->Fill(met, totalWeight);
      h_btagWeightDown->SetBinError(h_btagWeightDown->FindBin(met), newerror);

      totalWeight = puWeightUp * btagWeight * leptonSF * photonSF * topPtReweighting;
      olderror = h_puWeightUp->GetBinError(h_puWeightUp->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_puWeightUp->Fill(met, totalWeight);
      h_puWeightUp->SetBinError(h_puWeightUp->FindBin(met), newerror);

      totalWeight = puWeightDown * btagWeight * leptonSF * photonSF * topPtReweighting;
      olderror = h_puWeightDown->GetBinError(h_puWeightDown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_puWeightDown->Fill(met, totalWeight);
      h_puWeightDown->SetBinError(h_puWeightDown->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSFup * photonSF * topPtReweighting;
      olderror = h_leptonSFup->GetBinError(h_leptonSFup->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_leptonSFup->Fill(met, totalWeight);
      h_leptonSFup->SetBinError(h_leptonSFup->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSFdown * photonSF * topPtReweighting;
      olderror = h_leptonSFdown->GetBinError(h_leptonSFdown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_leptonSFdown->Fill(met, totalWeight);
      h_leptonSFdown->SetBinError(h_leptonSFdown->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFup * topPtReweighting;
      olderror = h_photonSFup->GetBinError(h_photonSFup->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_photonSFup->Fill(met, totalWeight);
      h_photonSFup->SetBinError(h_photonSFup->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFdown * topPtReweighting;
      olderror = h_photonSFdown->GetBinError(h_photonSFdown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_photonSFdown->Fill(met, totalWeight);
      h_photonSFdown->SetBinError(h_photonSFdown->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting * topPtReweighting;
      olderror = h_topPtUp->GetBinError(h_topPtUp->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_topPtUp->Fill(met, totalWeight);
      h_topPtUp->SetBinError(h_topPtUp->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      olderror = h_topPtDown->GetBinError(h_topPtDown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_topPtDown->Fill(met, totalWeight);
      h_topPtDown->SetBinError(h_topPtDown->FindBin(met), newerror);
    }

    for(int i = 0; i < tree_JECup->GetEntries(); i++) {
      tree_JECup->GetEntry(i);

      if(!inControlRegion(ngamma, nfake)) continue;

      if(!channel.Contains("b")) {
	btagWeight = 1.;
	btagWeightErr = 0.;
	btagWeightUp = 1.;
	btagWeightDown = 1.;
      }
      
      if(btagWeightUp < 0) btagWeightUp = 0.;
      if(btagWeightDown < 0) btagWeightDown = 0;
      
      if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
      
      if(btagWeight != btagWeight) continue;

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h->GetBinError(h->FindBin(met));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_JECup->Fill(met, totalWeight);
      h_JECup->SetBinError(h->FindBin(met), newerror);

    }

    for(int i = 0; i < tree_JECdown->GetEntries(); i++) {
      tree_JECdown->GetEntry(i);

      if(!inControlRegion(ngamma, nfake)) continue;

      if(!channel.Contains("b")) {
	btagWeight = 1.;
	btagWeightErr = 0.;
	btagWeightUp = 1.;
	btagWeightDown = 1.;
      }
      
      if(btagWeightUp < 0) btagWeightUp = 0.;
      if(btagWeightDown < 0) btagWeightDown = 0;
      
      if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
      
      if(btagWeight != btagWeight) continue;

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h->GetBinError(h->FindBin(met));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_JECdown->Fill(met, totalWeight);
      h_JECdown->SetBinError(h->FindBin(met), newerror);

    }

    double contamination = 0;

    for(int i = 0; i < tree_contam->GetEntries(); i++) {
      tree_contam->GetEntry(i);

      if(!inControlRegion(ngamma, nfake)) continue;

      if(!channel.Contains("b")) btagWeight = 1.;
      
      if(btagWeight != btagWeight) continue;

      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, chan, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;

      contamination += totalWeight;
    }

    h_acc->SetBinContent(h_acc->FindBin(index1, index2), h->Integral() / (0.438/3.) / 15000.);
    h_contamination->SetBinContent(h_contamination->FindBin(index1, index2), contamination / h->Integral());

    double xsec = h_xsec->GetBinContent(h_xsec->FindBin(index1, index2));
    
    h->Scale(xsec * 19712. / 15000.);
    h_btagWeightUp->Scale(xsec * 19712. / 15000.);
    h_btagWeightDown->Scale(xsec * 19712. / 15000.);
    h_puWeightUp->Scale(xsec * 19712. / 15000.);
    h_puWeightDown->Scale(xsec * 19712. / 15000.);
    h_topPtUp->Scale(xsec * 19712. / 15000.);
    h_topPtDown->Scale(xsec * 19712. / 15000.);
    h_JECup->Scale(xsec * 19712. / 15000.);
    h_JECdown->Scale(xsec * 19712. / 15000.);
    h_leptonSFup->Scale(xsec * 19712. / 15000.);
    h_leptonSFdown->Scale(xsec * 19712. / 15000.);
    h_photonSFup->Scale(xsec * 19712. / 15000.);
    h_photonSFdown->Scale(xsec * 19712. / 15000.);

    fSignalOut->cd();

    if(channel.Contains("ele")) {
      fSignalOut->cd("ele");
    }
    else {
      fSignalOut->cd("muon");
    }

    h->Write("signal"+code_t);

    for(int j = 0; j < h->GetNbinsX(); j++) {
      TH1D * h_flux_up = (TH1D*)h->Clone("clone_signal_"+code_t+"_flux_up");
      TH1D * h_flux_down = (TH1D*)h->Clone("clone_signal"+code_t+"_flux_down");
      
      Double_t centralValue = h->GetBinContent(j+1);
      Double_t statError = h->GetBinError(j+1);
      
      if(statError > 0.) h_flux_up->SetBinContent(j+1, centralValue + statError);
      if(centralValue > statError && statError > 0.) h_flux_down->SetBinContent(j+1, centralValue - statError);
      
      h_flux_up->Write("signal"+code_t+"_signal_stat_bin"+Form("%d", j+1)+"Up");
      h_flux_down->Write("signal"+code_t+"_signal_stat_bin"+Form("%d", j+1)+"Down");
    }
    
    h_btagWeightUp->Write("signal"+code_t+"_btagWeightUp");
    h_btagWeightDown->Write("signal"+code_t+"_btagWeightDown");
    h_puWeightUp->Write("signal"+code_t+"_puWeightUp");
    h_puWeightDown->Write("signal"+code_t+"_puWeightDown");
    h_topPtUp->Write("signal"+code_t+"_topPtUp");
    h_topPtDown->Write("signal"+code_t+"_topPtDown");
    h_JECup->Write("signal"+code_t+"_JECUp");
    h_JECdown->Write("signal"+code_t+"_JECDown");
    h_leptonSFup->Write("signal"+code_t+"_leptonSFUp");
    h_leptonSFdown->Write("signal"+code_t+"_leptonSFDown");
    h_photonSFup->Write("signal"+code_t+"_photonSFUp");
    h_photonSFdown->Write("signal"+code_t+"_photonSFDown");

    f->Close();

  }

  // draw acc and etc
  TCanvas * can_acc = new TCanvas("canvas_acc", "Plot", 10, 10, 2000, 2000);
  fillPotHoles(h_acc);
  h_acc->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_acc->GetXaxis()->SetRangeUser(0, 1600);
  h_acc->GetXaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_acc->GetYaxis()->SetTitleOffset(1.3);
  h_acc->GetYaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetRangeUser(0, 1600);
  h_acc->GetZaxis()->SetLabelSize(0.02);
  h_acc->Draw("colz");
  can_acc->SaveAs("acceptance_"+channel+".pdf");
  
  h_contamination->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_contamination->GetXaxis()->SetRangeUser(0, 1600);
  h_contamination->GetXaxis()->SetLabelSize(0.03);
  h_contamination->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_contamination->GetYaxis()->SetTitleOffset(1.3);
  h_contamination->GetYaxis()->SetLabelSize(0.03);
  h_contamination->GetYaxis()->SetRangeUser(0, 1600);
  h_contamination->GetZaxis()->SetLabelSize(0.02);
  h_contamination->Draw("colz");
  can_acc->SaveAs("contamination_"+channel+".pdf");
  
  delete can_acc;

  fSignalOut->Close();

  f_xsec->Close();
  
}

