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
#include <fstream>

#include <map>
#include <vector>
#include <stdio.h>

#include "rootRoutines.h"
#include "Binning.h"

using namespace std;

const int nChannels = 4;

TString channels[nChannels] = {"ele_bjj", "muon_bjj",
			       "ele_jjj", "muon_jjj"};

TString channelLabels[nChannels] = {"XYZ ele", "XYZ muon",
				    "XYZ ele (no b-tag)", "XYZ muon (no b-tag)"};

enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};
TString crNames[kNumControlRegions] = {"SR1", "SR2", "CR1", "CR2", "CR2a", "CR0", "SigmaPlot", "Any"};
TString crLabels[kNumControlRegions] = {"SR1", "SR2", "CR1", "CR2", "CR2a", "CR0", "SigmaPlot", "Pre-selection"};

enum pdfCorrelatesWith {kGG, kQQ, kQG, kNpdfCorrelations};
enum scaleCorrelatesWith {kTTbar, kV, kVV, kNscaleCorrelations};

class PlotMaker : public TObject {

  ClassDef(PlotMaker, 1);

 public:
  PlotMaker(int chanNo, int cr, bool useQCD, TString metType_, vector<Float_t> sf_qcd_, vector<Float_t> sfError_qcd_);
  virtual ~PlotMaker();

  void BookMCLayer(vector<TString> newNames, int color, TString limitName, TString legendEntry, int pdfCorr, int scaleCorr, 
		   vector<Float_t> scale, vector<Float_t> scaleErr) { 
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

    limitNames.push_back(limitName);

    pdfCorrelations.push_back(pdfCorr);
    scaleCorrelations.push_back(scaleCorr);

    fitScales.push_back(scale);
    fitScaleErrors.push_back(scaleErr);
  };

  void BookMCLayer(vector<TString> newNames, int color, TString limitName, TString legendEntry, int pdfCorr, int scaleCorr) {
    vector<Float_t> emptyVector;
    BookMCLayer(newNames, color, limitName, legendEntry, pdfCorr, scaleCorr, emptyVector, emptyVector);
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

  void SetDoRebinMET(bool v) { doRebinMET = v; };
  
  void RebinMET() {

    if(variables.size() == 0) return;

    data = (TH1D*)data->Rebin(nMetBins_2g, "data_reb", xbins_met_2g);

    siga = (TH1D*)siga->Rebin(nMetBins_2g, "siga_reb", xbins_met_2g);
    sigb = (TH1D*)sigb->Rebin(nMetBins_2g, "sigb_reb", xbins_met_2g);

    qcd = (TH1D*)qcd->Rebin(nMetBins_2g, "qcd_reb", xbins_met_2g);
    qcd_defUp = (TH1D*)qcd_defUp->Rebin(nMetBins_2g, "qcd_defUp_reb", xbins_met_2g);
    qcd_defDown = (TH1D*)qcd_defDown->Rebin(nMetBins_2g, "qcd_defDown_reb", xbins_met_2g);
  
    for(unsigned int i = 0; i < mc.size(); i++) {
      mc[i] = (TH1D*)mc[i]->Rebin(nMetBins_2g, "mc_reb", xbins_met_2g);

      mc_btagWeightUp[i] = (TH1D*)mc_btagWeightUp[i]->Rebin(nMetBins_2g, "mc_btagWeightUp_reb", xbins_met_2g);
      mc_btagWeightDown[i] = (TH1D*)mc_btagWeightDown[i]->Rebin(nMetBins_2g, "mc_btagWeightDown_reb", xbins_met_2g);

      mc_puWeightUp[i] = (TH1D*)mc_puWeightUp[i]->Rebin(nMetBins_2g, "mc_puWeightUp_reb", xbins_met_2g);
      mc_puWeightDown[i] = (TH1D*)mc_puWeightDown[i]->Rebin(nMetBins_2g, "mc_puWeightDown_reb", xbins_met_2g);

      mc_scaleUp[i] = (TH1D*)mc_scaleUp[i]->Rebin(nMetBins_2g, "mc_scaleUp_reb", xbins_met_2g);
      mc_scaleDown[i] = (TH1D*)mc_scaleDown[i]->Rebin(nMetBins_2g, "mc_scaleDown_reb", xbins_met_2g);

      mc_pdfUp[i] = (TH1D*)mc_pdfUp[i]->Rebin(nMetBins_2g, "mc_pdfUp_reb", xbins_met_2g);
      mc_pdfDown[i] = (TH1D*)mc_pdfDown[i]->Rebin(nMetBins_2g, "mc_pdfDown_reb", xbins_met_2g);

      mc_topPtUp[i] = (TH1D*)mc_topPtUp[i]->Rebin(nMetBins_2g, "mc_topPtUp_reb", xbins_met_2g);
      mc_topPtDown[i] = (TH1D*)mc_topPtDown[i]->Rebin(nMetBins_2g, "mc_topPtDown_reb", xbins_met_2g);

      mc_JECUp[i] = (TH1D*)mc_JECUp[i]->Rebin(nMetBins_2g, "mc_JECUp_reb", xbins_met_2g);
      mc_JECDown[i] = (TH1D*)mc_JECDown[i]->Rebin(nMetBins_2g, "mc_JECDown_reb", xbins_met_2g);

      mc_leptonSFUp[i] = (TH1D*)mc_leptonSFUp[i]->Rebin(nMetBins_2g, "mc_leptonSFUp_reb", xbins_met_2g);
      mc_leptonSFDown[i] = (TH1D*)mc_leptonSFDown[i]->Rebin(nMetBins_2g, "mc_leptonSFDown_reb", xbins_met_2g);

      mc_photonSFUp[i] = (TH1D*)mc_photonSFUp[i]->Rebin(nMetBins_2g, "mc_photonSFUp_reb", xbins_met_2g);
      mc_photonSFDown[i] = (TH1D*)mc_photonSFDown[i]->Rebin(nMetBins_2g, "mc_photonSFDown_reb", xbins_met_2g);
    }

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

  void ScaleByFit(unsigned int n, vector<TH1D*>& h, float sf, float sfError);

  void CalculateQCDNormalization();
  void ScaleQCD();

  void DetermineAxisRanges(unsigned int n);
  void DetermineLegendRanges(unsigned int n);

  void CreatePlot(unsigned int n);
  void CreatePlots() {
    MakeCanvas();
    CalculateQCDNormalization();
    for(unsigned int i = 0; i < variables.size(); i++) CreatePlot(i);
  };
  
  void METDifference();
  void SaveLimitOutputs();
  void SaveLimitOutputs_KillZeroBin();

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

  vector<TString> limitNames;

  vector<int> pdfCorrelations;
  vector<int> scaleCorrelations;

  vector< vector<Float_t> > fitScales;
  vector< vector<Float_t> > fitScaleErrors;

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
  TString metType;

  vector<Float_t> sf_qcd;
  vector<Float_t> sfError_qcd;

  bool doRebinMET;

  TLegend * leg;
  TLegend * legDrawSignal;
  TLegend * ratioLeg;

};

PlotMaker::PlotMaker(int chanNo, int cr, bool useQCD, TString metType_, vector<Float_t> sf_qcd_, vector<Float_t> sfError_qcd_) {

  channel = channels[chanNo];
  channelLabel = channelLabels[chanNo];
  controlRegion = cr;
  needsQCD = useQCD;
  metType = metType_;

  doRebinMET = false;

  sf_qcd.clear();
  sfError_qcd.clear();

  for(unsigned int i = 0; i < sf_qcd_.size(); i++) sf_qcd.push_back(sf_qcd_[i]);
  for(unsigned int i = 0; i < sfError_qcd_.size(); i++) sfError_qcd.push_back(sfError_qcd_[i]);

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

  limitNames.clear();

  pdfCorrelations.clear();
  scaleCorrelations.clear();

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

  sf_qcd.clear();
  sfError_qcd.clear();

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

  limitNames.clear();

  pdfCorrelations.clear();
  scaleCorrelations.clear();

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
  qcd_defUp = (TH1D*)input->Get(variables[n]+"_qcd_relIso_10_"+channel);
  qcd_defDown = (TH1D*)input->Get(variables[n]+"_qcd_relIso_m10_"+channel);
  
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

      mc_btagWeightUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_btagWeightUp"));
      mc_btagWeightDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_btagWeightDown"));
      
      mc_puWeightUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_puWeightUp"));
      mc_puWeightDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_puWeightDown"));
      
      mc_scaleUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_scaleUp"));
      mc_scaleDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_scaleDown"));
      
      mc_pdfUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_pdfUp"));
      mc_pdfDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_pdfDown"));
      
      mc_topPtUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_topPtUp"));
      mc_topPtDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_topPtDown"));
      
      mc_JECUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_JECUp"));
      mc_JECDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_JECDown"));
      
      mc_leptonSFUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_leptonSFUp"));
      mc_leptonSFDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_leptonSFDown"));
      
      mc_photonSFUp[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_photonSFUp"));
      mc_photonSFDown[i]->Add((TH1D*)input->Get(variables[n]+"_"+layerNames[i][j]+"_"+channel+"_photonSFDown"));

    }
  }

  for(unsigned int i = 0; i < mc.size(); i++) {
    if(fitScales[i].size() > 0) {
      ScaleByFit(i, mc, fitScales[i][0], fitScaleErrors[i][0]);
      ScaleByFit(i, mc_btagWeightUp, fitScales[i][1], fitScaleErrors[i][1]);
      ScaleByFit(i, mc_btagWeightDown, fitScales[i][2], fitScaleErrors[i][2]);
      ScaleByFit(i, mc_puWeightUp, fitScales[i][3], fitScaleErrors[i][3]);
      ScaleByFit(i, mc_puWeightDown, fitScales[i][4], fitScaleErrors[i][4]);
      ScaleByFit(i, mc_scaleUp, fitScales[i][0], fitScaleErrors[i][0]);
      ScaleByFit(i, mc_scaleDown, fitScales[i][0], fitScaleErrors[i][0]);
      ScaleByFit(i, mc_pdfUp, fitScales[i][0], fitScaleErrors[i][0]);
      ScaleByFit(i, mc_pdfDown, fitScales[i][0], fitScaleErrors[i][0]);
      ScaleByFit(i, mc_topPtUp, fitScales[i][5], fitScaleErrors[i][5]);
      ScaleByFit(i, mc_topPtDown, fitScales[i][6], fitScaleErrors[i][6]);
      ScaleByFit(i, mc_JECUp, fitScales[i][7], fitScaleErrors[i][7]);
      ScaleByFit(i, mc_JECDown, fitScales[i][8], fitScaleErrors[i][8]);
      ScaleByFit(i, mc_leptonSFUp, fitScales[i][9], fitScaleErrors[i][9]);
      ScaleByFit(i, mc_leptonSFDown, fitScales[i][10], fitScaleErrors[i][10]);
      ScaleByFit(i, mc_photonSFUp, fitScales[i][11], fitScaleErrors[i][11]);
      ScaleByFit(i, mc_photonSFDown, fitScales[i][12], fitScaleErrors[i][12]);
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

  if(n == 0) {
    TFile * fDurp = new TFile("durp.root", "UPDATE");
    TString name = "_"+channel+"_"+crNames[controlRegion];
    bkg->Write("bkg"+name);
    errors_stat->Write("stat"+name);
    errors_sys->Write("sys"+name);
    fDurp->Close();
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
  for(unsigned int i = 1; i < mc.size(); i++) {
    if(layerLegends[i] == layerLegends[i-1]) continue;
    leg->AddEntry(mc[i], layerLegends[i], "F");
  }
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
  for(unsigned int i = 1; i < mc.size(); i++) {
    if(layerLegends[i] == layerLegends[i-1]) continue;
    legDrawSignal->AddEntry(mc[i], layerLegends[i], "F");
  }
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
  reqText->AddText(channelLabel.ReplaceAll("XYZ", crLabels[controlRegion]));

  lumiHeader = new TPaveText(0.1, 0.901, 0.9, 0.94, "NDC");
  lumiHeader->SetFillColor(0);
  lumiHeader->SetFillStyle(0);
  lumiHeader->SetLineColor(0);
  lumiHeader->AddText("CMS Preliminary 2015     #sqrt{s} = 8 TeV     #intL = 19.7 fb^{-1}");
  //aps15 lumiHeader->AddText("Work in Progress     #sqrt{s} = 8 TeV     #intL = 19.7 fb^{-1}");

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

  DetermineAxisRanges(n);
  //DetermineLegendRanges(n);

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

void PlotMaker::ScaleByFit(unsigned int n, vector<TH1D*>& h, float sf, float sfError) {

  if(n >= h.size()) return;

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

    if(oldcontent == 0) continue;
    
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
    if(variables[i] == metType) {
      met_index = i;
      foundMET = true;
      break;
    }
  }

  if(!foundMET) {
    cout << endl << endl << "Can't normalize QCD in " << metType << " if you don't plot " << metType << "!" << endl << endl;
    return;
  }

  GetHistograms(met_index);

  if(sf_qcd.size() > 0) {
    qcdScale = sf_qcd[0];
    qcdScaleError = sfError_qcd[0];
    qcdScale_defUp = sf_qcd[17];
    qcdScale_defDown = sf_qcd[18];
    return;
  }

  const int endBin = data->GetXaxis()->FindBin(20) - 1;

  double n_data = data->Integral(1, endBin);
  double n_qcd = qcd->Integral(1, endBin);

  double n_qcd_defUp = qcd_defUp->Integral(1, endBin);
  double n_qcd_defDown = qcd_defDown->Integral(1, endBin);

  // With the MC subtraction, you can actually have 0 < n_qcd < 1
  if(n_qcd < 1.e-6) {
    qcdScale = 1.e-6;
    qcdScale_defUp = 1.e-6;
    qcdScale_defDown = 1.e-6;
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
    qcdScale = 1.e-6;
    qcdScale_defUp = 1.e-6;
    qcdScale_defDown = 1.e-6;
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

  if(doRebinMET) RebinMET();

  ScaleQCD();
  if(n == 0) {
    SaveLimitOutputs();
    // KillZeroBins doesn't quite work yet -- it gives a weird bin ~ -600 to 20
    //if(controlRegion != kSR1) SaveLimitOutputs();
    //else SaveLimitOutputs_KillZeroBin();
  }

  StackHistograms(n);
  if(n == 0) METDifference();

  CalculateRatio(n);
  if(n == 0)  MakeLegends();

  if(divideByBinWidth[n]) DivideWidth();

  SetStyles(n);

  padhi->cd();

  bkg->Draw("hist");
  for(unsigned int i = 0; i < mc.size(); i++) {
    if(!needsQCD && i == 0) continue;
    if(i > 0 && layerLegends[i] == layerLegends[i-1]) continue;
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
  h->Divide(bkg);

  TString outName = "met_differences.root";

  TFile * output = new TFile(outName, "UPDATE");

  h->Write();
  output->Close();
}

void PlotMaker::SaveLimitOutputs() {

  TString outName = "limitInputs_";
  if(channel.Contains("ele")) outName += channel(4, 3);
  if(channel.Contains("muon")) outName += channel(5, 3);
  outName += ".root";

  TFile * fLimits = new TFile(outName, "UPDATE");
  fLimits->cd();
  if(channel.Contains("ele")) {
    if(!(fLimits->GetDirectory("ele_"+crNames[controlRegion]))) fLimits->mkdir("ele_"+crNames[controlRegion]);
    fLimits->cd("ele_"+crNames[controlRegion]);
  }
  else {
    if(!(fLimits->GetDirectory("muon_"+crNames[controlRegion]))) fLimits->mkdir("muon_"+crNames[controlRegion]);
    fLimits->cd("muon_"+crNames[controlRegion]);
  }

  data->Write("data_obs");
  qcd->Write("qcd");

  if(channel.Contains("ele")) {
    qcd_defUp->Write("qcd_ele_qcdDefUp");
    qcd_defDown->Write("qcd_ele_qcdDefDown");
  }
  else {
    qcd_defUp->Write("qcd_muon_qcdDefUp");
    qcd_defDown->Write("qcd_muon_qcdDefDown");
  }

  for(int j = 0; j < qcd->GetNbinsX(); j++) {
    TH1D * h_flux_up = (TH1D*)qcd->Clone("clone_qcd_flux_up");
    TH1D * h_flux_down = (TH1D*)qcd->Clone("clone_qcd_flux_down");
    
    Double_t centralValue = qcd->GetBinContent(j+1);
    Double_t statError = qcd->GetBinError(j+1);
    
    if(statError > 0.) h_flux_up->SetBinContent(j+1, centralValue + statError);
    if(centralValue > statError && statError > 0.) h_flux_down->SetBinContent(j+1, centralValue - statError);
    
    // ttjets_ + ttjets_ele_SR2_stat_binX Up/Down
    TString statName = "qcd_qcd_";
    if(channel.Contains("ele")) statName += "ele";
    if(channel.Contains("muon")) statName += "muon";
    statName += "_"+crNames[controlRegion]+"_stat_bin"+Form("%d", j+1);
    
    h_flux_up->Write(statName + "Up");
    h_flux_down->Write(statName + "Down");
  }

  for(unsigned int i = 0; i < mc.size(); i++) {
    mc[i]->Write(limitNames[i]);

    for(int j = 0; j < mc[i]->GetNbinsX(); j++) {
      TH1D * h_flux_up = (TH1D*)mc[i]->Clone("clone_"+limitNames[i]+"_flux_up");
      TH1D * h_flux_down = (TH1D*)mc[i]->Clone("clone_"+limitNames[i]+"_flux_down");

      Double_t centralValue = mc[i]->GetBinContent(j+1);
      Double_t statError = mc[i]->GetBinError(j+1);

      if(statError > 0.) h_flux_up->SetBinContent(j+1, centralValue + statError);
      if(centralValue > statError && statError > 0.) h_flux_down->SetBinContent(j+1, centralValue - statError);

      // ttjets_ + ttjets_ele_SR2_stat_binX Up/Down
      TString statName = limitNames[i]+"_"+limitNames[i]+"_";
      if(channel.Contains("ele")) statName += "ele";
      if(channel.Contains("muon")) statName += "muon";
      statName += "_"+crNames[controlRegion]+"_stat_bin"+Form("%d", j+1);

      h_flux_up->Write(statName + "Up");
      h_flux_down->Write(statName + "Down");
    }
      
  }

  for(unsigned int i = 0; i < mc_btagWeightUp.size(); i++) mc_btagWeightUp[i]->Write(limitNames[i]+"_btagWeightUp");
  for(unsigned int i = 0; i < mc_btagWeightDown.size(); i++) mc_btagWeightDown[i]->Write(limitNames[i]+"_btagWeightDown");

  for(unsigned int i = 0; i < mc_puWeightUp.size(); i++) mc_puWeightUp[i]->Write(limitNames[i]+"_puWeightUp");
  for(unsigned int i = 0; i < mc_puWeightDown.size(); i++) mc_puWeightDown[i]->Write(limitNames[i]+"_puWeightDown");

  for(unsigned int i = 0; i < mc_topPtUp.size(); i++) mc_topPtUp[i]->Write(limitNames[i]+"_topPtUp");
  for(unsigned int i = 0; i < mc_topPtDown.size(); i++) mc_topPtDown[i]->Write(limitNames[i]+"_topPtDown");

  for(unsigned int i = 0; i < mc_scaleUp.size(); i++) {
    if(scaleCorrelations[i] == kTTbar) mc_scaleUp[i]->Write(limitNames[i]+"_scale_ttUp");
    if(scaleCorrelations[i] == kV) mc_scaleUp[i]->Write(limitNames[i]+"_scale_VUp");
    if(scaleCorrelations[i] == kVV) mc_scaleUp[i]->Write(limitNames[i]+"_scale_VVUp");
  }
  for(unsigned int i = 0; i < mc_scaleDown.size(); i++) {
    if(scaleCorrelations[i] == kTTbar) mc_scaleDown[i]->Write(limitNames[i]+"_scale_ttDown");
    if(scaleCorrelations[i] == kV) mc_scaleDown[i]->Write(limitNames[i]+"_scale_VDown");
    if(scaleCorrelations[i] == kVV) mc_scaleDown[i]->Write(limitNames[i]+"_scale_VVDown");
  }

  for(unsigned int i = 0; i < mc_pdfUp.size(); i++) {
    if(pdfCorrelations[i] == kGG) mc_pdfUp[i]->Write(limitNames[i]+"_pdf_ggUp");
    if(pdfCorrelations[i] == kQQ) mc_pdfUp[i]->Write(limitNames[i]+"_pdf_qqUp");
    if(pdfCorrelations[i] == kQG) mc_pdfUp[i]->Write(limitNames[i]+"_pdf_qgUp");
  }
  for(unsigned int i = 0; i < mc_pdfDown.size(); i++) {
    if(pdfCorrelations[i] == kGG) mc_pdfDown[i]->Write(limitNames[i]+"_pdf_ggDown");
    if(pdfCorrelations[i] == kQQ) mc_pdfDown[i]->Write(limitNames[i]+"_pdf_qqDown");
    if(pdfCorrelations[i] == kQG) mc_pdfDown[i]->Write(limitNames[i]+"_pdf_qgDown");
  }
  
  for(unsigned int i = 0; i < mc_JECUp.size(); i++) mc_JECUp[i]->Write(limitNames[i]+"_JECUp");
  for(unsigned int i = 0; i < mc_JECDown.size(); i++) mc_JECDown[i]->Write(limitNames[i]+"_JECDown");
  
  for(unsigned int i = 0; i < mc_leptonSFUp.size(); i++) {
    if(channel.Contains("ele")) mc_leptonSFUp[i]->Write(limitNames[i]+"_eleSFUp");
    if(channel.Contains("muon")) mc_leptonSFUp[i]->Write(limitNames[i]+"_muonSFUp");
  }
  for(unsigned int i = 0; i < mc_leptonSFDown.size(); i++) {
    if(channel.Contains("ele")) mc_leptonSFDown[i]->Write(limitNames[i]+"_eleSFDown");
    if(channel.Contains("muon")) mc_leptonSFDown[i]->Write(limitNames[i]+"_muonSFDown");
  }
    
  for(unsigned int i = 0; i < mc_photonSFUp.size(); i++) mc_photonSFUp[i]->Write(limitNames[i]+"_photonSFUp");
  for(unsigned int i = 0; i < mc_photonSFDown.size(); i++) mc_photonSFDown[i]->Write(limitNames[i]+"_photonSFDown");

  fLimits->Close();
 
}

void PlotMaker::SaveLimitOutputs_KillZeroBin() {

  TString outName = "limitInputs_";
  if(channel.Contains("ele")) outName += channel(4, 3);
  if(channel.Contains("muon")) outName += channel(5, 3);
  outName += ".root";

  TFile * fLimits = new TFile(outName, "UPDATE");
  fLimits->cd();
  if(channel.Contains("ele")) {
    if(!(fLimits->GetDirectory("ele_"+crNames[controlRegion]))) fLimits->mkdir("ele_"+crNames[controlRegion]);
    fLimits->cd("ele_"+crNames[controlRegion]);
  }
  else {
    if(!(fLimits->GetDirectory("muon_"+crNames[controlRegion]))) fLimits->mkdir("muon_"+crNames[controlRegion]);
    fLimits->cd("muon_"+crNames[controlRegion]);
  }

  TH1D * data_cut = (TH1D*)data->Rebin(nMetBins_2g_cut, "data_cut", xbins_met_2g_cut); data_cut->Write("data_obs");
  TH1D * qcd_cut = (TH1D*)qcd->Rebin(nMetBins_2g_cut, "qcd_cut", xbins_met_2g_cut); qcd_cut->Write("qcd");

  TH1D * qcd_defUp_cut = (TH1D*)qcd_defUp->Rebin(nMetBins_2g_cut, "qcd_defUp_cut", xbins_met_2g_cut);
  TH1D * qcd_defDown_cut = (TH1D*)qcd_defDown->Rebin(nMetBins_2g_cut, "qcd_defDown_cut", xbins_met_2g_cut);

  if(channel.Contains("ele")) {
    qcd_defUp_cut->Write("qcd_ele_qcdDefUp");
    qcd_defDown_cut->Write("qcd_ele_qcdDefDown");
  }
  else {
    qcd_defUp_cut->Write("qcd_muon_qcdDefUp");
    qcd_defDown_cut->Write("qcd_muon_qcdDefDown");
  }

  for(int j = 0; j < qcd_cut->GetNbinsX(); j++) {
    TH1D * h_flux_up = (TH1D*)qcd_cut->Clone("clone_qcd_flux_up");
    TH1D * h_flux_down = (TH1D*)qcd_cut->Clone("clone_qcd_flux_down");
    
    Double_t centralValue = qcd_cut->GetBinContent(j+1);
    Double_t statError = qcd_cut->GetBinError(j+1);
    
    if(statError > 0.) h_flux_up->SetBinContent(j+1, centralValue + statError);
    if(centralValue > statError && statError > 0.) h_flux_down->SetBinContent(j+1, centralValue - statError);
    
    // ttjets_ + ttjets_ele_SR2_stat_binX Up/Down
    TString statName = "qcd_qcd_";
    if(channel.Contains("ele")) statName += "ele";
    if(channel.Contains("muon")) statName += "muon";
    statName += "_"+crNames[controlRegion]+"_stat_bin"+Form("%d", j+1);
    
    h_flux_up->Write(statName + "Up");
    h_flux_down->Write(statName + "Down");
  }

  for(unsigned int i = 0; i < mc.size(); i++) {
    TH1D * mc_cut = (TH1D*)mc[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_cut", xbins_met_2g_cut);
    mc_cut->Write(limitNames[i]);

    for(int j = 0; j < mc_cut->GetNbinsX(); j++) {
      TH1D * h_flux_up = (TH1D*)mc_cut->Clone("clone_"+limitNames[i]+"_flux_up");
      TH1D * h_flux_down = (TH1D*)mc_cut->Clone("clone_"+limitNames[i]+"_flux_down");

      Double_t centralValue = mc_cut->GetBinContent(j+1);
      Double_t statError = mc_cut->GetBinError(j+1);

      if(statError > 0.) h_flux_up->SetBinContent(j+1, centralValue + statError);
      if(centralValue > statError && statError > 0.) h_flux_down->SetBinContent(j+1, centralValue - statError);

      // ttjets_ + ttjets_ele_SR2_stat_binX Up/Down
      TString statName = limitNames[i]+"_"+limitNames[i]+"_";
      if(channel.Contains("ele")) statName += "ele";
      if(channel.Contains("muon")) statName += "muon";
      statName += "_"+crNames[controlRegion]+"_stat_bin"+Form("%d", j+1);

      h_flux_up->Write(statName + "Up");
      h_flux_down->Write(statName + "Down");
    }
      
  }

  for(unsigned int i = 0; i < mc_btagWeightUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_btagWeightUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_btagWeightUp_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_btagWeightUp");
  }
  for(unsigned int i = 0; i < mc_btagWeightDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_btagWeightDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_btagWeightDown_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_btagWeightDown");
  }

  for(unsigned int i = 0; i < mc_puWeightUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_puWeightUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_puWeightUp_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_puWeightUp");
  }
  for(unsigned int i = 0; i < mc_puWeightDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_puWeightDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_puWeightDown_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_puWeightDown");
  }

  for(unsigned int i = 0; i < mc_topPtUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_topPtUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_topPtUp_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_topPtUp");
  }
  for(unsigned int i = 0; i < mc_topPtDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_topPtDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_topPtDown_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_topPtDown");
  }

  for(unsigned int i = 0; i < mc_scaleUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_scaleUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_scaleUp_cut", xbins_met_2g_cut);
    if(scaleCorrelations[i] == kTTbar) h_cut->Write(limitNames[i]+"_scale_ttUp");
    if(scaleCorrelations[i] == kV) h_cut->Write(limitNames[i]+"_scale_VUp");
    if(scaleCorrelations[i] == kVV) h_cut->Write(limitNames[i]+"_scale_VVUp");
  }
  for(unsigned int i = 0; i < mc_scaleDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_scaleDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_scaleDown_cut", xbins_met_2g_cut);
    if(scaleCorrelations[i] == kTTbar) h_cut->Write(limitNames[i]+"_scale_ttDown");
    if(scaleCorrelations[i] == kV) h_cut->Write(limitNames[i]+"_scale_VDown");
    if(scaleCorrelations[i] == kVV) h_cut->Write(limitNames[i]+"_scale_VVDown");
  }

  for(unsigned int i = 0; i < mc_pdfUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_pdfUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_pdfUp_cut", xbins_met_2g_cut);
    if(pdfCorrelations[i] == kGG) h_cut->Write(limitNames[i]+"_pdf_ggUp");
    if(pdfCorrelations[i] == kQQ) h_cut->Write(limitNames[i]+"_pdf_qqUp");
    if(pdfCorrelations[i] == kQG) h_cut->Write(limitNames[i]+"_pdf_qgUp");
  }
  for(unsigned int i = 0; i < mc_pdfDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_pdfDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_pdfDown_cut", xbins_met_2g_cut);
    if(pdfCorrelations[i] == kGG) h_cut->Write(limitNames[i]+"_pdf_ggDown");
    if(pdfCorrelations[i] == kQQ) h_cut->Write(limitNames[i]+"_pdf_qqDown");
    if(pdfCorrelations[i] == kQG) h_cut->Write(limitNames[i]+"_pdf_qgDown");
  }

  for(unsigned int i = 0; i < mc_JECUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_JECUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_JECUp_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_JECUp");
  }
  for(unsigned int i = 0; i < mc_JECDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_JECDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_JECDown_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_JECDown");
  }

  for(unsigned int i = 0; i < mc_leptonSFUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_leptonSFUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_leptonSFUp_cut", xbins_met_2g_cut);
    if(channel.Contains("ele")) h_cut->Write(limitNames[i]+"_eleSFUp");
    if(channel.Contains("muon")) h_cut->Write(limitNames[i]+"_muonSFUp");
  }
  for(unsigned int i = 0; i < mc_leptonSFDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_leptonSFDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_leptonSFDown_cut", xbins_met_2g_cut);
    if(channel.Contains("ele")) h_cut->Write(limitNames[i]+"_eleSFDown");
    if(channel.Contains("muon")) h_cut->Write(limitNames[i]+"_muonSFDown");
  }

  for(unsigned int i = 0; i < mc_photonSFUp.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_photonSFUp[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_photonSFUp_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_photonSFUp");
  }
  for(unsigned int i = 0; i < mc_photonSFDown.size(); i++) {
    TH1D * h_cut = (TH1D*)mc_photonSFDown[i]->Rebin(nMetBins_2g_cut, limitNames[i]+"_photonSFDown_cut", xbins_met_2g_cut);
    h_cut->Write(limitNames[i]+"_photonSFDown");
  }

  fLimits->Close();
 
}

void PlotMaker::DetermineAxisRanges(unsigned int n) {

  double padhi_max = -1.;
  double padhi_min = 1.e10;

  double multiply_up = 20.0;
  double multiply_down = 0.66;

  double multiply_up_linear = 1.3;

  for(Int_t ibin = 0; ibin < data->GetNbinsX(); ibin++) {
    if(xMaximums[n] > xMinimums[n] && ibin > data->FindBin(xMaximums[n])) break;
    double value_up = data->GetBinContent(ibin+1) + data->GetBinError(ibin+1) * multiply_up;
    double value_down = (data->GetBinContent(ibin+1) - data->GetBinError(ibin+1)) * multiply_down;

    if(value_up > padhi_max) padhi_max = value_up;
    if(value_down < padhi_min && value_down > 0.) padhi_min = value_down;
  }

  for(Int_t ibin = 0; ibin < mc.back()->GetNbinsX(); ibin++) {
    if(xMaximums[n] > xMinimums[n] && ibin > mc.back()->FindBin(xMaximums[n])) break;
    double value_up = mc.back()->GetBinContent(ibin+1) + mc.back()->GetBinError(ibin+1) * multiply_up;
    double value_down = (mc.back()->GetBinContent(ibin+1) - mc.back()->GetBinError(ibin+1)) * multiply_down;

    if(value_up > padhi_max) padhi_max = value_up;
    if(value_down < padhi_min && value_down > 0.) padhi_min = value_down;
  }

  for(Int_t ibin = 0; ibin < siga->GetNbinsX(); ibin++) {
    if(xMaximums[n] > xMinimums[n] && ibin > siga->FindBin(xMaximums[n])) break;
    double value_up = siga->GetBinContent(ibin+1) * multiply_up;
    double value_down = siga->GetBinContent(ibin+1) * multiply_down;

    if(value_up > padhi_max) padhi_max = value_up;
    if(value_down < padhi_min && value_down > 0.) padhi_min = value_down;
  }

  for(Int_t ibin = 0; ibin < sigb->GetNbinsX(); ibin++) {
    if(xMaximums[n] > xMinimums[n] && ibin > sigb->FindBin(xMaximums[n])) break;
    double value_up = sigb->GetBinContent(ibin+1) * multiply_up;
    double value_down = sigb->GetBinContent(ibin+1) * multiply_down;

    if(value_up > padhi_max) padhi_max = value_up;
    if(value_down < padhi_min && value_down > 0.) padhi_min = value_down;
  }

  double padlo_max = -1.;
  double padlo_min = 0.;

  for(Int_t ibin = 0; ibin < ratio->GetNbinsX(); ibin++) {
    if(xMaximums[n] > xMinimums[n] && ibin > ratio->FindBin(xMaximums[n])) break;
    double value = ratio->GetBinContent(ibin+1) + ratio->GetBinError(ibin+1) * multiply_up_linear;

    if(value > padlo_max) padlo_max = value;
  }

  for(Int_t ibin = 0; ibin < ratio_sys->GetNbinsX(); ibin++) {
    if(xMaximums[n] > xMinimums[n] && ibin > ratio_sys->FindBin(xMaximums[n])) break;
    double value = ratio_sys->GetBinContent(ibin+1) + ratio_sys->GetBinError(ibin+1) * multiply_up_linear;

    if(value > padlo_max) padlo_max = value;
  }

  yMinimums[n] = padhi_min;
  yMaximums[n] = padhi_max;

  ratioMinimums[n] = padlo_min;
  ratioMaximums[n] = padlo_max;

}

void PlotMaker::DetermineLegendRanges(unsigned int n) {

  if(!doDrawLegend[n] && !doDrawPrelim[n]) return;

  if(xMaximums[n] > xMinimums[n]) errors_sys->GetXaxis()->SetRangeUser(xMinimums[n], xMaximums[n]);

  padhi->cd();

  while(true) {
    bool thisOverlaps = false;

    for(Int_t ibin = 0; ibin < errors_sys->GetNbinsX(); ibin++) {
      
      errors_sys->GetYaxis()->SetRangeUser(yMinimums[n], yMaximums[n]);
      errors_sys->Draw("e2");
      leg->Draw("same");
      reqText->Draw("same");

      double binCenter = errors_sys->GetXaxis()->GetBinCenter(ibin+1);
      double sys_max = errors_sys->GetBinContent(ibin+1) + errors_sys->GetBinError(ibin+1);
      double siga_max = siga->GetBinContent(ibin+1);
      double sigb_max = sigb->GetBinContent(ibin+1);

      if(reqText->IsInside(binCenter, sys_max) ||
	 (doDrawSignal[n] && reqText->IsInside(binCenter, siga_max)) ||
	 (doDrawSignal[n] && reqText->IsInside(binCenter, sigb_max)) ||
	 leg->IsInside(binCenter, sys_max) ||
	 (doDrawSignal[n] && leg->IsInside(binCenter, siga_max)) ||
	 (doDrawSignal[n] && leg->IsInside(binCenter, sigb_max))) {
	
	yMaximums[n] = yMaximums[n] * 5.0;
	thisOverlaps = true;
	break;
      }

    }

    if(!thisOverlaps) break;
  }

}
