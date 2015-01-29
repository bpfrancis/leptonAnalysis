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

using namespace std;

const int nChannels = 4;
TString channels[nChannels] = {"ele_jjj", "ele_bjj",
			       "muon_jjj", "muon_bjj"};

TString channelLabels[nChannels] = {"XYZ e (no b-tag)", "XYZ e",
				    "XYZ #mu (no b-tag)", "XYZ #mu"};

enum controlRegions {kSR1, kSR2, kCR1, kCR2, kNumControlRegions};
TString crNames[kNumControlRegions] = {"SR1", "SR2", "CR1", "CR2"};

class PlotMaker : public TObject {

  ClassDef(PlotMaker, 1);

 public:
  PlotMaker(int chanNo, int cr, bool useQCD);
  ~PlotMaker();

  void BookMCLayer(vector<TString> newNames, int color, TString legendEntry) { 
    TH1D * h;
    mc.push_back(h);

    layerNames.push_back(newNames);
    layerColors.push_back(color);
    layerLegends.push_back(legendEntry);
  };

  void BookPlot(TString variable, bool divideByWidth,
		TString xaxisTitle, TString yaxisTitle,
		Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax,
		Float_t ratiomin, Float_t ratiomax,
		bool drawSignal, bool drawLegend, bool drawPrelim);

  // done once for all variables

  void MakeCanvas() {
    can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
    padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
    padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
    padhi->SetLogy(false);
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

  void CreatePlot(unsigned int n);
  void CreatePlots() {
    MakeCanvas();
    MakeLegends();
    for(unsigned int i = 0; i < variables.size(); i++) CreatePlot(i);
  };
  

 private:
  TFile * input;

  // only one copy of each below for all variables -- for each variable, copy over histograms
  
  TH1D * data;
  TH1D * qcd;
  TH1D * siga;
  TH1D * sigb;

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

  delete can;

  input->Close();

}

void PlotMaker::GetHistograms(unsigned int n) {

  data = (TH1D*)input->Get(variables[n]+"_gg_"+channel);
  qcd = (TH1D*)input->Get(variables[n]+"_qcd_"+channel);
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

  if(divideByBinWidth[n]) {
    data = (TH1D*)DivideByBinWidth(data);
    qcd = (TH1D*)DivideByBinWidth(qcd);
    siga = (TH1D*)DivideByBinWidth(siga);
    sigb = (TH1D*)DivideByBinWidth(sigb);
    
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

  legDrawSignal = (TLegend*)leg->Clone("legDrawSignal");
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
    ratio->GetXaxis()->SetRangeUser(xmin, xmax);
  }

  bkg->GetYaxis()->SetTitle(yaxisTitles[n]);
  ratio->GetYaxis()->SetTitle("Data / Background");
  ratio->GetYaxis()->SetLabelFont(63);
  ratio->GetYaxis()->SetLabelSize(48);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetYaxis()->SetNdivisions(508);

  bkg->GetYaxis()->SetRangeUser(yMinimums[n], yMaxmums[n]);
  ratio->GetYaxis()->SetRangeUser(ratioMinimums[n], ratioMaximums[n]);

  oneLine = new TLine(xmin, 1, xmax, 1);
  oneLine->SetLineStyle(2);

  siga->SetLineColor(kMagenta);
  siga->SetLineWidth(3);
  
  sigb->SetLineColor(kBlue);
  sigb->SetLineWidth(3);

}

void PlotMaker::CreatePlot(unsigned int n) {

  GetHistograms(n);
  StackHistograms(n);
  CalculateRatio(n);
  SetStyles(n);

  padhi->cd();

  bkg->Draw("hist");
  if(needsQCD) mc[0]->Draw("same hist");
  for(unsigned int i = 0; i < mc.size(); i++) mc[i]->Draw("same hist");
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
  legRatio->Draw("same");

  oneLine->Draw();  

  can->SaveAs(variables[n]+"_"+crNames[controlRegion]+".pdf");

}
