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

using namespace std;

void eventCounts() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TString backgrounds[10] = {"qcd", "ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};
  TString bkgLabels[10] = {"QCD", "t#bar{t} + Jets", "W + Jets", "Z/#gamma* + Jets", "Single top", 
			   "VV, V#gamma", "do not draw", "t#bar{t} + V", "do not draw", "t#bar{t} + #gamma"};
  int bkgColors[10] = {kSpring-6, kGray, kOrange-3, kYellow, kRed, kViolet-2, kAzure-2, kCyan, kOrange-5, 8};

  TString channels[6] = {"ele_Any", "muon_Any",
			 "ele_SR1", "muon_SR1",
			 "ele_SR2", "muon_SR2"};

  TFile * input = new TFile("../limitInputs_bjj.root", "READ");

  double value, error;

  vector<TH1D*> nphoEleBkgs, nphoMuonBkgs;
  for(int i = 0; i < 10; i++) {
    TH1D * h_npho = new TH1D("npho_ele_"+backgrounds[i], "npho_ele_"+backgrounds[i], 3, 0, 3);
    h_npho->Sumw2();
    nphoEleBkgs.push_back(h_npho);

    h_npho = new TH1D("npho_muon_"+backgrounds[i], "npho_muon_"+backgrounds[i], 3, 0, 3);
    h_npho->Sumw2();
    nphoMuonBkgs.push_back(h_npho);
  }

  TH1D * nphoEleData = new TH1D("npho_ele_data", "npho_ele_data", 3, 0, 3);
  nphoEleData->Sumw2();

  TH1D * nphoMuonData = new TH1D("npho_muon_data", "npho_muon_data", 3, 0, 3);
  nphoMuonData->Sumw2();

  TH1D * nphoEleSiga = new TH1D("npho_ele_siga", "npho_ele_siga", 3, 0, 3);
  nphoEleSiga->Sumw2();

  TH1D * nphoEleSigb = new TH1D("npho_ele_sigb", "npho_ele_sigb", 3, 0, 3);
  nphoEleSigb->Sumw2();
  
  TH1D * nphoMuonSiga = new TH1D("npho_muon_siga", "npho_muon_siga", 3, 0, 3);
  nphoMuonSiga->Sumw2();

  TH1D * nphoMuonSigb = new TH1D("npho_muon_sigb", "npho_muon_sigb", 3, 0, 3);
  nphoMuonSigb->Sumw2();

  for(int i = 0; i < 6; i++) {

    TH1D * h = (TH1D*)input->Get(channels[i]+"/data_obs");
    value = h->IntegralAndError(0, -1, error);
    cout << "data_obs " << channels[i] << " = " << value << " pm " << error << endl;

    if(channels[i] == "ele_Any") {
      int bin = nphoEleData->GetXaxis()->FindBin(0.);
      nphoEleData->SetBinContent(bin, value);
      nphoEleData->SetBinError(bin, error);
    }
    if(channels[i] == "muon_Any") {
      int bin = nphoMuonData->GetXaxis()->FindBin(0.);
      nphoMuonData->SetBinContent(bin, value);
      nphoMuonData->SetBinError(bin, error);
    }
    if(channels[i] == "ele_SR1") {
      int bin = nphoEleData->GetXaxis()->FindBin(1.);
      nphoEleData->SetBinContent(bin, value);
      nphoEleData->SetBinError(bin, error);
      
      double oldcontent = nphoEleData->GetBinContent(bin-1);
      double olderror = nphoEleData->GetBinError(bin-1);
      nphoEleData->SetBinContent(bin-1, oldcontent - value);
      nphoEleData->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "muon_SR1") {
      int bin = nphoMuonData->GetXaxis()->FindBin(1.);
      nphoMuonData->SetBinContent(bin, value);
      nphoMuonData->SetBinError(bin, error);
      
      double oldcontent = nphoMuonData->GetBinContent(bin-1);
      double olderror = nphoMuonData->GetBinError(bin-1);
      nphoMuonData->SetBinContent(bin-1, oldcontent - value);
      nphoMuonData->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "ele_SR2") {
      int bin = nphoEleData->GetXaxis()->FindBin(2.);
      nphoEleData->SetBinContent(bin, value);
      nphoEleData->SetBinError(bin, error);
      
      double oldcontent = nphoEleData->GetBinContent(bin-2);
      double olderror = nphoEleData->GetBinError(bin-2);
      nphoEleData->SetBinContent(bin-2, oldcontent - value);
      nphoEleData->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "muon_SR2") {
      int bin = nphoMuonData->GetXaxis()->FindBin(2.);
      nphoMuonData->SetBinContent(bin, value);
      nphoMuonData->SetBinError(bin, error);
      
      double oldcontent = nphoMuonData->GetBinContent(bin-2);
      double olderror = nphoMuonData->GetBinError(bin-2);
      nphoMuonData->SetBinContent(bin-2, oldcontent - value);
      nphoMuonData->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
    }

    double total = 0;
    double totalError = 0;

    for(int j = 0; j < 10; j++) {

      h = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]);
      value = h->IntegralAndError(0, -1, error);
      cout << backgrounds[j] << " " << channels[i] << " = " << value << " pm " << error << endl;

      if(i > 1 && j == 0) continue; // don't add QCD to total for SR1/2

      if(value > 0) {

	if(channels[i] == "ele_Any") {
	  int bin = nphoEleBkgs[j]->GetXaxis()->FindBin(0.);
	  nphoEleBkgs[j]->SetBinContent(bin, value);
	  nphoEleBkgs[j]->SetBinError(bin, error);
	}
	if(channels[i] == "muon_Any") {
	  int bin = nphoMuonBkgs[j]->GetXaxis()->FindBin(0.);
	  nphoMuonBkgs[j]->SetBinContent(bin, value);
	  nphoMuonBkgs[j]->SetBinError(bin, error);
	}
	if(channels[i] == "ele_SR1") {
	  int bin = nphoEleBkgs[j]->GetXaxis()->FindBin(1.);
	  nphoEleBkgs[j]->SetBinContent(bin, value);
	  nphoEleBkgs[j]->SetBinError(bin, error);
      
	  double oldcontent = nphoEleBkgs[j]->GetBinContent(bin-1);
	  double olderror = nphoEleBkgs[j]->GetBinError(bin-1);
	  nphoEleBkgs[j]->SetBinContent(bin-1, oldcontent - value);
	  nphoEleBkgs[j]->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
	}
	if(channels[i] == "muon_SR1") {
	  int bin = nphoMuonBkgs[j]->GetXaxis()->FindBin(1.);
	  nphoMuonBkgs[j]->SetBinContent(bin, value);
	  nphoMuonBkgs[j]->SetBinError(bin, error);
      
	  double oldcontent = nphoMuonBkgs[j]->GetBinContent(bin-1);
	  double olderror = nphoMuonBkgs[j]->GetBinError(bin-1);
	  nphoMuonBkgs[j]->SetBinContent(bin-1, oldcontent - value);
	  nphoMuonBkgs[j]->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
	}
	if(channels[i] == "ele_SR2") {
	  int bin = nphoEleBkgs[j]->GetXaxis()->FindBin(2.);
	  nphoEleBkgs[j]->SetBinContent(bin, value);
	  nphoEleBkgs[j]->SetBinError(bin, error);
      
	  double oldcontent = nphoEleBkgs[j]->GetBinContent(bin-2);
	  double olderror = nphoEleBkgs[j]->GetBinError(bin-2);
	  nphoEleBkgs[j]->SetBinContent(bin-2, oldcontent - value);
	  nphoEleBkgs[j]->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
	}
	if(channels[i] == "muon_SR2") {
	  int bin = nphoMuonBkgs[j]->GetXaxis()->FindBin(2.);
	  nphoMuonBkgs[j]->SetBinContent(bin, value);
	  nphoMuonBkgs[j]->SetBinError(bin, error);
      
	  double oldcontent = nphoMuonBkgs[j]->GetBinContent(bin-2);
	  double olderror = nphoMuonBkgs[j]->GetBinError(bin-2);
	  nphoMuonBkgs[j]->SetBinContent(bin-2, oldcontent - value);
	  nphoMuonBkgs[j]->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
	}

	total += value;
	totalError = sqrt(totalError*totalError + error*error);
      }
    }

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175");
    value = h->IntegralAndError(0, -1, error);
    cout << "signal 460/175 " << channels[i] << " = " << value << " pm " << error << endl;

    if(channels[i] == "ele_Any") {
      int bin = nphoEleSiga->GetXaxis()->FindBin(0.);
      nphoEleSiga->SetBinContent(bin, value);
      nphoEleSiga->SetBinError(bin, error);
    }
    if(channels[i] == "muon_Any") {
      int bin = nphoMuonSiga->GetXaxis()->FindBin(0.);
      nphoMuonSiga->SetBinContent(bin, value);
      nphoMuonSiga->SetBinError(bin, error);
    }
    if(channels[i] == "ele_SR1") {
      int bin = nphoEleSiga->GetXaxis()->FindBin(1.);
      nphoEleSiga->SetBinContent(bin, value);
      nphoEleSiga->SetBinError(bin, error);
      
      double oldcontent = nphoEleSiga->GetBinContent(bin-1);
      double olderror = nphoEleSiga->GetBinError(bin-1);
      nphoEleSiga->SetBinContent(bin-1, oldcontent - value);
      nphoEleSiga->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "muon_SR1") {
      int bin = nphoMuonSiga->GetXaxis()->FindBin(1.);
      nphoMuonSiga->SetBinContent(bin, value);
      nphoMuonSiga->SetBinError(bin, error);
      
      double oldcontent = nphoMuonSiga->GetBinContent(bin-1);
      double olderror = nphoMuonSiga->GetBinError(bin-1);
      nphoMuonSiga->SetBinContent(bin-1, oldcontent - value);
      nphoMuonSiga->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "ele_SR2") {
      int bin = nphoEleSiga->GetXaxis()->FindBin(2.);
      nphoEleSiga->SetBinContent(bin, value);
      nphoEleSiga->SetBinError(bin, error);
      
      double oldcontent = nphoEleSiga->GetBinContent(bin-2);
      double olderror = nphoEleSiga->GetBinError(bin-2);
      nphoEleSiga->SetBinContent(bin-2, oldcontent - value);
      nphoEleSiga->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "muon_SR2") {
      int bin = nphoMuonSiga->GetXaxis()->FindBin(2.);
      nphoMuonSiga->SetBinContent(bin, value);
      nphoMuonSiga->SetBinError(bin, error);
      
      double oldcontent = nphoMuonSiga->GetBinContent(bin-2);
      double olderror = nphoMuonSiga->GetBinError(bin-2);
      nphoMuonSiga->SetBinContent(bin-2, oldcontent - value);
      nphoMuonSiga->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
    }

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325");
    value = h->IntegralAndError(0, -1, error);
    cout << "signal 560/325 " << channels[i] << " = " << value << " pm " << error << endl;

    if(channels[i] == "ele_Any") {
      int bin = nphoEleSigb->GetXaxis()->FindBin(0.);
      nphoEleSigb->SetBinContent(bin, value);
      nphoEleSigb->SetBinError(bin, error);
    }
    if(channels[i] == "muon_Any") {
      int bin = nphoMuonSigb->GetXaxis()->FindBin(0.);
      nphoMuonSigb->SetBinContent(bin, value);
      nphoMuonSigb->SetBinError(bin, error);
    }
    if(channels[i] == "ele_SR1") {
      int bin = nphoEleSigb->GetXaxis()->FindBin(1.);
      nphoEleSigb->SetBinContent(bin, value);
      nphoEleSigb->SetBinError(bin, error);
      
      double oldcontent = nphoEleSigb->GetBinContent(bin-1);
      double olderror = nphoEleSigb->GetBinError(bin-1);
      nphoEleSigb->SetBinContent(bin-1, oldcontent - value);
      nphoEleSigb->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "muon_SR1") {
      int bin = nphoMuonSigb->GetXaxis()->FindBin(1.);
      nphoMuonSigb->SetBinContent(bin, value);
      nphoMuonSigb->SetBinError(bin, error);
      
      double oldcontent = nphoMuonSigb->GetBinContent(bin-1);
      double olderror = nphoMuonSigb->GetBinError(bin-1);
      nphoMuonSigb->SetBinContent(bin-1, oldcontent - value);
      nphoMuonSigb->SetBinError(bin-1, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "ele_SR2") {
      int bin = nphoEleSigb->GetXaxis()->FindBin(2.);
      nphoEleSigb->SetBinContent(bin, value);
      nphoEleSigb->SetBinError(bin, error);
      
      double oldcontent = nphoEleSigb->GetBinContent(bin-2);
      double olderror = nphoEleSigb->GetBinError(bin-2);
      nphoEleSigb->SetBinContent(bin-2, oldcontent - value);
      nphoEleSigb->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
    }
    if(channels[i] == "muon_SR2") {
      int bin = nphoMuonSigb->GetXaxis()->FindBin(2.);
      nphoMuonSigb->SetBinContent(bin, value);
      nphoMuonSigb->SetBinError(bin, error);
      
      double oldcontent = nphoMuonSigb->GetBinContent(bin-2);
      double olderror = nphoMuonSigb->GetBinError(bin-2);
      nphoMuonSigb->SetBinContent(bin-2, oldcontent - value);
      nphoMuonSigb->SetBinError(bin-2, sqrt(olderror*olderror - error*error));
    }

    cout << "total background = " << total << " pm " << totalError << endl;

    cout << endl << endl;
  }

  for(int i = 0; i < 10; i++) {
    for(int j = i+1; j < 10; j++) {
      nphoEleBkgs[i]->Add(nphoEleBkgs[j]);
      nphoMuonBkgs[i]->Add(nphoMuonBkgs[j]);
    }

    nphoEleBkgs[i]->SetFillColor(bkgColors[i]);
    nphoEleBkgs[i]->SetMarkerSize(0);
    nphoEleBkgs[i]->SetLineColor(1);
    nphoMuonBkgs[i]->SetFillColor(bkgColors[i]);
    nphoMuonBkgs[i]->SetMarkerSize(0);
    nphoMuonBkgs[i]->SetLineColor(1);
  }

  TH1D * ratioEle = (TH1D*)nphoEleData->Clone("ratio_ele");
  ratioEle->Reset();
  ratioEle->SetTitle("Data / Background");
  for(int i = 0; i < ratioEle->GetNbinsX(); i++) {
    if(nphoEleBkgs[0]->GetBinContent(i+1) == 0.) continue;
    ratioEle->SetBinContent(i+1, nphoEleData->GetBinContent(i+1) / nphoEleBkgs[0]->GetBinContent(i+1));
    ratioEle->SetBinError(i+1, nphoEleData->GetBinError(i+1) / nphoEleBkgs[0]->GetBinContent(i+1));
  }

  TH1D * ratioEle_stat = (TH1D*)nphoEleBkgs[0]->Clone("ratio_ele_stat");
  for(int i = 0; i < ratioEle_stat->GetNbinsX(); i++) {
    ratioEle_stat->SetBinContent(i+1, 1.);
    if(nphoEleBkgs[0]->GetBinContent(i+1) == 0.) ratioEle_stat->SetBinError(i+1, 0.);
    else ratioEle_stat->SetBinError(i+1, nphoEleBkgs[0]->GetBinError(i+1) / nphoEleBkgs[0]->GetBinContent(i+1));
  }

  TH1D * ratioMuon = (TH1D*)nphoMuonData->Clone("ratio_muon");
  ratioMuon->Reset();
  ratioMuon->SetTitle("Data / Background");
  for(int i = 0; i < ratioMuon->GetNbinsX(); i++) {
    if(nphoMuonBkgs[0]->GetBinContent(i+1) == 0.) continue;
    ratioMuon->SetBinContent(i+1, nphoMuonData->GetBinContent(i+1) / nphoMuonBkgs[0]->GetBinContent(i+1));
    ratioMuon->SetBinError(i+1, nphoMuonData->GetBinError(i+1) / nphoMuonBkgs[0]->GetBinContent(i+1));
  }

  TH1D * ratioMuon_stat = (TH1D*)nphoMuonBkgs[0]->Clone("ratio_muon_stat");
  for(int i = 0; i < ratioMuon_stat->GetNbinsX(); i++) {
    ratioMuon_stat->SetBinContent(i+1, 1.);
    if(nphoMuonBkgs[0]->GetBinContent(i+1) == 0.) ratioMuon_stat->SetBinError(i+1, 0.);
    else ratioMuon_stat->SetBinError(i+1, nphoMuonBkgs[0]->GetBinError(i+1) / nphoMuonBkgs[0]->GetBinContent(i+1));
  }

  nphoEleData->SetMarkerStyle(20);
  nphoEleData->SetMarkerSize(1.5);

  nphoMuonData->SetMarkerStyle(20);
  nphoMuonData->SetMarkerSize(1.5);

  nphoEleSiga->SetLineColor(kMagenta);
  nphoEleSiga->SetLineWidth(3);

  nphoMuonSiga->SetLineColor(kMagenta);
  nphoMuonSiga->SetLineWidth(3);

  nphoEleSigb->SetLineColor(kBlue);
  nphoEleSigb->SetLineWidth(3);

  nphoMuonSigb->SetLineColor(kBlue);
  nphoMuonSigb->SetLineWidth(3);
  
  TH1D * bkgErrors_ele = (TH1D*)nphoEleBkgs[0]->Clone("bkgErrors_ele");
  bkgErrors_ele->SetFillColor(kOrange+10);
  bkgErrors_ele->SetFillStyle(3154);
  bkgErrors_ele->SetMarkerSize(0);

  TH1D * bkgErrors_muon = (TH1D*)nphoMuonBkgs[0]->Clone("bkgErrors_muon");
  bkgErrors_muon->SetFillColor(kOrange+10);
  bkgErrors_muon->SetFillStyle(3154);
  bkgErrors_muon->SetMarkerSize(0);

  ratioEle_stat->SetFillStyle(1001);
  ratioEle_stat->SetFillColor(kGray+1);
  ratioEle_stat->SetLineColor(kGray+1);
  ratioEle_stat->SetMarkerColor(kGray+1);

  ratioMuon_stat->SetFillStyle(1001);
  ratioMuon_stat->SetFillColor(kGray+1);
  ratioMuon_stat->SetLineColor(kGray+1);
  ratioMuon_stat->SetMarkerColor(kGray+1);

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
  TPad * padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  TPad * padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
  padhi->SetLogy(true);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetBottomMargin(0);
  
  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);
  
  padhi->Draw();
  padlo->Draw();
  
  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, NULL, "brNDC");
  leg->SetNColumns(2);
  leg->AddEntry(nphoEleData, "Data", "LP");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(bkgErrors_ele, "Errors (stat. only)", "F");
  leg->AddEntry((TObject*)0, "", "");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  for(int i = 0; i < 10; i++) {
    if(bkgLabels[i] == "do not draw") continue;
    leg->AddEntry(nphoEleBkgs[i], bkgLabels[i], "F");
  }

  leg->AddEntry(nphoEleSiga, "GGM (460_175)", "L");
  leg->AddEntry(nphoEleSigb, "GGM (560_325)", "L");

  TPaveText * reqText_ele = new TPaveText(0.45, 0.47, 0.85, 0.57, "NDC");
  reqText_ele->SetFillColor(0);
  reqText_ele->SetFillStyle(0);
  reqText_ele->SetLineColor(0);
  reqText_ele->SetBorderSize(0);
  reqText_ele->AddText("ele");

  TPaveText * reqText_muon = new TPaveText(0.45, 0.47, 0.85, 0.57, "NDC");
  reqText_muon->SetFillColor(0);
  reqText_muon->SetFillStyle(0);
  reqText_muon->SetLineColor(0);
  reqText_muon->SetBorderSize(0);
  reqText_muon->AddText("muon");

  TPaveText * lumiHeader = new TPaveText(0.1, 0.901, 0.9, 0.94, "NDC");
  lumiHeader->SetFillColor(0);
  lumiHeader->SetFillStyle(0);
  lumiHeader->SetLineColor(0);
  lumiHeader->SetBorderSize(0);
  lumiHeader->AddText("CMS Preliminary 2015     #sqrt{s} = 8 TeV     #intL = 19.7 fb^{-1}");

  TLine * oneLine = new TLine(0, 1, 3, 1);
  oneLine->SetLineStyle(2);

  padhi->cd();

  nphoEleBkgs[0]->GetYaxis()->SetRangeUser(2, 5e6);
  nphoMuonBkgs[0]->GetYaxis()->SetRangeUser(2, 5e6);

  ratioEle->GetYaxis()->SetRangeUser(0.75, 2.15);
  ratioEle->GetYaxis()->SetTitle("Data / Background");
  ratioEle->GetXaxis()->SetTitle("N_{#gamma}");
  ratioEle->GetXaxis()->SetLabelFont(63);
  ratioEle->GetXaxis()->SetLabelSize(48);
  ratioEle->GetXaxis()->SetTitleSize(0.12);
  ratioEle->GetXaxis()->SetTitleOffset(0.6);

  ratioMuon->GetYaxis()->SetRangeUser(0.75, 2.15);
  ratioMuon->GetYaxis()->SetTitle("Data / Background");
  ratioMuon->GetXaxis()->SetTitle("N_{#gamma}");
  ratioMuon->GetXaxis()->SetLabelFont(63);
  ratioMuon->GetXaxis()->SetLabelSize(48);
  ratioMuon->GetXaxis()->SetTitleSize(0.12);
  ratioMuon->GetXaxis()->SetTitleOffset(0.6);

  for(int i = 0; i < 10; i++) {
    if(bkgLabels[i] == "do not draw") continue;

    if(i == 0) nphoEleBkgs[i]->Draw("hist");
    else nphoEleBkgs[i]->Draw("hist same");
  }

  nphoEleData->Draw("e1 same");
  bkgErrors_ele->Draw("e2 same");
  nphoEleBkgs[0]->Draw("axis same");
  nphoEleSiga->Draw("hist same");
  nphoEleSigb->Draw("hist same");
  
  lumiHeader->Draw();
  leg->Draw();
  reqText_ele->Draw();

  padlo->cd();
  
  ratioEle->Draw("e1");
  ratioEle_stat->Draw("e2 same");
  ratioEle->Draw("e1 same");
  ratioEle->Draw("axis same");
  oneLine->Draw();

  can->SaveAs("Nphotons_ele.pdf");

  padhi->cd();

  for(int i = 0; i < 10; i++) {
    if(bkgLabels[i] == "do not draw") continue;

    if(i == 0) nphoMuonBkgs[i]->Draw("hist");
    else nphoMuonBkgs[i]->Draw("hist same");
  }

  nphoMuonData->Draw("e1 same");
  bkgErrors_muon->Draw("e2 same");
  nphoMuonBkgs[0]->Draw("axis same");
  nphoMuonSiga->Draw("hist same");
  nphoMuonSigb->Draw("hist same");
  
  lumiHeader->Draw();
  leg->Draw();
  reqText_muon->Draw();

  padlo->cd();
  
  ratioMuon->Draw("e1");
  ratioMuon_stat->Draw("e2 same");
  ratioMuon->Draw("e1 same");
  ratioMuon->Draw("axis same");
  oneLine->Draw();

  can->SaveAs("Nphotons_muon.pdf");



  input->Close();

}
