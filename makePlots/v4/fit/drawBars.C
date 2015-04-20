#include <iostream>
#include <fstream>

using namespace std;

void drawBars() {

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TH1D * central_ttjets_ele = new TH1D("central_ttjets_ele", "central_ttjets_ele", 12, 0, 12);
  TH1D * central_ttjets_muon = new TH1D("central_ttjets_muon", "central_ttjets_muon", 12, 0, 12);

  TH1D * ttjets_ele = new TH1D("ttjets_ele", "ttjets_ele", 12, 0, 12);
  TH1D * ttjets_muon = new TH1D("ttjets_muon", "ttjets_muon", 12, 0, 12);

  double topSF, topSFerror, wjetsSF, wjetsSFerror, qcdSF, qcdSFerror;
  string systematic;

  ifstream eleFile;
  eleFile.open("fitResults_ele_bjj.txt");

  // skip the first line
  getline(eleFile, systematic);

  // get the central value (systematic is empty)
  eleFile >> topSF >> topSFerror >> wjetsSF >> wjetsSFerror >> qcdSF >> qcdSFerror;
  
  for(int i = 0; i < 12; i++) {
    central_ttjets_ele->SetBinContent(i+1, topSF);
    central_ttjets_ele->SetBinError(i+1, topSFerror);
  }

  for(int i = 0; i < 12; i++) {
    eleFile >> systematic >> topSF >> topSFerror >> wjetsSF >> wjetsSFerror >> qcdSF >> qcdSFerror;
    ttjets_ele->SetBinContent(i+1, topSF);
    ttjets_ele->SetBinError(i+1, topSFerror);
  }

  eleFile.close();

  ifstream muonFile;
  muonFile.open("fitResults_muon_bjj.txt");

  // skip the first line
  getline(muonFile, systematic);

  // get the central value (systematic is empty)
  muonFile >> topSF >> topSFerror >> wjetsSF >> wjetsSFerror >> qcdSF >> qcdSFerror;
  
  for(int i = 0; i < 12; i++) {
    central_ttjets_muon->SetBinContent(i+1, topSF);
    central_ttjets_muon->SetBinError(i+1, topSFerror);
  }

  for(int i = 0; i < 12; i++) {
    muonFile >> systematic >> topSF >> topSFerror >> wjetsSF >> wjetsSFerror >> qcdSF >> qcdSFerror;
    ttjets_muon->SetBinContent(i+1, topSF);
    ttjets_muon->SetBinError(i+1, topSFerror);
  }

  muonFile.close();
  
  TString names[12] = {"btagWeightUp", "btagWeightDown",
		       "pileupUp", "pileupDown",
		       "topPtUp", "topPtDown",
		       "JECUp", "JECDown",
		       "leptonSFUp", "leptonSFDown",
		       "photonSFUp", "photonSFDown"};

  for(int i = 0; i < 12; i++) {
    ttjets_ele->GetXaxis()->SetBinLabel(i+1, names[i]);
    ttjets_muon->GetXaxis()->SetBinLabel(i+1, names[i]);
  }

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);

  central_ttjets_ele->SetFillColor(kGray);
  central_ttjets_muon->SetFillColor(kGray);

  TLine * line_ttjets_ele = new TLine(0, central_ttjets_ele->GetBinContent(1), 16, central_ttjets_ele->GetBinContent(1)); line_ttjets_ele->SetLineStyle(2);
  TLine * line_ttjets_muon = new TLine(0, central_ttjets_muon->GetBinContent(1), 16, central_ttjets_muon->GetBinContent(1)); line_ttjets_muon->SetLineStyle(2);

  ttjets_ele->SetLineWidth(2);
  ttjets_muon->SetLineWidth(2);

  ttjets_ele->Draw("axis");
  central_ttjets_ele->Draw("e2 same");
  line_ttjets_ele->Draw("same");
  ttjets_ele->Draw("e1 same");
  ttjets_ele->Draw("axis same");

  can->SaveAs("ttjets_ele_bars.png");

  ttjets_muon->Draw("axis");
  central_ttjets_muon->Draw("e2 same");
  line_ttjets_muon->Draw("same");
  ttjets_muon->Draw("e1 same");
  ttjets_muon->Draw("axis same");

  can->SaveAs("ttjets_muon_bars.png");

}

