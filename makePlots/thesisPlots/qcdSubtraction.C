#include "../rootRoutines.h"

using namespace std;

void qcdSubtraction() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TFile * fEle = new TFile("qcdHistograms_ele_bjj_Any.root", "READ");
  TFile * fMuon = new TFile("qcdHistograms_muon_bjj_Any.root", "READ");

  TString backgrounds[23] = {"ttJetsHadronic", "ttJetsSemiLep", "ttJetsFullLep",
			     "W3JetsToLNu", "W4JetsToLNu",
			     "dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
			     "TBar_s", "TBar_t", "TBar_tW",
			     "T_s", "T_t", "T_tW",
			     "WW", "WZ", "ZZ",
			     "WGToLNuG", "ZGToLLG",
			     "TTWJets", "TTZJets",
			     "TTGamma"};


  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
  can->SetLogy(true);

  TH1D * qcd_ele = (TH1D*)fEle->Get("pfMET_t01_qcd_ele_bjj");
  TH1D * qcd_muon = (TH1D*)fMuon->Get("pfMET_t01_qcd_muon_bjj");

  TH1D * ttbar_ele = (TH1D*)fEle->Get("pfMET_t01_qcd_ttJetsHadronic_ele_bjj");
  TH1D * ttbar_muon = (TH1D*)fMuon->Get("pfMET_t01_qcd_ttJetsHadronic_muon_bjj");

  ttbar_ele->Add((TH1D*)fEle->Get("pfMET_t01_qcd_ttJetsSemiLep_ele_bjj"));
  ttbar_muon->Add((TH1D*)fMuon->Get("pfMET_t01_qcd_ttJetsSemiLep_muon_bjj"));

  ttbar_ele->Add((TH1D*)fEle->Get("pfMET_t01_qcd_ttJetsFullLep_ele_bjj"));
  ttbar_muon->Add((TH1D*)fMuon->Get("pfMET_t01_qcd_ttJetsFullLep_muon_bjj"));

  TH1D * wjets_ele = (TH1D*)fEle->Get("pfMET_t01_qcd_W3JetsToLNu_ele_bjj");
  TH1D * wjets_muon = (TH1D*)fMuon->Get("pfMET_t01_qcd_W3JetsToLNu_muon_bjj");

  wjets_ele->Add((TH1D*)fEle->Get("pfMET_t01_qcd_W4JetsToLNu_ele_bjj"));
  wjets_muon->Add((TH1D*)fMuon->Get("pfMET_t01_qcd_W4JetsToLNu_muon_bjj"));

  TH1D * other_ele = (TH1D*)fEle->Get("pfMET_t01_qcd_dy1JetsToLL_ele_bjj");
  TH1D * other_muon = (TH1D*)fMuon->Get("pfMET_t01_qcd_dy1JetsToLL_muon_bjj");

  for(int i = 6; i < 23; i++) {
    other_ele->Add((TH1D*)fEle->Get("pfMET_t01_qcd_"+backgrounds[i]+"_ele_bjj"));
    other_muon->Add((TH1D*)fMuon->Get("pfMET_t01_qcd_"+backgrounds[i]+"_muon_bjj"));
  }

  ttbar_ele->Add(wjets_ele);
  ttbar_ele->Add(other_ele);
  ttbar_muon->Add(wjets_muon);
  ttbar_muon->Add(other_muon);

  wjets_ele->Add(other_ele);
  wjets_muon->Add(other_muon);

  qcd_ele->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");
  qcd_ele->GetXaxis()->SetRangeUser(0, 300);
  qcd_ele->GetYaxis()->SetRangeUser(1e-2, 2e5);
  qcd_ele->SetMarkerStyle(20);
  qcd_ele->SetMarkerSize(1.5);

  qcd_muon->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");
  qcd_muon->GetXaxis()->SetRangeUser(0, 300);
  qcd_muon->GetYaxis()->SetRangeUser(1e-2, 2e5);
  qcd_muon->SetMarkerStyle(20);
  qcd_muon->SetMarkerSize(1.5);

  ttbar_ele->SetFillColor(kGray);
  ttbar_ele->SetMarkerSize(0);
  ttbar_ele->SetLineColor(1);

  ttbar_muon->SetFillColor(kGray);
  ttbar_muon->SetMarkerSize(0);
  ttbar_muon->SetLineColor(1);

  wjets_ele->SetFillColor(kOrange-3);
  wjets_ele->SetMarkerSize(0);
  wjets_ele->SetLineColor(1);

  wjets_muon->SetFillColor(kOrange-3);
  wjets_muon->SetMarkerSize(0);
  wjets_muon->SetLineColor(1);

  other_ele->SetFillColor(8);
  other_ele->SetMarkerSize(0);
  other_ele->SetLineColor(1);

  other_muon->SetFillColor(8);
  other_muon->SetMarkerSize(0);
  other_muon->SetLineColor(1);

  TLegend * leg = new TLegend(0.45, 0.65, 0.85, 0.85, NULL, "brNDC");
  leg->AddEntry(qcd_ele, "QCD Data", "LP");
  leg->AddEntry((TObject*)0, "QCD Selection on MC:", "");
  leg->AddEntry(ttbar_ele, "t#bar{t} + Jets", "F");
  leg->AddEntry(wjets_ele, "W + Jets", "F");
  leg->AddEntry(other_ele, "Other MC", "F");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  TPaveText * reqText_ele = new TPaveText(0.45, 0.55, 0.85, 0.62, "NDC");
  reqText_ele->SetFillColor(0);
  reqText_ele->SetFillStyle(0);
  reqText_ele->SetLineColor(0);
  reqText_ele->SetBorderSize(0);
  reqText_ele->AddText("ele");

  TPaveText * reqText_muon = new TPaveText(0.45, 0.55, 0.85, 0.62, "NDC");
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

  qcd_ele->Draw("e1");
  ttbar_ele->Draw("hist same");
  wjets_ele->Draw("hist same");
  other_ele->Draw("hist same");
  qcd_ele->Draw("e1 same");
  qcd_ele->Draw("axis same");

  leg->Draw();
  reqText_ele->Draw();
  lumiHeader->Draw();

  can->SaveAs("qcdSubtraction_ele.pdf");

  qcd_muon->Draw("e1");
  ttbar_muon->Draw("hist same");
  wjets_muon->Draw("hist same");
  other_muon->Draw("hist same");
  qcd_muon->Draw("e1 same");
  qcd_muon->Draw("axis same");

  leg->Draw();
  reqText_muon->Draw();
  lumiHeader->Draw();

  can->SaveAs("qcdSubtraction_muon.pdf");
    

}
