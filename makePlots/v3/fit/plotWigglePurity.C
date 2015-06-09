#include <vector>

using namespace std;

void compareWiggle(TString channel) {

  TString cr = "SR1";
  TString version = "chHadIso";

  TString namesShort[9] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};

  TFile * fBefore = new TFile("prewigglePurity_"+channel+"_"+cr+"_"+version+".root", "READ");
  TFile * fAfter = new TFile("wigglePurity_"+channel+"_"+cr+"_"+version+".root", "READ");

  TH1D * before = (TH1D*)fBefore->Get("pfMET_t01_ttjets_"+channel+"_"+cr);
  for(int i = 1; i < 9; i++) before->Add((TH1D*)fBefore->Get("pfMET_t01_"+namesShort[i]+"_"+channel+"_"+cr));

  TH1D * after = (TH1D*)fAfter->Get("pfMET_t01_ttjets_"+channel+"_"+cr);
  for(int i = 1; i < 9; i++) after->Add((TH1D*)fAfter->Get("pfMET_t01_"+namesShort[i]+"_"+channel+"_"+cr));

  const int nMetBins_1g = 10;
  Double_t xbins_met_1g[nMetBins_1g+1] = {0, 10, 20, 30, 40, 50, 75, 100, 150, 300, 800};

  TH1D * beforeReb = (TH1D*)before->Rebin(nMetBins_1g, "beforeReb", xbins_met_1g);
  TH1D * afterReb = (TH1D*)after->Rebin(nMetBins_1g, "afterReb", xbins_met_1g);

  TH1D * ratio = (TH1D*)afterReb->Clone("ratio");
  ratio->Reset();
  ratio->SetTitle("Scaled / Nominal");
  for(int i = 0; i < ratio->GetNbinsX(); i++) {
    if(beforeReb->GetBinContent(i+1) == 0.) continue;
    ratio->SetBinContent(i+1, afterReb->GetBinContent(i+1) / beforeReb->GetBinContent(i+1));
    ratio->SetBinError(i+1, afterReb->GetBinError(i+1) / beforeReb->GetBinContent(i+1));
  }

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);

  padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
  padhi->SetLogy(true);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetBottomMargin(0);
  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);
  padhi->Draw();
  padlo->Draw();

  afterReb->SetLineColor(kRed);

  ratio->GetXaxis()->SetTitle("MET");
  ratio->GetXaxis()->SetLabelFont(63);
  ratio->GetXaxis()->SetLabelSize(48);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleOffset(0.6);

  ratio->GetYaxis()->SetTitle("After / Before");
  ratio->GetYaxis()->SetLabelFont(63);
  ratio->GetYaxis()->SetLabelSize(48);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetYaxis()->SetNdivisions(508);

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, "", "brNDC");
  leg->AddEntry(beforeReb, "MC Background", "LP");
  leg->AddEntry(afterReb, "Photon purity adjusted", "LP");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  padhi->cd();
  afterReb->Draw("hist");
  beforeReb->Draw("hist same");
  leg->Draw("same");

  padlo->cd();
  ratio->Draw("e1");

  can->SaveAs("plots/compareWiggle_"+channel+"_"+cr+".pdf");
  delete can;

}

void compareWiggle(TString channel, TString background) {

  TString cr = "SR1";
  TString version = "chHadIso";

  TFile * fBefore = new TFile("prewigglePurity_"+channel+"_"+cr+"_"+version+".root", "READ");
  TFile * fAfter = new TFile("wigglePurity_"+channel+"_"+cr+"_"+version+".root", "READ");

  TH1D * before = (TH1D*)fBefore->Get("pfMET_t01_"+background+"_"+channel+"_"+cr);
  TH1D * after = (TH1D*)fAfter->Get("pfMET_t01_"+background+"_"+channel+"_"+cr);

  const int nMetBins_1g = 10;
  Double_t xbins_met_1g[nMetBins_1g+1] = {0, 10, 20, 30, 40, 50, 75, 100, 150, 300, 800};

  TH1D * beforeReb = (TH1D*)before->Rebin(nMetBins_1g, "beforeReb", xbins_met_1g);
  TH1D * afterReb = (TH1D*)after->Rebin(nMetBins_1g, "afterReb", xbins_met_1g);

  TH1D * ratio = (TH1D*)afterReb->Clone("ratio");
  ratio->Reset();
  ratio->SetTitle("Scaled / Nominal");
  for(int i = 0; i < ratio->GetNbinsX(); i++) {
    if(beforeReb->GetBinContent(i+1) == 0.) continue;
    ratio->SetBinContent(i+1, afterReb->GetBinContent(i+1) / beforeReb->GetBinContent(i+1));
    ratio->SetBinError(i+1, afterReb->GetBinError(i+1) / beforeReb->GetBinContent(i+1));
  }

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);

  padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
  padhi->SetLogy(true);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetBottomMargin(0);
  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);
  padhi->Draw();
  padlo->Draw();

  afterReb->SetLineColor(kRed);

  ratio->GetXaxis()->SetTitle("MET");
  ratio->GetXaxis()->SetLabelFont(63);
  ratio->GetXaxis()->SetLabelSize(48);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleOffset(0.6);

  ratio->GetYaxis()->SetTitle("After / Before");
  ratio->GetYaxis()->SetLabelFont(63);
  ratio->GetYaxis()->SetLabelSize(48);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetYaxis()->SetNdivisions(508);

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, "", "brNDC");
  leg->AddEntry(beforeReb, "MC Background", "LP");
  leg->AddEntry(afterReb, "Photon purity adjusted", "LP");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  padhi->cd();
  afterReb->Draw("hist");
  beforeReb->Draw("hist same");
  leg->Draw("same");

  padlo->cd();
  ratio->Draw("e1");

  can->SaveAs("plots/compareWiggle_"+channel+"_"+cr+"_"+background+".pdf");
  delete can;

}

void compareShapes(TString channel) {

  TString cr = "SR1";
  TString version = "chHadIso";

  TString namesFull[23] = {"ttJetsSemiLep", "ttJetsFullLep", "ttJetsHadronic",
			   "W3JetsToLNu", "W4JetsToLNu",
			   "dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
			   "TBar_s", "TBar_t", "TBar_tW", "T_s", "T_t", "T_tW",
			   "WW", "WZ", "ZZ",
			   "WGToLNuG", "ZGToLLG",
			   "TTWJets", "TTZJets",
			   "TTGamma"};

  TString namesShort[9] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};

  TFile * fBefore = new TFile("../fitTemplates.root", "READ");

  TH1D * prompt = (TH1D*)fBefore->Get("pfMET_t01_ttJetsSemiLep_"+channel+"_"+cr+"_matchPhoton");
  prompt->Add((TH1D*)fBefore->Get("pfMET_t01_ttJetsSemiLep_"+channel+"_"+cr+"_matchElectron"));

  for(int i = 1; i < 23; i++) {
    prompt->Add((TH1D*)fBefore->Get("pfMET_t01_"+namesFull[i]+"_"+channel+"_"+cr+"_matchPhoton"));
    prompt->Add((TH1D*)fBefore->Get("pfMET_t01_"+namesFull[i]+"_"+channel+"_"+cr+"_matchElectron"));
  }

  TH1D * nonprompt = (TH1D*)fBefore->Get("pfMET_t01_ttJetsSemiLep_"+channel+"_"+cr+"_matchJet");
  for(int i = 1; i < 23; i++) nonprompt->Add((TH1D*)fBefore->Get("pfMET_t01_"+namesFull[i]+"_"+channel+"_"+cr+"_matchJet"));

  const int nMetBins_1g = 10;
  Double_t xbins_met_1g[nMetBins_1g+1] = {0, 10, 20, 30, 40, 50, 75, 100, 150, 300, 800};

  TH1D * promptReb = (TH1D*)prompt->Rebin(nMetBins_1g, "promptReb", xbins_met_1g);
  TH1D * nonpromptReb = (TH1D*)nonprompt->Rebin(nMetBins_1g, "nonpromptReb", xbins_met_1g);

  promptReb->Scale(1./promptReb->Integral());
  nonpromptReb->Scale(1./nonpromptReb->Integral());

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
  can->SetLogy(true);

  nonpromptReb->SetLineColor(kRed);

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, "", "brNDC");
  leg->AddEntry(promptReb, "Prompt MC Background", "LP");
  leg->AddEntry(nonpromptReb, "Non-prompt MC Background", "LP");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  nonpromptReb->Draw("hist");
  promptReb->Draw("hist same");
  leg->Draw("same");

  can->SaveAs("plots/compareShapes_"+channel+"_"+cr+".pdf");
  delete can;

}

void plotWigglePurity() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  compareWiggle("ele_bjj");
  compareShapes("ele_bjj");

  compareWiggle("muon_bjj");
  compareShapes("muon_bjj");

  TString namesShort[9] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};
  for(int i = 0; i < 9; i++) {
    compareWiggle("ele_bjj", namesShort[i]);
    compareWiggle("muon_bjj", namesShort[i]);
  }

}
