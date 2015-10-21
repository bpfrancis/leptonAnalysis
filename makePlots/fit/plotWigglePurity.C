#include <vector>

using namespace std;

void compareWiggle(TString folder, TString metCutName, TString version) {

  TString names[9] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};

  TFile * input = new TFile("limitInputsScaled_bjj_"+version+metCutName+".root", "READ");

  TH1D * before = (TH1D*)input->Get(folder+"/ttjets_noScaling");
  TH1D * after = (TH1D*)input->Get(folder+"/ttjets");

  for(int i = 1; i < 9; i++) {
    if(names[i] == "ttgamma") {
      before->Add((TH1D*)input->Get(folder+"/"+names[i]+"_noScaling"));
      after->Add((TH1D*)input->Get(folder+"/"+names[i]));
    }
    else {
      before->Add((TH1D*)input->Get(folder+"/"+names[i]));
      after->Add((TH1D*)input->Get(folder+"/"+names[i]));
    }
  }

  TH1D * ratio = (TH1D*)after->Clone("ratio");
  ratio->Reset();
  ratio->SetTitle("Scaled / Nominal");
  for(int i = 0; i < ratio->GetNbinsX(); i++) {
    if(before->GetBinContent(i+1) == 0.) continue;
    ratio->SetBinContent(i+1, after->GetBinContent(i+1) / before->GetBinContent(i+1));
    ratio->SetBinError(i+1, after->GetBinError(i+1) / before->GetBinContent(i+1));
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

  after->SetLineColor(kRed);

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
  leg->AddEntry(before, "MC Background", "LP");
  leg->AddEntry(after, "Photon purity adjusted", "LP");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  padhi->cd();
  after->Draw("hist");
  before->Draw("hist same");
  leg->Draw("same");

  padlo->cd();
  ratio->Draw("e1");

  can->SaveAs("plots/compareWiggle_"+version+metCutName+"_"+folder+".pdf");
  delete can;

}

void plotWigglePurity() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  compareWiggle("ele_SR1", "", "chHadIso");
  compareWiggle("ele_SR1", "_metCut_50", "chHadIso");

  compareWiggle("muon_SR1", "", "chHadIso");
  compareWiggle("muon_SR1", "_metCut_50", "chHadIso");

  compareWiggle("ele_SR2", "", "chHadIso");
  compareWiggle("ele_SR2", "_metCut_50", "chHadIso");

  compareWiggle("muon_SR2", "", "chHadIso");
  compareWiggle("muon_SR2", "_metCut_50", "chHadIso");

}
