void plotDibosonSplit() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  //gStyle->SetOptTitle(0);
  
  TFile * _file0 = new TFile("ele_sr1.root", "READ");
  TFile * _file1 = new TFile("muon_sr1.root", "READ");
  
  TH1D * ele = (TH1D*)_file0->Get("wg");
  TH1D * muon = (TH1D*)_file1->Get("wg");

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
  
  ele->SetMarkerSize(0);
  ele->SetFillColor(kOrange+10);
  ele->SetFillStyle(3154);
  ele->SetTitle("W#gamma SR1");
  ele->GetXaxis()->SetTitle("MET (GeV)");
  //ele->GetYaxis()->SetTitle("Events");
  ele->Draw("e2");
  
  muon->Draw("e1 same");

  TLegend * leg = new TLegend(0.5, 0.7, 0.85, 0.85, NULL, "brNDC");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);
  leg->AddEntry(ele, "ele #oplus stat. #oplus syst.", "F");
  leg->AddEntry(muon, "muon #oplus stat. #oplus syst.", "LP");
  leg->Draw("same");

  can->SaveAs("split.png");
  
}
