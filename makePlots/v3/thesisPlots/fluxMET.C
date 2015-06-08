void go(TString channel) {

  TFile * input = new TFile("../limitInputs_bjj.root", "READ");

  TString names[9] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma"};

  TH1D * ttgamma = (TH1D*)input->Get(channel+"_SR1/ttgamma");
  ttgamma->Add((TH1D*)input->Get(channel+"_SR2/ttgamma"));

  TH1D * ttjets = (TH1D*)input->Get(channel+"_SR1/ttjets");
  ttjets->Add((TH1D*)input->Get(channel+"_SR2/ttjets"));

  double total_ttx = ttgamma->Integral() + ttjets->Integral();

  TH1D * bkg = (TH1D*)ttjets->Clone("bkg");

  TH1D * bkg_ttgammaOnly = (TH1D*)ttgamma->Clone("bkg_ttgammaOnly");
  bkg_ttgammaOnly->Scale(total_ttx / ttgamma->Integral());

  TH1D * bkg_ttjetsOnly = (TH1D*)ttjets->Clone("bkg_ttjetsOnly");
  bkg_ttjetsOnly->Scale(total_ttx / ttjets->Integral());

  for(int i = 1; i < 9; i++) {
    TH1D * h = (TH1D*)input->Get(channel+"_SR1/"+names[i]);
    h->Add((TH1D*)input->Get(channel+"_SR2/"+names[i]));

    bkg->Add(h);
    if(names[i] != "ttgamma") {
      bkg_ttgammaOnly->Add(h);
      bkg_ttjetsOnly->Add(h);
    }

  }

  TH1D * bkg_syst = (TH1D*)bkg->Clone("bkg_syst");
  for(int i = 0; i < bkg_syst->GetNbinsX(); i++) {
    double a = ttgamma->GetBinContent(i+1);
    double b = ttjets->GetBinContent(i+1);

    bkg_syst->SetBinError(i+1, sqrt(a*a + b*b));
  }

  TCanvas * can = new TCanvas("can", "can", 10, 10, 800, 800);
  can->SetLogy(true);

  bkg->Scale(1./bkg->Integral());
  bkg_syst->Scale(1./bkg_syst->Integral());
  bkg_ttgammaOnly->Scale(1./bkg_ttgammaOnly->Integral());
  bkg_ttjetsOnly->Scale(1./bkg_ttjetsOnly->Integral());
  
  bkg_syst->SetFillColor(kOrange+10);
  bkg_syst->SetFillStyle(3154);
  bkg_syst->SetMarkerSize(0);

  bkg->SetLineColor(kBlack);
  bkg->SetLineWidth(3);
  bkg->GetXaxis()->SetTitle("MET");
  bkg->GetYaxis()->SetRangeUser(9.e-4, 5.e-1);

  bkg_ttgammaOnly->SetLineColor(8);
  bkg_ttgammaOnly->SetLineWidth(3);
  bkg_ttgammaOnly->SetMarkerSize(0);
  
  bkg_ttjetsOnly->SetLineColor(kBlue);
  bkg_ttjetsOnly->SetLineWidth(3);
  bkg_ttjetsOnly->SetMarkerSize(0);

  bkg->Draw("hist");
  bkg_syst->Draw("e2 same");
  bkg_ttgammaOnly->Draw("e1 same");
  bkg_ttjetsOnly->Draw("e1 same");
  bkg->Draw("hist same");
  bkg->Draw("axis same");

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, channel+" (SR1 + SR2)", "brNDC");
  leg->AddEntry(bkg, "Background", "LP");
  leg->AddEntry(bkg_syst, "tt#gamma/ttjets #pm 100%", "F");
  leg->AddEntry(bkg_ttgammaOnly, "Bkg: tt#gamma Only", "LP");
  leg->AddEntry(bkg_ttjetsOnly, "Bkg: ttjets Only", "LP");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);
  
  TPaveText * lumiHeader = new TPaveText(0.1, 0.901, 0.9, 0.94, "NDC");
  lumiHeader->SetFillColor(0);
  lumiHeader->SetFillStyle(0);
  lumiHeader->SetLineColor(0);
  lumiHeader->SetBorderSize(0);
  lumiHeader->AddText("CMS Preliminary 2015     #sqrt{s} = 8 TeV     #intL = 19.7 fb^{-1}");
  
  leg->Draw("same");
  lumiHeader->Draw("same");

  can->SaveAs("fluxMET_"+channel+".pdf");
  delete can;

}

void fluxMET() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  go("ele");
  go("muon");

}
