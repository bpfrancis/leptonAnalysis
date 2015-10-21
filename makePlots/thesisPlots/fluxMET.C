void DivideByBinWidth(TH1D*& h) {

  for(Int_t i = 0; i < h->GetNbinsX(); i++) {
    Double_t val = h->GetBinContent(i+1);
    Double_t err = h->GetBinError(i+1);
    Double_t width = h->GetBinWidth(i+1);
    
    h->SetBinContent(i+1, val / width);
    h->SetBinError(i+1, err / width);
  }

}

void go(TString channel) {

  TFile * input = new TFile("../limitInputs_bjj.root", "READ");

  TString names[9] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma"};

  TH1D * ttgamma = (TH1D*)input->Get(channel+"_SR1/ttgamma");
  TH1D * ttjets = (TH1D*)input->Get(channel+"_SR1/ttjets");

  double total_ttx = ttgamma->Integral() + ttjets->Integral();

  TH1D * bkg = (TH1D*)ttjets->Clone("bkg");

  TH1D * bkg_ttgammaOnly = (TH1D*)ttgamma->Clone("bkg_ttgammaOnly");
  bkg_ttgammaOnly->Scale(total_ttx / ttgamma->Integral());

  TH1D * bkg_ttjetsOnly = (TH1D*)ttjets->Clone("bkg_ttjetsOnly");
  bkg_ttjetsOnly->Scale(total_ttx / ttjets->Integral());

  for(int i = 1; i < 9; i++) {
    if(names[i] == "ttgamma") continue;
    
    TH1D * h = (TH1D*)input->Get(channel+"_SR1/"+names[i]);
    bkg->Add(h);
    bkg_ttgammaOnly->Add(h);
    bkg_ttjetsOnly->Add(h);
  }
  bkg->Add(ttgamma);

  TH1D * bkg_syst = (TH1D*)bkg->Clone("bkg_syst");
  for(int i = 0; i < bkg_syst->GetNbinsX(); i++) {

    double flux_ttg = ttgamma->GetBinContent(i+1);
    double flux_ttj = ttjets->GetBinContent(i+1);

    if(channel.Contains("ele")) {
      flux_ttg *= 0.11;
      flux_ttj *= 0.13;
    }
    else {
      flux_ttg *= 0.16;
      flux_ttj *= 0.17;
    }
 
    bkg_syst->SetBinError(i+1, sqrt(flux_ttg*flux_ttg + flux_ttj*flux_ttj));
  }

  TH1D * ratio_ttgammaOnly = (TH1D*)bkg_ttgammaOnly->Clone("ratio_ttgOnly");
  for(int i = 0; i < ratio_ttgammaOnly->GetNbinsX(); i++) {
    ratio_ttgammaOnly->SetBinContent(i+1, bkg_ttgammaOnly->GetBinContent(i+1) / bkg->GetBinContent(i+1));
    ratio_ttgammaOnly->SetBinError(i+1, bkg_ttgammaOnly->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  TH1D * ratio_ttjetsOnly = (TH1D*)bkg_ttjetsOnly->Clone("ratio_ttjOnly");
  for(int i = 0; i < ratio_ttjetsOnly->GetNbinsX(); i++) {
    ratio_ttjetsOnly->SetBinContent(i+1, bkg_ttjetsOnly->GetBinContent(i+1) / bkg->GetBinContent(i+1));
    ratio_ttjetsOnly->SetBinError(i+1, bkg_ttjetsOnly->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  TH1D * ratio_syst = (TH1D*)bkg_syst->Clone("ratio_syst");
  for(int i = 0; i < ratio_syst->GetNbinsX(); i++) {
    ratio_syst->SetBinContent(i+1, 1.);
    ratio_syst->SetBinError(i+1, bkg_syst->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  TCanvas * can = new TCanvas("can", "can", 10, 10, 800, 800);
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

  DivideByBinWidth(bkg);
  DivideByBinWidth(bkg_syst);
  DivideByBinWidth(bkg_ttgammaOnly);
  DivideByBinWidth(bkg_ttjetsOnly);

  padhi->cd();

  bkg_syst->SetFillColor(kOrange+10);
  bkg_syst->SetFillStyle(3154);
  bkg_syst->SetMarkerSize(0);

  bkg->SetLineColor(kBlack);
  //bkg->SetLineWidth(3);
  bkg->GetXaxis()->SetTitle("MET");
  bkg->GetYaxis()->SetTitle("Events / GeV");
  //bkg->GetYaxis()->SetRangeUser(9.e-4, 5.e-1);

  bkg_ttgammaOnly->SetLineColor(8);
  //bkg_ttgammaOnly->SetLineWidth(3);
  bkg_ttgammaOnly->SetMarkerSize(0);
  
  bkg_ttjetsOnly->SetLineColor(kBlue);
  //bkg_ttjetsOnly->SetLineWidth(3);
  bkg_ttjetsOnly->SetMarkerSize(0);

  bkg->Draw("hist");
  bkg_syst->Draw("e2 same");
  bkg_ttgammaOnly->Draw("e1 same");
  bkg_ttjetsOnly->Draw("e1 same");
  bkg->Draw("hist same");
  bkg->Draw("axis same");

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, channel+" (SR1)", "brNDC");
  leg->AddEntry(bkg, "Background", "LP");
  leg->AddEntry(bkg_syst, "Purity fit systematic", "F");
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

  padlo->cd();

  ratio_syst->GetXaxis()->SetTitle("MET (GeV)");
  ratio_syst->GetYaxis()->SetRangeUser(0.5, 1.5);
  ratio_syst->SetFillStyle(1001);
  ratio_syst->SetFillColor(kGray);
  ratio_syst->SetLineColor(kGray);
  ratio_syst->SetMarkerColor(kGray);
  ratio_syst->Draw("e2");

  TLine * oneLine = new TLine(0, 1, 800, 1);
  oneLine->SetLineStyle(2);
  oneLine->Draw("same");
  
  ratio_ttgammaOnly->SetLineColor(8);
  ratio_ttgammaOnly->SetMarkerColor(8);
  ratio_ttgammaOnly->Draw("e1 same");

  ratio_ttjetsOnly->SetLineColor(kBlue);
  ratio_ttgammaOnly->SetMarkerColor(kBlue);
  ratio_ttjetsOnly->Draw("e1 same");

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
