void plot(TH1D * ele, TH1D * muon, TString name) {

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);

  TH1D * ele_bars = (TH1D*)ele->Clone("ele_bars_"+name);
  ele_bars->SetFillColor(kOrange+10);
  ele_bars->SetFillStyle(3154);
  ele_bars->SetMarkerSize(0);

  ele->SetLineColor(kBlack);
  ele->GetXaxis()->SetRangeUser(0, 400);

  muon->SetLineColor(kRed);
  muon->SetLineWidth(2);
  muon->SetMarkerColor(kRed);

  TLegend * leg = new TLegend(0.7, 0.6, 0.85, 0.85, "Lepton P_{T}", "brNDC");
  leg->AddEntry(ele, "ele", "LP");
  leg->AddEntry(ele_bars, "ele stat.", "F");
  leg->AddEntry(muon, "muon", "LP");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  ele->Draw("hist");
  ele_bars->Draw("e2 same");
  muon->Draw("e1 same");
  leg->Draw("same");

  can->SaveAs(name+".pdf");

  delete can;

}

void plotRainbow(TH1D * ele_ttjets, TH1D * ele_wjets, TH1D * ele_zjets, TH1D * ele_singleTop, TH1D * ele_diboson, TH1D * ele_ttw, TH1D * ele_ttgamma,
		 TH1D * muon, TString name) {

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
  can->SetLogy(true);

  ele_ttjets->Add(ele_wjets); // just mc[0] didn't have everything else added

  ele_ttjets->SetFillColor(kGray);
  ele_ttjets->SetMarkerSize(0);
  ele_ttjets->SetLineColor(1);
  
  ele_wjets->SetFillColor(kOrange-3);
  ele_wjets->SetMarkerSize(0);
  ele_wjets->SetLineColor(1);

  ele_zjets->SetFillColor(kYellow);
  ele_zjets->SetMarkerSize(0);
  ele_zjets->SetLineColor(1);

  ele_singleTop->SetFillColor(kRed);
  ele_singleTop->SetMarkerSize(0);
  ele_singleTop->SetLineColor(1);

  ele_diboson->SetFillColor(kViolet-2);
  ele_diboson->SetMarkerSize(0);
  ele_diboson->SetLineColor(1);

  ele_ttw->SetFillColor(kCyan);
  ele_ttw->SetMarkerSize(0);
  ele_ttw->SetLineColor(1);

  ele_ttgamma->SetFillColor(8);
  ele_ttgamma->SetMarkerSize(0);
  ele_ttgamma->SetLineColor(1);

  TH1D * ele_bars = (TH1D*)ele_ttjets->Clone("ele_bars_"+name);
  ele_bars->SetFillColor(kOrange+10);
  ele_bars->SetFillStyle(3154);
  ele_bars->SetMarkerSize(0);
  
  muon->SetLineColor(kBlack);
  muon->SetMarkerColor(kBlack);
  muon->SetMarkerStyle(20);
  muon->SetMarkerSize(1.5);
  
  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, "MC Backgrounds", "brNDC");
  leg->SetNColumns(2);
  leg->AddEntry(muon, "Muon Total", "LP");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(ele_bars, "Electron Total", "F");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(ele_ttjets, "t#bar{t} + Jets", "F");
  leg->AddEntry(ele_wjets, "W + Jets", "F");
  leg->AddEntry(ele_zjets, "Z/#gamma* + Jets", "F");
  leg->AddEntry(ele_singleTop, "Single Top", "F");
  leg->AddEntry(ele_diboson, "VV, V#gamma", "F");
  leg->AddEntry(ele_ttw, "t#bar{t} + W/Z", "F");
  leg->AddEntry(ele_ttgamma, "t#bar{t} + #gamma", "F");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  ele_ttjets->GetXaxis()->SetRangeUser(0., 399.);
  ele_ttjets->GetXaxis()->SetTitle("Lepton P_{T}");

  ele_ttjets->Draw("hist");
  ele_wjets->Draw("hist same");
  ele_zjets->Draw("hist same");
  ele_singleTop->Draw("hist same");
  ele_diboson->Draw("hist same");
  ele_ttw->Draw("hist same");
  ele_ttgamma->Draw("hist same");

  ele_bars->Draw("e2 same");
  muon->Draw("e1 same");
  ele_ttjets->Draw("axis same");
  leg->Draw("same");

  can->SaveAs(name+".pdf");

  delete can;

}

void plotLeptonPtComparison() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TString channel = "bjj_SR2";

  TFile * input = new TFile("../durp_leptonPt.root", "READ");

  TH1D * bkg_ele = (TH1D*)input->Get("bkg_ele_"+channel);
  TH1D * bkg_muon = (TH1D*)input->Get("bkg_muon_"+channel);

  TH1D * ttjets_ele = (TH1D*)input->Get("ttjets_ele_"+channel);
  TH1D * wjets_ele = (TH1D*)input->Get("wjets_ele_"+channel);
  TH1D * zjets_ele = (TH1D*)input->Get("zjets_ele_"+channel);
  TH1D * singleTop_ele = (TH1D*)input->Get("singleTop_ele_"+channel);
  TH1D * diboson_ele = (TH1D*)input->Get("diboson_ele_"+channel);
  TH1D * ttw_ele = (TH1D*)input->Get("ttW_ele_"+channel);
  TH1D * ttgamma_ele = (TH1D*)input->Get("ttgamma_ele_"+channel);

  TH1D * ttgamma_muon = (TH1D*)input->Get("ttgamma_muon_"+channel);

  plot(bkg_ele, bkg_muon, "leptonPt_SR2_totalBkg");
  plot(ttgamma_ele, ttgamma_muon, "leptonPt_SR2_ttgamma");
  plotRainbow(ttjets_ele, wjets_ele, zjets_ele, singleTop_ele, diboson_ele, ttw_ele, ttgamma_ele,
	      bkg_muon, "leptonPt_SR2_tasteTheRainbow");
  
  input->Close();

}

