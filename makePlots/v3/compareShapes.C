void compareShapes() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  overlaySignalRegions("ele");
  overlaySignalRegions("muon");

  TString vars[6] = {"pfMET", "HT", "HT_jets", "dPhi_leadPhoton_l", "dPhi_met_l", "dPhi_met_ht"};

  for(int i = 0; i < 6; i++) {
    overlayChannels(vars[i], "CR0");
    overlayChannels(vars[i], "CR1");
    overlayChannels(vars[i], "SR1");
    overlayChannels(vars[i], "SR2");
  }

  overlayLeptonPt("CR0");
  overlayLeptonPt("CR1");
  overlayLeptonPt("SR1");
  overlayLeptonPt("SR2");

}

void overlaySignalRegions(TString channel) {

  TString names[9] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "qcd"};

  TFile * input = new TFile("limitInputs_bjj.root", "READ");

  TH1D * h_sr1 = (TH1D*)input->Get(channel+"_SR1/ttjets");
  TH1D * h_sr2 = (TH1D*)input->Get(channel+"_SR2/ttjets");

  for(int i = 1; i < 9; i++) {

    if(names[i] != "ttgamma") continue;

    TH1D * h1 = (TH1D*)input->Get(channel+"_SR1/"+names[i]);
    TH1D * h2 = (TH1D*)input->Get(channel+"_SR2/"+names[i]);

    if(h1->Integral() > 0. && h2->Integral() > 0.) {
      h_sr1->Add(h1);
      h_sr2->Add(h2);
    }
    
  }

  TH1D * h_sr1_sys = (TH1D*)h_sr1->Clone("sr1_sys");

  h_sr1->Scale(1./h_sr1->Integral());
  h_sr1_sys->Scale(1./h_sr1_sys->Integral());
  h_sr2->Scale(1./h_sr2->Integral());

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);
  padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetBottomMargin(0);

  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);

  padhi->Draw();
  padlo->Draw();

  padhi->cd();

  h_sr1_sys->SetFillColor(kOrange+10);
  h_sr1_sys->SetFillStyle(3154);
  h_sr1_sys->SetMarkerSize(0);

  h_sr2->SetMarkerStyle(20);
  h_sr2->SetMarkerSize(1.5);

  h_sr1_sys->GetYaxis()->SetRangeUser(0., 0.5);

  h_sr1_sys->Draw("e2");
  h_sr1->Draw("hist same");
  h_sr2->Draw("e1 same");
  
  TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, channel.Data(), "brNDC");
  leg->AddEntry(h_sr1_sys, "SR1", "LPF");
  leg->AddEntry(h_sr2, "SR2", "LP");
  leg->Draw("same");

  padlo->cd();

  TH1D * ratio = (TH1D*)h_sr2->Clone("ratio");
  ratio->Divide(h_sr1);

  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");

  TH1D * ratio_sys = (TH1D*)h_sr1_sys->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    if(h_sr1->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, h_sr1_sys->GetBinError(i+1) / h_sr1->GetBinContent(i+1));
    ratio_sys->SetBinContent(i+1, 1.);
  }

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("plot_"+channel+".png");

  input->Close();

}

void overlayChannels(TString variable, TString mode) {

  if(mode == "CR0" && variable.Contains("Photon")) return;

  TString names[24] = {"ttJetsSemiLep", "ttJetsFullLep", "ttJetsHadronic",
		       "W1JetsToLNu", "W2JetsToLNu", "W3JetsToLNu", "W4JetsToLNu",
		       "dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
		       "TBar_s", "TBar_t", "TBar_tW", "T_s", "T_t", "T_tW",
		       "WW", "WZ", "ZZ",
		       "TTWJets", "TTZJets",
		       "ttA_2to5",
		       "qcd"};

  TFile * inputEle = new TFile("histograms_ele_bjj_"+mode+".root", "READ");
  TFile * inputMuon = new TFile("histograms_muon_bjj_"+mode+".root", "READ");

  TH1D * h_ele = (TH1D*)inputEle->Get(variable+"_ttJetsSemiLep_ele_bjj");
  TH1D * h_muon = (TH1D*)inputMuon->Get(variable+"_ttJetsSemiLep_muon_bjj");

  for(int i = 1; i < 24; i++) {

    if(!(names[i].Contains("ttJets")) || !(names[i].Contains("ttA"))) continue;

    TH1D * he = (TH1D*)inputEle->Get(variable+"_"+names[i]+"_ele_bjj");
    TH1D * hm = (TH1D*)inputMuon->Get(variable+"_"+names[i]+"_muon_bjj");

    double n_ele, n_muon;
    double err_ele, err_muon;
    
    n_ele = he->IntegralAndError(1, -1, err_ele);
    n_muon = hm->IntegralAndError(1, -1, err_muon);

    if(err_ele > n_ele || err_muon > n_muon) continue;

    if(n_ele > 0. && n_muon > 0.) {
      h_ele->Add(he);
      h_ele->Add(hm);
    }

  }

  if(variable.Contains("dPhi")) {
    h_ele->Rebin(3);
    h_muon->Rebin(3);
  }

  TH1D * h_ele_sys = (TH1D*)h_ele->Clone("ele_sys");

  h_ele->Scale(1./h_ele->Integral());
  h_ele_sys->Scale(1./h_ele_sys->Integral());
  h_muon->Scale(1./h_muon->Integral());

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);
  padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetBottomMargin(0);
  //if(variable == "pfMET" || variable == "HT" || variable == "HT_jets") padhi->SetLogy(true);

  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);

  padhi->Draw();
  padlo->Draw();

  padhi->cd();

  h_ele_sys->SetFillColor(kOrange+10);
  h_ele_sys->SetFillStyle(3154);
  h_ele_sys->SetMarkerSize(0);

  h_muon->SetMarkerStyle(20);
  h_muon->SetMarkerSize(1.5);

  if(variable.Contains("dPhi")) h_ele_sys->GetYaxis()->SetRangeUser(0., 0.15);
  if(variable.Contains("HT") && mode == "SR2") h_ele_sys->GetYaxis()->SetRangeUser(0., 0.7);
  else h_ele_sys->GetYaxis()->SetRangeUser(0., 0.5);  

  h_ele_sys->Draw("e2");
  h_ele->Draw("hist same");
  h_muon->Draw("e1 same");
  
  TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, mode.Data(), "brNDC");
  leg->AddEntry(h_ele_sys, "ele", "LPF");
  leg->AddEntry(h_muon, "muon", "LP");
  leg->Draw("same");

  padlo->cd();

  TH1D * ratio = (TH1D*)h_muon->Clone("ratio");
  ratio->Divide(h_ele);

  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");

  TH1D * ratio_sys = (TH1D*)h_ele_sys->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    if(h_ele->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, h_ele_sys->GetBinError(i+1) / h_ele->GetBinContent(i+1));
    ratio_sys->SetBinContent(i+1, 1.);
  }

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("plot_"+variable+"_"+mode+".png");

  inputEle->Close();
  inputMuon->Close();

}

void overlayLeptonPt(TString mode) {

  TString names[24] = {"ttJetsSemiLep", "ttJetsFullLep", "ttJetsHadronic",
		       "W1JetsToLNu", "W2JetsToLNu", "W3JetsToLNu", "W4JetsToLNu",
		       "dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
		       "TBar_s", "TBar_t", "TBar_tW", "T_s", "T_t", "T_tW",
		       "WW", "WZ", "ZZ",
		       "TTWJets", "TTZJets",
		       "ttA_2to5",
		       "qcd"};

  TFile * inputEle = new TFile("histograms_ele_bjj_"+mode+".root", "READ");
  TFile * inputMuon = new TFile("histograms_muon_bjj_"+mode+".root", "READ");

  TH1D * h_ele = (TH1D*)inputEle->Get("ele_pt_ttJetsSemiLep_ele_bjj");
  TH1D * h_muon = (TH1D*)inputMuon->Get("muon_pt_ttJetsSemiLep_muon_bjj");

  for(int i = 1; i < 24; i++) {

    if(!(names[i].Contains("ttJets")) || !(names[i].Contains("ttA"))) continue;

    TH1D * he = (TH1D*)inputEle->Get("ele_pt_"+names[i]+"_ele_bjj");
    TH1D * hm = (TH1D*)inputMuon->Get("muon_pt_"+names[i]+"_muon_bjj");

    double n_ele, n_muon;
    double err_ele, err_muon;
    
    n_ele = he->IntegralAndError(1, -1, err_ele);
    n_muon = hm->IntegralAndError(1, -1, err_muon);

    if(err_ele > n_ele || err_muon > n_muon) continue;

    if(n_ele > 0. && n_muon > 0.) {
      h_ele->Add(he);
      h_ele->Add(hm);
    }

  }

  TH1D * h_ele_sys = (TH1D*)h_ele->Clone("ele_sys");

  h_ele->Scale(1./h_ele->Integral());
  h_ele_sys->Scale(1./h_ele_sys->Integral());
  h_muon->Scale(1./h_muon->Integral());

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);
  padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetBottomMargin(0);
  padhi->SetLogy(true);

  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);

  padhi->Draw();
  padlo->Draw();

  padhi->cd();

  h_ele_sys->SetFillColor(kOrange+10);
  h_ele_sys->SetFillStyle(3154);
  h_ele_sys->SetMarkerSize(0);

  h_muon->SetMarkerStyle(20);
  h_muon->SetMarkerSize(1.5);

  h_ele_sys->GetXaxis()->SetRangeUser(0, 1000);
  //h_ele_sys->GetYaxis()->SetRangeUser(0., 0.5);  

  h_ele_sys->Draw("e2");
  h_ele->Draw("hist same");
  h_muon->Draw("e1 same");
  
  TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, mode.Data(), "brNDC");
  leg->AddEntry(h_ele_sys, "ele", "LPF");
  leg->AddEntry(h_muon, "muon", "LP");
  leg->Draw("same");

  padlo->cd();

  TH1D * ratio = (TH1D*)h_muon->Clone("ratio");
  ratio->Divide(h_ele);

  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");
  ratio->GetXaxis()->SetRangeUser(0, 1000);

  TH1D * ratio_sys = (TH1D*)h_ele_sys->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    if(h_ele->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, h_ele_sys->GetBinError(i+1) / h_ele->GetBinContent(i+1));
    ratio_sys->SetBinContent(i+1, 1.);
  }

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("plot_lepton_pt_"+mode+".png");

  inputEle->Close();
  inputMuon->Close();

}
