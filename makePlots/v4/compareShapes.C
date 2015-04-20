void compareShapes() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  /*
  overlayCR1SR1("ele");
  overlayCR1SR1("muon");

  overlayCR2SR2("ele");
  overlayCR2SR2("muon");
  */

  overlayThing("ele");
  overlayThing("muon");
}

void overlayCR1SR1(TString channel) {

  TString names[10] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

  TFile * input = new TFile("limitInputs_bjj.root", "READ");

  TH1D * h_cr1 = (TH1D*)input->Get(channel+"_CR1/ttjets");
  TH1D * h_sr1 = (TH1D*)input->Get(channel+"_SR1/ttjets");

  for(int i = 1; i < 10; i++) {

    if(channel == "ele") {
      // the e-->g fake rate SFs are measured in SR1 so this is borked
      //if(names[i] == "zjets" || names[i] == "vgamma") continue;
    }
    if(channel == "muon") {
      //if(names[i] == "zjets" || names[i] == "vgamma") continue;
      //if(names[i] == "diboson" || names[i] == "ttW" || names[i] == "ttZ" || names[i] == "vgamma") continue;
    }

    TH1D * hc = (TH1D*)input->Get(channel+"_CR1/"+names[i]);
    TH1D * h1 = (TH1D*)input->Get(channel+"_SR1/"+names[i]);

    double nc, n1, ec, e1;
    nc = hc->IntegralAndError(0, -1, ec);
    n1 = h1->IntegralAndError(0, -1, e1);

    if(nc < 1.e-6 || n1 < 1.e-6 || ec > nc || e1 > n1) {cout << "skipping negligible " << names[i] << endl; continue;}

    h_cr1->Add(hc);
    h_sr1->Add(h1);
    
  }

  h_cr1->Scale(1./h_cr1->Integral());
  h_sr1->Scale(1./h_sr1->Integral());

  TH1D * h_cr1_sys = (TH1D*)h_cr1->Clone("cr1_sys");
  TH1D * h_sr1_sys = (TH1D*)h_sr1->Clone("sr1_sys");

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

  h_cr1_sys->SetFillColor(kOrange+10);
  h_cr1_sys->SetFillStyle(3154);
  h_cr1_sys->SetMarkerSize(0);
  //h_cr1_sys->GetYaxis()->SetRangeUser(0., 0.5);

  h_sr1->SetMarkerStyle(20);
  h_sr1->SetMarkerSize(1.5);

  h_cr1_sys->Draw("e2");
  h_cr1->Draw("hist same");
  h_sr1->Draw("e1 same");
  
  TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, channel.Data(), "brNDC");
  leg->AddEntry(h_cr1_sys, "CR1", "LPF");
  leg->AddEntry(h_sr1, "SR1", "LP");
  leg->Draw("same");

  padlo->cd();

  TH1D * ratio = (TH1D*)h_sr1->Clone("ratio");
  ratio->Divide(h_cr1);

  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetTitle("MET");

  TH1D * ratio_sys = (TH1D*)h_cr1_sys->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    if(h_cr1->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, h_cr1_sys->GetBinError(i+1) / h_cr1->GetBinContent(i+1));
    ratio_sys->SetBinContent(i+1, 1.);
  }

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("cr1_to_sr1_"+channel+".png");

  TFile * output = new TFile("extraErrors.root", "UPDATE");
  ratio->Write("cr1_to_sr1_"+channel);
  output->Close();

  input->Close();

}

void overlayCR2SR2(TString channel) {

  TString names[10] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

  TFile * input = new TFile("limitInputs_bjj.root", "READ");

  TH1D * h_cr2 = (TH1D*)input->Get(channel+"_CR2/ttjets");
  TH1D * h_sr2 = (TH1D*)input->Get(channel+"_SR2/ttjets");

  for(int i = 1; i < 10; i++) {

//if(names[i] != "ttjets" && names[i] != "ttgamma") continue;

    TH1D * hc = (TH1D*)input->Get(channel+"_CR2/"+names[i]);
    TH1D * h2 = (TH1D*)input->Get(channel+"_SR2/"+names[i]);

    double nc, n2, ec, e2;
    nc = hc->IntegralAndError(0, -1, ec);
    n2 = h2->IntegralAndError(0, -1, e2);

    if(nc < 1.e-6 || n2 < 1.e-6) continue;
    if(ec > nc || e2 > n2) continue;

    h_cr2->Add(hc);
    h_sr2->Add(h2);
    
  }

  h_cr2->Scale(1./h_cr2->Integral());
  h_sr2->Scale(1./h_sr2->Integral());

  TH1D * h_cr2_sys = (TH1D*)h_cr2->Clone("cr2_sys");
  TH1D * h_sr2_sys = (TH1D*)h_sr2->Clone("sr2_sys");

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

  h_cr2_sys->SetFillColor(kOrange+10);
  h_cr2_sys->SetFillStyle(3154);
  h_cr2_sys->SetMarkerSize(0);
  //h_cr2_sys->GetYaxis()->SetRangeUser(0., 0.5);

  h_sr2->SetMarkerStyle(20);
  h_sr2->SetMarkerSize(1.5);

  h_cr2_sys->Draw("e2");
  h_cr2->Draw("hist same");
  h_sr2->Draw("e1 same");
  
  TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, channel.Data(), "brNDC");
  leg->AddEntry(h_cr2_sys, "CR2", "LPF");
  leg->AddEntry(h_sr2, "SR2", "LP");
  leg->Draw("same");

  padlo->cd();

  TH1D * ratio = (TH1D*)h_sr2->Clone("ratio");
  ratio->Divide(h_cr2);

  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetTitle("MET");

  TH1D * ratio_sys = (TH1D*)h_cr2_sys->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    if(h_cr2->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, h_cr2_sys->GetBinError(i+1) / h_cr2->GetBinContent(i+1));
    ratio_sys->SetBinContent(i+1, 1.);
  }

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("cr2_to_sr2_"+channel+".png");

  TFile * output = new TFile("extraErrors.root", "UPDATE");
  ratio->Write("cr2_to_sr2_"+channel);
  output->Close();

  input->Close();

}

void overlaySRSR(TString channel) {

  TString names[10] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

  TFile * input = new TFile("limitInputs_bjj.root", "READ");

  TH1D * h_cr1 = (TH1D*)input->Get(channel+"_CR1/ttjets");
  TH1D * h_sr1 = (TH1D*)input->Get(channel+"_SR1/ttjets");
  TH1D * h_sr2 = (TH1D*)input->Get(channel+"_SR2/ttjets");

  for(int i = 1; i < 10; i++) {

    if(names[i] != "ttgamma") continue;

    TH1D * hc = (TH1D*)input->Get(channel+"_CR1/"+names[i]);
    TH1D * h1 = (TH1D*)input->Get(channel+"_SR1/"+names[i]);
    TH1D * h2 = (TH1D*)input->Get(channel+"_SR2/"+names[i]);

    double nc, n1, n2, ec, e1, e2;
    nc = hc->IntegralAndError(0, -1, ec);
    n1 = h1->IntegralAndError(0, -1, e1);
    n2 = h2->IntegralAndError(0, -1, e2);

    if(nc < 1.e-6 || n1 < 1.e-6 || n2 < 1.e-6) continue;
    if(ec > nc || e1 > n1 || e2 > n2) continue;

    h_cr1->Add(hc);
    h_sr1->Add(h1);
    h_sr2->Add(h2);
    
  }

  TH1D * h_cr1_sys = (TH1D*)h_cr1->Clone("cr1_sys");
  TH1D * h_sr1_sys = (TH1D*)h_sr1->Clone("sr1_sys");

  h_cr1->Scale(1./h_cr1->Integral());
  h_cr1_sys->Scale(1./h_cr1_sys->Integral());
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
  ratio->GetXaxis()->SetTitle("MET");

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

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("sr1_to_sr2_"+channel+".png");

  TFile * output = new TFile("extraErrors.root", "UPDATE");
  ratio->Write("sr1_to_sr2_"+channel);
  output->Close();

  input->Close();

}

void overlayChannels(TString variable, TString mode) {

  if(mode == "CR0" && variable.Contains("Photon")) return;

  TString names[24] = {"ttJetsSemiLep", "ttJetsFullLep", "ttJetsHadronic",
		       "W3JetsToLNu", "W4JetsToLNu",
		       "dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
		       "TBar_s", "TBar_t", "TBar_tW", "T_s", "T_t", "T_tW",
		       "WW", "WZ", "ZZ",
		       "WGToLNuG", "ZGToLLG",
		       "TTWJets", "TTZJets",
		       "TTGamma",
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
      h_muon->Add(hm);
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
  //if(variable == "pfMET_t01" || variable == "HT" || variable == "HT_jets") padhi->SetLogy(true);

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
  ratio->GetXaxis()->SetTitle(variable);

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

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
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
		       "W3JetsToLNu", "W4JetsToLNu",
		       "dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
		       "TBar_s", "TBar_t", "TBar_tW", "T_s", "T_t", "T_tW",
		       "WW", "WZ", "ZZ",
		       "WGToLNuG", "ZGToLLG",
		       "TTWJets", "TTZJets",
		       "TTGamma",
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

  ratio->GetXaxis()->SetTitle("Lepton PT");
  ratio->GetXaxis()->SetRangeUser(0, 1000);

  ratio->SetMarkerStyle(1);
  ratio->SetMarkerSize(1);

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

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
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

void overlaySignalRegions(TString variable, TString channel) {

  TString names[24] = {"ttJetsSemiLep", "ttJetsFullLep", "ttJetsHadronic",
		       "W3JetsToLNu", "W4JetsToLNu",
		       "dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
		       "TBar_s", "TBar_t", "TBar_tW", "T_s", "T_t", "T_tW",
		       "WW", "WZ", "ZZ",
		       "WGToLNuG", "ZGToLLG",
		       "TTWJets", "TTZJets",
		       "TTGamma",
		       "qcd"};

  TFile * inputSR1 = new TFile("histograms_"+channel+"_bjj_SR1.root", "READ");
  TFile * inputSR2 = new TFile("histograms_"+channel+"_bjj_SR2.root", "READ");

  TH1D * h_sr1 = (TH1D*)inputSR1->Get(variable+"_ttJetsSemiLep_"+channel+"_bjj");
  TH1D * h_sr2 = (TH1D*)inputSR2->Get(variable+"_ttJetsSemiLep_"+channel+"_bjj");

  for(int i = 1; i < 24; i++) {

    if(!(names[i].Contains("ttJets")) || !(names[i].Contains("ttA"))) continue;

    TH1D * h1 = (TH1D*)inputSR1->Get(variable+"_"+names[i]+"_"+channel+"_bjj");
    TH1D * h2 = (TH1D*)inputSR2->Get(variable+"_"+names[i]+"_"+channel+"_bjj");

    double n_1, n_2;
    double err_1, err_2;
    
    n_1 = h1->IntegralAndError(1, -1, err_1);
    n_2 = h2->IntegralAndError(1, -1, err_2);

    if(err_1 > n_1 || err_2 > n_2) continue;

    if(n_1 > 0. && n_2 > 0.) {
      h_sr1->Add(h1);
      h_sr2->Add(h2);
    }

  }

  if(variable.Contains("dPhi")) {
    h_sr1->Rebin(3);
    h_sr2->Rebin(3);
  }

  TH1D * h_sr1_sys = (TH1D*)h_sr1->Clone("1_sys");

  h_sr1->Scale(1./h_sr1->Integral());
  h_sr1_sys->Scale(1./h_sr1_sys->Integral());
  h_sr2->Scale(1./h_sr2->Integral());

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);
  padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetBottomMargin(0);
  //if(variable == "pfMET_t01" || variable == "HT" || variable == "HT_jets") padhi->SetLogy(true);

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

  if(variable.Contains("dPhi")) h_sr1_sys->GetYaxis()->SetRangeUser(0., 0.15);
  if(variable.Contains("HT")) h_sr1_sys->GetYaxis()->SetRangeUser(0., 0.7);
  else h_sr1_sys->GetYaxis()->SetRangeUser(0., 0.5);  

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
  ratio->GetXaxis()->SetTitle(variable);

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

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("plot_"+variable+"_"+channel+".png");

  inputSR1->Close();
  inputSR2->Close();

}

void overlayThing(TString channel) {

  TString names[10] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

  TFile * input = new TFile("limitInputs_bjj.root", "READ");

  TH1D * h_mc = (TH1D*)input->Get(channel+"_SR1/ttjets");
  TH1D * h_data = (TH1D*)input->Get(channel+"_CR1/data_obs");

  for(int i = 1; i < 10; i++) {
    TH1D * h = (TH1D*)input->Get(channel+"_SR1/"+names[i]);

    double n, err;
    n = h->IntegralAndError(0, -1, err);

    if(n < 1.e-6) continue;
    if(err > n) continue;

    h_mc->Add(h);
  }

  h_mc->Scale(1./h_mc->Integral());
  h_data->Scale(1./h_data->Integral());

  TH1D * h_mc_sys = (TH1D*)h_mc->Clone("mc_sys");
  TH1D * h_data_sys = (TH1D*)h_data->Clone("data_sys");

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

  h_mc_sys->SetFillColor(kOrange+10);
  h_mc_sys->SetFillStyle(3154);
  h_mc_sys->SetMarkerSize(0);
  h_mc_sys->GetYaxis()->SetRangeUser(0., 0.5);

  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(1.5);

  h_mc_sys->Draw("e2");
  h_mc->Draw("hist same");
  h_data->Draw("e1 same");
  
  TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, channel.Data(), "brNDC");
  leg->AddEntry(h_mc_sys, "MC in SR1", "LPF");
  leg->AddEntry(h_data, "Data in CR1", "LP");
  leg->Draw("same");

  padlo->cd();

  TH1D * ratio = (TH1D*)h_data->Clone("ratio");
  ratio->Divide(h_mc);

  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetTitle("MET");

  TH1D * ratio_sys = (TH1D*)h_mc_sys->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    if(h_mc->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, h_mc_sys->GetBinError(i+1) / h_mc->GetBinContent(i+1));
    ratio_sys->SetBinContent(i+1, 1.);
  }

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("dataCR1_to_mcSR1_"+channel+".png");

  TFile * output = new TFile("extraErrors.root", "UPDATE");
  ratio->Write("dataCR1_to_mcSR1_"+channel);
  output->Close();

  input->Close();

}

void overlayClosure1(TString channel) {

  TString names[10] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

  TFile * input = new TFile("limitInputs_bjj.root", "READ");

  TH1D * h_mc = (TH1D*)input->Get(channel+"_SR1/ttjets");
  TH1D * h_data = (TH1D*)input->Get(channel+"_CR1/data_obs");

  for(int i = 1; i < 10; i++) {
    TH1D * h = (TH1D*)input->Get(channel+"_SR1/"+names[i]);

    double n, err;
    n = h->IntegralAndError(0, -1, err);

    if(n < 1.e-6) continue;
    if(err > n) continue;

    h_mc->Add(h);
  }

  h_mc->Scale(1./h_mc->Integral());
  h_data->Scale(1./h_data->Integral());

  TH1D * h_mc_sys = (TH1D*)h_mc->Clone("mc_sys");
  TH1D * h_data_sys = (TH1D*)h_data->Clone("data_sys");

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

  h_mc_sys->SetFillColor(kOrange+10);
  h_mc_sys->SetFillStyle(3154);
  h_mc_sys->SetMarkerSize(0);
  h_mc_sys->GetYaxis()->SetRangeUser(0., 0.5);

  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(1.5);

  h_mc_sys->Draw("e2");
  h_mc->Draw("hist same");
  h_data->Draw("e1 same");
  
  TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, channel.Data(), "brNDC");
  leg->AddEntry(h_mc_sys, "MC in SR1", "LPF");
  leg->AddEntry(h_data, "Data in CR1", "LP");
  leg->Draw("same");

  padlo->cd();

  TH1D * ratio = (TH1D*)h_data->Clone("ratio");
  ratio->Divide(h_mc);

  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetTitle("MET");

  TH1D * ratio_sys = (TH1D*)h_mc_sys->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    if(h_mc->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, h_mc_sys->GetBinError(i+1) / h_mc->GetBinContent(i+1));
    ratio_sys->SetBinContent(i+1, 1.);
  }

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  TLine * oneLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1, 
			      ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()) + ratio->GetXaxis()->GetBinWidth(ratio->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");
  oneLine->Draw();
  
  can->SaveAs("dataCR1_to_mcSR1_"+channel+".png");

  TFile * output = new TFile("extraErrors.root", "UPDATE");
  ratio->Write("dataCR1_to_mcSR1_"+channel);
  output->Close();

  input->Close();

}
