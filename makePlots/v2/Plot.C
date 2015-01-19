void Plot() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TString variable = "dR_leadPhoton_l";

  TFile * f_data_gamma = new TFile("data_ele_bjj_gamma.root", "READ");
  TFile * f_data_fake = new TFile("data_ele_bjj_fake.root", "READ");

  TFile * f_ttA_gamma = new TFile("ttA_2to5_ele_bjj_gamma.root", "READ");
  TFile * f_ttA_fake = new TFile("ttA_2to5_ele_bjj_fake.root", "READ");

  TFile * f_ttJetsSemiLep_gamma = new TFile("ttJetsSemiLep_ele_bjj_gamma.root", "READ");
  TFile * f_ttJetsFullLep_gamma = new TFile("ttJetsFullLep_ele_bjj_gamma.root", "READ");
  TFile * f_ttJetsHadronic_gamma = new TFile("ttJetsHadronic_ele_bjj_gamma.root", "READ");

  TFile * f_ttJetsSemiLep_fake = new TFile("ttJetsSemiLep_ele_bjj_fake.root", "READ");
  TFile * f_ttJetsFullLep_fake = new TFile("ttJetsFullLep_ele_bjj_fake.root", "READ");
  TFile * f_ttJetsHadronic_fake = new TFile("ttJetsHadronic_ele_bjj_fake.root", "READ");
 
  TH1D * h_data_gamma = (TH1D*)f_data_gamma->Get(variable);
  TH1D * h_data_fake = (TH1D*)f_data_fake->Get(variable);

  TH1D * h_ttA_gamma = (TH1D*)f_ttA_gamma->Get(variable);
  TH1D * h_ttA_fake = (TH1D*)f_ttA_fake->Get(variable);

  TH1D * h_ttJets_gamma = (TH1D*)f_ttJetsSemiLep_gamma->Get(variable);
  h_ttJets_gamma->Add((TH1D*)f_ttJetsFullLep_gamma->Get(variable));
  h_ttJets_gamma->Add((TH1D*)f_ttJetsHadronic_gamma->Get(variable));
  h_ttJets_gamma->Add(h_ttA_gamma);

  TH1D * h_ttJets_fake = (TH1D*)f_ttJetsSemiLep_fake->Get(variable);
  h_ttJets_fake->Add((TH1D*)f_ttJetsFullLep_fake->Get(variable));
  h_ttJets_fake->Add((TH1D*)f_ttJetsHadronic_fake->Get(variable));
  h_ttJets_fake->Add(h_ttA_fake);

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 800, 800);
  //can->SetLogy(true);

  h_data_gamma->Rebin(5);
  h_data_fake->Rebin(5);

  h_ttJets_gamma->Rebin(5);
  h_ttJets_fake->Rebin(5);
  h_ttA_gamma->Rebin(5);
  h_ttA_fake->Rebin(5);

  h_data_gamma->Scale(1./h_data_gamma->Integral());
  h_data_fake->Scale(1./h_data_fake->Integral());

  h_ttJets_gamma->SetLineColor(kGray);
  h_ttJets_gamma->Scale(1./h_ttJets_gamma->Integral());
  h_ttJets_fake->SetLineColor(kGray);
  h_ttJets_fake->Scale(1./h_ttJets_fake->Integral());
  
  h_ttA_gamma->SetLineColor(8);
  h_ttA_gamma->Scale(1./h_ttA_gamma->Integral());
  h_ttA_fake->SetLineColor(8);
  h_ttA_fake->Scale(1./h_ttA_fake->Integral());


  h_data_gamma->Draw("e1");
  h_ttJets_gamma->Draw("hist same");
  h_ttA_gamma->Draw("hist same");
  can->SaveAs("gamma.png");

  h_data_fake->Draw("e1");
  h_ttA_fake->Draw("hist same");
  h_ttJets_fake->Draw("hist same");
  can->SaveAs("fake.png");

  h_ttJets_gamma->SetLineColor(kBlue);
  h_ttJets_fake->SetLineColor(kRed);
  h_ttJets_gamma->Draw("hist");
  h_ttJets_fake->Draw("hist same");
  can->SaveAs("ttJets.png");

}
