void Go(TString variable, bool logy, int rebin) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

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

  TCanvas * can = new TCanvas("canvas_"+variable, "Plot", 10, 10, 800, 800);
  can->SetLogy(logy);

  if(rebin > 0) {
    h_data_gamma->Rebin(rebin);
    h_data_fake->Rebin(rebin);
    
    h_ttJets_gamma->Rebin(rebin);
    h_ttJets_fake->Rebin(rebin);
    h_ttA_gamma->Rebin(rebin);
    h_ttA_fake->Rebin(rebin);
  }

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
  can->SaveAs("plots/gamma_"+variable+".png");

  h_data_fake->Draw("e1");
  h_ttA_fake->Draw("hist same");
  h_ttJets_fake->Draw("hist same");
  can->SaveAs("plots/fake_"+variable+".png");

  h_ttJets_gamma->SetLineColor(kBlue);
  h_ttJets_fake->SetLineColor(kRed);
  h_ttJets_gamma->Draw("e1");
  h_ttJets_fake->Draw("e1 same");
  can->SaveAs("plots/ttJets_"+variable+".png");

  cout << "KS ttJets_" << variable << " gamma/fake: " << h_ttJets_gamma->KolmogorovTest(h_ttJets_fake) << endl;

}

void Plot() {

  // Go(TString variable, bool logy)
  Go("pfMET", true, 0);
  Go("HT_jets", true, 0);
  Go("leadPhotonEta", false, 5);
  Go("leadPhotonEt", true, 2);
  Go("ele_pt", true, 2);
  Go("ele_eta", false, 5);


}
