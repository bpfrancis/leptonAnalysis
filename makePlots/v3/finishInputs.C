const int nBackgrounds = 9;
TString names[nBackgrounds] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "qcd"};

void finishInputs() {

  TFile * fInputs = new TFile("limitInputs_bjj.root", "UPDATE");
  TFile * fDifferences = new TFile("met_differences.root", "READ");

  TH1D * diff_ele_cr1 = (TH1D*)fDifferences->Get("ele_bjj_CR1");
  TH1D * diff_muon_cr1 = (TH1D*)fDifferences->Get("muon_bjj_CR1");

  for(int i = 0; i < nBackgrounds; i++) {

    TH1D * h_ele_sr1 = (TH1D*)fInputs->Get("ele_SR1/"+names[i]);
    TH1D * h_ele_sr1_up = (TH1D*)h_ele_sr1->Clone(names[i]+"_extraSystematicUp");
    TH1D * h_ele_sr1_down = (TH1D*)h_ele_sr1->Clone(names[i]+"_extraSystematicDown");
    
    TH1D * h_muon_sr1 = (TH1D*)fInputs->Get("muon_SR1/"+names[i]);
    TH1D * h_muon_sr1_up = (TH1D*)h_ele_sr1->Clone(names[i]+"_extraSystematicUp");
    TH1D * h_muon_sr1_down = (TH1D*)h_ele_sr1->Clone(names[i]+"_extraSystematicDown");
    
    TH1D * h_ele_sr2 = (TH1D*)fInputs->Get("ele_SR2/"+names[i]);
    TH1D * h_ele_sr2_up = (TH1D*)h_ele_sr2->Clone(names[i]+"_extraSystematicUp");
    TH1D * h_ele_sr2_down = (TH1D*)h_ele_sr2->Clone(names[i]+"_extraSystematicDown");

    TH1D * h_muon_sr2 = (TH1D*)fInputs->Get("muon_SR2/"+names[i]);
    TH1D * h_muon_sr2_up = (TH1D*)h_ele_sr2->Clone(names[i]+"_extraSystematicUp");
    TH1D * h_muon_sr2_down = (TH1D*)h_ele_sr2->Clone(names[i]+"_extraSystematicDown");
    
    for(int j = 0; j < diff_ele_cr1->GetNbinsX(); j++) {
      double ratio = diff_ele_cr1->GetBinContent(j+1);
      
      h_ele_sr1_up->SetBinContent(j+1, h_ele_sr1_up->GetBinContent(j+1) * (1. + ratio));
      h_ele_sr1_down->SetBinContent(j+1, h_ele_sr1_down->GetBinContent(j+1) * (1. - ratio));
      
      h_muon_sr1_up->SetBinContent(j+1, h_muon_sr1_up->GetBinContent(j+1) * (1. + ratio));
      h_muon_sr1_down->SetBinContent(j+1, h_muon_sr1_down->GetBinContent(j+1) * (1. - ratio));
      
      h_ele_sr2_up->SetBinContent(j+1, h_ele_sr2_up->GetBinContent(j+1) * (1. + ratio));
      h_ele_sr2_down->SetBinContent(j+1, h_ele_sr2_down->GetBinContent(j+1) * (1. - ratio));
      
      h_muon_sr2_up->SetBinContent(j+1, h_muon_sr2_up->GetBinContent(j+1) * (1. + ratio));
      h_muon_sr2_down->SetBinContent(j+1, h_muon_sr2_down->GetBinContent(j+1) * (1. - ratio));
      
    }
    
    fInputs->cd("ele_SR1");
    h_ele_sr1_up->Write();
    h_ele_sr1_down->Write();
    
    fInputs->cd("muon_SR1");
    h_muon_sr1_up->Write();
    h_muon_sr1_down->Write();

    fInputs->cd("ele_SR2");
    h_ele_sr2_up->Write();
    h_ele_sr2_down->Write();
    
    fInputs->cd("muon_SR2");
    h_muon_sr2_up->Write();
    h_muon_sr2_down->Write();
    
  }

  fInputs->Close();
  fDifferences->Close();
}
