const int nBackgrounds = 9;
TString names[nBackgrounds] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "qcd"};

void finishInputs() {

  TFile * fInputs = new TFile("limitInputs_bjj.root", "UPDATE");
  TFile * fDifferences = new TFile("met_differences.root", "READ");

  TH1D * diff_ele_cr1 = (TH1D*)fDifferences->Get("ele_bjj_CR1");
  TH1D * diff_muon_cr1 = (TH1D*)fDifferences->Get("muon_bjj_CR1");

  for(int i = 0; i < nBackgrounds; i++) {

    TH1D * h_ele_sr1 = (TH1D*)fInputs->Get("ele_SR1/"+names[i]);
    TH1D * h_ele_sr1_up = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_extraSystematicUp");
    TH1D * h_ele_sr1_down = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_extraSystematicDown");
    
    TH1D * h_muon_sr1 = (TH1D*)fInputs->Get("muon_SR1/"+names[i]);
    TH1D * h_muon_sr1_up = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_extraSystematicUp");
    TH1D * h_muon_sr1_down = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_extraSystematicDown");
    
    TH1D * h_ele_sr2 = (TH1D*)fInputs->Get("ele_SR2/"+names[i]);
    TH1D * h_ele_sr2_up = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_extraSystematicUp");
    TH1D * h_ele_sr2_down = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_extraSystematicDown");

    TH1D * h_muon_sr2 = (TH1D*)fInputs->Get("muon_SR2/"+names[i]);
    TH1D * h_muon_sr2_up = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_extraSystematicUp");
    TH1D * h_muon_sr2_down = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_extraSystematicDown");
    
    for(int j = 0; j < diff_ele_cr1->GetNbinsX(); j++) {
      double ratio_ele = diff_ele_cr1->GetBinContent(j+1);
      double ratio_muon = diff_muon_cr1->GetBinContent(j+1);
      
      // sr1 skips 0-20 bin
      // so sr1(1) gets cr1(2) for its systematic
      if(j > 0) {
	h_ele_sr1_up->SetBinContent(j, h_ele_sr1_up->GetBinContent(j) * (1. + ratio_ele));
	h_ele_sr1_down->SetBinContent(j, h_ele_sr1_down->GetBinContent(j) * (1. - ratio_ele));
      
	h_muon_sr1_up->SetBinContent(j, h_muon_sr1_up->GetBinContent(j) * (1. + ratio_muon));
	h_muon_sr1_down->SetBinContent(j, h_muon_sr1_down->GetBinContent(j) * (1. - ratio_muon));
      }

      h_ele_sr2_up->SetBinContent(j+1, h_ele_sr2_up->GetBinContent(j+1) * (1. + ratio_ele));
      h_ele_sr2_down->SetBinContent(j+1, h_ele_sr2_down->GetBinContent(j+1) * (1. - ratio_ele));
      
      h_muon_sr2_up->SetBinContent(j+1, h_muon_sr2_up->GetBinContent(j+1) * (1. + ratio_muon));
      h_muon_sr2_down->SetBinContent(j+1, h_muon_sr2_down->GetBinContent(j+1) * (1. - ratio_muon));
      
    }
    
    fInputs->cd("ele_SR1");
    h_ele_sr1_up->Write(names[i]+"_extraSystematicUp");
    h_ele_sr1_down->Write(names[i]+"_extraSystematicDown");
    
    fInputs->cd("muon_SR1");
    h_muon_sr1_up->Write(names[i]+"_extraSystematicUp");
    h_muon_sr1_down->Write(names[i]+"_extraSystematicDown");

    fInputs->cd("ele_SR2");
    h_ele_sr2_up->Write(names[i]+"_extraSystematicUp");
    h_ele_sr2_down->Write(names[i]+"_extraSystematicDown");
    
    fInputs->cd("muon_SR2");
    h_muon_sr2_up->Write(names[i]+"_extraSystematicUp");
    h_muon_sr2_down->Write(names[i]+"_extraSystematicDown");
    
  }

  fInputs->Close();
  fDifferences->Close();
}
