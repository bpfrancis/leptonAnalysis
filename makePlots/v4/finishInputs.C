const int nBackgrounds = 10;
TString names[nBackgrounds] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

void finishInputs() {

  TFile * fInputs = new TFile("limitInputs_bjj.root", "UPDATE");
  TFile * fDifferences = new TFile("met_differences.root", "READ");
  TFile * fOthers = new TFile("extraErrors.root", "READ");

  TH1D * diff_ele_cr1 = (TH1D*)fDifferences->Get("ele_bjj_CR1");
  TH1D * diff_muon_cr1 = (TH1D*)fDifferences->Get("muon_bjj_CR1");

  TH1D * diff_ele_cr2 = (TH1D*)fDifferences->Get("ele_bjj_CR2");
  TH1D * diff_muon_cr2 = (TH1D*)fDifferences->Get("muon_bjj_CR2");

  TH1D * convert_ele_sr1 = (TH1D*)fOthers->Get("cr1_to_sr1_ele");
  TH1D * convert_muon_sr1 = (TH1D*)fOthers->Get("cr1_to_sr1_muon");

  TH1D * convert_ele_sr2 = (TH1D*)fOthers->Get("cr2_to_sr2_ele");
  TH1D * convert_muon_sr2 = (TH1D*)fOthers->Get("cr2_to_sr2_muon");

  for(int i = 0; i < nBackgrounds; i++) {

    // 2 systematics for each channel
    // differenceCR1_ele
    // convertCR1_ele
    TH1D * h_ele_sr1 = (TH1D*)fInputs->Get("ele_SR1/"+names[i]);
    TH1D * h_ele_sr1_up = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_differenceCR1_eleUp");
    TH1D * h_ele_sr1_down = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_differenceCR1_eleDown");
    TH1D * h_ele_sr1_up2 = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_convertCR1_eleUp");
    TH1D * h_ele_sr1_down2 = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_convertCR1_eleDown");
    
    TH1D * h_muon_sr1 = (TH1D*)fInputs->Get("muon_SR1/"+names[i]);
    TH1D * h_muon_sr1_up = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_differenceCR1_muonUp");
    TH1D * h_muon_sr1_down = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_differenceCR1_muonDown");
    TH1D * h_muon_sr1_up2 = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_convertCR1_muonUp");
    TH1D * h_muon_sr1_down2 = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_convertCR1_muonDown");
    
    TH1D * h_ele_sr2 = (TH1D*)fInputs->Get("ele_SR2/"+names[i]);
    TH1D * h_ele_sr2_up = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_differenceCR2_eleUp");
    TH1D * h_ele_sr2_down = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_differenceCR2_eleDown");
    TH1D * h_ele_sr2_up2 = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_convertCR2_eleUp");
    TH1D * h_ele_sr2_down2 = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_convertCR2_eleDown");

    TH1D * h_muon_sr2 = (TH1D*)fInputs->Get("muon_SR2/"+names[i]);
    TH1D * h_muon_sr2_up = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_differenceCR2_muonUp");
    TH1D * h_muon_sr2_down = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_differenceCR2_muonDown");
    TH1D * h_muon_sr2_up2 = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_convertCR2_muonUp");
    TH1D * h_muon_sr2_down2 = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_convertCR2_muonDown");
    
    for(int j = 0; j < diff_ele_cr1->GetNbinsX(); j++) {
      // (data - bkg) in cr1
      h_ele_sr1_up->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. + diff_ele_cr1->GetBinContent(j+1)));
      h_ele_sr1_down->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. - diff_ele_cr1->GetBinContent(j+1)));

      h_muon_sr1_up->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. + diff_muon_cr1->GetBinContent(j+1)));
      h_muon_sr1_down->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. - diff_muon_cr1->GetBinContent(j+1)));

      // (bkg sr1 / bkg sr2)
      h_ele_sr1_up2->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * convert_ele_sr1->GetBinContent(j+1));
      h_ele_sr1_down2->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (2. - convert_ele_sr1->GetBinContent(j+1)));

      h_muon_sr1_up2->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * convert_muon_sr1->GetBinContent(j+1));
      h_muon_sr1_down2->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (2. - convert_muon_sr1->GetBinContent(j+1)));
    }

    for(int j = 0; j < diff_ele_cr2->GetNbinsX(); j++) {
      // (data - bkg) in cr2
      h_ele_sr2_up->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. + diff_ele_cr2->GetBinContent(j+1)));
      h_ele_sr2_down->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. - diff_ele_cr2->GetBinContent(j+1)));

      h_muon_sr2_up->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. + diff_muon_cr2->GetBinContent(j+1)));
      h_muon_sr2_down->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. - diff_muon_cr2->GetBinContent(j+1)));

      // (bkg sr2 / bkg sr2)
      h_ele_sr2_up2->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * convert_ele_sr2->GetBinContent(j+1));
      h_ele_sr2_down2->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (2. - convert_ele_sr2->GetBinContent(j+1)));

      h_muon_sr2_up2->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * convert_muon_sr2->GetBinContent(j+1));
      h_muon_sr2_down2->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (2. - convert_muon_sr2->GetBinContent(j+1)));
    }

    fInputs->cd("ele_SR1");
    h_ele_sr1_up->Write(names[i]+"_differenceCR1_eleUp");
    h_ele_sr1_down->Write(names[i]+"_differenceCR1_eleDown");
    h_ele_sr1_up2->Write(names[i]+"_convertCR1_eleUp");
    h_ele_sr1_down2->Write(names[i]+"_convertCR1_eleDown");

    fInputs->cd("muon_SR1");
    h_muon_sr1_up->Write(names[i]+"_differenceCR1_muonUp");
    h_muon_sr1_down->Write(names[i]+"_differenceCR1_muonDown");
    h_muon_sr1_up2->Write(names[i]+"_convertCR1_muonUp");
    h_muon_sr1_down2->Write(names[i]+"_convertCR1_muonDown");

    fInputs->cd("ele_SR2");
    h_ele_sr2_up->Write(names[i]+"_differenceCR2_eleUp");
    h_ele_sr2_down->Write(names[i]+"_differenceCR2_eleDown");
    h_ele_sr2_up2->Write(names[i]+"_convertCR2_eleUp");
    h_ele_sr2_down2->Write(names[i]+"_convertCR2_eleDown");

    fInputs->cd("muon_SR2");
    h_muon_sr2_up->Write(names[i]+"_differenceCR2_muonUp");
    h_muon_sr2_down->Write(names[i]+"_differenceCR2_muonDown");
    h_muon_sr2_up2->Write(names[i]+"_convertCR2_muonUp");
    h_muon_sr2_down2->Write(names[i]+"_convertCR2_muonDown");
    
  } // for each background

  fInputs->Close();
  fDifferences->Close();
  fOthers->Close();
}
