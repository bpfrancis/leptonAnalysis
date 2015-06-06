const int nBackgrounds = 10;
TString names[nBackgrounds] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

void finishInputs() {

  TFile * fInputs = new TFile("limitInputs_bjj.root", "UPDATE");

  // create difference fraction in CR1
  TH1D * h_diffFraction_ele_cr1 = (TH1D*)fInputs->Get("ele_CR1/data_obs");
  TH1D * h_diffFraction_muon_cr1 = (TH1D*)fInputs->Get("muon_CR1/data_obs");

  TH1D * h_bkg_ele_cr1 = (TH1D*)fInputs->Get("ele_CR1/ttjets");
  TH1D * h_bkg_muon_cr1 = (TH1D*)fInputs->Get("muon_CR1/ttjets");

  for(int i = 1; i < nBackgrounds; i++) {
    h_bkg_ele_cr1->Add((TH1D*)fInputs->Get("ele_CR1/"+names[i]));
    if(names[i] != "qcd") h_bkg_muon_cr1->Add((TH1D*)fInputs->Get("muon_CR1/"+names[i]));
  }

  h_diffFraction_ele_cr1->Add(h_bkg_ele_cr1, -1.);
  h_diffFraction_ele_cr1->Divide(h_bkg_ele_cr1);

  h_diffFraction_muon_cr1->Add(h_bkg_muon_cr1, -1.);
  h_diffFraction_muon_cr1->Divide(h_bkg_muon_cr1);

  // create the sig(chi)/chi 
  //for chi = sr1/cr1
  //and chi = sr2/sr1
  // already have h_bkg_*_cr1 available
  TH1D * h_convFactor_ele_sr1 = (TH1D*)fInputs->Get("ele_SR1/ttjets");
  TH1D * h_convFactor_muon_sr1 = (TH1D*)fInputs->Get("muon_SR1/ttjets");

  TH1D * h_systB_ele = (TH1D*)fInputs->Get("ele_SR2/ttjets");
  TH1D * h_systB_muon = (TH1D*)fInputs->Get("muon_SR2/ttjets");

  for(int i = 1; i < nBackgrounds; i++) {
    if(names[i] == "qcd") continue;

    h_convFactor_ele_sr1->Add((TH1D*)fInputs->Get("ele_SR1/"+names[i]));
    h_convFactor_muon_sr1->Add((TH1D*)fInputs->Get("muon_SR1/"+names[i]));

    h_systB_ele->Add((TH1D*)fInputs->Get("ele_SR2/"+names[i]));
    h_systB_muon->Add((TH1D*)fInputs->Get("muon_SR2/"+names[i]));
  }

  h_systB_ele->Divide(h_convFactor_ele_sr1);
  h_systB_muon->Divide(h_convFactor_muon_sr1);

  h_convFactor_ele_sr1->Divide(h_bkg_ele_cr1);
  h_convFactor_muon_sr1->Divide(h_bkg_muon_cr1);

  // now set the convFactors to sig(it)/it
  for(int i = 0; i < h_convFactor_ele_sr1->GetNbinsX(); i++) {
    double val = h_convFactor_ele_sr1->GetBinContent(i+1);
    double err = h_convFactor_ele_sr1->GetBinError(i+1);
    h_convFactor_ele_sr1->SetBinContent(i+1, err/val);
  }

  for(int i = 0; i < h_convFactor_muon_sr1->GetNbinsX(); i++) {
    double val = h_convFactor_muon_sr1->GetBinContent(i+1);
    double err = h_convFactor_muon_sr1->GetBinError(i+1);
    h_convFactor_muon_sr1->SetBinContent(i+1, err/val);
  }

  for(int i = 0; i < h_systB_ele->GetNbinsX(); i++) {
    double val = h_systB_ele->GetBinContent(i+1);
    double err = h_systB_ele->GetBinError(i+1);
    h_systB_ele->SetBinContent(i+1, err/val);
  }

  for(int i = 0; i < h_systB_muon->GetNbinsX(); i++) {
    double val = h_systB_muon->GetBinContent(i+1);
    double err = h_systB_muon->GetBinError(i+1);
    h_systB_muon->SetBinContent(i+1, err/val);
  }

  // now we have three things:
  // h_diffFraction_ele_cr1,
  // h_convFactor_ele_sr1, and
  // h_systB_ele
  // the third is done and only for SR2, the first two need to be added in quadrature

  TH1D * h_systA_ele = (TH1D*)h_diffFraction_ele_cr1->Clone("systA_ele");
  for(int i = 0; i < h_systA_ele->GetNbinsX(); i++) {
    double x = h_diffFraction_ele_cr1->GetBinContent(i+1);
    double y = h_convFactor_ele_sr1->GetBinContent(i+1);

    h_systA_ele->SetBinContent(i+1, sqrt(x*x + y*y));
  }

  TH1D * h_systA_muon = (TH1D*)h_diffFraction_muon_cr1->Clone("systA_muon");
  for(int i = 0; i < h_systA_muon->GetNbinsX(); i++) {
    double x = h_diffFraction_muon_cr1->GetBinContent(i+1);
    double y = h_convFactor_muon_sr1->GetBinContent(i+1);

    h_systA_muon->SetBinContent(i+1, sqrt(x*x + y*y));
  }

  // so we have our systematics:
  // h_systA_* for SR1 and SR2
  // h_systB_* for just SR2
  // now store these fluctuations for all backgrounds

  for(int i = 0; i < nBackgrounds; i++) {

    TH1D * h_ele_sr1 = (TH1D*)fInputs->Get("ele_SR1/"+names[i]);
    TH1D * h_ele_sr1_up = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_systA_eleUp");
    TH1D * h_ele_sr1_down = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_systA_eleDown");

    TH1D * h_muon_sr1 = (TH1D*)fInputs->Get("muon_SR1/"+names[i]);
    TH1D * h_muon_sr1_up = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_systA_muonUp");
    TH1D * h_muon_sr1_down = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_systA_muonDown");
    
    TH1D * h_ele_sr2 = (TH1D*)fInputs->Get("ele_SR2/"+names[i]);
    TH1D * h_ele_sr2_up = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systA_eleUp");
    TH1D * h_ele_sr2_down = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systA_eleDown");
    TH1D * h_ele_sr2_Bup = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systB_eleUp");
    TH1D * h_ele_sr2_Bdown = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systB_eleDown");

    TH1D * h_muon_sr2 = (TH1D*)fInputs->Get("muon_SR2/"+names[i]);
    TH1D * h_muon_sr2_up = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systA_muonUp");
    TH1D * h_muon_sr2_down = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systA_muonDown");
    TH1D * h_muon_sr2_Bup = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systB_muonUp");
    TH1D * h_muon_sr2_Bdown = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systB_muonDown");

    for(int j = 0; j < h_systA_ele->GetNbinsX(); j++) {

      // systA for both sr1 and sr2
      h_ele_sr1_up->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. + h_systA_ele->GetBinContent(j+1)));
      h_ele_sr1_down->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. - h_systA_ele->GetBinContent(j+1)));

      h_muon_sr1_up->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. + h_systA_muon->GetBinContent(j+1)));
      h_muon_sr1_down->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. - h_systA_muon->GetBinContent(j+1)));

      h_ele_sr2_up->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. + h_systA_ele->GetBinContent(j+1)));
      h_ele_sr2_down->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. - h_systA_ele->GetBinContent(j+1)));

      h_muon_sr2_up->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. + h_systA_muon->GetBinContent(j+1)));
      h_muon_sr2_down->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. - h_systA_muon->GetBinContent(j+1)));

      // systB for just sr2
      h_ele_sr2_Bup->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. + h_systB_ele->GetBinContent(j+1)));
      h_ele_sr2_Bdown->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. - h_systB_ele->GetBinContent(j+1)));

      h_muon_sr2_Bup->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. + h_systB_muon->GetBinContent(j+1)));
      h_muon_sr2_Bdown->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. - h_systB_muon->GetBinContent(j+1)));

    }

    fInputs->cd("ele_SR1");
    h_ele_sr1_up->Write(names[i]+"_userSystA_eleUp");
    h_ele_sr1_down->Write(names[i]+"_userSystA_eleDown");

    fInputs->cd("muon_SR1");
    h_muon_sr1_up->Write(names[i]+"_userSystA_muonUp");
    h_muon_sr1_down->Write(names[i]+"_userSystA_muonDown");

    fInputs->cd("ele_SR2");
    h_ele_sr2_up->Write(names[i]+"_userSystA_eleUp");
    h_ele_sr2_down->Write(names[i]+"_userSystA_eleDown");
    h_ele_sr2_Bup->Write(names[i]+"_userSystB_eleUp");
    h_ele_sr2_Bdown->Write(names[i]+"_userSystB_eleDown");

    fInputs->cd("muon_SR2");
    h_muon_sr2_up->Write(names[i]+"_userSystA_muonUp");
    h_muon_sr2_down->Write(names[i]+"_userSystA_muonDown");
    h_muon_sr2_Bup->Write(names[i]+"_userSystB_muonUp");
    h_muon_sr2_Bdown->Write(names[i]+"_userSystB_muonDown");
    
  } // for each background

  fInputs->Close();
}
