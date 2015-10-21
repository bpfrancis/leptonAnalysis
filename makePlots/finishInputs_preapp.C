#include "Binning.h"

const int nBackgrounds = 10;
TString names[nBackgrounds] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

void finish_SR1() {

  TFile * fInputs = new TFile("limitInputs_bjj.root", "UPDATE");

  // create difference fraction in CR1
  TH1D * h_systA_ele = (TH1D*)fInputs->Get("ele_CR1/data_obs");
  TH1D * h_systA_muon = (TH1D*)fInputs->Get("muon_CR1/data_obs");

  TH1D * h_bkg_ele_cr1 = (TH1D*)fInputs->Get("ele_CR1/ttjets");
  TH1D * h_bkg_muon_cr1 = (TH1D*)fInputs->Get("muon_CR1/ttjets");

  for(int i = 1; i < nBackgrounds; i++) {
    h_bkg_ele_cr1->Add((TH1D*)fInputs->Get("ele_CR1/"+names[i]));
    if(names[i] != "qcd") h_bkg_muon_cr1->Add((TH1D*)fInputs->Get("muon_CR1/"+names[i]));
  }

 h_systA_ele->Add(h_bkg_ele_cr1, -1.);
 h_systA_ele->Divide(h_bkg_ele_cr1);

 h_systA_muon->Add(h_bkg_muon_cr1, -1.);
 h_systA_muon->Divide(h_bkg_muon_cr1);

  // create the sig(chi)/chi for systB/C
  //for chi = sr1/cr1 (systB)
  //and chi = sr2/sr1 (systC)
  // already have h_bkg_*_cr1 available

  TH1D * h_systB_ele = (TH1D*)fInputs->Get("ele_SR1/ttjets");
  TH1D * h_systB_muon = (TH1D*)fInputs->Get("muon_SR1/ttjets");

  for(int i = 1; i < nBackgrounds; i++) {
    if(names[i] == "qcd") continue;

    h_systB_ele->Add((TH1D*)fInputs->Get("ele_SR1/"+names[i]));
    h_systB_muon->Add((TH1D*)fInputs->Get("muon_SR1/"+names[i]));
  }

  h_systB_ele->Divide(h_bkg_ele_cr1);
  h_systB_muon->Divide(h_bkg_muon_cr1);

  // now set the convFactors to sig(it)/it
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

  // so we have our systematics:
  // h_systA_* and h_systB_* for SR1 and SR2
  // h_systC_* for just SR2
  // now store these fluctuations for all backgrounds

  for(int i = 0; i < nBackgrounds; i++) {

    // SR1
    TH1D * h_ele_sr1 = (TH1D*)fInputs->Get("ele_SR1/"+names[i]);
    TH1D * h_ele_sr1_upA = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_systA_eleUp");
    TH1D * h_ele_sr1_downA = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_systA_eleDown");
    TH1D * h_ele_sr1_upB = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_systB_eleUp");
    TH1D * h_ele_sr1_downB = (TH1D*)h_ele_sr1->Clone("ele_SR1_"+names[i]+"_systB_eleDown");

    TH1D * h_muon_sr1 = (TH1D*)fInputs->Get("muon_SR1/"+names[i]);
    TH1D * h_muon_sr1_upA = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_systA_muonUp");
    TH1D * h_muon_sr1_downA = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_systA_muonDown");
    TH1D * h_muon_sr1_upB = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_systB_muonUp");
    TH1D * h_muon_sr1_downB = (TH1D*)h_muon_sr1->Clone("muon_SR1_"+names[i]+"_systB_muonDown");
    
    for(int j = 0; j < h_systA_ele->GetNbinsX(); j++) {

      // systA for both sr1 and sr2
      h_ele_sr1_upA->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. + h_systA_ele->GetBinContent(j+1)));
      h_ele_sr1_downA->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. - h_systA_ele->GetBinContent(j+1)));

      h_muon_sr1_upA->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. + h_systA_muon->GetBinContent(j+1)));
      h_muon_sr1_downA->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. - h_systA_muon->GetBinContent(j+1)));

      // systB for both sr1 and sr2
      h_ele_sr1_upB->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. + h_systB_ele->GetBinContent(j+1)));
      h_ele_sr1_downB->SetBinContent(j+1, h_ele_sr1->GetBinContent(j+1) * (1. - h_systB_ele->GetBinContent(j+1)));

      h_muon_sr1_upB->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. + h_systB_muon->GetBinContent(j+1)));
      h_muon_sr1_downB->SetBinContent(j+1, h_muon_sr1->GetBinContent(j+1) * (1. - h_systB_muon->GetBinContent(j+1)));

    }

    fInputs->cd("ele_SR1");
    h_ele_sr1_upA->Write(names[i]+"_userSystA_eleUp");
    h_ele_sr1_downA->Write(names[i]+"_userSystA_eleDown");
    h_ele_sr1_upB->Write(names[i]+"_userSystB_eleUp");
    h_ele_sr1_downB->Write(names[i]+"_userSystB_eleDown");

    fInputs->cd("muon_SR1");
    h_muon_sr1_upA->Write(names[i]+"_userSystA_muonUp");
    h_muon_sr1_downA->Write(names[i]+"_userSystA_muonDown");
    h_muon_sr1_upB->Write(names[i]+"_userSystB_muonUp");
    h_muon_sr1_downB->Write(names[i]+"_userSystB_muonDown");
    
  } // for each background

  fInputs->Close();

}
  
void finish_SR2() {

  TFile * fInputs = new TFile("limitInputs_bjj.root", "UPDATE");

  // create difference fraction in CR1
  TH1D * h_systA_ele = (TH1D*)fInputs->Get("ele_CR1/data_obs");
  TH1D * h_systA_muon = (TH1D*)fInputs->Get("muon_CR1/data_obs");

  TH1D * h_bkg_ele_cr1 = (TH1D*)fInputs->Get("ele_CR1/ttjets");
  TH1D * h_bkg_muon_cr1 = (TH1D*)fInputs->Get("muon_CR1/ttjets");

  for(int i = 1; i < nBackgrounds; i++) {
    h_bkg_ele_cr1->Add((TH1D*)fInputs->Get("ele_CR1/"+names[i]));
    if(names[i] != "qcd") h_bkg_muon_cr1->Add((TH1D*)fInputs->Get("muon_CR1/"+names[i]));
  }

  h_systA_ele = (TH1D*)h_systA_ele->Rebin(nMetBins_2g, "systA_ele_rebin", xbins_met_2g);
  h_systA_muon = (TH1D*)h_systA_muon->Rebin(nMetBins_2g, "systA_muon_rebin", xbins_met_2g);
  
  h_bkg_ele_cr1 = (TH1D*)h_bkg_ele_cr1->Rebin(nMetBins_2g, "ele_cr1_rebin", xbins_met_2g);
  h_bkg_muon_cr1 = (TH1D*)h_bkg_muon_cr1->Rebin(nMetBins_2g, "muon_cr1_rebin", xbins_met_2g);
  
 h_systA_ele->Add(h_bkg_ele_cr1, -1.);
 h_systA_ele->Divide(h_bkg_ele_cr1);

 h_systA_muon->Add(h_bkg_muon_cr1, -1.);
 h_systA_muon->Divide(h_bkg_muon_cr1);

  // create the sig(chi)/chi for systB/C
  //for chi = sr1/cr1 (systB)
  //and chi = sr2/sr1 (systC)
  // already have h_bkg_*_cr1 available

  TH1D * h_systB_ele = (TH1D*)fInputs->Get("ele_SR1/ttjets");
  TH1D * h_systB_muon = (TH1D*)fInputs->Get("muon_SR1/ttjets");

  TH1D * h_systC_ele = (TH1D*)fInputs->Get("ele_SR2/ttjets");
  TH1D * h_systC_muon = (TH1D*)fInputs->Get("muon_SR2/ttjets");

  for(int i = 1; i < nBackgrounds; i++) {
    if(names[i] == "qcd") continue;

    h_systB_ele->Add((TH1D*)fInputs->Get("ele_SR1/"+names[i]));
    h_systB_muon->Add((TH1D*)fInputs->Get("muon_SR1/"+names[i]));

    h_systC_ele->Add((TH1D*)fInputs->Get("ele_SR2/"+names[i]));
    h_systC_muon->Add((TH1D*)fInputs->Get("muon_SR2/"+names[i]));
  }

  h_systB_ele = (TH1D*)h_systB_ele->Rebin(nMetBins_2g, "systB_ele_rebin", xbins_met_2g);
  h_systB_muon = (TH1D*)h_systB_muon->Rebin(nMetBins_2g, "systB_muon_rebin", xbins_met_2g);

  h_systC_ele = (TH1D*)h_systC_ele->Rebin(nMetBins_2g, "systC_ele_rebin", xbins_met_2g);
  h_systC_muon = (TH1D*)h_systC_muon->Rebin(nMetBins_2g, "systC_muon_rebin", xbins_met_2g);

  h_systC_ele->Divide(h_systB_ele);
  h_systC_muon->Divide(h_systB_muon);

  h_systB_ele->Divide(h_bkg_ele_cr1);
  h_systB_muon->Divide(h_bkg_muon_cr1);

  // now set the convFactors to sig(it)/it
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

  for(int i = 0; i < h_systC_ele->GetNbinsX(); i++) {
    double val = h_systC_ele->GetBinContent(i+1);
    double err = h_systC_ele->GetBinError(i+1);
    h_systC_ele->SetBinContent(i+1, err/val);
  }

  for(int i = 0; i < h_systC_muon->GetNbinsX(); i++) {
    double val = h_systC_muon->GetBinContent(i+1);
    double err = h_systC_muon->GetBinError(i+1);
    h_systC_muon->SetBinContent(i+1, err/val);
  }

  // so we have our systematics:
  // h_systA_* and h_systB_* for SR1 and SR2
  // h_systC_* for just SR2
  // now store these fluctuations for all backgrounds

  for(int i = 0; i < nBackgrounds; i++) {

    // SR2
    TH1D * h_ele_sr2 = (TH1D*)fInputs->Get("ele_SR2/"+names[i]);
    TH1D * h_ele_sr2_upA = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systA_eleUp");
    TH1D * h_ele_sr2_downA = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systA_eleDown");
    TH1D * h_ele_sr2_upB = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systB_eleUp");
    TH1D * h_ele_sr2_downB = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systB_eleDown");
    TH1D * h_ele_sr2_upC = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systC_eleUp");
    TH1D * h_ele_sr2_downC = (TH1D*)h_ele_sr2->Clone("ele_SR2_"+names[i]+"_systC_eleDown");

    TH1D * h_muon_sr2 = (TH1D*)fInputs->Get("muon_SR2/"+names[i]);
    TH1D * h_muon_sr2_upA = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systA_muonUp");
    TH1D * h_muon_sr2_downA = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systA_muonDown");
    TH1D * h_muon_sr2_upB = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systB_muonUp");
    TH1D * h_muon_sr2_downB = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systB_muonDown");
    TH1D * h_muon_sr2_upC = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systC_muonUp");
    TH1D * h_muon_sr2_downC = (TH1D*)h_muon_sr2->Clone("muon_SR2_"+names[i]+"_systC_muonDown");
    
    for(int j = 0; j < h_systA_ele->GetNbinsX(); j++) {

      // systA for both sr1 and sr2
      h_ele_sr2_upA->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. + h_systA_ele->GetBinContent(j+1)));
      h_ele_sr2_downA->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. - h_systA_ele->GetBinContent(j+1)));

      h_muon_sr2_upA->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. + h_systA_muon->GetBinContent(j+1)));
      h_muon_sr2_downA->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. - h_systA_muon->GetBinContent(j+1)));

      // systB for both sr1 and sr2
      h_ele_sr2_upB->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. + h_systB_ele->GetBinContent(j+1)));
      h_ele_sr2_downB->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. - h_systB_ele->GetBinContent(j+1)));

      h_muon_sr2_upB->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. + h_systB_muon->GetBinContent(j+1)));
      h_muon_sr2_downB->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. - h_systB_muon->GetBinContent(j+1)));

      // systC for just sr2
      h_ele_sr2_upC->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. + h_systC_ele->GetBinContent(j+1)));
      h_ele_sr2_downC->SetBinContent(j+1, h_ele_sr2->GetBinContent(j+1) * (1. - h_systC_ele->GetBinContent(j+1)));

      h_muon_sr2_upC->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. + h_systC_muon->GetBinContent(j+1)));
      h_muon_sr2_downC->SetBinContent(j+1, h_muon_sr2->GetBinContent(j+1) * (1. - h_systC_muon->GetBinContent(j+1)));

    }

    fInputs->cd("ele_SR2");
    h_ele_sr2_upA->Write(names[i]+"_userSystA_eleUp");
    h_ele_sr2_downA->Write(names[i]+"_userSystA_eleDown");
    h_ele_sr2_upB->Write(names[i]+"_userSystB_eleUp");
    h_ele_sr2_downB->Write(names[i]+"_userSystB_eleDown");
    h_ele_sr2_upC->Write(names[i]+"_userSystC_eleUp");
    h_ele_sr2_downC->Write(names[i]+"_userSystC_eleDown");

    fInputs->cd("muon_SR2");
    h_muon_sr2_upA->Write(names[i]+"_userSystA_muonUp");
    h_muon_sr2_downA->Write(names[i]+"_userSystA_muonDown");
    h_muon_sr2_upB->Write(names[i]+"_userSystB_muonUp");
    h_muon_sr2_downB->Write(names[i]+"_userSystB_muonDown");
    h_muon_sr2_upC->Write(names[i]+"_userSystC_muonUp");
    h_muon_sr2_downC->Write(names[i]+"_userSystC_muonDown");
    
  } // for each background

  fInputs->Close();
}

void finishInputs_preapp() {

  finish_SR1();
  finish_SR2();

}
