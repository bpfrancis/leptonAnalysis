void go(TH2D * h, TString xaxisLabel, TString yaxisLabel, TString headerLabel, bool useLogx, TString name) {

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);
  can->SetLogx(useLogx);

  TPaveText * header = new TPaveText(0.1, 0.901, 0.9, 0.94, "NDC");
  header->SetFillColor(0);
  header->SetFillStyle(0);
  header->SetLineColor(0);
  header->AddText(headerLabel);

  h->GetXaxis()->SetTitle(xaxisLabel);
  h->GetXaxis()->SetLabelFont(63);
  h->GetXaxis()->SetLabelSize(48);
  
  h->GetYaxis()->SetTitle(yaxisLabel);
  h->GetYaxis()->SetLabelFont(63);
  h->GetYaxis()->SetLabelSize(48);
  
  if(xaxisLabel.Contains("P_{T}")) {
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.6);
  }
  else {
    h->GetYaxis()->SetTitleSize(0.035);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.6);
  }

  h->GetZaxis()->SetLabelSize(0.02);
  h->Draw("colz");
  header->Draw("same");
  can->SaveAs(name);

  delete can;
}
  

void plotLeptonPhotonSFs() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TFile * fLepton = new TFile("/eos/uscms/store/user/bfrancis/data/lepton_SF_8TeV_53x_baseline.root", "READ");
  TFile * fPhoton = new TFile("/eos/uscms/store/user/bfrancis/data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root", "READ");

  TH2D * sf_electron = (TH2D*)fLepton->Get("TightEleIdIsoSF");
  TH2D * sf_SingleElectronTrigger = (TH2D*)fLepton->Get("TightEleTriggerSF");
  TH2D * sf_muon = (TH2D*)fLepton->Get("mu_pt_eta_full_id_iso_hlt_8TeV");

  TH2D * sf_photon_id = (TH2D*)fPhoton->Get("PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
  TH2D * sf_photon_veto = (TH2D*)fPhoton->Get("PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");

  //go(sf_electron, "|#eta|", "P_{T} (GeV/c)", "Electron ID #times Iso", false, "sfElectron.pdf");
  //go(sf_SingleElectronTrigger, "|#eta|", "P_{T} (GeV/c)", "SingleElectronTrigger", false, "sfSingleElectronTrigger.pdf");
  go(sf_muon, "P_{T} (GeV/c)", "|#eta|", "Muon ID #times Iso #times Trigger", false, "sfMuon.pdf");

  sf_electron->Multiply(sf_SingleElectronTrigger);
  go(sf_electron, "|#eta|", "P_{T} (GeV/c)", "Electron ID #times Iso #times Trigger", false, "sfElectronWithTrigger.pdf");

  sf_photon_id->GetZaxis()->SetRangeUser(0.95, 1.05);
  sf_photon_veto->GetZaxis()->SetRangeUser(0.95, 1.05);

  go(sf_photon_id, "P_{T} (GeV/c)", "|#eta|", "Photon ID", false, "sfPhotonID.pdf");
  go(sf_photon_veto, "P_{T} (GeV/c)", "|#eta|", "Photon ConvSafeVeto", false, "sfPhotonVeto.pdf");

  fLepton->Close();
  fPhoton->Close();

}
