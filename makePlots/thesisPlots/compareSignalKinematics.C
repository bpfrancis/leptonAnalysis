void go(TString varName, TString channel, TString controlRegion, int nbins, double xlo, double xhi, bool sameBinoMass) {

  TFile * input_A = new TFile("durp/signal_contamination_mst_460_m1_175.root", "READ");
  TFile * input_B = new TFile("durp/signal_contamination_mst_910_m1_575.root", "READ");
  TFile * input_C = new TFile("durp/signal_contamination_mst_810_m1_200.root", "READ");
  TFile * input_D = new TFile("durp/signal_contamination_mst_1410_m1_200.root", "READ");
  TFile * input_E = new TFile("durp/signal_contamination_mst_910_m1_175.root", "READ");

  TTree * tree_A = (TTree*)input_A->Get(channel+"_signalTree");
  TTree * tree_B = (TTree*)input_B->Get(channel+"_signalTree");
  TTree * tree_C = (TTree*)input_C->Get(channel+"_signalTree");
  TTree * tree_D = (TTree*)input_D->Get(channel+"_signalTree");
  TTree * tree_E = (TTree*)input_E->Get(channel+"_signalTree");

  Float_t met, nphotons, var;

  //tree_A->SetBranchAddress("pfMET_t01", &met);
  tree_A->SetBranchAddress("Ngamma", &nphotons);
  tree_A->SetBranchAddress(varName, &var);

  //tree_B->SetBranchAddress("pfMET_t01", &met);
  tree_B->SetBranchAddress("Ngamma", &nphotons);
  tree_B->SetBranchAddress(varName, &var);

  //tree_C->SetBranchAddress("pfMET_t01", &met);
  tree_C->SetBranchAddress("Ngamma", &nphotons);
  tree_C->SetBranchAddress(varName, &var);
  
  //tree_D->SetBranchAddress("pfMET_t01", &met);
  tree_D->SetBranchAddress("Ngamma", &nphotons);
  tree_D->SetBranchAddress(varName, &var);

  //tree_E->SetBranchAddress("pfMET_t01", &met);
  tree_E->SetBranchAddress("Ngamma", &nphotons);
  tree_E->SetBranchAddress(varName, &var);

  TH1D * h_A = new TH1D("h_a", "h_a", nbins, xlo, xhi);
  TH1D * h_B = new TH1D("h_b", "h_b", nbins, xlo, xhi);
  TH1D * h_C = new TH1D("h_c", "h_c", nbins, xlo, xhi);
  TH1D * h_D = new TH1D("h_d", "h_d", nbins, xlo, xhi);
  TH1D * h_E = new TH1D("h_e", "h_e", nbins, xlo, xhi);

  for(int i = 0; i < tree_A->GetEntries(); i++) {
    tree_A->GetEntry(i);
    if(nphotons != 1 && controlRegion == "SR1") continue;
    if(nphotons < 2 && controlRegion == "SR2") continue;

    h_A->Fill(var);
  }

  for(int i = 0; i < tree_B->GetEntries(); i++) {
    tree_B->GetEntry(i);
    if(nphotons != 1 && controlRegion == "SR1") continue;
    if(nphotons < 2 && controlRegion == "SR2") continue;

    h_B->Fill(var);
  }

  for(int i = 0; i < tree_C->GetEntries(); i++) {
    tree_C->GetEntry(i);
    if(nphotons != 1 && controlRegion == "SR1") continue;
    if(nphotons < 2 && controlRegion == "SR2") continue;

    h_C->Fill(var);
  }

  for(int i = 0; i < tree_D->GetEntries(); i++) {
    tree_D->GetEntry(i);
    if(nphotons != 1 && controlRegion == "SR1") continue;
    if(nphotons < 2 && controlRegion == "SR2") continue;

    h_D->Fill(var);
  }
  
  for(int i = 0; i < tree_E->GetEntries(); i++) {
    tree_E->GetEntry(i);
    if(nphotons != 1 && controlRegion == "SR1") continue;
    if(nphotons < 2 && controlRegion == "SR2") continue;

    h_E->Fill(var);
  }

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);

  h_A->Scale(1./h_A->Integral());
  h_B->Scale(1./h_B->Integral());
  h_C->Scale(1./h_C->Integral());
  h_D->Scale(1./h_D->Integral());
  h_E->Scale(1./h_E->Integral());

  h_A->GetXaxis()->SetTitle(varName);

  h_B->SetLineColor(kRed);
  h_C->SetLineColor(kBlue);
  h_D->SetLineColor(8);
  h_E->SetLineColor(kViolet);

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, "#Delta", "brNDC");
  leg->AddEntry(h_a, "A 285", "L");
  
  if(sameBinoMass) {
    leg->AddEntry(h_e, "E 735", "L");
    h_A->Draw("hist");
    h_E->Draw("hist same");
    leg->Draw("same");
  }

  else {
    leg->AddEntry(h_b, "B 335", "L");
    leg->AddEntry(h_c, "C 610", "L");
    leg->AddEntry(h_d, "D 1210", "L");
    leg->AddEntry(h_e, "E 735", "L");
    
    h_A->Draw("hist");
    h_B->Draw("hist same");
    h_C->Draw("hist same");
    h_D->Draw("hist same");
    h_E->Draw("hist same");
    leg->Draw("same");
  }

  can->SaveAs("durp/plots/"+varName+"_"+channel+"_"+controlRegion+".pdf");

  delete can;

  input_A->Close();
  input_B->Close();
  input_C->Close();
  input_D->Close();
  input_E->Close();
  
}

void compareSignalKinematics() {

  go("leadPhotonEt", "muon_bjj", "Any", 50, 0., 1000., false);
  go("trailPhotonEt", "muon_bjj", "Any", 50, 0., 1000., false);
  go("btag1_pt", "muon_bjj", "Any", 50, 0., 1000., false);
  go("pfMET_t01", "muon_bjj", "Any", 50, 0., 1000., false);
  go("Njets", "muon_bjj", "Any", 10, 3., 13., false);

}
