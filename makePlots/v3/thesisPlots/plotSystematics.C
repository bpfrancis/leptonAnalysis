void plotSystematics() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  go("ele_SR1");
  go("ele_SR2");
  go("muon_SR1");
  go("muon_SR2");

}

void go(TString channel) {

  TString names[10] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma", "vgamma", "qcd"};

  TString systematics[19] = {"btagWeight", "puWeight", "JEC", "eleSF", "muonSF", "photonSF", "topPt",
			     "userSystA_ele", "userSystA_muon",
			     "userSystB_ele", "userSystB_muon",
			     "scale_tt", "scale_V", "scale_VV",
			     "pdf_gg", "pdf_qq", "pdf_gq",
			     "ele_qcdDef", "muon_qcdDef"};

  TFile * input = new TFile("../limitInputs_bjj.root", "READ");

  TH1D * h_central;
  TH1D * h_up;
  TH1D * h_down;

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);

  cout << endl << endl << "Channel: " << channel << endl;
    
  for(int i = 0; i < 19; i++) {

    h_central = (TH1D*)input->Get(channel+"/ttjets");
    h_up = (TH1D*)h_central->Clone(systematics[i]+"Up");
    h_down = (TH1D*)h_central->Clone(systematics[i]+"Down");

    h_central->Reset();
    h_up->Reset();
    h_down->Reset();

    for(int j = 0; j < 10; j++) {

      if(names[j] == "qcd") continue;

      TH1D * h = (TH1D*)input->Get(channel+"/"+names[j]);
    
      h_central->Add(h);
      
      TH1D * hup = (TH1D*)input->Get(channel+"/"+names[j]+"_"+systematics[i]+"Up");
      TH1D * hdown = (TH1D*)input->Get(channel+"/"+names[j]+"_"+systematics[i]+"Down");
      
      if(!hup || !hdown) {
	h_up->Add(h);
	h_down->Add(h);
      }
      else {
	h_up->Add(hup);
	h_down->Add(hdown);
      }

    }
    
    double n_central = h_central->Integral();
    double n_up = h_up->Integral();
    double n_down = h_down->Integral();

    double diff_up = fabs(n_up - n_central) / n_central;
    double diff_down = fabs(n_down - n_central) / n_central;

    if(diff_up > diff_down) cout << systematics[i] << " -- " << diff_up * 100. << endl;
    else cout << systematics[i] << " -- " << diff_down * 100. << endl;

    h_central->Scale(1./h_central->Integral());
    h_up->Scale(1./h_up->Integral());
    h_down->Scale(1./h_down->Integral());
    
    h_up->SetLineColor(kRed);
    h_down->SetLineColor(kBlue);
    
    TLegend * leg = new TLegend(0.55, 0.65, 0.85, 0.85, systematics[i].Data(), "brNDC");
    leg->AddEntry(h_up, "Up", "LP");
    leg->AddEntry(h_central, "Central", "LP");
    leg->AddEntry(h_down, "Down", "LP");
    
    h_up->Draw("hist");
    h_central->Draw("hist same");
    h_down->Draw("hist same");
    leg->Draw("same");
    
    can->SaveAs("systematicPlots/"+channel+"_"+systematics[i]+".pdf");
    
  }

  cout << "----------------------" << endl << endl;


}

