void qcdSubtraction() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TFile * fEle = new TFile("qcdHistograms_ele_bjj_Any.root", "READ");
  TFile * fMuon = new TFile("qcdHistograms_muon_bjj_Any.root", "READ");

  TFile * fCentral = new TFile("../");

  //durp;
}
