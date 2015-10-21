#include <vector>

using namespace std;

void qcdContamination() {

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TString channel = "ele_bjj_eQCDTree";

  TString prefix = "/eos/uscms/store/user/MYUSERNAME/stopBino_inputs/signal_contamination_";

  const int nBkgs = 23;

  TString backgrounds[nBkgs] = {"ttJetsSemiLep", "ttJetsFullLep", "ttJetsHadronic",
				"TBar_s", "TBar_tW", "TBar_t", "T_s", "T_t", "T_tW",
				"W3JetsToLNu", "W4JetsToLNu",
				"dy1JetsToLL", "dy2JetsToLL", "dy3JetsToLL", "dy4JetsToLL",
				"ZGToLLG", "WGToLNuG",
				"WW", "WZ", "ZZ",
				"TTZJets", "TTWJets",
				"TTGamma"};

  double xsec[nBkgs] = {245.8 * 0.438, 245.8 * 0.105, 245.8 * 0.457,
			1.76, 11.1, 30.7, 3.79, 56.4, 11.1,
			650. * 3. * 12234.4 / 37509., 264. * 3. * 12234.4 / 37509.,
			666.7 * 3. * 1177.3 / 3503.71, 215.1 * 3. * 1177.3 / 3503.71, 66.07 * 3. * 1177.3 / 3503.71, 27.38 * 3. * 1177.3 / 3503.71,
			159.1, 553.9,
			57.1097, 32.3161, 8.25561,
			0.2057, 0.232,
			0.033 * 9 + 0.148 * 12 + 0.8};

  TFile * fData = channel.Contains("ele") ? new TFile("/eos/uscms/store/user/MYUSERNAME/stopBino_inputs/SingleElectron.root") : 
    new TFile("/eos/uscms/store/user/MYUSERNAME/stopBino_inputs/SingleMu.root");
  TTree * tData = (TTree*)fData->Get(channel);

  const int nMetBins = 10;
  Double_t xbins_met[nMetBins+1] = {0, 10, 20, 30, 40, 50, 75, 100, 150, 300, 650};

  TH1D * hData = new TH1D("hData", "hData", nMetBins, xbins_met); hData->Sumw2();
  TH1D * hBkg = new TH1D("hBkg", "hBkg", nMetBins, xbins_met); hBkg->Sumw2();

  Float_t met;
  tData->SetBranchAddress("pfMET_t01", &met);

  for(int i = 0; i < tData->GetEntries(); i++) {
    tData->GetEntry(i);
    hData->Fill(met);
  }

  TFile * fBkg;

  for(int i = 0; i < nBkgs; i++) {

    if(!(backgrounds[i].Contains("ttJets")) && !(backgrounds[i].Contains("W3")) && !(backgrounds[i].Contains("W4"))) continue;

    cout << "Start working on " << backgrounds[i] << endl;

    fBkg = new TFile(prefix+backgrounds[i]+".root", "READ");
    TTree * tBkg = (TTree*)fBkg->Get(channel);

    TH1D * h_ngen = (TH1D*)fBkg->Get("nEvents_"+backgrounds[i]);

    double ngen = h_ngen->Integral();

    double weight = 19712 * xsec[i] / ngen;

    tBkg->SetBranchAddress("pfMET_t01", &met);

    for(int j = 0; j < tBkg->GetEntries(); j++) {
      tBkg->GetEntry(j);
      hBkg->Fill(met, weight);
    }

    fBkg->Close();

  }

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 800, 800);
  can->SetLogy(true);

  hData->GetXaxis()->SetTitle("MET (GeV)");
  if(channel.Contains("ele")) hData->GetYaxis()->SetRangeUser(0.1, 2.e4);
  //else {}

  hData->SetLineWidth(2);

  hData->Draw("e1");
  
  hBkg->SetFillColor(kGray);
  hBkg->Draw("hist same");
  hData->Draw("e1 same");
  hData->Draw("axis same");

  TLegend * leg = new TLegend(0.45, 0.6, 0.85, 0.85, channel.Contains("ele") ? "Electron Channel" : "Muon Channel", "brNDC");
  leg->AddEntry(hData, "QCD Data", "LP");
  leg->AddEntry(hBkg, "MC Contamination", "F");
  leg->Draw("same");

  can->SaveAs("qcdContamination_"+channel+".png");
  
}
