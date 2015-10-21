#include <iostream>
#include <iomanip>

using namespace std;

void cutflow_signal() {

  const int nCuts = 14;
  TString cutNames[nCuts] = {
    "All events",
    "JSON",
    "MET filters",
    "nPV $\\geq$ 1",
    "== 1 tight lepton",
    "==0 loose leptons",
    "HLT",
    "nJets $\\geq$ 3",
    "nBtags $\\geq$ 1",
    "N($\\gamma$, fake) $\\geq$ 1",
    "SR1",
    "CR1",
    "SR2",
    "CR2"};
  
  TFile * inputA = new TFile("cutflow_mst_460_m1_175.root", "READ");
  TFile * inputB = new TFile("cutflow_mst_560_m1_375.root", "READ");

  TH1D * h_photon_A = (TH1D*)inputA->Get("cutflow_photons");
  TH1D * h_photon_B = (TH1D*)inputB->Get("cutflow_photons");
  
  TH1D * h_ele_A = (TH1D*)inputA->Get("cutflow_ele");
  TH1D * h_ele_B = (TH1D*)inputB->Get("cutflow_ele");
  TH1D * h_muon_A = (TH1D*)inputA->Get("cutflow_muon");
  TH1D * h_muon_B = (TH1D*)inputB->Get("cutflow_muon");

  TFile * fAcceptanceA = new TFile ("signal_contamination_mst_460_m1_175.root", "READ");
  TFile * fAcceptanceB = new TFile ("signal_contamination_mst_560_m1_375.root", "READ");

  TTree * signalTree_ele_A = (TTree*)fAcceptanceA->Get("ele_bjj_signalTree");
  TTree * signalTree_ele_B = (TTree*)fAcceptanceB->Get("ele_bjj_signalTree");
  TTree * signalTree_muon_A = (TTree*)fAcceptanceA->Get("muon_bjj_signalTree");
  TTree * signalTree_muon_B = (TTree*)fAcceptanceB->Get("muon_bjj_signalTree");

  TTree * fakeTree_ele_A = (TTree*)fAcceptanceA->Get("ele_bjj_fakeTree");
  TTree * fakeTree_ele_B = (TTree*)fAcceptanceB->Get("ele_bjj_fakeTree");
  TTree * fakeTree_muon_A = (TTree*)fAcceptanceA->Get("muon_bjj_fakeTree");
  TTree * fakeTree_muon_B = (TTree*)fAcceptanceB->Get("muon_bjj_fakeTree");

  h_ele_A->SetBinContent(11,
			 (Double_t)signalTree_ele_A->Draw("Njets", "Ngamma == 1", "goff"));
  h_ele_B->SetBinContent(11,
			 (Double_t)signalTree_ele_B->Draw("Njets", "Ngamma == 1", "goff"));
  h_muon_A->SetBinContent(11,
			  (Double_t)signalTree_muon_A->Draw("Njets", "Ngamma == 1", "goff"));
  h_muon_B->SetBinContent(11,
			  (Double_t)signalTree_muon_B->Draw("Njets", "Ngamma == 1", "goff"));

  h_ele_A->SetBinContent(12,
			 (Double_t)fakeTree_ele_A->Draw("Njets", "Ngamma == 0 && Nfake == 1", "goff"));
  h_ele_B->SetBinContent(12,
			 (Double_t)fakeTree_ele_B->Draw("Njets", "Ngamma == 0 && Nfake == 1", "goff"));
  h_muon_A->SetBinContent(12,
			  (Double_t)fakeTree_muon_A->Draw("Njets", "Ngamma == 0 && Nfake == 1", "goff"));
  h_muon_B->SetBinContent(12,
			  (Double_t)fakeTree_muon_B->Draw("Njets", "Ngamma == 0 && Nfake == 1", "goff"));

  h_ele_A->SetBinContent(13,
			 (Double_t)signalTree_ele_A->Draw("Njets", "Ngamma >= 2", "goff"));
  h_ele_B->SetBinContent(13,
			 (Double_t)signalTree_ele_B->Draw("Njets", "Ngamma >= 2", "goff"));
  h_muon_A->SetBinContent(13,
			  (Double_t)signalTree_muon_A->Draw("Njets", "Ngamma >= 2", "goff"));
  h_muon_B->SetBinContent(13,
			  (Double_t)signalTree_muon_B->Draw("Njets", "Ngamma >= 2", "goff"));

  h_ele_A->SetBinContent(14,
			 (Double_t)fakeTree_ele_A->Draw("Njets", "Ngamma == 0 && Nfake >= 2", "goff"));
  h_ele_B->SetBinContent(14,
			 (Double_t)fakeTree_ele_B->Draw("Njets", "Ngamma == 0 && Nfake >= 2", "goff"));
  h_muon_A->SetBinContent(14,
			  (Double_t)fakeTree_muon_A->Draw("Njets", "Ngamma == 0 && Nfake >= 2", "goff"));
  h_muon_B->SetBinContent(14,
			  (Double_t)fakeTree_muon_B->Draw("Njets", "Ngamma == 0 && Nfake >= 2", "goff"));

  fAcceptanceA->Close();
  fAcceptanceB->Close();

  int thisEleValue, lastEleValue;
  int thisMuonValue, lastMuonValue;
  
  cout << "Table A (460_175)" << endl << endl;

  cout << "\\hline \\hline" << endl;
  cout << "Cut & \\multicolumn{2}{|c|}{Number of events (\\% of previous)} \\\\" << endl;
  cout << " & ele & muon \\\\" << endl;
  cout << "\\hline" << endl;
  for(int i = 0; i < nCuts; i++) {
    thisEleValue = h_ele_A->GetBinContent(i+1);
    thisMuonValue = h_muon_A->GetBinContent(i+1);

    cout << cutNames[i] << " & ";
    if(i == 0) {
      cout << thisEleValue << " & ";
      cout << thisMuonValue << " \\\\" << endl;
    }
    else {
      cout << thisEleValue << " (" << setprecision(0) << fixed << 100. * (double)thisEleValue / (double)lastEleValue << ") & ";
      cout << thisMuonValue << " (" << 100. * (double)thisMuonValue / (double)lastMuonValue << ") \\\\" << endl;
    }

    lastEleValue = thisEleValue;
    lastMuonValue = thisMuonValue;
  }
  cout << "\\hline" << endl;

  cout << endl << endl << endl;

  cout << "Table B (560_375)" << endl << endl;

  cout << "\\hline \\hline" << endl;
  cout << "Cut & \\multicolumn{2}{|c|}{Number of events (\\% of previous)} \\\\" << endl;
  cout << " & ele & muon \\\\" << endl;
  cout << "\\hline" << endl;
  for(int i = 0; i < nCuts; i++) {
    thisEleValue = h_ele_B->GetBinContent(i+1);
    thisMuonValue = h_muon_B->GetBinContent(i+1);

    cout << cutNames[i] << " & ";
    if(i == 0) {
      cout << thisEleValue << " & ";
      cout << thisMuonValue << " \\\\" << endl;
    }
    else {
      cout << thisEleValue << " (" << setprecision(0) << fixed << 100. * (double)thisEleValue / (double)lastEleValue << ") & ";
      cout << thisMuonValue << " (" << 100. * (double)thisMuonValue / (double)lastMuonValue << ") \\\\" << endl;
    }

    lastEleValue = thisEleValue;
    lastMuonValue = thisMuonValue;
  }
  cout << "\\hline" << endl;

  cout << endl << endl << endl;

  const int nPhotonCuts = 11;
  TString photonCutNames[nPhotonCuts] = {
    "All candidates",
    "$|\\eta| < 1.4442$",
    "$E_{T} > 20$ GeV",
    "$H/E < 0.05$",
    "Conv.-safe ele veto",
    "Neutral had. iso",
    "Photon iso",
    "Charged had. iso",
    "$\\sigma_{i\\eta i\\eta} < 0.012$",
    "$\\Delta R(\\gamma, \\mu) \\geq 0.7$",
    "$\\Delta R(\\gamma, e) \\geq 0.7$"};

  int thisValueA, thisValueB;
  int lastValueA, lastValueB;
  
  cout << "Photon Table" << endl << endl;

  cout << "\\hline \\hline" << endl;
  cout << "Cut & \\multicolumn{2}{|c|}{Number of photons (\\% of previous)} \\\\" << endl;
  cout << " & 460\\_175 & 560\\_375 \\\\" << endl;
  cout << "\\hline" << endl;
  for(int i = 0; i < nPhotonCuts; i++) {
    thisValueA = h_photon_A->GetBinContent(i+1);
    thisValueB = h_photon_B->GetBinContent(i+1);
    
    cout << photonCutNames[i] << " & ";
    if(i == 0) {
      cout << thisValueA << " & ";
      cout << thisValueB << " \\\\" << endl;
    }
    else {
      cout << thisValueA << " (" << setprecision(0) << fixed << 100. * (double)thisValueA / (double)lastValueA << ") & ";
      cout << thisValueB << " (" << 100. * (double)thisValueB / (double)lastValueB << ") \\\\" << endl;
    }

    lastValueA = thisValueA;
    lastValueB = thisValueB;
  }
  cout << "\\hline" << endl;
  
}
  
