#include <vector>

const int nMetBins = 5;
Double_t xbins[nMetBins+1] = {50, 75, 100, 150, 300, 800};

void CutSR1LowMET() {

  rebinStuff("ele_SR1");
  rebinStuff("muon_SR1");

  copyStuff("ele_SR2");
  copyStuff("muon_SR2");
  
}

void rebinStuff(TString channel) {

  vector<TH1D*> rebinnedHistos;
  vector<TString> nameList;
  
  TFile * input = new TFile("/uscms_data/d2/bfrancis/stopBino/limits/CMSSW_7_1_5/src/LimitSetting/condor_submit/limitInputs_bjj.root", "READ");
  input->cd(channel);
  
  TIter next(gDirectory->GetListOfKeys());
  TKey * key;

  while((key = (TKey*)next())) {
    TClass * cl = gROOT->GetClass(key->GetClassName());
    if(!cl->InheritsFrom("TH1D")) continue;

    TH1D * h = (TH1D*)key->ReadObj();
    TString hName = key->GetName();

    rebinnedHistos.push_back((TH1D*)h->Rebin(nMetBins, hName+"_reb", xbins));
    nameList.push_back(hName);
    
  }

  TFile * output = new TFile("limitInputs_bjj_reb.root", "UPDATE");

  output->mkdir(channel);
  output->cd(channel);

  for(unsigned int i = 0; i < rebinnedHistos.size(); i++) {
    rebinnedHistos[i]->Write(nameList[i]);
  }

  output->Close();
  input->Close();
  
}

void copyStuff(TString channel) {

  vector<TH1D*> clonedHistos;
  vector<TString> nameList;
  
  TFile * input = new TFile("/uscms_data/d2/bfrancis/stopBino/limits/CMSSW_7_1_5/src/LimitSetting/condor_submit/limitInputs_bjj.root", "READ");
  input->cd(channel);
  
  TIter next(gDirectory->GetListOfKeys());
  TKey * key;

  while((key = (TKey*)next())) {
    TClass * cl = gROOT->GetClass(key->GetClassName());
    if(!cl->InheritsFrom("TH1D")) continue;

    TH1D * h = (TH1D*)key->ReadObj();
    TString hName = key->GetName();

    clonedHistos.push_back((TH1D*)h->Clone(hName+"_clone"));
    nameList.push_back(hName);
    
  }

  TFile * output = new TFile("limitInputs_bjj_reb.root", "UPDATE");

  output->mkdir(channel);
  output->cd(channel);

  for(unsigned int i = 0; i < clonedHistos.size(); i++) {
    clonedHistos[i]->Write(nameList[i]);
  }

  output->Close();
  input->Close();
  
}
