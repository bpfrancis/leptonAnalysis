#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"

#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>

using namespace std;

class HistogramMaker : public TObject {

  ClassDef(HistogramMaker, 1);

 public:
  // HistogramMaker("ele_bjj", "noSigmaIetaIeta", 1, 1, -1.)
  HistogramMaker(TString chan, TString treeType, int nRequiredPhotons, double cutMET, double cutSigmaIetaIeta);
  ~HistogramMaker();

  void BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi, bool doAbsValue = false);
  void BookHistogram(TString variable, Int_t nBins, Double_t* customBins, bool doAbsValue = false);

  void BookHistogram2D(TString var_x, TString var_y, Int_t nBins_x, Float_t xlo, Float_t xhi, Int_t nBins_y, Float_t ylo, Float_t yhi, Float_t zlo = 0.0, Float_t zhi = -1.0);

  void LoadLeptonSFs(TString fileName);
  void LoadPhotonSFs(TString fileName);

  void HistogramData(TString input);
  void HistogramQCD(TString input);
  
  void HistogramMCBackground(TString fileName, TString scanName,
			     Double_t xsec, bool removeTTA, bool reweightTop);
  void HistogramSignal(TString input, TString output);

  void ResetHistograms() {
    for(unsigned int i = 0; i < variables.size(); i++) {
      histograms[i]->Reset();
      histograms_btagUp[i]->Reset();
      histograms_btagDown[i]->Reset();
      histograms_pileupUp[i]->Reset();
      histograms_pileupDown[i]->Reset();
      histograms_topPtUp[i]->Reset();
      histograms_topPtDown[i]->Reset();
      histograms_leptonSFUp[i]->Reset();
      histograms_leptonSFDown[i]->Reset();
      histograms_photonSFUp[i]->Reset();
      histograms_photonSFDown[i]->Reset();
      histograms_JECUp[i]->Reset();
      histograms_JECDown[i]->Reset();
      histograms_mcQCD[i]->Reset();
    }
    for(unsigned int i = 0; i < variables_2d.size(); i++) histograms_2d[i]->Reset();
  }
  
  void GetLeptonSF(Float_t lepton_pt, Float_t lepton_eta, TString chan, Float_t& central, Float_t& up, Float_t& down);
  void GetLeptonSF(vector<Float_t> vars, TString chan, Float_t& central, Float_t& up, Float_t& down) { 
    if(chan.Contains("ele")) GetLeptonSF(varNumber["ele_pt"], varNumber["ele_eta"], chan, central, up, down);
    else GetLeptonSF(varNumber["muon_pt"], varNumber["muon_eta"], chan, central, up, down);
  };
  
  void GetPhotonSF(Float_t lead_photon_et, Float_t lead_photon_eta, Float_t trail_photon_et, Float_t trail_photon_eta, Float_t nphotons, 
		   Float_t& central, Float_t& up, Float_t& down);
  void GetPhotonSF(vector<Float_t> vars, Float_t& central, Float_t& up, Float_t& down) { 
    if(varNumber["Nphotons"] == 0) {
      central = 1.;
      up = 1.;
      down = 1.;
    }
    else if(varNumber["Nphotons"] == 1 && vars.size() >= 21) GetPhotonSF(varNumber["leadPhotonEt"], varNumber["leadPhotonEta"], -1., -1., varNumber["Nphotons"], central, up, down); 
    else if(vars.size() >= 27) GetPhotonSF(varNumber["leadPhotonEt"], varNumber["leadPhotonEta"], varNumber["trailPhotonEt"], varNumber["trailPhotonEta"], varNumber["Nphotons"], central, up, down);
  };

  void SetSigmaIetaIetaCut(double cut) { sigmaIetaIetaCut = cut; }

 private:

  vector<TString> variables;
  vector<pair<TString, TString> > variables_2d;

  map<TString, unsigned int> varNumber;

  vector<TH1D*> histograms;
  vector<TH1D*> histograms_btagUp;
  vector<TH1D*> histograms_btagDown;
  vector<TH1D*> histograms_pileupUp;
  vector<TH1D*> histograms_pileupDown;
  vector<TH1D*> histograms_topPtUp;
  vector<TH1D*> histograms_topPtDown;
  vector<TH1D*> histograms_leptonSFUp;
  vector<TH1D*> histograms_leptonSFDown;
  vector<TH1D*> histograms_photonSFUp;
  vector<TH1D*> histograms_photonSFDown;
  vector<TH1D*> histograms_JECUp;
  vector<TH1D*> histograms_JECDown;
  vector<TH1D*> histograms_mcQCD;

  vector<TH2D*> histograms_2d;

  TFile * fLeptonSF;
  TH2D * sf_muon;
  TH2D * sf_electron;
  TH2D * sf_SingleElectronTrigger;

  TFile * fPhotonSF;
  TH2D * sf_photon_id;
  TH2D * sf_photon_veto;

  Int_t intLumi;

  TString channel;
  TString photonType;
  int photonReq;
  double metCut;
  double sigmaIetaIetaCut;
  vector<bool> absValue;
};

HistogramMaker::HistogramMaker(TString chan, TString treeType, int nRequiredPhotons, double cutMET, double cutSigmaIetaIeta) {
  
  channel = chan;
  photonType = treeType;
  photonReq = nRequiredPhotons;
  metCut = cutMET;
  sigmaIetaIetaCut = cutSigmaIetaIeta;
  intLumi = 19712;

  variables.clear();
  variables_2d.clear();
  absValue.clear();

  histograms.clear();
  histograms_btagUp.clear();
  histograms_btagDown.clear();
  histograms_pileupUp.clear();
  histograms_pileupDown.clear();
  histograms_topPtUp.clear();
  histograms_topPtDown.clear();
  histograms_leptonSFUp.clear();
  histograms_leptonSFDown.clear();
  histograms_photonSFUp.clear();
  histograms_photonSFDown.clear();
  histograms_JECUp.clear();
  histograms_JECDown.clear();
  histograms_mcQCD.clear();

}

HistogramMaker::~HistogramMaker() {

  variables.clear();
  variables_2d.clear();
  absValue.clear();

  histograms.clear();
  histograms_btagUp.clear();
  histograms_btagDown.clear();
  histograms_pileupUp.clear();
  histograms_pileupDown.clear();
  histograms_topPtUp.clear();
  histograms_topPtDown.clear();
  histograms_leptonSFUp.clear();
  histograms_leptonSFDown.clear();
  histograms_photonSFUp.clear();
  histograms_photonSFDown.clear();
  histograms_JECUp.clear();
  histograms_JECDown.clear();
  histograms_mcQCD.clear();

  fLeptonSF->Close();
  fPhotonSF->Close();

}

void HistogramMaker::HistogramData(TString input) {

  ResetHistograms();

  vector<Float_t> vars;
  vars.resize(variables.size());

  TFile * f = new TFile(input, "READ");
  TTree * tree = (TTree*)f->Get(channel+"_"+photonType+"Tree");

  for(unsigned int i = 0; i < variables.size(); i++) tree->SetBranchAddress(variables[i], &(vars[i]));

  for(int nEvent = 0; nEvent < tree->GetEntries(); nEvent++) {
    tree->GetEntry(nEvent);

    if(metCut > 0. && varNumber["pfMET"] >= metCut) continue;
    
    for(unsigned int i = 0; i < variables.size(); i++) {
      if(absValue[i]) vars[i] = fabs(vars[i]);
      if(variables[i] != "Nphotons" && (int)varNumber["Nphotons"] != photonReq && photonReq >= 0) continue;

      if(variables[i] != "leadSigmaIetaIeta" && varNumber["Nphotons"] >= 1) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["leadSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["leadSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      if(variables[i] != "trailSigmaIetaIeta" && varNumber["Nphotons"] >= 2) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["trailSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["trailSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      histograms[i]->Fill(vars[i]);

      for(unsigned int j = 0; j < variables_2d.size(); j++) {
	if(variables[i] == variables_2d[j].first) {
	  for(unsigned int k = 0; k < variables_2d.size(); k++) {
	    if(variables[k] == variables_2d[j].second) {
	      histograms_2d[j]->Fill(vars[i], vars[k]);
	    }
	  }
	}
      }

    }

  }

  delete tree;

  TString outName = "data_"+channel+"_";
  outName += (sigmaIetaIetaCut >= 0) ? "gamma" : "fake";
  outName += ".root";

  TFile * out = new TFile(outName, "RECREATE");

  TH1D * hOut;
  for(unsigned int i = 0; i < variables.size(); i++) {
    hOut = (TH1D*)histograms[i]->Clone(variables[i]);
    hOut->Write();
  }

  TH2D * hOut_2d;
  for(unsigned int i = 0; i < variables_2d.size(); i++) {
    hOut_2d = (TH2D*)histograms_2d[i]->Clone(variables_2d[i].second+"_"+variables_2d[i].first);
    hOut_2d->Write();
  }

  out->Close();

  f->Close();
}
  
void HistogramMaker::HistogramQCD(TString input) {

  ResetHistograms();

  vector<Float_t> vars;
  vars.resize(variables.size());

  TFile * f = new TFile(input, "READ");

  TString treeName = channel+"_";
  treeName += (channel.Contains("ele")) ? "eQCD" : "muQCD";
  if(!(photonType.Contains("signal"))) treeName += photonType;
  treeName += "Tree";

  TTree * tree = (TTree*)f->Get(treeName);

  for(unsigned int i = 0; i < variables.size(); i++) tree->SetBranchAddress(variables[i], &(vars[i]));

  for(int nEvent = 0; nEvent < tree->GetEntries(); nEvent++) {
    tree->GetEntry(nEvent);

    if(metCut > 0. && vars[1] >= metCut) continue;
    
    for(unsigned int i = 0; i < variables.size(); i++) {
      if(variables[i] != "Nphotons" && (int)varNumber["Nphotons"] != photonReq && photonReq >= 0) continue;
      if(absValue[i]) vars[i] = fabs(vars[i]);

      if(variables[i] != "leadSigmaIetaIeta" && varNumber["Nphotons"] >= 1) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["leadSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["leadSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      if(variables[i] != "trailSigmaIetaIeta" && varNumber["Nphotons"] >= 2) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["trailSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["trailSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      histograms[i]->Fill(vars[i]);

      for(unsigned int j = 0; j < variables_2d.size(); j++) {
	if(variables[i] == variables_2d[j].first) {
	  for(unsigned int k = 0; k < variables_2d.size(); k++) {
	    if(variables[k] == variables_2d[j].second) {
	      histograms_2d[j]->Fill(vars[i], vars[k]);
	    }
	  }
	}
      }

    }

  }

  delete tree;

  TString outName = "qcd_"+channel+"_";
  outName += (sigmaIetaIetaCut >= 0) ? "gamma" : "fake";
  outName += ".root";

  TFile * out = new TFile(outName, "RECREATE");

  TH1D * hOut;

  for(unsigned int i = 0; i < variables.size(); i++) {
    hOut = (TH1D*)histograms[i]->Clone(variables[i]);
    hOut->Write();
  }

  TH2D * hOut_2d;
  for(unsigned int i = 0; i < variables_2d.size(); i++) {
    hOut_2d = (TH2D*)histograms_2d[i]->Clone(variables_2d[i].second+"_"+variables_2d[i].first);
    hOut_2d->Write();
  }

  out->Close();

  f->Close();
}

void HistogramMaker::HistogramMCBackground(TString fileName, TString scanName,
					   Double_t xsec, bool removeTTA, bool reweightTop) {

  ResetHistograms();

  vector<Float_t> vars;
  vars.resize(variables.size());
  
  TFile * f = new TFile(fileName, "READ");
  TTree * tree = (TTree*)f->Get(channel+"_"+photonType+"Tree");
  TTree * tree_JECUp = (TTree*)f->Get(channel+"_"+photonType+"Tree_JECup");
  TTree * tree_JECDown = (TTree*)f->Get(channel+"_"+photonType+"Tree_JECdown");

  TString qcdTreeName = channel+"_";
  qcdTreeName += (channel.Contains("ele")) ? "eQCD" : "muQCD";
  if(!(photonType.Contains("signal"))) qcdTreeName += photonType;
  qcdTreeName += "Tree";

  TTree * tree_qcd = (TTree*)f->Get(qcdTreeName);

  for(unsigned int i = 0; i < variables.size(); i++) {
    tree->SetBranchAddress(variables[i], &(vars[i]));
    tree_JECUp->SetBranchAddress(variables[i], &(vars[i]));
    tree_JECDown->SetBranchAddress(variables[i], &(vars[i]));
    tree_qcd->SetBranchAddress(variables[i], &(vars[i]));
  }

  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr;
  Float_t puWeightUp, puWeightDown, btagWeightUp, btagWeightDown;
  Float_t overlaps_ttA;
  Float_t topPtReweighting;

  tree->SetBranchAddress("pileupWeight", &puWeight);
  tree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree->SetBranchAddress("btagWeight", &btagWeight);
  tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  tree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  tree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  tree->SetBranchAddress("pileupWeightUp", &puWeightUp);
  tree->SetBranchAddress("pileupWeightDown", &puWeightDown);
  if(removeTTA) tree->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
  if(reweightTop) tree->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  tree_JECUp->SetBranchAddress("pileupWeight", &puWeight);
  tree_JECUp->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree_JECUp->SetBranchAddress("btagWeight", &btagWeight);
  tree_JECUp->SetBranchAddress("btagWeightErr", &btagWeightErr);
  if(removeTTA) tree_JECUp->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
  if(reweightTop) tree_JECUp->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  tree_JECDown->SetBranchAddress("pileupWeight", &puWeight);
  tree_JECDown->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree_JECDown->SetBranchAddress("btagWeight", &btagWeight);
  tree_JECDown->SetBranchAddress("btagWeightErr", &btagWeightErr);
  if(removeTTA) tree_JECDown->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
  if(reweightTop) tree_JECDown->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  tree_qcd->SetBranchAddress("pileupWeight", &puWeight);
  tree_qcd->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree_qcd->SetBranchAddress("btagWeight", &btagWeight);
  tree_qcd->SetBranchAddress("btagWeightErr", &btagWeightErr);
  if(removeTTA) tree_qcd->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
  if(reweightTop) tree_qcd->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  Float_t leptonSF, leptonSFUp, leptonSFDown;
  Float_t photonSF, photonSFUp, photonSFDown;

  for(int nEvent = 0; nEvent < tree->GetEntries(); nEvent++) {
    tree->GetEntry(nEvent);

    if(!(channel.Contains("b"))) {
      btagWeight = 1.;
      btagWeightErr = 0.;
      btagWeightUp = 1.;
      btagWeightDown = 1.;
    }

    if(btagWeightUp < 0) btagWeightUp = 0.;
    if(btagWeightDown < 0) btagWeightDown = 0;
    if(btagWeight != btagWeight) continue;
    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    if(topPtReweighting < 0) topPtReweighting = 1.;

    if(removeTTA && overlaps_ttA > 0.001) continue;

    if(metCut > 0. && vars[1] >= metCut) continue;
    
    GetLeptonSF(vars, channel, leptonSF, leptonSFUp, leptonSFDown);
    GetPhotonSF(vars, photonSF, photonSFUp, photonSFDown);

    Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
    Float_t addError2_puOnly = btagWeight*btagWeight*puWeightErr*puWeightErr;
    Float_t addError2_btagOnly =  puWeight*puWeight*btagWeightErr*btagWeightErr;

    for(unsigned int i = 0; i < variables.size(); i++) {
      if(variables[i] != "Nphotons" && (int)varNumber["Nphotons"] != photonReq && photonReq >= 0) continue;
      if(absValue[i]) vars[i] = fabs(vars[i]);

      if(variables[i] != "leadSigmaIetaIeta" && varNumber["Nphotons"] >= 1) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["leadSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["leadSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      if(variables[i] != "trailSigmaIetaIeta" && varNumber["Nphotons"] >= 2) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["trailSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["trailSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting;

      Float_t oldError = histograms[i]->GetBinError(histograms[i]->FindBin(vars[i]));
      Float_t newerror = sqrt(oldError*oldError + addError2);
      histograms[i]->Fill(vars[i], totalWeight);
      histograms[i]->SetBinError(histograms[i]->FindBin(vars[i]), newerror);

      for(unsigned int j = 0; j < variables_2d.size(); j++) {
	if(variables[i] == variables_2d[j].first) {
	  for(unsigned int k = 0; k < variables_2d.size(); k++) {
	    if(variables[k] == variables_2d[j].second) {
	      histograms_2d[j]->Fill(vars[i], vars[k], totalWeight);
	    }
	  }
	}
      }

      totalWeight = puWeight * btagWeightUp * leptonSF * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_btagUp[i]->GetBinError(histograms_btagUp[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2_puOnly);
      histograms_btagUp[i]->Fill(vars[i], totalWeight);
      histograms_btagUp[i]->SetBinError(histograms_btagUp[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeight * btagWeightDown * leptonSF * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_btagDown[i]->GetBinError(histograms_btagDown[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2_puOnly);
      histograms_btagDown[i]->Fill(vars[i], totalWeight);
      histograms_btagDown[i]->SetBinError(histograms_btagDown[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeightUp * btagWeight * leptonSF * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_pileupUp[i]->GetBinError(histograms_pileupUp[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2_btagOnly);
      histograms_pileupUp[i]->Fill(vars[i], totalWeight);
      histograms_pileupUp[i]->SetBinError(histograms_pileupUp[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeightDown * btagWeight * leptonSF * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_pileupDown[i]->GetBinError(histograms_pileupDown[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2_btagOnly);
      histograms_pileupDown[i]->Fill(vars[i], totalWeight);
      histograms_pileupDown[i]->SetBinError(histograms_pileupDown[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeight * btagWeight * leptonSFUp * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_leptonSFUp[i]->GetBinError(histograms_leptonSFUp[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2);
      histograms_leptonSFUp[i]->Fill(vars[i], totalWeight);
      histograms_leptonSFUp[i]->SetBinError(histograms_leptonSFUp[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeight * btagWeight * leptonSFDown * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_leptonSFDown[i]->GetBinError(histograms_leptonSFDown[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2);
      histograms_leptonSFDown[i]->Fill(vars[i], totalWeight);
      histograms_leptonSFDown[i]->SetBinError(histograms_leptonSFDown[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFUp;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_photonSFUp[i]->GetBinError(histograms_photonSFUp[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2);
      histograms_photonSFUp[i]->Fill(vars[i], totalWeight);
      histograms_photonSFUp[i]->SetBinError(histograms_photonSFUp[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFDown;
      if(reweightTop) totalWeight *= topPtReweighting;
      oldError = histograms_photonSFDown[i]->GetBinError(histograms_photonSFDown[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2);
      histograms_photonSFDown[i]->Fill(vars[i], totalWeight);
      histograms_photonSFDown[i]->SetBinError(histograms_photonSFDown[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTop) totalWeight *= topPtReweighting * topPtReweighting;
      oldError = histograms_topPtUp[i]->GetBinError(histograms_topPtUp[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2);
      histograms_topPtUp[i]->Fill(vars[i], totalWeight);
      histograms_topPtUp[i]->SetBinError(histograms_topPtUp[i]->FindBin(vars[i]), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      oldError = histograms_topPtDown[i]->GetBinError(histograms_topPtDown[i]->FindBin(vars[i]));
      newerror = sqrt(oldError*oldError + addError2);
      histograms_topPtDown[i]->Fill(vars[i], totalWeight);
      histograms_topPtDown[i]->SetBinError(histograms_topPtDown[i]->FindBin(vars[i]), newerror);

    }

  }

  delete tree;

  for(int nEvent = 0; nEvent < tree_JECUp->GetEntries(); nEvent++) {
    tree_JECUp->GetEntry(nEvent);
    
    if(!(channel.Contains("b"))) {
      btagWeight = 1.;
      btagWeightErr = 0.;
      btagWeightUp = 1.;
      btagWeightDown = 1.;
    }

    if(btagWeightUp < 0) btagWeightUp = 0.;
    if(btagWeightDown < 0) btagWeightDown = 0;
    if(btagWeight != btagWeight) continue;
    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    if(topPtReweighting < 0) topPtReweighting = 1.;

    if(removeTTA && overlaps_ttA > 0.001) continue;

    if(metCut > 0. && vars[1] >= metCut) continue;
    
    GetLeptonSF(vars, channel, leptonSF, leptonSFUp, leptonSFDown);
    GetPhotonSF(vars, photonSF, photonSFUp, photonSFDown);

    double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
    if(reweightTop) totalWeight *= topPtReweighting;

    for(unsigned int i = 0; i < variables.size(); i++) {
      if(variables[i] != "Nphotons" && (int)varNumber["Nphotons"] != photonReq && photonReq >= 0) continue;
      if(absValue[i]) vars[i] = fabs(vars[i]);
      if(variables[i] != "leadSigmaIetaIeta" && varNumber["Nphotons"] >= 1) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["leadSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["leadSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }
      if(variables[i] != "trailSigmaIetaIeta" && varNumber["Nphotons"] >= 2) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["trailSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["trailSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      histograms_JECUp[i]->Fill(vars[i], totalWeight);
    }

  }

  delete tree_JECUp;

  for(int nEvent = 0; nEvent < tree_JECDown->GetEntries(); nEvent++) {
    tree_JECDown->GetEntry(nEvent);
    
    if(!(channel.Contains("b"))) {
      btagWeight = 1.;
      btagWeightErr = 0.;
      btagWeightUp = 1.;
      btagWeightDown = 1.;
    }

    if(btagWeightUp < 0) btagWeightUp = 0.;
    if(btagWeightDown < 0) btagWeightDown = 0;
    if(btagWeight != btagWeight) continue;
    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    if(topPtReweighting < 0) topPtReweighting = 1.;

    if(removeTTA && overlaps_ttA > 0.001) continue;

    if(metCut > 0. && vars[1] >= metCut) continue;
    
    GetLeptonSF(vars, channel, leptonSF, leptonSFUp, leptonSFDown);
    GetPhotonSF(vars, photonSF, photonSFUp, photonSFDown);

    double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
    if(reweightTop) totalWeight *= topPtReweighting;

    for(unsigned int i = 0; i < variables.size(); i++) {
      if(variables[i] != "Nphotons" && (int)varNumber["Nphotons"] != photonReq && photonReq >= 0) continue;
      if(absValue[i]) vars[i] = fabs(vars[i]);
      if(variables[i] != "leadSigmaIetaIeta" && varNumber["Nphotons"] >= 1) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["leadSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["leadSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }
      if(variables[i] != "trailSigmaIetaIeta" && varNumber["Nphotons"] >= 2) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["trailSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["trailSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      histograms_JECDown[i]->Fill(vars[i], totalWeight);
    }

  }

  delete tree_JECDown;

  for(int nEvent = 0; nEvent < tree_qcd->GetEntries(); nEvent++) {
    tree_qcd->GetEntry(nEvent);
    
    if(!(channel.Contains("b"))) {
      btagWeight = 1.;
      btagWeightErr = 0.;
      btagWeightUp = 1.;
      btagWeightDown = 1.;
    }

    if(btagWeightUp < 0) btagWeightUp = 0.;
    if(btagWeightDown < 0) btagWeightDown = 0;
    if(btagWeight != btagWeight) continue;
    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    if(topPtReweighting < 0) topPtReweighting = 1.;

    if(removeTTA && overlaps_ttA > 0.001) continue;

    if(metCut > 0. && vars[1] >= metCut) continue;
    
    GetLeptonSF(vars, channel, leptonSF, leptonSFUp, leptonSFDown);
    GetPhotonSF(vars, photonSF, photonSFUp, photonSFDown);

    double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
    if(reweightTop) totalWeight *= topPtReweighting;

    Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;

    for(unsigned int i = 0; i < variables.size(); i++) {
      if(variables[i] != "Nphotons" && (int)varNumber["Nphotons"] != photonReq && photonReq >= 0) continue;
      if(absValue[i]) vars[i] = fabs(vars[i]);

      if(variables[i] != "leadSigmaIetaIeta" && varNumber["Nphotons"] >= 1) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["leadSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["leadSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      if(variables[i] != "trailSigmaIetaIeta" && varNumber["Nphotons"] >= 2) {
	// negative cut means invert (ie require fake), so reject passing photons
	if(sigmaIetaIetaCut < 0 && varNumber["trailSigmaIetaIeta"] < -1. * sigmaIetaIetaCut) continue;
	if(sigmaIetaIetaCut >= 0 && varNumber["trailSigmaIetaIeta"] >= sigmaIetaIetaCut) continue;
      }

      Float_t oldError = histograms_mcQCD[i]->GetBinError(histograms_mcQCD[i]->FindBin(vars[i]));
      Float_t newerror = sqrt(oldError*oldError + addError2);
      histograms_mcQCD[i]->Fill(vars[i], totalWeight);
      histograms_mcQCD[i]->SetBinError(histograms_mcQCD[i]->FindBin(vars[i]), newerror);
    }

  }

  delete tree_qcd;

  // remember to scale by xsec

  TH1D * h_nGen = (TH1D*)f->Get("nEvents_"+scanName);
  Double_t mcNGen = h_nGen->Integral();

  Double_t xsecScaling = intLumi * xsec / mcNGen;
  for(unsigned int i = 0; i < variables.size(); i++) {
    histograms[i]->Scale(xsecScaling);
    histograms_btagUp[i]->Scale(xsecScaling);
    histograms_btagDown[i]->Scale(xsecScaling);
    histograms_pileupUp[i]->Scale(xsecScaling);
    histograms_pileupDown[i]->Scale(xsecScaling);
    histograms_topPtUp[i]->Scale(xsecScaling);
    histograms_topPtDown[i]->Scale(xsecScaling);
    histograms_leptonSFUp[i]->Scale(xsecScaling);
    histograms_leptonSFDown[i]->Scale(xsecScaling);
    histograms_photonSFUp[i]->Scale(xsecScaling);
    histograms_photonSFDown[i]->Scale(xsecScaling);
    histograms_JECUp[i]->Scale(xsecScaling);
    histograms_JECDown[i]->Scale(xsecScaling);
    histograms_mcQCD[i]->Scale(xsecScaling);
  }

  for(unsigned int i = 0; i < variables_2d.size(); i++) histograms_2d[i]->Scale(xsecScaling);

  TString outName = scanName+"_"+channel+"_";
  outName += (sigmaIetaIetaCut >= 0) ? "gamma" : "fake";
  outName += ".root";

  TFile * out = new TFile(outName, "RECREATE");

  TH1D * hOut;
  for(unsigned int i = 0; i < variables.size(); i++) {
    hOut = (TH1D*)histograms[i]->Clone(variables[i]);
    hOut->Write();

    hOut = (TH1D*)histograms_btagUp[i]->Clone(variables[i]+"_btagUp");
    hOut->Write();

    hOut = (TH1D*)histograms_btagDown[i]->Clone(variables[i]+"_btagDown");
    hOut->Write();

    hOut = (TH1D*)histograms_pileupUp[i]->Clone(variables[i]+"_pileupUp");
    hOut->Write();

    hOut = (TH1D*)histograms_pileupDown[i]->Clone(variables[i]+"_pileupDown");
    hOut->Write();

    hOut = (TH1D*)histograms_leptonSFUp[i]->Clone(variables[i]+"_leptonSFUp");
    hOut->Write();

    hOut = (TH1D*)histograms_leptonSFDown[i]->Clone(variables[i]+"_leptonSFDown");
    hOut->Write();

    hOut = (TH1D*)histograms_photonSFUp[i]->Clone(variables[i]+"_photonSFUp");
    hOut->Write();

    hOut = (TH1D*)histograms_photonSFDown[i]->Clone(variables[i]+"_photonSFDown");
    hOut->Write();

    hOut = (TH1D*)histograms_topPtUp[i]->Clone(variables[i]+"_topPtUp");
    hOut->Write();

    hOut = (TH1D*)histograms_topPtDown[i]->Clone(variables[i]+"_topPtDown");
    hOut->Write();

    hOut = (TH1D*)histograms_JECUp[i]->Clone(variables[i]+"_JECUp");
    hOut->Write();

    hOut = (TH1D*)histograms_JECDown[i]->Clone(variables[i]+"_JECDown");
    hOut->Write();

    hOut = (TH1D*)histograms_mcQCD[i]->Clone(variables[i]+"_mcQCD");
    hOut->Write();

  }

  TH2D * hOut_2d;
  for(unsigned int i = 0; i < variables_2d.size(); i++) {
    hOut_2d = (TH2D*)histograms_2d[i]->Clone(variables_2d[i].second+"_"+variables_2d[i].first);
    hOut_2d->Write();
  }

  out->Close();

  f->Close();

};

void HistogramMaker::HistogramSignal(TString input, TString output) {};

void HistogramMaker::LoadLeptonSFs(TString fileName) {
  
  fLeptonSF = new TFile(fileName, "READ");
  sf_electron = (TH2D*)fLeptonSF->Get("TightEleIdIsoSF");
  sf_SingleElectronTrigger = (TH2D*)fLeptonSF->Get("TightEleTriggerSF");

  sf_muon = (TH2D*)fLeptonSF->Get("mu_pt_eta_full_id_iso_hlt_8TeV");

}

void HistogramMaker::LoadPhotonSFs(TString fileName) {
  
  fPhotonSF = new TFile(fileName, "READ");
  sf_photon_id = (TH2D*)fPhotonSF->Get("PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
  sf_photon_veto = (TH2D*)fPhotonSF->Get("PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");

}

void HistogramMaker::BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi, bool doAbsValue) {
  
  variables.push_back(variable);
  varNumber[variable] = variables.size();

  absValue.push_back(doAbsValue);
  if(doAbsValue) {
    xlo = 0.;
    nBins *= 2;
  }

  TH1D * h = new TH1D("h_"+variable, "h_"+variable, nBins, xlo, xhi);
  h->Sumw2();
  histograms.push_back(h);

  TH1D * h_btagUp = new TH1D("h_btagUp_"+variable, "h_btagUp_"+variable, nBins, xlo, xhi);
  h_btagUp->Sumw2();
  histograms_btagUp.push_back(h_btagUp);

  TH1D * h_btagDown = new TH1D("h_btagDown_"+variable, "h_btagDown_"+variable, nBins, xlo, xhi);
  h_btagDown->Sumw2();
  histograms_btagDown.push_back(h_btagDown);

  TH1D * h_pileupUp = new TH1D("h_pileupUp_"+variable, "h_pileupUp_"+variable, nBins, xlo, xhi);
  h_pileupUp->Sumw2();
  histograms_pileupUp.push_back(h_pileupUp);

  TH1D * h_pileupDown = new TH1D("h_pileupDown_"+variable, "h_pileupDown_"+variable, nBins, xlo, xhi);
  h_pileupDown->Sumw2();
  histograms_pileupDown.push_back(h_pileupDown);

  TH1D * h_topPtUp = new TH1D("h_topPtUp_"+variable, "h_topPtUp_"+variable, nBins, xlo, xhi);
  h_topPtUp->Sumw2();
  histograms_topPtUp.push_back(h_topPtUp);

  TH1D * h_topPtDown = new TH1D("h_topPtDown_"+variable, "h_topPtDown_"+variable, nBins, xlo, xhi);
  h_topPtDown->Sumw2();
  histograms_topPtDown.push_back(h_topPtDown);

  TH1D * h_leptonSFUp = new TH1D("h_leptonSFUp_"+variable, "h_leptonSFUp_"+variable, nBins, xlo, xhi);
  h_leptonSFUp->Sumw2();
  histograms_leptonSFUp.push_back(h_leptonSFUp);

  TH1D * h_leptonSFDown = new TH1D("h_leptonSFDown_"+variable, "h_leptonSFDown_"+variable, nBins, xlo, xhi);
  h_leptonSFDown->Sumw2();
  histograms_leptonSFDown.push_back(h_leptonSFDown);

  TH1D * h_photonSFUp = new TH1D("h_photonSFUp_"+variable, "h_photonSFUp_"+variable, nBins, xlo, xhi);
  h_photonSFUp->Sumw2();
  histograms_photonSFUp.push_back(h_photonSFUp);

  TH1D * h_photonSFDown = new TH1D("h_photonSFDown_"+variable, "h_photonSFDown_"+variable, nBins, xlo, xhi);
  h_photonSFDown->Sumw2();
  histograms_photonSFDown.push_back(h_photonSFDown);

  TH1D * h_JECUp = new TH1D("h_JECUp_"+variable, "h_JECUp_"+variable, nBins, xlo, xhi);
  h_JECUp->Sumw2();
  histograms_JECUp.push_back(h_JECUp);

  TH1D * h_JECDown = new TH1D("h_JECDown_"+variable, "h_JECDown_"+variable, nBins, xlo, xhi);
  h_JECDown->Sumw2();
  histograms_JECDown.push_back(h_JECDown);

  TH1D * h_mcQCD = new TH1D("h_mcQCD_"+variable, "h_mcQCD_"+variable, nBins, xlo, xhi);
  h_mcQCD->Sumw2();
  histograms_mcQCD.push_back(h_mcQCD);

}

void HistogramMaker::BookHistogram(TString variable, Int_t nBins, Double_t* customBins, bool doAbsValue) {

  variables.push_back(variable);
  varNumber[variable] = variables.size();

  absValue.push_back(doAbsValue);

  TH1D * h = new TH1D("h_"+variable, "h_"+variable, nBins, customBins);
  h->Sumw2();
  histograms.push_back(h);
  
  TH1D * h_btagUp = new TH1D("h_btagUp_"+variable, "h_btagUp_"+variable, nBins, customBins);
  h_btagUp->Sumw2();
  histograms_btagUp.push_back(h_btagUp);

  TH1D * h_btagDown = new TH1D("h_btagDown_"+variable, "h_btagDown_"+variable, nBins, customBins);
  h_btagDown->Sumw2();
  histograms_btagDown.push_back(h_btagDown);

  TH1D * h_pileupUp = new TH1D("h_pileupUp_"+variable, "h_pileupUp_"+variable, nBins, customBins);
  h_pileupUp->Sumw2();
  histograms_pileupUp.push_back(h_pileupUp);

  TH1D * h_pileupDown = new TH1D("h_pileupDown_"+variable, "h_pileupDown_"+variable, nBins, customBins);
  h_pileupDown->Sumw2();
  histograms_pileupDown.push_back(h_pileupDown);

  TH1D * h_topPtUp = new TH1D("h_topPtUp_"+variable, "h_topPtUp_"+variable, nBins, customBins);
  h_topPtUp->Sumw2();
  histograms_topPtUp.push_back(h_topPtUp);

  TH1D * h_topPtDown = new TH1D("h_topPtDown_"+variable, "h_topPtDown_"+variable, nBins, customBins);
  h_topPtDown->Sumw2();
  histograms_topPtDown.push_back(h_topPtDown);

  TH1D * h_leptonSFUp = new TH1D("h_leptonSFUp_"+variable, "h_leptonSFUp_"+variable, nBins, customBins);
  h_leptonSFUp->Sumw2();
  histograms_leptonSFUp.push_back(h_leptonSFUp);

  TH1D * h_leptonSFDown = new TH1D("h_leptonSFDown_"+variable, "h_leptonSFDown_"+variable, nBins, customBins);
  h_leptonSFDown->Sumw2();
  histograms_leptonSFDown.push_back(h_leptonSFDown);

  TH1D * h_photonSFUp = new TH1D("h_photonSFUp_"+variable, "h_photonSFUp_"+variable, nBins, customBins);
  h_photonSFUp->Sumw2();
  histograms_photonSFUp.push_back(h_photonSFUp);

  TH1D * h_photonSFDown = new TH1D("h_photonSFDown_"+variable, "h_photonSFDown_"+variable, nBins, customBins);
  h_photonSFDown->Sumw2();
  histograms_photonSFDown.push_back(h_photonSFDown);

  TH1D * h_JECUp = new TH1D("h_JECUp_"+variable, "h_JECUp_"+variable, nBins, customBins);
  h_JECUp->Sumw2();
  histograms_JECUp.push_back(h_JECUp);

  TH1D * h_JECDown = new TH1D("h_JECDown_"+variable, "h_JECDown_"+variable, nBins, customBins);
  h_JECDown->Sumw2();
  histograms_JECDown.push_back(h_JECDown);

  TH1D * h_mcQCD = new TH1D("h_mcQCD_"+variable, "h_mcQCD_"+variable, nBins, customBins);
  h_mcQCD->Sumw2();
  histograms_mcQCD.push_back(h_mcQCD);

}

void HistogramMaker::BookHistogram2D(TString var_x, TString var_y, Int_t nBins_x, Float_t xlo, Float_t xhi, Int_t nBins_y, Float_t ylo, Float_t yhi, Float_t zlo, Float_t zhi) {
  
  variables_2d.push_back(make_pair(var_x, var_y));

  TString hName = "h_"+var_y+"_vs_"+var_x;

  TH2D * h = new TH2D(hName, hName, nBins_x, xlo, xhi, nBins_y, ylo, yhi);
  h->Sumw2();
  if(zhi > zlo) h->GetZaxis()->SetRangeUser(zlo, zhi);
  histograms_2d.push_back(h);

}

void HistogramMaker::GetLeptonSF(Float_t lepton_pt, Float_t lepton_eta, TString chan, Float_t& central, Float_t& up, Float_t& down) {

  Float_t pt, eta, error;

  if(chan.Contains("ele")) {
    pt = min(lepton_pt, (float)199.);
    pt = max(pt, (float)15.);
    eta = min(fabs(lepton_eta), (double)2.39);

    Float_t id_val = sf_electron->GetBinContent(sf_electron->FindBin(eta, pt));
    Float_t id_error = sf_electron->GetBinError(sf_electron->FindBin(eta, pt));

    Float_t trigger_val = sf_SingleElectronTrigger->GetBinContent(sf_SingleElectronTrigger->FindBin(eta, pt));
    Float_t trigger_error = sf_SingleElectronTrigger->GetBinError(sf_SingleElectronTrigger->FindBin(eta, pt));

    central = id_val * trigger_val;
    error = central * sqrt(id_error*id_error/(id_val*id_val) + trigger_error*trigger_error/(trigger_val*trigger_val));

    up = central + error;
    down = central - error;
  }

  else {
    pt = min(lepton_pt, (float)499.);
    pt = max(pt, (float)10.);
    eta = min(fabs(lepton_eta), (double)2.09);

    central = sf_muon->GetBinContent(sf_muon->FindBin(pt, eta));
    error = sf_muon->GetBinError(sf_muon->FindBin(pt, eta));

    up = error;
    down = 2. * central - error;
  }

  return;  
}

void HistogramMaker::GetPhotonSF(Float_t lead_photon_et, Float_t lead_photon_eta, Float_t trail_photon_et, Float_t trail_photon_eta, Float_t nphotons, 
			    Float_t& central, Float_t& up, Float_t& down) {

  if(nphotons == 0) {
    central = 1.;
    up = 1.;
    down = 1.;
    return;
  }

  Float_t et, eta, error;

  if(nphotons == 1) {
    et = min(lead_photon_et, (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(lead_photon_eta), (double)1.44441);

    Float_t id_val = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    central = id_val * veto_val;
    error = central * sqrt(id_error*id_error/(id_val*id_val) + veto_error*veto_error/(veto_val*veto_val));
  }

  else if(nphotons >= 2) {
    // lead photon
    et = min(lead_photon_et, (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(lead_photon_eta), (double)1.44441);

    Float_t id_val_lead = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error_lead = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val_lead = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error_lead = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    // trail photon
    et = min(trail_photon_et, (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(trail_photon_eta), (double)1.44441);

    Float_t id_val_trail = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error_trail = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val_trail = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error_trail = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    central = id_val_lead * veto_val_lead * id_val_trail * veto_val_trail;
    error = central * sqrt(id_error_lead*id_error_lead/(id_val_lead*id_val_lead) +
			   veto_error_lead*veto_error_lead/(veto_val_lead*veto_val_lead) +
			   id_error_trail*id_error_trail/(id_val_trail*id_val_trail) +
			   veto_error_trail*veto_error_trail/(veto_val_trail*veto_val_trail));
  }

  up = central + error;
  down = central - error;

  return;
}
