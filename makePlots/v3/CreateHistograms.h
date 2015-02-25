#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>
#include <vector>
#include <stdio.h>
#include <stdarg.h>
#include <exception>

#include "rootRoutines.h"
#include "Binning.h"

using namespace std;

const int nChannels = 4;

TString channels[nChannels] = {"ele_jjj", "muon_jjj",
                               "ele_bjj", "muon_bjj"};

unsigned int nBtagReq[nChannels] = {0, 0,
                                    1, 1};

TString qcdChannels_noSigmaIetaIeta[nChannels] = {"ele_jjj_eQCDnoSigmaIetaIetaTree", "muon_jjj_muQCDnoSigmaIetaIetaTree",
                                                  "ele_jjj_veto_eQCDnoSigmaIetaIetaTree", "muon_jjj_veto_muQCDnoSigmaIetaIetaTree"};

enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kNumControlRegions};

TString crNames[kNumControlRegions] = {"SR1", "SR2", "CR1", "CR2", "CR2a", "CR0"};

class HistogramMaker : public TObject {
  
  ClassDef(HistogramMaker, 1);

 public:
  HistogramMaker(int chanNo, bool blind, int cRegion, Float_t cutOnMet);
  ~HistogramMaker();
  
  Float_t getValue(TString name) {
    if(varMap.find(name) != varMap.end()) return varMap.find(name)->second;
    return -1.e6;
  };

  Float_t getValue(unsigned int i) {
    if(i > variables.size()) return -1.e-6;
    return getValue(variables[i]);
  };

  Int_t getIntegerValue(TString name) { return (Int_t)getValue(name); };

  bool inControlRegion() {
    switch(controlRegion) {
    case kSR1:
      return getIntegerValue("Ngamma") == 1;
    case kSR2:
      return getIntegerValue("Ngamma") >= 2;
    case kCR1:
      return (getIntegerValue("Ngamma") == 0 && getIntegerValue("Nfake") == 1);
    case kCR2:
      return (getIntegerValue("Ngamma") == 0 && getIntegerValue("Nfake") >= 2);
    case kCR2a:
      return (getIntegerValue("Ngamma") == 1 && getIntegerValue("Nfake") == 1) || (getIntegerValue("Ngamma") == 0 && getIntegerValue("Nfake") >= 2);
    case kCR0:
      return getIntegerValue("Ngamma") == 0;
    default:
      return false;
    }
  };

  bool isBlindedRegion() { return (blinded && (controlRegion == kSR2 || controlRegion == kSR1)); };

  bool passMetCut() {
    if(metCut > 0.) return (getValue("pfMET") >= metCut);
    return true;
  };
    
  bool checkBtagging() {
    if(nBtagReq[channel] == 0) {
      btagWeight = 1.;
      btagWeightErr = 0.;
      btagWeightUp = 1.;
      btagWeightDown = 1.;
    }

    if(btagWeightUp < 0) btagWeightUp = 0.;
    if(btagWeightDown < 0) btagWeightDown = 0;

    if(btagWeightErr > 20. || btagWeightErr != btagWeightErr) btagWeightErr = btagWeight;
    
    return (btagWeight == btagWeight);
  };

  bool LoadMCBackground(TString fileName, TString scanName,
			Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfErrorUp, Double_t pdfErrorDown,
			bool removeTTA, bool reweightTop,
			Double_t fitScaling = -1.0, Double_t fitScalingError = 0.0);
  
  void SetTrees(TTree * gg, TTree * qcd, TTree * sig_a, TTree * sig_b);
  
  void BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi);
  void BookHistogram(TString variable, Int_t nBins, Double_t* customBins);

  void BookHistogram2D(TString var_x, TString var_y, Int_t nBins_x, Float_t xlo, Float_t xhi, Int_t nBins_y, Float_t ylo, Float_t yhi, Float_t zlo = 0.0, Float_t zhi = -1.0);

  void FillData();
  void FillQCD();
  void FillMCBackgrounds();
  void FillSignal();

  void FillHistograms() {
    cout << endl << "Filling Data..." << endl;
    FillData();
    cout << "Filling QCD..." << endl;
    FillQCD();
    cout << "Filling MC Backgrounds..." << endl;
    FillMCBackgrounds();
    cout << "Filling Signal..." << endl;
    FillSignal();
    cout << "Done filling!" << endl << endl;
  }

  void SubtractMCFromQCD();
  void NormalizeQCD();

  void SaveOutput();

  void GetLeptonSF(Float_t& central, Float_t& up, Float_t& down);
  void GetPhotonSF(Float_t& central, Float_t& up, Float_t& down);
  
  void GetLeptonSF(Float_t lepton_pt, Float_t lepton_eta, Float_t& central, Float_t& up, Float_t& down);
  void GetPhotonSF(Float_t lead_photon_et, Float_t lead_photon_eta, Float_t trail_photon_et, Float_t trail_photon_eta, Float_t nphotons, 
		   Float_t& central, Float_t& up, Float_t& down);

  void LoadLeptonSFs(TString fileName) {
    fLeptonSF = new TFile(fileName, "READ");
    sf_electron = (TH2D*)fLeptonSF->Get("TightEleIdIsoSF");
    sf_SingleElectronTrigger = (TH2D*)fLeptonSF->Get("TightEleTriggerSF");
    sf_muon = (TH2D*)fLeptonSF->Get("mu_pt_eta_full_id_iso_hlt_8TeV");
  };
  void LoadPhotonSFs(TString fileName) {
    fPhotonSF = new TFile(fileName, "READ");
    sf_photon_id = (TH2D*)fPhotonSF->Get("PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
    sf_photon_veto = (TH2D*)fPhotonSF->Get("PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
  };

  void CreateDatacards();
  
 private:
  
  // vectors of objects to keep open
  vector<TFile*> mcFiles;
  vector<TTree*> mcTrees;
  vector<TTree*> mcTrees_JECup;
  vector<TTree*> mcTrees_JECdown;
  vector<TTree*> mcQCDTrees;

  // internal histogram naming
  vector<TString> mcNames;

  // xsec scaling
  vector<Double_t> mcNGen;
  vector<Double_t> crossSections;

  // xsec systematics
  vector<Double_t> scaleErrUp;
  vector<Double_t> scaleErrDown;
  vector<Double_t> pdfErrUp;
  vector<Double_t> pdfErrDown;

  // reweighting flags
  vector<bool> removeTTAoverlap;
  vector<bool> reweightTopPt;

  // additional scalings
  vector<Double_t> fitScale;
  vector<Double_t> fitScaleError;

  // mcHistograms[process][variable]
  vector< vector<TH1D*> > mcHistograms;
  vector< vector<TH1D*> > mcHistograms_btagWeightUp;
  vector< vector<TH1D*> > mcHistograms_btagWeightDown;
  vector< vector<TH1D*> > mcHistograms_puWeightUp;
  vector< vector<TH1D*> > mcHistograms_puWeightDown;
  vector< vector<TH1D*> > mcHistograms_scaleUp;
  vector< vector<TH1D*> > mcHistograms_scaleDown;
  vector< vector<TH1D*> > mcHistograms_pdfUp;
  vector< vector<TH1D*> > mcHistograms_pdfDown;
  vector< vector<TH1D*> > mcHistograms_topPtUp;
  vector< vector<TH1D*> > mcHistograms_topPtDown;
  vector< vector<TH1D*> > mcHistograms_JECup;
  vector< vector<TH1D*> > mcHistograms_JECdown;
  vector< vector<TH1D*> > mcHistograms_leptonSFup;
  vector< vector<TH1D*> > mcHistograms_leptonSFdown;
  vector< vector<TH1D*> > mcHistograms_photonSFup;
  vector< vector<TH1D*> > mcHistograms_photonSFdown;
  vector< vector<TH2D*> > mcHistograms_2d;

  vector< vector<TH1D*> > mcQCDHistograms;
  vector< vector<TH2D*> > mcQCDHistograms_2d;

  vector< vector<TH1D*> > mcQCDHistograms_relIso_5;
  vector< vector<TH1D*> > mcQCDHistograms_relIso_10;

  // h_xyz[variable]
  vector<TH1D*> h_gg;
  vector<TH2D*> h_gg_2d;

  vector<TH1D*> h_qcd;
  vector<TH2D*> h_qcd_2d;

  vector<TH1D*> h_qcd_relIso_5;
  vector<TH1D*> h_qcd_relIso_10;

  vector<TH1D*> h_siga;
  vector<TH2D*> h_siga_2d;

  vector<TH1D*> h_sigb;
  vector<TH2D*> h_sigb_2d;

  TTree * ggTree;
  TTree * qcdTree;
  
  TTree * sigaTree;
  TTree * sigbTree;

  vector<TString> variables;
  vector<pair<TString, TString> > variables_2d;
  map<TString, Float_t> varMap;

  TFile * fLeptonSF;
  TH2D * sf_muon;
  TH2D * sf_electron;
  TH2D * sf_SingleElectronTrigger;

  TFile * fPhotonSF;
  TH2D * sf_photon_id;
  TH2D * sf_photon_veto;

  Int_t intLumi_int;

  int channel;
  bool blinded;

  int controlRegion;

  Float_t metCut;

  TString req;

  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr;
  Float_t puWeightUp, puWeightDown, btagWeightUp, btagWeightDown;
  Float_t overlaps_ttA;
  Float_t topPtReweighting;

  Float_t relIso;

};

HistogramMaker::HistogramMaker(int chanNo, bool blind, int cRegion, Float_t cutOnMet) :
channel(chanNo),
  blinded(blind),
  controlRegion(cRegion),
  metCut(cutOnMet)
{
  req = channels[chanNo];

  intLumi_int = 19712;

  variables.clear();
  variables_2d.clear();
  varMap.clear();

  h_gg.clear();
  h_gg_2d.clear();

  h_qcd.clear();
  h_qcd_2d.clear();
  
  h_qcd_relIso_5.clear();
  h_qcd_relIso_10.clear();

  h_siga.clear();
  h_siga_2d.clear();

  h_sigb.clear();
  h_sigb_2d.clear();

  mcFiles.clear();
  mcTrees.clear();
  mcTrees_JECup.clear();
  mcTrees_JECdown.clear();
  mcQCDTrees.clear();
  mcNGen.clear();
  crossSections.clear();
  scaleErrUp.clear();
  scaleErrDown.clear();
  pdfErrUp.clear();
  pdfErrDown.clear();

  mcNames.clear();

  removeTTAoverlap.clear();
  reweightTopPt.clear();

  fitScale.clear();
  fitScaleError.clear();

  mcHistograms.clear();
  mcHistograms_btagWeightUp.clear();
  mcHistograms_btagWeightDown.clear();
  mcHistograms_puWeightUp.clear();
  mcHistograms_puWeightDown.clear();
  mcHistograms_scaleUp.clear();
  mcHistograms_scaleDown.clear();
  mcHistograms_pdfUp.clear();
  mcHistograms_pdfDown.clear();
  mcHistograms_topPtUp.clear();
  mcHistograms_topPtDown.clear();
  mcHistograms_JECup.clear();
  mcHistograms_JECdown.clear();
  mcHistograms_leptonSFup.clear();
  mcHistograms_leptonSFdown.clear();
  mcHistograms_photonSFup.clear();
  mcHistograms_photonSFdown.clear();
  mcHistograms_2d.clear();

  mcQCDHistograms.clear();
  mcQCDHistograms_2d.clear();

  mcQCDHistograms_relIso_5.clear();
  mcQCDHistograms_relIso_10.clear();

}

HistogramMaker::~HistogramMaker() { 

  variables.clear();
  variables_2d.clear();
  varMap.clear();

  delete ggTree;
  delete qcdTree;

  mcHistograms.clear();
  mcHistograms_btagWeightUp.clear();
  mcHistograms_btagWeightDown.clear();
  mcHistograms_puWeightUp.clear();
  mcHistograms_puWeightDown.clear();
  mcHistograms_scaleUp.clear();
  mcHistograms_scaleDown.clear();
  mcHistograms_pdfUp.clear();
  mcHistograms_pdfDown.clear();
  mcHistograms_topPtUp.clear();
  mcHistograms_topPtDown.clear();
  mcHistograms_JECup.clear();
  mcHistograms_JECdown.clear();
  mcHistograms_leptonSFup.clear();
  mcHistograms_leptonSFdown.clear();
  mcHistograms_photonSFup.clear();
  mcHistograms_photonSFdown.clear();
  mcHistograms_2d.clear();

  mcQCDHistograms.clear();
  mcQCDHistograms_2d.clear();

  mcQCDHistograms_relIso_5.clear();
  mcQCDHistograms_relIso_10.clear();

  mcTrees.clear();
  mcTrees_JECup.clear();
  mcTrees_JECdown.clear();
  mcQCDTrees.clear();
  mcFiles.clear();
  mcNGen.clear();
  crossSections.clear();
  scaleErrUp.clear();
  scaleErrDown.clear();
  pdfErrUp.clear();
  pdfErrDown.clear();

  mcNames.clear();

  removeTTAoverlap.clear();
  reweightTopPt.clear();

  fitScale.clear();
  fitScaleError.clear();

  fLeptonSF->Close();
  fPhotonSF->Close();

  delete sigaTree;
  delete sigbTree;
    
  h_gg.clear();
  h_gg_2d.clear();

  h_qcd.clear();
  h_qcd_2d.clear();

  h_qcd_relIso_5.clear();
  h_qcd_relIso_10.clear();

  h_siga.clear();
  h_siga_2d.clear();
    
  h_sigb.clear();
  h_sigb_2d.clear();

}

void HistogramMaker::SetTrees(TTree * gg, TTree * qcd, TTree * sig_a, TTree * sig_b) {

  ggTree = gg;
  qcdTree = qcd;
  
  sigaTree = sig_a;
  sigbTree = sig_b;

}

bool HistogramMaker::LoadMCBackground(TString fileName, TString scanName,
				      Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfErrorUp, Double_t pdfErrorDown,
				      bool removeTTA, bool reweightTop,
				      Double_t fitScaling, Double_t fitScalingError) {

  mcFiles.push_back(new TFile(fileName, "READ"));
  if(!mcFiles.back()) {
    cout << "Could not load TFile " << fileName << endl;
    return false;
  }

  TString signalString = channels[channel]+"_noSigmaIetaIetaTree";
  TString qcdString = qcdChannels_noSigmaIetaIeta[channel];

  mcTrees.push_back((TTree*)mcFiles.back()->Get(signalString));
  if(!mcTrees.back()) {
    cout << "Could not load TTree " << signalString << " from TFile " << fileName << endl;
    return false;
  }

  mcTrees_JECup.push_back((TTree*)mcFiles.back()->Get(signalString+"_JECup"));
  if(!mcTrees_JECup.back()) {
    cout << "Could not load TTree " << signalString << "_JECup from TFile " << fileName << endl;
    return false;
  }
  
  mcTrees_JECdown.push_back((TTree*)mcFiles.back()->Get(signalString+"_JECdown"));
  if(!mcTrees_JECdown.back()) {
    cout << "Could not load TTree " << signalString << "_JECdown from TFile " << fileName << endl;
    return false;
  }

  mcQCDTrees.push_back((TTree*)mcFiles.back()->Get(qcdString));
  if(!mcQCDTrees.back()) {
    cout << "Could not load TTree " << qcdString << " from TFile " << fileName << endl;
    return false;
  }

  TH1D * h_nGen = (TH1D*)mcFiles.back()->Get("nEvents_"+scanName);
  if(!h_nGen) {
    cout << "Could not load histogram " << "nEvents_" << scanName << " from TFile " << fileName << endl;
    return false;
  }
  mcNGen.push_back(h_nGen->Integral());

  crossSections.push_back(xsec);
  scaleErrUp.push_back(scaleErrorUp);
  scaleErrDown.push_back(scaleErrorDown);
  pdfErrUp.push_back(pdfErrorUp);
  pdfErrDown.push_back(pdfErrorDown);  

  mcNames.push_back(scanName);

  removeTTAoverlap.push_back(removeTTA);
  reweightTopPt.push_back(reweightTop);

  fitScale.push_back(fitScaling);
  fitScaleError.push_back(fitScalingError);

  mcHistograms.resize(mcHistograms.size() + 1);
  mcHistograms_btagWeightUp.resize(mcHistograms_btagWeightUp.size() + 1);
  mcHistograms_btagWeightDown.resize(mcHistograms_btagWeightDown.size() + 1);
  mcHistograms_puWeightUp.resize(mcHistograms_puWeightUp.size() + 1);
  mcHistograms_puWeightDown.resize(mcHistograms_puWeightDown.size() + 1);
  mcHistograms_scaleUp.resize(mcHistograms_scaleUp.size() + 1);
  mcHistograms_scaleDown.resize(mcHistograms_scaleDown.size() + 1);
  mcHistograms_pdfUp.resize(mcHistograms_pdfUp.size() + 1);
  mcHistograms_pdfDown.resize(mcHistograms_pdfDown.size() + 1);
  mcHistograms_topPtUp.resize(mcHistograms_topPtUp.size() + 1);
  mcHistograms_topPtDown.resize(mcHistograms_topPtDown.size() + 1);
  mcHistograms_JECup.resize(mcHistograms_JECup.size() + 1);
  mcHistograms_JECdown.resize(mcHistograms_JECdown.size() + 1);
  mcHistograms_leptonSFup.resize(mcHistograms_leptonSFup.size() + 1);
  mcHistograms_leptonSFdown.resize(mcHistograms_leptonSFdown.size() + 1);
  mcHistograms_photonSFup.resize(mcHistograms_photonSFup.size() + 1);
  mcHistograms_photonSFdown.resize(mcHistograms_photonSFdown.size() + 1);
  mcHistograms_2d.resize(mcHistograms_2d.size() + 1);

  mcQCDHistograms.resize(mcQCDHistograms.size() + 1);
  mcQCDHistograms_2d.resize(mcQCDHistograms_2d.size() + 1);

  mcQCDHistograms_relIso_5.resize(mcQCDHistograms_relIso_5.size() + 1);
  mcQCDHistograms_relIso_10.resize(mcQCDHistograms_relIso_10.size() + 1);

  return true;
}

void HistogramMaker::BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi) {
  
  variables.push_back(variable);
  varMap[variable] = 0.;

  TH1D * gg = new TH1D(variable+"_gg_"+req, variable, nBins, xlo, xhi);
  gg->Sumw2();
  h_gg.push_back(gg);

  TH1D * qcd = new TH1D(variable+"_qcd_"+req, variable, nBins, xlo, xhi);
  qcd->Sumw2();
  h_qcd.push_back(qcd);

  h_qcd_relIso_5.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_5_"+req));
  h_qcd_relIso_10.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_10_"+req));

  TH1D * h_bkg;
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_"+mcNames[i]+"_"+req, variable, nBins, xlo, xhi);
    h_bkg->Sumw2();
    mcHistograms[i].push_back(h_bkg);

    TH1D * h_bkg_btagWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightUp");
    mcHistograms_btagWeightUp[i].push_back(h_bkg_btagWeightUp);

    TH1D * h_bkg_btagWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightDown");
    mcHistograms_btagWeightDown[i].push_back(h_bkg_btagWeightDown);

    TH1D * h_bkg_puWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_puWeightUp");
    mcHistograms_puWeightUp[i].push_back(h_bkg_puWeightUp);

    TH1D * h_bkg_puWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_puWeightDown");
    mcHistograms_puWeightDown[i].push_back(h_bkg_puWeightDown);

    TH1D * h_bkg_scaleUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleUp");
    mcHistograms_scaleUp[i].push_back(h_bkg_scaleUp);

    TH1D * h_bkg_scaleDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleDown");
    mcHistograms_scaleDown[i].push_back(h_bkg_scaleDown);

    TH1D * h_bkg_pdfUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfUp");
    mcHistograms_pdfUp[i].push_back(h_bkg_pdfUp);

    TH1D * h_bkg_pdfDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfDown");
    mcHistograms_pdfDown[i].push_back(h_bkg_pdfDown);

    TH1D * h_bkg_topPtUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtUp");
    mcHistograms_topPtUp[i].push_back(h_bkg_topPtUp);

    TH1D * h_bkg_topPtDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtDown");
    mcHistograms_topPtDown[i].push_back(h_bkg_topPtDown);

    TH1D * h_bkg_JECup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECUp");
    mcHistograms_JECup[i].push_back(h_bkg_JECup);

    TH1D * h_bkg_JECdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECDown");
    mcHistograms_JECdown[i].push_back(h_bkg_JECdown);

    TH1D * h_bkg_leptonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFUp");
    mcHistograms_leptonSFup[i].push_back(h_bkg_leptonSFup);

    TH1D * h_bkg_leptonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFDown");
    mcHistograms_leptonSFdown[i].push_back(h_bkg_leptonSFdown);

    TH1D * h_bkg_photonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFUp");
    mcHistograms_photonSFup[i].push_back(h_bkg_photonSFup);

    TH1D * h_bkg_photonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFDown");
    mcHistograms_photonSFdown[i].push_back(h_bkg_photonSFdown);
  }

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_qcd_"+mcNames[i]+"_"+req, variable, nBins, xlo, xhi);
    h_bkg->Sumw2();
    mcQCDHistograms[i].push_back(h_bkg);

    mcQCDHistograms_relIso_5[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_5_"+mcNames[i]+"_"+req));
    mcQCDHistograms_relIso_10[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_10_"+mcNames[i]+"_"+req));
  }

  TH1D * sig_a = new TH1D(variable+"_a_"+req, variable, nBins, xlo, xhi);
  sig_a->Sumw2();
  h_siga.push_back(sig_a);

  TH1D * sig_b = new TH1D(variable+"_b_"+req, variable, nBins, xlo, xhi);
  sig_b->Sumw2();
  h_sigb.push_back(sig_b);

}

void HistogramMaker::BookHistogram(TString variable, Int_t nBins, Double_t* customBins) {

  variables.push_back(variable);
  varMap[variable] = 0.;

  TH1D * gg = new TH1D(variable+"_gg_"+req, variable, nBins, customBins);
  gg->Sumw2();
  h_gg.push_back(gg);
  
  TH1D * qcd = new TH1D(variable+"_qcd_"+req, variable, nBins, customBins);
  qcd->Sumw2();
  h_qcd.push_back(qcd);

  h_qcd_relIso_5.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_5_"+req));
  h_qcd_relIso_10.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_10_"+req));

  TH1D * h_bkg;
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_"+mcNames[i]+"_"+req, variable, nBins, customBins);
    h_bkg->Sumw2();
    mcHistograms[i].push_back(h_bkg);

    TH1D * h_bkg_btagWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightUp");
    mcHistograms_btagWeightUp[i].push_back(h_bkg_btagWeightUp);

    TH1D * h_bkg_btagWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_btagWeightDown");
    mcHistograms_btagWeightDown[i].push_back(h_bkg_btagWeightDown);

    TH1D * h_bkg_puWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_puWeightUp");
    mcHistograms_puWeightUp[i].push_back(h_bkg_puWeightUp);

    TH1D * h_bkg_puWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_puWeightDown");
    mcHistograms_puWeightDown[i].push_back(h_bkg_puWeightDown);

    TH1D * h_bkg_scaleUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleUp");
    mcHistograms_scaleUp[i].push_back(h_bkg_scaleUp);

    TH1D * h_bkg_scaleDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_scaleDown");
    mcHistograms_scaleDown[i].push_back(h_bkg_scaleDown);

    TH1D * h_bkg_pdfUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfUp");
    mcHistograms_pdfUp[i].push_back(h_bkg_pdfUp);

    TH1D * h_bkg_pdfDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_pdfDown");
    mcHistograms_pdfDown[i].push_back(h_bkg_pdfDown);
    
    TH1D * h_bkg_topPtUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtUp");
    mcHistograms_topPtUp[i].push_back(h_bkg_topPtUp);

    TH1D * h_bkg_topPtDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_topPtDown");
    mcHistograms_topPtDown[i].push_back(h_bkg_topPtDown);

    TH1D * h_bkg_JECup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECUp");
    mcHistograms_JECup[i].push_back(h_bkg_JECup);

    TH1D * h_bkg_JECdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_JECDown");
    mcHistograms_JECdown[i].push_back(h_bkg_JECdown);

    TH1D * h_bkg_leptonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFUp");
    mcHistograms_leptonSFup[i].push_back(h_bkg_leptonSFup);

    TH1D * h_bkg_leptonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_leptonSFDown");
    mcHistograms_leptonSFdown[i].push_back(h_bkg_leptonSFdown);

    TH1D * h_bkg_photonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFUp");
    mcHistograms_photonSFup[i].push_back(h_bkg_photonSFup);

    TH1D * h_bkg_photonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+req+"_photonSFDown");
    mcHistograms_photonSFdown[i].push_back(h_bkg_photonSFdown);
  }
  
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH1D(variable+"_qcd_"+mcNames[i]+"_"+req, variable, nBins, customBins);
    h_bkg->Sumw2();
    mcQCDHistograms[i].push_back(h_bkg);

    mcQCDHistograms_relIso_5[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_5_"+mcNames[i]+"_"+req));
    mcQCDHistograms_relIso_10[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_10_"+mcNames[i]+"_"+req));
  }

  TH1D * sig_a = new TH1D(variable+"_a_"+req, variable, nBins, customBins);
  sig_a->Sumw2();
  h_siga.push_back(sig_a);

  TH1D * sig_b = new TH1D(variable+"_b_"+req, variable, nBins, customBins);
  sig_b->Sumw2();
  h_sigb.push_back(sig_b);
  
}

void HistogramMaker::BookHistogram2D(TString var_x, TString var_y, Int_t nBins_x, Float_t xlo, Float_t xhi, Int_t nBins_y, Float_t ylo, Float_t yhi, Float_t zlo, Float_t zhi) {
  
  variables_2d.push_back(make_pair(var_x, var_y));

  TH2D * gg = new TH2D(var_y+"_vs_"+var_x+"_gg_"+req, var_y+"_vs_"+var_x, nBins_x, xlo, xhi, nBins_y, ylo, yhi);
  gg->Sumw2();
  if(zhi > zlo) gg->GetZaxis()->SetRangeUser(zlo, zhi);
  h_gg_2d.push_back(gg);

  TH2D * qcd = new TH2D(var_y+"_vs_"+var_x+"_qcd_"+req, var_y+"_vs_"+var_x, nBins_x, xlo, xhi, nBins_y, ylo, yhi);
  qcd->Sumw2();
  if(zhi > zlo) qcd->GetZaxis()->SetRangeUser(zlo, zhi);
  h_qcd_2d.push_back(qcd);
  
  TH2D * h_bkg;
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH2D(var_y+"_vs_"+var_x+"_"+mcNames[i]+"_"+req, var_y+"_vs_"+var_x, nBins_x, xlo, xhi, nBins_y, ylo, yhi);
    h_bkg->Sumw2();
    if(zhi > zlo) h_bkg->GetZaxis()->SetRangeUser(zlo, zhi);
    mcHistograms_2d[i].push_back(h_bkg);
  }

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    h_bkg = new TH2D(var_y+"_vs_"+var_x+"_qcd_"+mcNames[i]+"_"+req, var_y+"_vs_"+var_x, nBins_x, xlo, xhi, nBins_y, ylo, yhi);
    h_bkg->Sumw2();
    if(zhi > zlo) h_bkg->GetZaxis()->SetRangeUser(zlo, zhi);
    mcQCDHistograms_2d[i].push_back(h_bkg);
  }

  TH2D * sig_a = new TH2D(var_y+"_vs_"+var_x+"_a_"+req, var_y+"_vs_"+var_x, nBins_x, xlo, xhi, nBins_y, ylo, yhi);
  sig_a->Sumw2();
  if(zhi > zlo) sig_a->GetZaxis()->SetRangeUser(zlo, zhi);
  h_siga_2d.push_back(sig_a);

  TH2D * sig_b = new TH2D(var_y+"_vs_"+var_x+"_b_"+req, var_y+"_vs_"+var_x, nBins_x, xlo, xhi, nBins_y, ylo, yhi);
  sig_b->Sumw2();
  if(zhi > zlo) sig_b->GetZaxis()->SetRangeUser(zlo, zhi);
  h_sigb_2d.push_back(sig_b);

}

void HistogramMaker::FillData() {
  
  for(unsigned int i = 0; i < variables.size(); i++) ggTree->SetBranchAddress(variables[i], &(varMap[variables[i]]));

  for(int i = 0; i < ggTree->GetEntries(); i++) {
    ggTree->GetEntry(i);

    if(!passMetCut()) continue;
    if(!inControlRegion()) continue;
    if(isBlindedRegion()) continue;

    for(unsigned int j = 0; j < variables.size(); j++) {

      for(unsigned int k = 0; k < variables_2d.size(); k++) {
	if(variables[j] == variables_2d[k].first) {
	  for(unsigned int m = 0; m < variables.size(); m++) {
	    if(variables[m] == variables_2d[k].second) {
	      h_gg_2d[k]->Fill(getValue(j), getValue(m));
	    }
	  }
	}
      }

      h_gg[j]->Fill(getValue(j));
    }

  }

  ggTree->ResetBranchAddresses();

}
  
void HistogramMaker::FillQCD() {

  for(unsigned int i = 0; i < variables.size(); i++) qcdTree->SetBranchAddress(variables[i], &(varMap[variables[i]]));

  if(req.Contains("ele")) qcdTree->SetBranchAddress("ele_relIso", &relIso);
  else qcdTree->SetBranchAddress("muon_relIso", &relIso);

  for(int i = 0; i < qcdTree->GetEntries(); i++) {
    qcdTree->GetEntry(i);
    
    if(!passMetCut()) continue;
    if(!inControlRegion()) continue;

    for(unsigned int j = 0; j < variables.size(); j++) {

      for(unsigned int k = 0; k < variables_2d.size(); k++) {
	if(variables[j] == variables_2d[k].first) {
	  for(unsigned int m = 0; m < variables.size(); m++) {
	    if(variables[m] == variables_2d[k].second) {
	      h_qcd_2d[k]->Fill(getValue(j), getValue(m));
	    }
	  }
	}
      }

      h_qcd[j]->Fill(getValue(j));

      if(relIso > 0.2 * 1.05) h_qcd_relIso_5[j]->Fill(getValue(j));
      if(relIso > 0.2 * 1.10) h_qcd_relIso_10[j]->Fill(getValue(j));

    }

  }

  qcdTree->ResetBranchAddresses();

}

void HistogramMaker::FillMCBackgrounds() {

  for(unsigned int i = 0; i < variables.size(); i++) {
    
    for(unsigned int j = 0; j < mcTrees.size(); j++) {
      mcTrees[j]->SetBranchAddress(variables[i], &(varMap[variables[i]]));
      mcTrees_JECup[j]->SetBranchAddress(variables[i], &(varMap[variables[i]]));
      mcTrees_JECdown[j]->SetBranchAddress(variables[i], &(varMap[variables[i]]));
      mcQCDTrees[j]->SetBranchAddress(variables[i], &(varMap[variables[i]]));
    }
    
  }
  
  for(unsigned int i = 0; i < mcTrees.size(); i++) {
    mcTrees[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees[i]->SetBranchAddress("pileupWeightErr", &puWeightErr);
    mcTrees[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcTrees[i]->SetBranchAddress("btagWeightErr", &btagWeightErr);
    mcTrees[i]->SetBranchAddress("btagWeightUp", &btagWeightUp);
    mcTrees[i]->SetBranchAddress("btagWeightDown", &btagWeightDown);
    mcTrees[i]->SetBranchAddress("pileupWeightUp", &puWeightUp);
    mcTrees[i]->SetBranchAddress("pileupWeightDown", &puWeightDown);

    mcTrees_JECup[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees_JECup[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcTrees_JECdown[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees_JECdown[i]->SetBranchAddress("btagWeight", &btagWeight);

    mcQCDTrees[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightErr", &puWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcQCDTrees[i]->SetBranchAddress("btagWeightErr", &btagWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeightUp", &btagWeightUp);
    mcQCDTrees[i]->SetBranchAddress("btagWeightDown", &btagWeightDown);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightUp", &puWeightUp);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightDown", &puWeightDown);

    if(req.Contains("ele")) mcQCDTrees[i]->SetBranchAddress("ele_relIso", &relIso);
    else mcQCDTrees[i]->SetBranchAddress("muon_relIso", &relIso);

    if(removeTTAoverlap[i]) {
      mcTrees[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
      mcTrees_JECup[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
      mcTrees_JECdown[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
      mcQCDTrees[i]->SetBranchAddress("overlaps_ttA", &overlaps_ttA);
    }

    if(reweightTopPt[i]) {
      mcTrees[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
      mcTrees_JECup[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
      mcTrees_JECdown[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
      mcQCDTrees[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    }

  }

  Float_t leptonSF, leptonSFup, leptonSFdown;
  Float_t photonSF, photonSFup, photonSFdown;

  for(unsigned int i = 0; i < mcTrees.size(); i++) {
    
    cout << " -- Starting on " << mcNames[i] << endl;

    for(int j = 0; j < mcTrees[i]->GetEntries(); j++) {
      mcTrees[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!passMetCut()) continue;
      if(!inControlRegion()) continue;
      
      if(removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      if(fitScale[i] > 0) addError2 = fitScale[i]*fitScale[i]*puWeight*puWeight*btagWeightErr*btagWeightErr + 
			    fitScale[i]*fitScale[i]*btagWeight*btagWeight*puWeightErr*puWeightErr +
			    puWeight*puWeight*btagWeight*btagWeight*fitScaleError[i]*fitScaleError[i];
      
      Float_t addError2_puOnly = btagWeight*btagWeight*puWeightErr*puWeightErr;
      if(fitScale[i] > 0) addError2_puOnly = fitScale[i]*fitScale[i]*btagWeight*btagWeight*puWeightErr*puWeightErr +
			    puWeight*puWeight*btagWeight*btagWeight*fitScaleError[i]*fitScaleError[i];
      
      for(unsigned int k = 0; k < variables.size(); k++) {
	
	double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	
	Float_t oldError = mcHistograms[i][k]->GetBinError(mcHistograms[i][k]->FindBin(getValue(k)));
	Float_t newerror = sqrt(oldError*oldError + addError2);
	mcHistograms[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms[i][k]->SetBinError(mcHistograms[i][k]->FindBin(getValue(k)), newerror);
					
	for(unsigned int m = 0; m < variables_2d.size(); m++) {
	  if(variables[k] == variables_2d[m].first) {
	    for(unsigned int n = 0; n < variables.size(); n++) {
	      if(variables[n] == variables_2d[m].second) {
		mcHistograms_2d[i][m]->Fill(getValue(k), getValue(n), totalWeight);
	      }
	    }
	  }
	}
	
	totalWeight = puWeight * btagWeightUp * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_btagWeightUp[i][k]->GetBinError(mcHistograms_btagWeightUp[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2_puOnly);
	mcHistograms_btagWeightUp[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_btagWeightUp[i][k]->SetBinError(mcHistograms_btagWeightUp[i][k]->FindBin(getValue(k)), newerror);
	
	totalWeight = puWeight * btagWeightDown * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_btagWeightDown[i][k]->GetBinError(mcHistograms_btagWeightDown[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2_puOnly);
	mcHistograms_btagWeightDown[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_btagWeightDown[i][k]->SetBinError(mcHistograms_btagWeightDown[i][k]->FindBin(getValue(k)), newerror);
	
	totalWeight = puWeightUp * btagWeight * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_puWeightUp[i][k]->GetBinError(mcHistograms_puWeightUp[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2_puOnly);
	mcHistograms_puWeightUp[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_puWeightUp[i][k]->SetBinError(mcHistograms_puWeightUp[i][k]->FindBin(getValue(k)), newerror);

	totalWeight = puWeightDown * btagWeight * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_puWeightDown[i][k]->GetBinError(mcHistograms_puWeightDown[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2_puOnly);
	mcHistograms_puWeightDown[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_puWeightDown[i][k]->SetBinError(mcHistograms_puWeightDown[i][k]->FindBin(getValue(k)), newerror);

	totalWeight = puWeight * btagWeight * leptonSFup * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_leptonSFup[i][k]->GetBinError(mcHistograms_leptonSFup[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_leptonSFup[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_leptonSFup[i][k]->SetBinError(mcHistograms[i][k]->FindBin(getValue(k)), newerror);

	totalWeight = puWeight * btagWeight * leptonSFdown * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_leptonSFdown[i][k]->GetBinError(mcHistograms_leptonSFdown[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_leptonSFdown[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_leptonSFdown[i][k]->SetBinError(mcHistograms[i][k]->FindBin(getValue(k)), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSFup;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_photonSFup[i][k]->GetBinError(mcHistograms_photonSFup[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_photonSFup[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_photonSFup[i][k]->SetBinError(mcHistograms[i][k]->FindBin(getValue(k)), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSFdown;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_photonSFdown[i][k]->GetBinError(mcHistograms_photonSFdown[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_photonSFdown[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_photonSFdown[i][k]->SetBinError(mcHistograms[i][k]->FindBin(getValue(k)), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSF;
	if(reweightTopPt[i]) totalWeight *= topPtReweighting * topPtReweighting;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_topPtUp[i][k]->GetBinError(mcHistograms_topPtUp[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_topPtUp[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_topPtUp[i][k]->SetBinError(mcHistograms[i][k]->FindBin(getValue(k)), newerror);

	totalWeight = puWeight * btagWeight * leptonSF * photonSF;
	if(fitScale[i] > 0) totalWeight *= fitScale[i];
	oldError = mcHistograms_topPtDown[i][k]->GetBinError(mcHistograms_topPtDown[i][k]->FindBin(getValue(k)));
	newerror = sqrt(oldError*oldError + addError2);
	mcHistograms_topPtDown[i][k]->Fill(getValue(k), totalWeight);
	mcHistograms_topPtDown[i][k]->SetBinError(mcHistograms[i][k]->FindBin(getValue(k)), newerror);
	
      }

    }

  }

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    for(unsigned int j = 0; j < mcHistograms[i].size(); j++) {
      for(int k = 0; k < mcHistograms[i][j]->GetNbinsX(); k++) {
	Double_t content = mcHistograms[i][j]->GetBinContent(k+1);
	Double_t error = mcHistograms[i][j]->GetBinError(k+1);
	
	mcHistograms_scaleUp[i][j]->SetBinContent(k+1, content);
	mcHistograms_scaleDown[i][j]->SetBinContent(k+1, content);
	mcHistograms_pdfUp[i][j]->SetBinContent(k+1, content);
	mcHistograms_pdfDown[i][j]->SetBinContent(k+1, content);
	
	mcHistograms_scaleUp[i][j]->SetBinError(k+1, error);
	mcHistograms_scaleDown[i][j]->SetBinError(k+1, error);
	mcHistograms_pdfUp[i][j]->SetBinError(k+1, error);
	mcHistograms_pdfDown[i][j]->SetBinError(k+1, error);
      }
    }
  }
  
  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    for(unsigned int j = 0; j < mcHistograms[i].size(); j++) {
      mcHistograms[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_btagWeightUp[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_btagWeightDown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_puWeightUp[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_puWeightDown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_scaleUp[i][j]->Scale(intLumi_int * (crossSections[i] + scaleErrUp[i]) / mcNGen[i]);
      mcHistograms_scaleDown[i][j]->Scale(intLumi_int * (crossSections[i] - scaleErrDown[i]) / mcNGen[i]);
      mcHistograms_pdfUp[i][j]->Scale(intLumi_int * (crossSections[i] + pdfErrUp[i]) / mcNGen[i]);
      mcHistograms_pdfDown[i][j]->Scale(intLumi_int * (crossSections[i] - pdfErrDown[i]) / mcNGen[i]);
      mcHistograms_topPtUp[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_topPtDown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_leptonSFup[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_leptonSFdown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_photonSFup[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      mcHistograms_photonSFdown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    }
  }

  for(unsigned int i = 0; i < mcHistograms_2d.size(); i++) {
    for(unsigned int j = 0; j < mcHistograms_2d[i].size(); j++) {
      mcHistograms_2d[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    }
  }
  
  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECup[i]->GetEntries(); j++) {
      mcTrees_JECup[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!passMetCut()) continue;
      if(!inControlRegion()) continue;
      
      if(removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      if(fitScale[i] > 0) totalWeight *= fitScale[i];
      
      for(unsigned int k = 0; k < variables.size(); k++) mcHistograms_JECup[i][k]->Fill(getValue(k), totalWeight);
      
    }
    
    for(unsigned int j = 0; j < variables.size(); j++) mcHistograms_JECup[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    
  }
  
  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECdown[i]->GetEntries(); j++) {
      mcTrees_JECdown[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!passMetCut()) continue;
      if(!inControlRegion()) continue;
      
      if(removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      if(fitScale[i] > 0) totalWeight *= fitScale[i];
      
      for(unsigned int k = 0; k < variables.size(); k++) mcHistograms_JECdown[i][k]->Fill(getValue(k), totalWeight);
      
    }
    
    for(unsigned int j = 0; j < variables.size(); j++) mcHistograms_JECdown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    
  }
  
  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) {

    for(int j = 0; j < mcQCDTrees[i]->GetEntries(); j++) {
      mcQCDTrees[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!passMetCut()) continue;
      if(!inControlRegion()) continue;

      if(removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;
      if(topPtReweighting < 0) topPtReweighting = 1.;

      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      if(fitScale[i] > 0) addError2 = fitScale[i]*fitScale[i]*puWeight*puWeight*btagWeightErr*btagWeightErr + 
			    fitScale[i]*fitScale[i]*btagWeight*btagWeight*puWeightErr*puWeightErr +
			    puWeight*puWeight*btagWeight*btagWeight*fitScaleError[i]*fitScaleError[i];
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      if(fitScale[i] > 0) totalWeight *= fitScale[i];

      for(unsigned int k = 0; k < variables.size(); k++) {
		
	Float_t oldError = mcQCDHistograms[i][k]->GetBinError(mcQCDHistograms[i][k]->FindBin(getValue(k)));
	Float_t newerror = sqrt(oldError*oldError + addError2);
	mcQCDHistograms[i][k]->Fill(getValue(k), totalWeight);
	mcQCDHistograms[i][k]->SetBinError(mcQCDHistograms[i][k]->FindBin(getValue(k)), newerror);
	
	if(relIso > 0.2 * 1.05) {
	  oldError = mcQCDHistograms_relIso_5[i][k]->GetBinError(mcQCDHistograms_relIso_5[i][k]->FindBin(getValue(k)));
	  newerror = sqrt(oldError*oldError + addError2);
	  mcQCDHistograms_relIso_5[i][k]->Fill(getValue(k), totalWeight);
	  mcQCDHistograms_relIso_5[i][k]->SetBinError(mcQCDHistograms_relIso_5[i][k]->FindBin(getValue(k)), newerror);
	}

	if(relIso > 0.2 * 1.10) {
	  oldError = mcQCDHistograms_relIso_10[i][k]->GetBinError(mcQCDHistograms_relIso_10[i][k]->FindBin(getValue(k)));
	  newerror = sqrt(oldError*oldError + addError2);
	  mcQCDHistograms_relIso_10[i][k]->Fill(getValue(k), totalWeight);
	  mcQCDHistograms_relIso_10[i][k]->SetBinError(mcQCDHistograms_relIso_10[i][k]->FindBin(getValue(k)), newerror);
	}
	
	for(unsigned int m = 0; m < variables_2d.size(); m++) {
	  if(variables[k] == variables_2d[m].first) {
	    for(unsigned int n = 0; n < variables.size(); n++) {
	      if(variables[n] == variables_2d[m].second) {
		mcQCDHistograms_2d[i][m]->Fill(getValue(k), getValue(n), totalWeight);
	      }
	    }
	  }
	}
	
      }
      
    }
    
    for(unsigned int j = 0; j < variables.size(); j++) mcQCDHistograms[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    for(unsigned int j = 0; j < variables_2d.size(); j++) mcQCDHistograms_2d[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);

    for(unsigned int j = 0; j < variables.size(); j++) mcQCDHistograms_relIso_5[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    for(unsigned int j = 0; j < variables.size(); j++) mcQCDHistograms_relIso_10[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    
  }

  for(unsigned int i = 0; i < mcTrees.size(); i++) mcTrees[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) mcTrees_JECup[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) mcTrees_JECdown[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) mcQCDTrees[i]->ResetBranchAddresses();

}

void HistogramMaker::FillSignal() {

  for(unsigned int i = 0; i < variables.size(); i++) {
    sigaTree->SetBranchAddress(variables[i], &(varMap[variables[i]]));
    sigbTree->SetBranchAddress(variables[i], &(varMap[variables[i]]));
  }

  sigaTree->SetBranchAddress("pileupWeight", &puWeight);
  sigaTree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigaTree->SetBranchAddress("btagWeight", &btagWeight);
  sigaTree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigaTree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigaTree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigaTree->SetBranchAddress("pileupWeightUp", &puWeightUp);
  sigaTree->SetBranchAddress("pileupWeightDown", &puWeightDown);
  sigaTree->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  sigbTree->SetBranchAddress("pileupWeight", &puWeight);
  sigbTree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  sigbTree->SetBranchAddress("btagWeight", &btagWeight);
  sigbTree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  sigbTree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  sigbTree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  sigbTree->SetBranchAddress("pileupWeightUp", &puWeightUp);
  sigbTree->SetBranchAddress("pileupWeightDown", &puWeightDown);
  sigbTree->SetBranchAddress("TopPtReweighting", &topPtReweighting);

  Float_t leptonSF, leptonSFup, leptonSFdown;
  Float_t photonSF, photonSFup, photonSFdown;

  for(int i = 0; i < sigaTree->GetEntries(); i++) {
    sigaTree->GetEntry(i);

    if(!checkBtagging()) continue;
    if(!passMetCut()) continue;
    if(!inControlRegion()) continue;
    
    if(removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;
    if(topPtReweighting < 0) topPtReweighting = 1.;
    
    GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(photonSF, photonSFup, photonSFdown);

    Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;

    if(topPtReweighting < 0) topPtReweighting = 1.;

    for(unsigned int j = 0; j < variables.size(); j++) {

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h_siga[j]->GetBinError(h_siga[j]->FindBin(getValue(j)));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_siga[j]->Fill(getValue(j), totalWeight);
      h_siga[j]->SetBinError(h_siga[j]->FindBin(getValue(j)), newerror);

      for(unsigned int k = 0; k < variables_2d.size(); k++) {
	if(variables[j] == variables_2d[k].first) {
	  for(unsigned int m = 0; m < variables.size(); m++) {
	    if(variables[m] == variables_2d[k].second) {
	      h_siga_2d[k]->Fill(getValue(j), getValue(m), totalWeight);
	    }
	  }
	}
      }

    }

  }
  for(unsigned int j = 0; j < varMap.size(); j++) h_siga[j]->Scale(intLumi_int * 0.147492 / 15000.);
  for(unsigned int j = 0; j < variables_2d.size(); j++) h_siga_2d[j]->Scale(intLumi_int * 0.147492 / 15000.);

  for(int i = 0; i < sigbTree->GetEntries(); i++) {
    sigbTree->GetEntry(i);

    if(!checkBtagging()) continue;
    if(!passMetCut()) continue;
    if(!inControlRegion()) continue;
    
    if(removeTTAoverlap[i] && overlaps_ttA > 0.001) continue;
    if(topPtReweighting < 0) topPtReweighting = 1.;
    
    GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
    GetPhotonSF(photonSF, photonSFup, photonSFdown);

    Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;

    if(topPtReweighting < 0) topPtReweighting = 1.;

    for(unsigned int j = 0; j < variables.size(); j++) {

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h_sigb[j]->GetBinError(h_sigb[j]->FindBin(getValue(j)));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_sigb[j]->Fill(getValue(j), totalWeight);
      h_sigb[j]->SetBinError(h_sigb[j]->FindBin(getValue(j)), newerror);

      for(unsigned int k = 0; k < variables_2d.size(); k++) {
	if(variables[j] == variables_2d[k].first) {
	  for(unsigned int m = 0; m < variables.size(); m++) {
	    if(variables[m] == variables_2d[k].second) {
	      h_sigb_2d[k]->Fill(getValue(j), getValue(m), totalWeight);
	    }
	  }
	}
      }

    }

  }

  for(unsigned int j = 0; j < variables.size(); j++) h_sigb[j]->Scale(intLumi_int * 0.0399591 / 15000.);
  for(unsigned int j = 0; j < variables_2d.size(); j++) h_sigb_2d[j]->Scale(intLumi_int * 0.0399591 / 15000.);

  sigaTree->ResetBranchAddresses();
  sigbTree->ResetBranchAddresses();

}

void HistogramMaker::SubtractMCFromQCD() {

  for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms[i].size(); j++) {
      h_qcd[j]->Add(mcQCDHistograms[i][j], -1.);
    }
  }

  for(unsigned int i = 0; i < mcQCDHistograms_relIso_5.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms_relIso_5[i].size(); j++) {
      h_qcd_relIso_5[j]->Add(mcQCDHistograms_relIso_5[i][j], -1.);
    }
  }

  for(unsigned int i = 0; i < mcQCDHistograms_relIso_10.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms_relIso_10[i].size(); j++) {
      h_qcd_relIso_10[j]->Add(mcQCDHistograms_relIso_10[i][j], -1.);
    }
  }

  for(unsigned int i = 0; i < mcQCDHistograms_2d.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms_2d[i].size(); j++) {
      h_qcd_2d[j]->Add(mcQCDHistograms_2d[i][j], -1.);
    }
  }

  for(unsigned int i = 0; i < h_qcd.size(); i++) {
    for(Int_t j = 0; j < h_qcd[i]->GetNbinsX(); j++) {
      if(h_qcd[i]->GetBinContent(j+1) < 0) {
	h_qcd[i]->SetBinContent(j+1, 0.);
	h_qcd[i]->SetBinError(j+1, 0.);
      }
    }
  }

  for(unsigned int i = 0; i < h_qcd_relIso_5.size(); i++) {
    for(Int_t j = 0; j < h_qcd_relIso_5[i]->GetNbinsX(); j++) {
      if(h_qcd_relIso_5[i]->GetBinContent(j+1) < 0) {
	h_qcd_relIso_5[i]->SetBinContent(j+1, 0.);
	h_qcd_relIso_5[i]->SetBinError(j+1, 0.);
      }
    }
  }

  for(unsigned int i = 0; i < h_qcd_relIso_10.size(); i++) {
    for(Int_t j = 0; j < h_qcd_relIso_10[i]->GetNbinsX(); j++) {
      if(h_qcd_relIso_10[i]->GetBinContent(j+1) < 0) {
	h_qcd_relIso_10[i]->SetBinContent(j+1, 0.);
	h_qcd_relIso_10[i]->SetBinError(j+1, 0.);
      }
    }
  }

  for(unsigned int i = 0; i < h_qcd_2d.size(); i++) {
    for(Int_t j = 0; j < h_qcd_2d[i]->GetNbinsX(); j++) {
      for(Int_t k = 0; k < h_qcd_2d[i]->GetNbinsY(); k++) {
	if(h_qcd_2d[i]->GetBinContent(j+1, k+1) < 0) {
	  h_qcd_2d[i]->SetBinContent(j+1, k+1, 0.);
	  h_qcd_2d[i]->SetBinError(j+1, k+1, 0.);
	}
      }
    }
  }

}

void HistogramMaker::NormalizeQCD() {

  unsigned int met_index = 0;
  for(unsigned int i = 0; i < variables.size(); i++) {
    if(variables[i] == "pfMET") {
      met_index = i;
      break;
    }
  }

  double n_qcd_before = h_qcd[met_index]->Integral();

  const int endBin = h_gg[met_index]->GetXaxis()->FindBin(20) - 1;

  double n_sig = h_gg[met_index]->Integral(0, endBin);
  double n_qcd = h_qcd[met_index]->Integral(0, endBin);

  double relIso_5_n_qcd = h_qcd_relIso_5[met_index]->Integral(0, endBin);
  double relIso_10_n_qcd = h_qcd_relIso_10[met_index]->Integral(0, endBin);

  if(n_qcd < 1) return;

  double n_mc = 0;
  for(unsigned int i = 0; i < mcHistograms.size(); i++) n_mc += mcHistograms[i][met_index]->Integral(0, mcHistograms[i][met_index]->GetXaxis()->FindBin(20) - 1);

  double sigma_sig = 0;
  double sigma_qcd = 0;
  double sigma_mc = 0;

  double relIso_5_sigma_qcd = 0;
  double relIso_10_sigma_qcd = 0;
  
  for(int i = 0; i < endBin; i++) {
    sigma_sig += h_gg[met_index]->GetBinError(i+1) * h_gg[met_index]->GetBinError(i+1);
    sigma_qcd += h_qcd[met_index]->GetBinError(i+1) * h_qcd[met_index]->GetBinError(i+1);

    relIso_5_sigma_qcd += h_qcd_relIso_5[met_index]->GetBinError(i+1) * h_qcd_relIso_5[met_index]->GetBinError(i+1);
    relIso_10_sigma_qcd += h_qcd_relIso_10[met_index]->GetBinError(i+1) * h_qcd_relIso_10[met_index]->GetBinError(i+1);

    for(unsigned int j = 0; j < mcHistograms.size(); j++) sigma_mc += mcHistograms[j][met_index]->GetBinError(i+1) * mcHistograms[j][met_index]->GetBinError(i+1);
  }

  sigma_sig = sqrt(sigma_sig);
  sigma_qcd = sqrt(sigma_qcd);
  sigma_mc = sqrt(sigma_mc);

  relIso_5_sigma_qcd = sqrt(relIso_5_sigma_qcd);
  relIso_10_sigma_qcd = sqrt(relIso_10_sigma_qcd);

  double scale = (n_sig - n_mc) / n_qcd;
  
  double relIso_5_scale = (n_sig - n_mc) / relIso_5_n_qcd;
  double relIso_10_scale = (n_sig - n_mc) / relIso_10_n_qcd;

  if(scale < 0) scale = 1.e-6;

  for(unsigned int i = 0; i < h_qcd.size(); i++) {

    for(int j = 0; j < h_qcd[i]->GetNbinsX(); j++) {
      
      double newError = sigma_sig*sigma_sig + sigma_mc*sigma_mc + scale*scale*sigma_qcd*sigma_qcd;
      newError *= h_qcd[i]->GetBinContent(j+1)*h_qcd[i]->GetBinContent(j+1) / (n_qcd*n_qcd);
      newError += h_qcd[i]->GetBinError(j+1)*h_qcd[i]->GetBinError(j+1) * scale*scale;
      newError = sqrt(newError);

      h_qcd[i]->SetBinContent(j+1, h_qcd[i]->GetBinContent(j+1) * scale);
      h_qcd[i]->SetBinError(j+1, newError);

      newError = sigma_sig*sigma_sig + sigma_mc*sigma_mc + relIso_5_scale*relIso_5_scale*relIso_5_sigma_qcd*relIso_5_sigma_qcd;
      newError *= h_qcd_relIso_5[i]->GetBinContent(j+1)*h_qcd_relIso_5[i]->GetBinContent(j+1) / (relIso_5_n_qcd*relIso_5_n_qcd);
      newError += h_qcd_relIso_5[i]->GetBinError(j+1)*h_qcd_relIso_5[i]->GetBinError(j+1) * relIso_5_scale*relIso_5_scale;
      newError = sqrt(newError);

      h_qcd_relIso_5[i]->SetBinContent(j+1, h_qcd_relIso_5[i]->GetBinContent(j+1) * relIso_5_scale);
      h_qcd_relIso_5[i]->SetBinError(j+1, newError);

      newError = sigma_sig*sigma_sig + sigma_mc*sigma_mc + relIso_10_scale*relIso_10_scale*relIso_10_sigma_qcd*relIso_10_sigma_qcd;
      newError *= h_qcd_relIso_5[i]->GetBinContent(j+1)*h_qcd_relIso_5[i]->GetBinContent(j+1) / (relIso_10_n_qcd*relIso_10_n_qcd);
      newError += h_qcd_relIso_5[i]->GetBinError(j+1)*h_qcd_relIso_5[i]->GetBinError(j+1) * relIso_10_scale*relIso_10_scale;
      newError = sqrt(newError);

      h_qcd_relIso_5[i]->SetBinContent(j+1, h_qcd_relIso_5[i]->GetBinContent(j+1) * relIso_10_scale);
      h_qcd_relIso_5[i]->SetBinError(j+1, newError);

    }

  }

  for(unsigned int i = 0; i < h_qcd_2d.size(); i++) h_qcd_2d[i]->Scale(scale);

  double n_qcd_after = h_qcd[met_index]->Integral();

  cout << endl << "NormalizeQCD(): Overall scaling factor of " << n_qcd_after / n_qcd_before << " applied." << endl << endl;

}

void HistogramMaker::SaveOutput() {

  TString outName = "histograms_"+req+"_"+crNames[controlRegion]+".root";

  TFile * fOut = new TFile(outName, "UPDATE");

  for(unsigned int i = 0; i < h_gg.size(); i++) h_gg[i]->Write();
  for(unsigned int i = 0; i < h_gg_2d.size(); i++) h_gg_2d[i]->Write();
  for(unsigned int i = 0; i < h_qcd.size(); i++) h_qcd[i]->Write();
  for(unsigned int i = 0; i < h_qcd_relIso_5.size(); i++) h_qcd_relIso_5[i]->Write();
  for(unsigned int i = 0; i < h_qcd_relIso_10.size(); i++) h_qcd_relIso_10[i]->Write();
  for(unsigned int i = 0; i < h_qcd_2d.size(); i++) h_qcd_2d[i]->Write();

  for(unsigned int i = 0; i < h_siga.size(); i++) h_siga[i]->Write();
  for(unsigned int i = 0; i < h_siga_2d.size(); i++) h_siga_2d[i]->Write();
  for(unsigned int i = 0; i < h_sigb.size(); i++) h_sigb[i]->Write();
  for(unsigned int i = 0; i < h_sigb_2d.size(); i++) h_sigb_2d[i]->Write();

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    for(unsigned int j = 0; j < mcHistograms[i].size(); j++) {
      mcHistograms[i][j]->Write();
      mcHistograms_btagWeightUp[i][j]->Write();
      mcHistograms_btagWeightDown[i][j]->Write();
      mcHistograms_puWeightUp[i][j]->Write();
      mcHistograms_puWeightDown[i][j]->Write();
      mcHistograms_scaleUp[i][j]->Write();
      mcHistograms_scaleDown[i][j]->Write();
      mcHistograms_pdfUp[i][j]->Write();
      mcHistograms_pdfDown[i][j]->Write();
      mcHistograms_topPtUp[i][j]->Write();
      mcHistograms_topPtDown[i][j]->Write();
      mcHistograms_JECup[i][j]->Write();
      mcHistograms_JECdown[i][j]->Write();
      mcHistograms_leptonSFup[i][j]->Write();
      mcHistograms_leptonSFdown[i][j]->Write();
      mcHistograms_photonSFup[i][j]->Write();
      mcHistograms_photonSFdown[i][j]->Write();
    }
  }

  for(unsigned int i = 0; i < mcHistograms_2d.size(); i++) {
    for(unsigned int j = 0; j < mcHistograms_2d[i].size(); j++) {
      mcHistograms_2d[i][j]->Write();
    }
  }

  fOut->Close();

}

void HistogramMaker::GetLeptonSF(Float_t& central, Float_t& up, Float_t& down) {

  Float_t lepton_pt = (req.Contains("ele")) ? getValue("ele_pt") : getValue("muon_pt");
  Float_t lepton_eta = (req.Contains("ele")) ? getValue("ele_eta") : getValue("muon_eta");

  Float_t pt, eta, error;

  if(req.Contains("ele")) {
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

void HistogramMaker::GetPhotonSF(Float_t& central, Float_t& up, Float_t& down) {

  if(getIntegerValue("Nphotons") == 0) {
    central = 1.;
    up = 1.;
    down = 1.;
    return;
  }

  Float_t et, eta;
  Float_t error = 0;

  if(getIntegerValue("Nphotons") == 1) {
    et = min(getValue("leadPhotonEt"), (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(getValue("leadPhotonEta")), (double)1.44441);

    Float_t id_val = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    central = id_val * veto_val;
    error = central * sqrt(id_error*id_error/(id_val*id_val) + veto_error*veto_error/(veto_val*veto_val));
  }

  else if(getIntegerValue("Nphotons") >= 2) {
    // lead photon
    et = min(getValue("leadPhotonEt"), (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(getValue("leadPhotonEta")), (double)1.44441);

    Float_t id_val_lead = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error_lead = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val_lead = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error_lead = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    // trail photon
    et = min(getValue("trailPhotonEt"), (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(getValue("trailPhotonEta")), (double)1.44441);

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

void HistogramMaker::GetLeptonSF(Float_t lepton_pt, Float_t lepton_eta, Float_t& central, Float_t& up, Float_t& down) {

  Float_t pt, eta, error;

  if(req.Contains("ele")) {
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

  Float_t et, eta;
  Float_t error = 0;

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

void HistogramMaker::CreateDatacards() {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  Double_t xbins[31];
  xbins[0] = 0;
  xbins[1] = 55;
  for(int i = 1; i < 29; i++) xbins[i+1] = (mst[i] + mst[i-1])/2.;
  xbins[30] = 6510;

  Double_t ybins[33];
  ybins[0] = 0;
  ybins[1] = 12.5;
  for(int i = 1; i < 31; i++) ybins[i+1] = (mBino[i] + mBino[i-1])/2.;
  ybins[32] = 2175;

  char code[100];
  int index1, index2;

  TFile * f_xsec = new TFile("../../data/stop-bino_xsecs.root", "READ");
  TH2D * h_xsec = (TH2D*)f_xsec->Get("real_xsec");
  TH2D * h_xsec_errors = (TH2D*)f_xsec->Get("real_errors");

  TH2D * h_acc = new TH2D("acc_"+req, "acc_"+req, 30, xbins, 32, ybins);
  TH2D * h_contamination = new TH2D("contamination_"+req, "contamination_"+req, 30, xbins, 32, ybins);

  TString outName = "limitInputs_";
  if(req.Contains("ele")) outName += req(4, 3);
  if(req.Contains("muon")) outName += req(5, 3);
  outName += ".root";

  TFile * fSignalOut = new TFile(outName, "UPDATE");
  if(req.Contains("ele")) {
    if(!(fSignalOut->GetDirectory("ele_"+crNames[controlRegion]))) fSignalOut->mkdir("ele_"+crNames[controlRegion]);
    fSignalOut->cd("ele_"+crNames[controlRegion]);
  }
  else {
    if(!(fSignalOut->GetDirectory("muon_"+crNames[controlRegion]))) fSignalOut->mkdir("muon_"+crNames[controlRegion]);
    fSignalOut->cd("muon_"+crNames[controlRegion]);
  }

  for(int imass = 0; imass < 899; imass++) {

    index1 = mst[int(imass)/31];
    index2 = mBino[int(imass)%31];

    if(index1 < index2) continue;

    sprintf(code, "_mst_%d_m1_%d", index1, index2);
    TString code_t = code;

    TFile * f = new TFile("/eos/uscms/store/user/bfrancis/inputs_v4/acceptance/signal_contamination"+code_t+".root", "READ");
    if(f->IsZombie()) {
      f->Close();
      continue;
    }
    
    TTree * tree = (TTree*)f->Get(req+"_noSigmaIetaIetaTree");
    TTree * tree_JECup = (TTree*)f->Get(req+"_noSigmaIetaIetaTree_JECup");
    TTree * tree_JECdown = (TTree*)f->Get(req+"_noSigmaIetaIetaTree_JECdown");

    TTree * tree_contam;
    if(req.Contains("ele")) tree_contam = (TTree*)f->Get("ele_jjj_veto_eQCDnoSigmaIetaIetaTree");
    else if(req.Contains("muon")) tree_contam = (TTree*)f->Get("muon_jjj_veto_muQCDnoSigmaIetaIetaTree");

    if(!tree || !tree_JECup || !tree_JECdown || !tree_contam) {
      f->Close();
      continue;
    }

    Float_t met, ngamma, nfake, nphotons;
    Float_t lepton_pt, lepton_eta;
    Float_t lead_photon_et, lead_photon_eta;
    Float_t trail_photon_et, trail_photon_eta;

    tree->SetBranchAddress("pfMET", &met);
    tree_JECup->SetBranchAddress("pfMET", &met);
    tree_JECdown->SetBranchAddress("pfMET", &met);
    tree_contam->SetBranchAddress("pfMET", &met);

    tree->SetBranchAddress("Ngamma", &ngamma);
    tree_JECup->SetBranchAddress("Ngamma", &ngamma);
    tree_JECdown->SetBranchAddress("Ngamma", &ngamma);
    tree_contam->SetBranchAddress("Ngamma", &ngamma);

    tree->SetBranchAddress("Nfake", &nfake);
    tree_JECup->SetBranchAddress("Nfake", &nfake);
    tree_JECdown->SetBranchAddress("Nfake", &nfake);
    tree_contam->SetBranchAddress("Nfake", &nfake);

    tree->SetBranchAddress("Nphotons", &nphotons);
    tree_JECup->SetBranchAddress("Nphotons", &nphotons);
    tree_JECdown->SetBranchAddress("Nphotons", &nphotons);
    tree_contam->SetBranchAddress("Nphotons", &nphotons);

    tree->SetBranchAddress("leadPhotonEt", &lead_photon_et);
    tree_JECup->SetBranchAddress("leadPhotonEt", &lead_photon_et);
    tree_JECdown->SetBranchAddress("leadPhotonEt", &lead_photon_et);
    tree_contam->SetBranchAddress("leadPhotonEt", &lead_photon_et);

    tree->SetBranchAddress("leadPhotonEta", &lead_photon_eta);
    tree_JECup->SetBranchAddress("leadPhotonEta", &lead_photon_eta);
    tree_JECdown->SetBranchAddress("leadPhotonEta", &lead_photon_eta);
    tree_contam->SetBranchAddress("leadPhotonEta", &lead_photon_eta);

    tree->SetBranchAddress("trailPhotonEta", &trail_photon_eta);
    tree_JECup->SetBranchAddress("trailPhotonEta", &trail_photon_eta);
    tree_JECdown->SetBranchAddress("trailPhotonEta", &trail_photon_eta);
    tree_contam->SetBranchAddress("trailPhotonEta", &trail_photon_eta);

    tree->SetBranchAddress("trailPhotonEt", &trail_photon_et);
    tree_JECup->SetBranchAddress("trailPhotonEt", &trail_photon_et);
    tree_JECdown->SetBranchAddress("trailPhotonEt", &trail_photon_et);
    tree_contam->SetBranchAddress("trailPhotonEt", &trail_photon_et);

    if(req.Contains("ele")) {
      tree->SetBranchAddress("ele_pt", &lepton_pt);
      tree_JECup->SetBranchAddress("ele_pt", &lepton_pt);
      tree_JECdown->SetBranchAddress("ele_pt", &lepton_pt);
      tree_contam->SetBranchAddress("ele_pt", &lepton_pt);

      tree->SetBranchAddress("ele_eta", &lepton_eta);
      tree_JECup->SetBranchAddress("ele_eta", &lepton_eta);
      tree_JECdown->SetBranchAddress("ele_eta", &lepton_eta);
      tree_contam->SetBranchAddress("ele_eta", &lepton_eta);
    }
    else if(req.Contains("muon")) {
      tree->SetBranchAddress("muon_pt", &lepton_pt);
      tree_JECup->SetBranchAddress("muon_pt", &lepton_pt);
      tree_JECdown->SetBranchAddress("muon_pt", &lepton_pt);
      tree_contam->SetBranchAddress("muon_pt", &lepton_pt);

      tree->SetBranchAddress("muon_eta", &lepton_eta);
      tree_JECup->SetBranchAddress("muon_eta", &lepton_eta);
      tree_JECdown->SetBranchAddress("muon_eta", &lepton_eta);
      tree_contam->SetBranchAddress("muon_eta", &lepton_eta);
    }

    tree->SetBranchAddress("pileupWeight", &puWeight);
    tree->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree->SetBranchAddress("btagWeight", &btagWeight);
    tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    
    tree_JECup->SetBranchAddress("pileupWeight", &puWeight);
    tree_JECup->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree_JECup->SetBranchAddress("btagWeight", &btagWeight);
    tree_JECup->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree_JECup->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree_JECup->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree_JECup->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree_JECup->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree_JECup->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    
    tree_JECdown->SetBranchAddress("pileupWeight", &puWeight);
    tree_JECdown->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree_JECdown->SetBranchAddress("btagWeight", &btagWeight);
    tree_JECdown->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree_JECdown->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree_JECdown->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree_JECdown->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree_JECdown->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree_JECdown->SetBranchAddress("TopPtReweighting", &topPtReweighting);

    tree_contam->SetBranchAddress("pileupWeight", &puWeight);
    tree_contam->SetBranchAddress("pileupWeightErr", &puWeightErr);
    tree_contam->SetBranchAddress("btagWeight", &btagWeight);
    tree_contam->SetBranchAddress("btagWeightErr", &btagWeightErr);
    tree_contam->SetBranchAddress("btagWeightUp", &btagWeightUp);
    tree_contam->SetBranchAddress("btagWeightDown", &btagWeightDown);
    tree_contam->SetBranchAddress("pileupWeightUp", &puWeightUp);
    tree_contam->SetBranchAddress("pileupWeightDown", &puWeightDown);
    tree_contam->SetBranchAddress("TopPtReweighting", &topPtReweighting);

    TH1D * h = new TH1D("signal"+code_t, "signal"+code_t, nMetBins, xbins_met); h->Sumw2();

    TH1D * h_btagWeightUp = new TH1D("signal"+code_t+"_btagWeightUp", "signal"+code_t+"_btagWeightUp", nMetBins, xbins_met); h_btagWeightUp->Sumw2();
    TH1D * h_btagWeightDown = new TH1D("signal"+code_t+"_btagWeightDown", "signal"+code_t+"_btagWeightDown", nMetBins, xbins_met); h_btagWeightDown->Sumw2();

    TH1D * h_puWeightUp = new TH1D("signal"+code_t+"_puWeightUp", "signal"+code_t+"_puWeightUp", nMetBins, xbins_met); h_puWeightUp->Sumw2();
    TH1D * h_puWeightDown = new TH1D("signal"+code_t+"_puWeightDown", "signal"+code_t+"_puWeightDown", nMetBins, xbins_met); h_puWeightDown->Sumw2();

    TH1D * h_topPtUp = new TH1D("signal"+code_t+"_topPtUp", "signal"+code_t+"_topPtUp", nMetBins, xbins_met); h_topPtUp->Sumw2();
    TH1D * h_topPtDown = new TH1D("signal"+code_t+"_topPtDown", "signal"+code_t+"_topPtDown", nMetBins, xbins_met); h_topPtDown->Sumw2();

    TH1D * h_JECup = new TH1D("signal"+code_t+"_JECUp", "signal"+code_t+"_JECUp", nMetBins, xbins_met); h_JECup->Sumw2();
    TH1D * h_JECdown = new TH1D("signal"+code_t+"_JECDown", "signal"+code_t+"_JECDown", nMetBins, xbins_met); h_JECdown->Sumw2();

    TH1D * h_leptonSFup = new TH1D("signal"+code_t+"_leptonSFUp", "signal"+code_t+"_leptonSFUp", nMetBins, xbins_met); h_leptonSFup->Sumw2();
    TH1D * h_leptonSFdown = new TH1D("signal"+code_t+"_leptonSFDown", "signal"+code_t+"_leptonSFDown", nMetBins, xbins_met); h_leptonSFdown->Sumw2();

    TH1D * h_photonSFup = new TH1D("signal"+code_t+"_photonSFUp", "signal"+code_t+"_photonSFUp", nMetBins, xbins_met); h_photonSFup->Sumw2();
    TH1D * h_photonSFdown = new TH1D("signal"+code_t+"_photonSFDown", "signal"+code_t+"_photonSFDown", nMetBins, xbins_met); h_photonSFdown->Sumw2();

    for(int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);

      if(controlRegion == kSR1 && !(ngamma == 1)) continue;
      if(controlRegion == kSR2 && !(ngamma >= 2)) continue;
      if(controlRegion == kCR1 && !(ngamma == 0 && nfake == 1)) continue;
      if(controlRegion == kCR2 && !(ngamma == 0 && nfake >= 2)) continue;
      if(controlRegion == kCR2 && !((ngamma == 1 && nfake == 1) || (ngamma == 0 && nfake >= 2))) continue;
      if(controlRegion == kCR0 && !(ngamma == 0)) continue;

      if(!checkBtagging()) continue;

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;

      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;

      if(totalWeight < 0) continue;

      Float_t olderror = h->GetBinError(h->FindBin(met));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h->Fill(met, totalWeight);
      h->SetBinError(h->FindBin(met), newerror);

      totalWeight = puWeight * btagWeightUp * leptonSF * photonSF * topPtReweighting;
      olderror = h_btagWeightUp->GetBinError(h_btagWeightUp->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_btagWeightUp->Fill(met, totalWeight);
      h_btagWeightUp->SetBinError(h_btagWeightUp->FindBin(met), newerror);

      totalWeight = puWeight * btagWeightDown * leptonSF * photonSF * topPtReweighting;
      olderror = h_btagWeightDown->GetBinError(h_btagWeightDown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_btagWeightDown->Fill(met, totalWeight);
      h_btagWeightDown->SetBinError(h_btagWeightDown->FindBin(met), newerror);

      totalWeight = puWeightUp * btagWeight * leptonSF * photonSF * topPtReweighting;
      olderror = h_puWeightUp->GetBinError(h_puWeightUp->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_puWeightUp->Fill(met, totalWeight);
      h_puWeightUp->SetBinError(h_puWeightUp->FindBin(met), newerror);

      totalWeight = puWeightDown * btagWeight * leptonSF * photonSF * topPtReweighting;
      olderror = h_puWeightDown->GetBinError(h_puWeightDown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_puWeightDown->Fill(met, totalWeight);
      h_puWeightDown->SetBinError(h_puWeightDown->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSFup * photonSF * topPtReweighting;
      olderror = h_leptonSFup->GetBinError(h_leptonSFup->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_leptonSFup->Fill(met, totalWeight);
      h_leptonSFup->SetBinError(h_leptonSFup->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSFdown * photonSF * topPtReweighting;
      olderror = h_leptonSFdown->GetBinError(h_leptonSFdown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_leptonSFdown->Fill(met, totalWeight);
      h_leptonSFdown->SetBinError(h_leptonSFdown->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFup * topPtReweighting;
      olderror = h_photonSFup->GetBinError(h_photonSFup->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_photonSFup->Fill(met, totalWeight);
      h_photonSFup->SetBinError(h_photonSFup->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFdown * topPtReweighting;
      olderror = h_photonSFdown->GetBinError(h_photonSFdown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_photonSFdown->Fill(met, totalWeight);
      h_photonSFdown->SetBinError(h_photonSFdown->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting * topPtReweighting;
      olderror = h_topPtUp->GetBinError(h_topPtUp->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_topPtUp->Fill(met, totalWeight);
      h_topPtUp->SetBinError(h_topPtUp->FindBin(met), newerror);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      olderror = h_topPtDown->GetBinError(h_topPtDown->FindBin(met));
      newerror = sqrt(olderror*olderror + addError2);
      h_topPtDown->Fill(met, totalWeight);
      h_topPtDown->SetBinError(h_topPtDown->FindBin(met), newerror);
    }

    for(int i = 0; i < tree_JECup->GetEntries(); i++) {
      tree_JECup->GetEntry(i);

      if(controlRegion == kSR1 && !(ngamma == 1)) continue;
      if(controlRegion == kSR2 && !(ngamma >= 2)) continue;
      if(controlRegion == kCR1 && !(ngamma == 0 && nfake == 1)) continue;
      if(controlRegion == kCR2 && !(ngamma == 0 && nfake >= 2)) continue;
      if(controlRegion == kCR2 && !((ngamma == 1 && nfake == 1) || (ngamma == 0 && nfake >= 2))) continue;
      if(controlRegion == kCR0 && !(ngamma == 0)) continue;

      if(!checkBtagging()) continue;

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h->GetBinError(h->FindBin(met));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_JECup->Fill(met, totalWeight);
      h_JECup->SetBinError(h->FindBin(met), newerror);

    }

    for(int i = 0; i < tree_JECdown->GetEntries(); i++) {
      tree_JECdown->GetEntry(i);

      if(controlRegion == kSR1 && !(ngamma == 1)) continue;
      if(controlRegion == kSR2 && !(ngamma >= 2)) continue;
      if(controlRegion == kCR1 && !(ngamma == 0 && nfake == 1)) continue;
      if(controlRegion == kCR2 && !(ngamma == 0 && nfake >= 2)) continue;
      if(controlRegion == kCR2 && !((ngamma == 1 && nfake == 1) || (ngamma == 0 && nfake >= 2))) continue;
      if(controlRegion == kCR0 && !(ngamma == 0)) continue;

      if(!checkBtagging()) continue;

      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      
      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;
      Float_t olderror = h->GetBinError(h->FindBin(met));
      Float_t newerror = sqrt(olderror*olderror + addError2);
      h_JECdown->Fill(met, totalWeight);
      h_JECdown->SetBinError(h->FindBin(met), newerror);

    }

    double contamination = 0;

    for(int i = 0; i < tree_contam->GetEntries(); i++) {
      tree_contam->GetEntry(i);

      if(controlRegion == kSR1 && !(ngamma == 1)) continue;
      if(controlRegion == kSR2 && !(ngamma >= 2)) continue;
      if(controlRegion == kCR1 && !(ngamma == 0 && nfake == 1)) continue;
      if(controlRegion == kCR2 && !(ngamma == 0 && nfake >= 2)) continue;
      if(controlRegion == kCR2 && !((ngamma == 1 && nfake == 1) || (ngamma == 0 && nfake >= 2))) continue;
      if(controlRegion == kCR0 && !(ngamma == 0)) continue;

      if(!checkBtagging()) continue;

      if(topPtReweighting < 0) topPtReweighting = 1.;
      
      Float_t leptonSF, leptonSFup, leptonSFdown;
      Float_t photonSF, photonSFup, photonSFdown;

      GetLeptonSF(lepton_pt, lepton_eta, leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(lead_photon_et, lead_photon_eta, trail_photon_et, trail_photon_eta, nphotons, 
		  photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF * topPtReweighting;

      contamination += totalWeight;
    }

    h_acc->SetBinContent(h_acc->FindBin(index1, index2), h->Integral() / (0.438/3.) / 15000.);
    h_contamination->SetBinContent(h_contamination->FindBin(index1, index2), contamination / h->Integral());

    double xsec = h_xsec->GetBinContent(h_xsec->FindBin(index1, index2));
    
    h->Scale(xsec * 19712. / 15000.);
    h_btagWeightUp->Scale(xsec * 19712. / 15000.);
    h_btagWeightDown->Scale(xsec * 19712. / 15000.);
    h_puWeightUp->Scale(xsec * 19712. / 15000.);
    h_puWeightDown->Scale(xsec * 19712. / 15000.);
    h_topPtUp->Scale(xsec * 19712. / 15000.);
    h_topPtDown->Scale(xsec * 19712. / 15000.);
    h_JECup->Scale(xsec * 19712. / 15000.);
    h_JECdown->Scale(xsec * 19712. / 15000.);
    h_leptonSFup->Scale(xsec * 19712. / 15000.);
    h_leptonSFdown->Scale(xsec * 19712. / 15000.);
    h_photonSFup->Scale(xsec * 19712. / 15000.);
    h_photonSFdown->Scale(xsec * 19712. / 15000.);

    fSignalOut->cd();

    if(req.Contains("ele")) {
      fSignalOut->cd("ele_"+crNames[controlRegion]);
    }
    else {
      fSignalOut->cd("muon_"+crNames[controlRegion]);
    }

    h->Write("signal"+code_t);

    for(int j = 0; j < h->GetNbinsX(); j++) {
      TH1D * h_flux_up = (TH1D*)h->Clone("clone_signal_"+code_t+"_flux_up");
      TH1D * h_flux_down = (TH1D*)h->Clone("clone_signal"+code_t+"_flux_down");
      
      Double_t centralValue = h->GetBinContent(j+1);
      Double_t statError = h->GetBinError(j+1);
      
      if(statError > 0.) h_flux_up->SetBinContent(j+1, centralValue + statError);
      if(centralValue > statError && statError > 0.) h_flux_down->SetBinContent(j+1, centralValue - statError);
      
      h_flux_up->Write("signal"+code_t+"_signal_stat_bin"+Form("%d", j+1)+"Up");
      h_flux_down->Write("signal"+code_t+"_signal_stat_bin"+Form("%d", j+1)+"Down");
    }
    
    h_btagWeightUp->Write("signal"+code_t+"_btagWeightUp");
    h_btagWeightDown->Write("signal"+code_t+"_btagWeightDown");
    h_puWeightUp->Write("signal"+code_t+"_puWeightUp");
    h_puWeightDown->Write("signal"+code_t+"_puWeightDown");
    h_topPtUp->Write("signal"+code_t+"_topPtUp");
    h_topPtDown->Write("signal"+code_t+"_topPtDown");
    h_JECup->Write("signal"+code_t+"_JECUp");
    h_JECdown->Write("signal"+code_t+"_JECDown");
    h_leptonSFup->Write("signal"+code_t+"_leptonSFUp");
    h_leptonSFdown->Write("signal"+code_t+"_leptonSFDown");
    h_photonSFup->Write("signal"+code_t+"_photonSFUp");
    h_photonSFdown->Write("signal"+code_t+"_photonSFDown");

    f->Close();

  }

  // draw acc and etc
  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);
  fillPotHoles(h_acc);
  h_acc->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_acc->GetXaxis()->SetRangeUser(0, 1600);
  h_acc->GetXaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_acc->GetYaxis()->SetTitleOffset(1.3);
  h_acc->GetYaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetRangeUser(0, 1600);
  h_acc->GetZaxis()->SetLabelSize(0.02);
  h_acc->Draw("colz");

  TString plotName = "acc/acceptance_"+req+"_";
  if(controlRegion == kSR1) plotName += "SR1";
  if(controlRegion == kSR2) plotName += "SR2";
  if(controlRegion == kCR1) plotName += "CR1";
  if(controlRegion == kCR2) plotName += "CR2";
  if(controlRegion == kCR2a) plotName += "CR2a";
  if(controlRegion == kCR0) plotName += "CR0";
  plotName += ".pdf";

  can->SaveAs(plotName);
  
  h_contamination->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_contamination->GetXaxis()->SetRangeUser(0, 1600);
  h_contamination->GetXaxis()->SetLabelSize(0.03);
  h_contamination->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_contamination->GetYaxis()->SetTitleOffset(1.3);
  h_contamination->GetYaxis()->SetLabelSize(0.03);
  h_contamination->GetYaxis()->SetRangeUser(0, 1600);
  h_contamination->GetZaxis()->SetLabelSize(0.02);
  h_contamination->Draw("colz");

  plotName = "acc/contamination_"+req+"_";
  if(controlRegion == kSR1) plotName += "SR1";
  if(controlRegion == kSR2) plotName += "SR2";
  if(controlRegion == kCR1) plotName += "CR1";
  if(controlRegion == kCR2) plotName += "CR2";
  if(controlRegion == kCR2a) plotName += "CR2a";
  if(controlRegion == kCR0) plotName += "CR0";
  plotName += ".pdf";

  can->SaveAs(plotName);

  delete can;

  fSignalOut->Close();

  f_xsec->Close();
  
}
