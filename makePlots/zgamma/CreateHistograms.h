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

TString channels[nChannels] = {"ele_bjj", "muon_bjj",
			       "ele_jjj", "muon_jjj"};

unsigned int nBtagReq[nChannels] = {1, 1,
				    0, 0};

TString qcdChannels[nChannels] = {"ele_bjj_eQCDTree", "muon_bjj_muQCDTree",
				  "ele_jjj_eQCDTree", "muon_jjj_muQCDTree"};
TString qcdChannels_fakePhotons[nChannels] = {"ele_bjj_eQCDfakeTree", "muon_bjj_muQCDfakeTree",
					      "ele_jjj_eQCDfakeTree", "muon_jjj_muQCDfakeTree"};

enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};

TString crNames[kNumControlRegions] = {"SR1", "SR2", "CR1", "CR2", "CR2a", "CR0", "SigmaPlot", "Any"};

enum photonModes {kSignal, kFake, kNumPhotonModes};

class HistogramMaker : public TObject {
  
  ClassDef(HistogramMaker, 1);

 public:
  HistogramMaker(int chanNo, bool blind, int cRegion, Float_t cutOnMet, int mode, TString metType_, bool useNormalTopReweighting_);
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

  bool onZMass() {
    Float_t zmass = getValue("z_mass");
    return zmass > 81.0 && zmass < 101.0;
  }

  bool hasGoodPhotons() {

    return true;

    if(getIntegerValue("Nphotons") == 0) return true;
    if(getIntegerValue("Nphotons") == 1) {
      //bool chHadIso = getValue("leadChargedHadronIso") < 2.6;
      bool nHadIso = getValue("leadNeutralHadronIso") < 3.5;
      bool photonIso = getValue("leadPhotonIso") < 1.3;
      //bool sIetaIeta = getValue("leadSigmaIetaIeta") < 0.012;
      return nHadIso && photonIso;
    }
    if(getIntegerValue("Nphotons") >= 2) {
      //bool chHadIso = getValue("leadChargedHadronIso") < 2.6 && getValue("trailChargedHadronIso") < 2.6;
      bool nHadIso = getValue("leadNeutralHadronIso") < 3.5 && getValue("trailNeutralHadronIso") < 3.5;
      bool photonIso = getValue("leadPhotonIso") < 1.3 && getValue("trailPhotonIso") < 1.3;
      //bool sIetaIeta = getValue("leadSigmaIetaIeta") < 0.012 && getValue("trailSigmaIetaIeta") < 0.012;
      return nHadIso && photonIso;
    }

    return false;
  };

  bool inControlRegion() {
    switch(controlRegion) {
    case kSR1:
      return getIntegerValue("Ngamma") == 1;
    case kSR2:
      return getIntegerValue("Ngamma") >= 2;
    case kCR1:
      return (getIntegerValue("Ngamma") == 0 && getIntegerValue("Nfake") == 1) && hasGoodPhotons();
    case kCR2:
      return (getIntegerValue("Ngamma") == 0 && getIntegerValue("Nfake") >= 2) && hasGoodPhotons();
    case kCR2a:
      return (getIntegerValue("Ngamma") + getIntegerValue("Nfake") >= 2) && hasGoodPhotons();
    case kCR0:
      return getIntegerValue("Ngamma") == 0 && hasGoodPhotons();
    case kSigmaPlot:
      return getIntegerValue("Nphotons") == 1 && hasGoodPhotons();
    case kAny:
      return true;
    default:
      return false;
    }
  };

  bool isBlindedRegion() { return (blinded && (controlRegion == kSR2 || controlRegion == kSR1)); };

  bool passMetCut() {
    if(metCut > 0.) return (getValue(metType) < metCut);
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
			bool remove_whizard = false, bool remove_madgraph = false, bool reweightTop = false,
			Double_t fitScaling = -1.0, Double_t fitScalingError = 0.0);
  
  void SetTrees(TTree * gg, TTree * qcd);
  
  void BookHistogram(TString variable, Int_t nBins, Float_t xlo, Float_t xhi);
  void BookHistogram(TString variable, Int_t nBins, Double_t* customBins);

  void FillData();
  void FillQCD();
  void FillMCBackgrounds();

  void FillHistograms() {
    cout << endl << "Filling Data..." << endl;
    FillData();
    cout << "Filling QCD..." << endl;
    FillQCD();
    cout << "Filling MC Backgrounds..." << endl;
    FillMCBackgrounds();
    cout << "Done filling!" << endl << endl;
  }

  void SubtractMCFromQCD();

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

  void FillWithError(TH1D*& h, Float_t x, Float_t weight, Float_t addError2) {
    Float_t oldError = h->GetBinError(h->FindBin(x));
    Float_t newError = sqrt(oldError*oldError + addError2);
    h->Fill(x, weight);
    h->SetBinError(h->FindBin(x), newError);
  };

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
  vector<bool> removeWhizardOverlap;
  vector<bool> removeMadgraphOverlap;
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

  vector< vector<TH1D*> > mcQCDHistograms;
  vector< vector<TH1D*> > mcQCDHistograms_relIso_10;
  vector< vector<TH1D*> > mcQCDHistograms_relIso_m10;

  // h_xyz[variable]
  vector<TH1D*> h_gg;

  vector<TH1D*> h_qcd;

  vector<TH1D*> h_qcd_relIso_10;
  vector<TH1D*> h_qcd_relIso_m10;

  TTree * ggTree;
  TTree * qcdTree;
  
  vector<TString> variables;
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
  int photonMode;
  TString metType;
  bool useNormalTopReweighting;

  TString req;

  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr;
  Float_t puWeightUp, puWeightDown, btagWeightUp, btagWeightDown;
  Float_t overlaps_whizard, overlaps_madgraph;
  Float_t topPtReweighting;

  Float_t relIso1, relIso2;

};

HistogramMaker::HistogramMaker(int chanNo, bool blind, int cRegion, Float_t cutOnMet, int mode, TString metType_, bool useNormalTopReweighting_) :
  channel(chanNo),
  blinded(blind),
  controlRegion(cRegion),
  metCut(cutOnMet),
  photonMode(mode),
  metType(metType_),
  useNormalTopReweighting(useNormalTopReweighting_)
{
  req = channels[chanNo];

  intLumi_int = 19712;

  variables.clear();
  varMap.clear();

  h_gg.clear();

  h_qcd.clear();
  h_qcd_relIso_10.clear();
  h_qcd_relIso_m10.clear();

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

  removeWhizardOverlap.clear();
  removeMadgraphOverlap.clear();
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

  mcQCDHistograms.clear();
  mcQCDHistograms_relIso_10.clear();
  mcQCDHistograms_relIso_m10.clear();

}

HistogramMaker::~HistogramMaker() { 

  variables.clear();
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

  mcQCDHistograms.clear();
  mcQCDHistograms_relIso_10.clear();
  mcQCDHistograms_relIso_m10.clear();

  mcTrees.clear();
  mcTrees_JECup.clear();
  mcTrees_JECdown.clear();
  mcFiles.clear();
  mcQCDTrees.clear();
  mcNGen.clear();
  crossSections.clear();
  scaleErrUp.clear();
  scaleErrDown.clear();
  pdfErrUp.clear();
  pdfErrDown.clear();

  mcNames.clear();

  removeWhizardOverlap.clear();
  removeMadgraphOverlap.clear();
  reweightTopPt.clear();

  fitScale.clear();
  fitScaleError.clear();

  fLeptonSF->Close();
  fPhotonSF->Close();

  h_gg.clear();

  h_qcd.clear();
  h_qcd_relIso_10.clear();
  h_qcd_relIso_m10.clear();

}

void HistogramMaker::SetTrees(TTree * gg, TTree * qcd) {
  ggTree = gg;
  qcdTree = qcd;
}

bool HistogramMaker::LoadMCBackground(TString fileName, TString scanName,
				      Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfErrorUp, Double_t pdfErrorDown,
				      bool remove_whizard, bool remove_madgraph, bool reweightTop,
				      Double_t fitScaling, Double_t fitScalingError) {

  mcFiles.push_back(new TFile(fileName, "READ"));
  if(!mcFiles.back()) {
    cout << "Could not load TFile " << fileName << endl;
    return false;
  }

  TString signalString, qcdString;
  if(controlRegion == kSR1 || controlRegion == kSR2 || controlRegion == kCR0 || controlRegion == kAny) {
    signalString = channels[channel]+"_signalTree";
    qcdString = qcdChannels[channel];
  }
  else if(photonMode == kFake) {
    signalString = channels[channel]+"_fakeTree";
    qcdString = qcdChannels_fakePhotons[channel];
  }
  else {
    cout << "Invalid photonMode!" << endl;
    return false;
  }

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

  removeWhizardOverlap.push_back(remove_whizard);
  removeMadgraphOverlap.push_back(remove_madgraph);
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

  mcQCDHistograms.resize(mcQCDHistograms.size() + 1);
  mcQCDHistograms_relIso_10.resize(mcQCDHistograms_relIso_10.size() + 1);
  mcQCDHistograms_relIso_m10.resize(mcQCDHistograms_relIso_m10.size() + 1);

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

  h_qcd_relIso_10.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_10_"+req));
  h_qcd_relIso_m10.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_m10_"+req));

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

    mcQCDHistograms_relIso_10[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_10_"+mcNames[i]+"_"+req));
    mcQCDHistograms_relIso_m10[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_m10_"+mcNames[i]+"_"+req));
  }

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

  h_qcd_relIso_10.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_10_"+req));
  h_qcd_relIso_m10.push_back((TH1D*)qcd->Clone(variable+"_qcd_relIso_m10_"+req));

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

    mcQCDHistograms_relIso_10[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_10_"+mcNames[i]+"_"+req));
    mcQCDHistograms_relIso_m10[i].push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_m10_"+mcNames[i]+"_"+req));
  }
  
}

void HistogramMaker::FillData() {
  
  for(unsigned int i = 0; i < variables.size(); i++) ggTree->SetBranchAddress(variables[i], &(varMap[variables[i]]));

  for(int i = 0; i < ggTree->GetEntries(); i++) {
    ggTree->GetEntry(i);

    if(!passMetCut()) continue;
    if(!inControlRegion()) continue;
    if(isBlindedRegion()) continue;

    bool isFromZ = onZMass();

    for(unsigned int j = 0; j < variables.size(); j++) {
      if(variables[j] != "z_mass" && !isFromZ) continue;
      h_gg[j]->Fill(getValue(j));
    }

  }

  ggTree->ResetBranchAddresses();

}

void HistogramMaker::FillQCD() {

  for(unsigned int i = 0; i < variables.size(); i++) qcdTree->SetBranchAddress(variables[i], &(varMap[variables[i]]));

  if(req.Contains("ele")) {
    qcdTree->SetBranchAddress("ele1_relIso", &relIso1);
    qcdTree->SetBranchAddress("ele2_relIso", &relIso2);
  }
  else {
    qcdTree->SetBranchAddress("muon1_relIso", &relIso1);
    qcdTree->SetBranchAddress("muon2_relIso", &relIso2);
  }

  for(int i = 0; i < qcdTree->GetEntries(); i++) {
    qcdTree->GetEntry(i);
    
    if(!passMetCut()) continue;
    if(!inControlRegion()) continue;

    bool isFromZ = onZMass();

    for(unsigned int j = 0; j < variables.size(); j++) {

      if(variables[j] != "z_mass" && !isFromZ) continue;

      Float_t val = getValue(j);

      if(relIso1 > 0.25 && relIso2 > 0.25) h_qcd[j]->Fill(val); // central
      if(relIso1 > 0.25 * 1.1 && relIso2 > 0.25 * 1.1) h_qcd_relIso_10[j]->Fill(val); // +10%
      if(relIso1 > 0.25 * 0.9 && relIso2 > 0.25 * 0.9) h_qcd_relIso_m10[j]->Fill(val); // -10%

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

    if(req.Contains("ele")) {
      mcQCDTrees[i]->SetBranchAddress("ele1_relIso", &relIso1);
      mcQCDTrees[i]->SetBranchAddress("ele2_relIso", &relIso2);
    }
    else {
      mcQCDTrees[i]->SetBranchAddress("muon1_relIso", &relIso1);
      mcQCDTrees[i]->SetBranchAddress("muon2_relIso", &relIso2);
    }

    if(removeWhizardOverlap[i]) {
      mcTrees[i]->SetBranchAddress("overlaps_whizard", &overlaps_whizard);
      mcTrees_JECup[i]->SetBranchAddress("overlaps_whizard", &overlaps_whizard);
      mcTrees_JECdown[i]->SetBranchAddress("overlaps_whizard", &overlaps_whizard);
      mcQCDTrees[i]->SetBranchAddress("overlaps_whizard", &overlaps_whizard);
    }

    if(removeMadgraphOverlap[i]) {
      mcTrees[i]->SetBranchAddress("overlaps_madgraph", &overlaps_madgraph);
      mcTrees_JECup[i]->SetBranchAddress("overlaps_madgraph", &overlaps_madgraph);
      mcTrees_JECdown[i]->SetBranchAddress("overlaps_madgraph", &overlaps_madgraph);
      mcQCDTrees[i]->SetBranchAddress("overlaps_madgraph", &overlaps_madgraph);
    }

    TString topPtReweightingName = (useNormalTopReweighting) ? "TopPtReweighting" : "TopPtReweighting_ttHbb";

    if(reweightTopPt[i]) {
      mcTrees[i]->SetBranchAddress(topPtReweightingName, &topPtReweighting);
      mcTrees_JECup[i]->SetBranchAddress(topPtReweightingName, &topPtReweighting);
      mcTrees_JECdown[i]->SetBranchAddress(topPtReweightingName, &topPtReweighting);
      mcQCDTrees[i]->SetBranchAddress(topPtReweightingName, &topPtReweighting);
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
      
      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      Float_t addError2 = puWeight*puWeight*btagWeightErr*btagWeightErr + btagWeight*btagWeight*puWeightErr*puWeightErr;
      if(fitScale[i] > 0) addError2 = fitScale[i]*fitScale[i]*puWeight*puWeight*btagWeightErr*btagWeightErr + 
			    fitScale[i]*fitScale[i]*btagWeight*btagWeight*puWeightErr*puWeightErr +
			    puWeight*puWeight*btagWeight*btagWeight*fitScaleError[i]*fitScaleError[i];
      
      Float_t addError2_puOnly = btagWeight*btagWeight*puWeightErr*puWeightErr;
      if(fitScale[i] > 0) addError2_puOnly = fitScale[i]*fitScale[i]*btagWeight*btagWeight*puWeightErr*puWeightErr +
			    puWeight*puWeight*btagWeight*btagWeight*fitScaleError[i]*fitScaleError[i];
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      if(fitScale[i] > 0) totalWeight *= fitScale[i];

      if(totalWeight < 1.e-6) continue;
      if(addError2 != addError2) continue;

      bool isFromZ = onZMass();

      for(unsigned int k = 0; k < variables.size(); k++) {
	
	if(variables[k] != "z_mass" && !isFromZ) continue;

	Float_t val = getValue(k);

	FillWithError(mcHistograms[i][k], val, totalWeight, addError2);
	
	mcHistograms_btagWeightUp[i][k]->Fill(val, totalWeight * btagWeightUp/btagWeight);
	mcHistograms_btagWeightDown[i][k]->Fill(val, totalWeight * btagWeightDown/btagWeight);
	
	mcHistograms_puWeightUp[i][k]->Fill(val, totalWeight * puWeightUp/puWeight);
	mcHistograms_puWeightDown[i][k]->Fill(val, totalWeight * puWeightDown/puWeight);

	mcHistograms_leptonSFup[i][k]->Fill(val, totalWeight * leptonSFup/leptonSF);
	mcHistograms_leptonSFdown[i][k]->Fill(val, totalWeight * leptonSFdown/leptonSF);

	mcHistograms_photonSFup[i][k]->Fill(val, totalWeight * photonSFup/photonSF);
	mcHistograms_photonSFdown[i][k]->Fill(val, totalWeight * photonSFdown/photonSF);

	if(reweightTopPt[i]) {
	  mcHistograms_topPtUp[i][k]->Fill(val, totalWeight * topPtReweighting);
	  mcHistograms_topPtDown[i][k]->Fill(val, totalWeight / topPtReweighting);
	}
	else {
	  mcHistograms_topPtUp[i][k]->Fill(val, totalWeight);
	  mcHistograms_topPtDown[i][k]->Fill(val, totalWeight);
	}
	
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

  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECup[i]->GetEntries(); j++) {
      mcTrees_JECup[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!passMetCut()) continue;
      if(!inControlRegion()) continue;
      
      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      if(fitScale[i] > 0) totalWeight *= fitScale[i];
      
      bool isFromZ = onZMass();

      for(unsigned int k = 0; k < variables.size(); k++) {
	if(variables[k] != "z_mass" && !isFromZ) continue;
	mcHistograms_JECup[i][k]->Fill(getValue(k), totalWeight);
      }
      
    }
    
    for(unsigned int j = 0; j < variables.size(); j++) mcHistograms_JECup[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    
  }
  
  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECdown[i]->GetEntries(); j++) {
      mcTrees_JECdown[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!passMetCut()) continue;
      if(!inControlRegion()) continue;
      
      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      if(fitScale[i] > 0) totalWeight *= fitScale[i];
      
      bool isFromZ = onZMass();

      for(unsigned int k = 0; k < variables.size(); k++) {
	if(variables[k] != "z_mass" && !isFromZ) continue;
	mcHistograms_JECdown[i][k]->Fill(getValue(k), totalWeight);
      }
      
    }
    
    for(unsigned int j = 0; j < variables.size(); j++) mcHistograms_JECdown[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    
  }

  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) {
    
    for(int j = 0; j < mcQCDTrees[i]->GetEntries(); j++) {
      mcQCDTrees[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!passMetCut()) continue;
      if(!inControlRegion()) continue;
      
      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      if(fitScale[i] > 0) totalWeight *= fitScale[i];
      
      Float_t addError2 = puWeightErr*puWeightErr/puWeight/puWeight + btagWeightErr*btagWeightErr/btagWeight/btagWeight;
      if(fitScale[i] > 0.) addError2 += fitScaleError[i]*fitScaleError[i]/fitScale[i]/fitScale[i];
      addError2 *= totalWeight*totalWeight;
      
      if(totalWeight < 1.e-6) continue;
      if(addError2 != addError2) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      bool isFromZ = onZMass();

      for(unsigned int k = 0; k < variables.size(); k++) {
	
	if(variables[k] != "z_mass" && !isFromZ) continue;

	Float_t val = getValue(k);
	
	if(relIso1 > 0.25 && relIso2 > 0.25) FillWithError(mcQCDHistograms[i][k], val, totalWeight, addError2); // central
	
	if(relIso1 > 0.25 * 1.1 && relIso2 > 0.25 * 1.1) mcQCDHistograms_relIso_10[i][k]->Fill(val, totalWeight); // +10%
	if(relIso1 > 0.25 * 0.9 && relIso2 > 0.25 * 0.9) mcQCDHistograms_relIso_m10[i][k]->Fill(val, totalWeight); // -10%
	
      }
      
    }

    for(unsigned int j = 0; j < variables.size(); j++) mcQCDHistograms[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    for(unsigned int j = 0; j < variables.size(); j++) mcQCDHistograms_relIso_10[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    for(unsigned int j = 0; j < variables.size(); j++) mcQCDHistograms_relIso_m10[i][j]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);

  }
  
  for(unsigned int i = 0; i < mcTrees.size(); i++) mcTrees[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) mcTrees_JECup[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) mcTrees_JECdown[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) mcQCDTrees[i]->ResetBranchAddresses();

}

void HistogramMaker::SubtractMCFromQCD() {

  for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms[i].size(); j++) {
      h_qcd[j]->Add(mcQCDHistograms[i][j], -1.);
    }
  }

  for(unsigned int i = 0; i < mcQCDHistograms_relIso_10.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms_relIso_10[i].size(); j++) {
      h_qcd_relIso_10[j]->Add(mcQCDHistograms_relIso_10[i][j], -1.);
    }
  }

  for(unsigned int i = 0; i < mcQCDHistograms_relIso_m10.size(); i++) {
    for(unsigned int j = 0; j < mcQCDHistograms_relIso_m10[i].size(); j++) {
      h_qcd_relIso_m10[j]->Add(mcQCDHistograms_relIso_m10[i][j], -1.);
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

  for(unsigned int i = 0; i < h_qcd_relIso_10.size(); i++) {
    for(Int_t j = 0; j < h_qcd_relIso_10[i]->GetNbinsX(); j++) {
      if(h_qcd_relIso_10[i]->GetBinContent(j+1) < 0) {
	h_qcd_relIso_10[i]->SetBinContent(j+1, 0.);
	h_qcd_relIso_10[i]->SetBinError(j+1, 0.);
      }
    }
  }

  for(unsigned int i = 0; i < h_qcd_relIso_m10.size(); i++) {
    for(Int_t j = 0; j < h_qcd_relIso_m10[i]->GetNbinsX(); j++) {
      if(h_qcd_relIso_m10[i]->GetBinContent(j+1) < 0) {
	h_qcd_relIso_m10[i]->SetBinContent(j+1, 0.);
	h_qcd_relIso_m10[i]->SetBinError(j+1, 0.);
      }
    }
  }

}

void HistogramMaker::SaveOutput() {

  TString outName = "histograms_"+req+"_"+crNames[controlRegion]+".root";

  TFile * fOut = new TFile(outName, "UPDATE");

  for(unsigned int i = 0; i < h_gg.size(); i++) h_gg[i]->Write();
  for(unsigned int i = 0; i < h_qcd.size(); i++) h_qcd[i]->Write();
  for(unsigned int i = 0; i < h_qcd_relIso_10.size(); i++) h_qcd_relIso_10[i]->Write();
  for(unsigned int i = 0; i < h_qcd_relIso_m10.size(); i++) h_qcd_relIso_m10[i]->Write();

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

  fOut->Close();

}

void HistogramMaker::GetLeptonSF(Float_t& central, Float_t& up, Float_t& down) {

  Float_t lead_pt, trail_pt, lead_eta, trail_eta;
  if(req.Contains("ele")) {
    lead_pt = getValue("ele1_pt");
    trail_pt = getValue("ele2_pt");
    lead_eta = getValue("ele1_eta");
    trail_eta = getValue("ele2_eta");
  }
  else {
    lead_pt = getValue("muon1_pt");
    trail_pt = getValue("muon2_pt");
    lead_eta = getValue("muon1_eta");
    trail_eta = getValue("muon2_eta");
  } 

  Float_t pt1, pt2;
  Float_t eta1, eta2;
  Float_t error;

  if(req.Contains("ele")) {
    pt1 = min(lead_pt, (float)199.);
    pt1 = max(pt1, (float)15.);
    eta1 = min(fabs(lead_eta), (double)2.39);

    pt2 = min(trail_pt, (float)199.);
    pt2 = max(pt2, (float)15.);
    eta2 = min(fabs(trail_eta), (double)2.39);

    Float_t id1_val = sf_electron->GetBinContent(sf_electron->FindBin(eta1, pt1));
    Float_t id1_error = sf_electron->GetBinError(sf_electron->FindBin(eta1, pt1));

    Float_t trigger1_val = sf_SingleElectronTrigger->GetBinContent(sf_SingleElectronTrigger->FindBin(eta1, pt1));
    Float_t trigger1_error = sf_SingleElectronTrigger->GetBinError(sf_SingleElectronTrigger->FindBin(eta1, pt1));

    Float_t id2_val = sf_electron->GetBinContent(sf_electron->FindBin(eta2, pt2));
    Float_t id2_error = sf_electron->GetBinError(sf_electron->FindBin(eta2, pt2));

    Float_t trigger2_val = sf_SingleElectronTrigger->GetBinContent(sf_SingleElectronTrigger->FindBin(eta2, pt2));
    Float_t trigger2_error = sf_SingleElectronTrigger->GetBinError(sf_SingleElectronTrigger->FindBin(eta2, pt2));

    Float_t x = 1. - trigger1_val;
    Float_t y = 1. - trigger2_val;

    Float_t trigger_val = 1. - x*y;
    Float_t trigger_error = x * y * sqrt(trigger1_error*trigger1_error/x/x + trigger2_error*trigger2_error/y/y);


    central = id1_val * id2_val * trigger_val;
    error = central * sqrt(id1_error*id1_error/id1_val/id1_val + id2_error*id2_error/id2_val/id2_val + trigger_error*trigger_error/trigger_val/trigger_val);

    up = central + error;
    down = central - error;
  }

  else {
    pt1 = min(lead_pt, (float)499.);
    pt1 = max(pt1, (float)10.);
    eta1 = min(fabs(lead_eta), (double)2.09);

    pt2 = min(trail_pt, (float)499.);
    pt2 = max(pt2, (float)10.);
    eta2 = min(fabs(trail_eta), (double)2.09);

    Float_t id1_val = sf_muon->GetBinContent(sf_muon->FindBin(pt1, eta1));
    Float_t id1_error = (sf_muon->GetBinError(sf_muon->FindBin(pt1, eta1)) * id1_val) - id1_val;

    Float_t id2_val = sf_muon->GetBinContent(sf_muon->FindBin(pt2, eta2));
    Float_t id2_error = (sf_muon->GetBinError(sf_muon->FindBin(pt2, eta2)) * id2_val) - id2_val;

    central = id1_val * id2_val;
    error = central * sqrt(id1_error*id1_error/id1_val/id1_val + id2_error*id2_error/id2_val/id2_val);

    up = central + error;
    down = central - error;

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

