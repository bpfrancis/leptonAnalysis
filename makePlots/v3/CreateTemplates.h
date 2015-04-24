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

enum controlRegions {kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny, kNumControlRegions};

TString crNames[kNumControlRegions] = {"SR1", "SR2", "CR1", "CR2", "CR2a", "CR0", "SigmaPlot", "Any"};

enum photonModes {kSignal, kFake, kNumPhotonModes};

class TemplateMaker : public TObject {
  
  ClassDef(TemplateMaker, 1);

 public:
  TemplateMaker(TString var, TString chan, int cRegion, Int_t n_Bins, Float_t x_lo, Float_t x_hi, Float_t cutOnMet);
  ~TemplateMaker();
  
  bool inControlRegion() {
    switch(controlRegion) {
    case kSR1:
      return Ngamma == 1;
    case kSR2:
      return Ngamma >= 2;
    case kCR1:
      return (Ngamma == 0 && Nfake == 1);
    case kCR2:
      return (Ngamma == 0 && Nfake >= 2);
    case kCR2a:
      return (Ngamma + Nfake >= 2);
    case kCR0:
      return Ngamma == 0;
    case kSigmaPlot:
      return Nphotons == 1;
    case kAny:
      return true;
    default:
      return false;
    }
  };

  bool checkBtagging() {
    if(!(channel.Contains("b"))) {
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
			bool remove_whizard = false, bool remove_madgraph = false, bool reweightTop = false);
  
  void SetTrees(TTree * gg, TTree * qcd);
  void SetAddresses();
  void BookTemplates();

  //durp
  void FillData();
  void FillQCD();
  void FillMCBackgrounds();

  void FillTemplates() {
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

  // mcHistograms[process]
  vector<TH1D*> mcHistograms;
  vector<TH1D*> mcHistograms_btagWeightUp;
  vector<TH1D*> mcHistograms_btagWeightDown;
  vector<TH1D*> mcHistograms_puWeightUp;
  vector<TH1D*> mcHistograms_puWeightDown;
  vector<TH1D*> mcHistograms_scaleUp;
  vector<TH1D*> mcHistograms_scaleDown;
  vector<TH1D*> mcHistograms_pdfUp;
  vector<TH1D*> mcHistograms_pdfDown;
  vector<TH1D*> mcHistograms_topPtUp;
  vector<TH1D*> mcHistograms_topPtDown;
  vector<TH1D*> mcHistograms_JECup;
  vector<TH1D*> mcHistograms_JECdown;
  vector<TH1D*> mcHistograms_leptonSFup;
  vector<TH1D*> mcHistograms_leptonSFdown;
  vector<TH1D*> mcHistograms_photonSFup;
  vector<TH1D*> mcHistograms_photonSFdown;

  vector<TH1D*> mcHistograms_matchPhoton;
  vector<TH1D*> mcHistograms_matchPhoton_btagWeightUp;
  vector<TH1D*> mcHistograms_matchPhoton_btagWeightDown;
  vector<TH1D*> mcHistograms_matchPhoton_puWeightUp;
  vector<TH1D*> mcHistograms_matchPhoton_puWeightDown;
  vector<TH1D*> mcHistograms_matchPhoton_scaleUp;
  vector<TH1D*> mcHistograms_matchPhoton_scaleDown;
  vector<TH1D*> mcHistograms_matchPhoton_pdfUp;
  vector<TH1D*> mcHistograms_matchPhoton_pdfDown;
  vector<TH1D*> mcHistograms_matchPhoton_topPtUp;
  vector<TH1D*> mcHistograms_matchPhoton_topPtDown;
  vector<TH1D*> mcHistograms_matchPhoton_JECup;
  vector<TH1D*> mcHistograms_matchPhoton_JECdown;
  vector<TH1D*> mcHistograms_matchPhoton_leptonSFup;
  vector<TH1D*> mcHistograms_matchPhoton_leptonSFdown;
  vector<TH1D*> mcHistograms_matchPhoton_photonSFup;
  vector<TH1D*> mcHistograms_matchPhoton_photonSFdown;

  vector<TH1D*> mcHistograms_matchElectron;
  vector<TH1D*> mcHistograms_matchElectron_btagWeightUp;
  vector<TH1D*> mcHistograms_matchElectron_btagWeightDown;
  vector<TH1D*> mcHistograms_matchElectron_puWeightUp;
  vector<TH1D*> mcHistograms_matchElectron_puWeightDown;
  vector<TH1D*> mcHistograms_matchElectron_scaleUp;
  vector<TH1D*> mcHistograms_matchElectron_scaleDown;
  vector<TH1D*> mcHistograms_matchElectron_pdfUp;
  vector<TH1D*> mcHistograms_matchElectron_pdfDown;
  vector<TH1D*> mcHistograms_matchElectron_topPtUp;
  vector<TH1D*> mcHistograms_matchElectron_topPtDown;
  vector<TH1D*> mcHistograms_matchElectron_JECup;
  vector<TH1D*> mcHistograms_matchElectron_JECdown;
  vector<TH1D*> mcHistograms_matchElectron_leptonSFup;
  vector<TH1D*> mcHistograms_matchElectron_leptonSFdown;
  vector<TH1D*> mcHistograms_matchElectron_photonSFup;
  vector<TH1D*> mcHistograms_matchElectron_photonSFdown;

  vector<TH1D*> mcHistograms_matchJet;
  vector<TH1D*> mcHistograms_matchJet_btagWeightUp;
  vector<TH1D*> mcHistograms_matchJet_btagWeightDown;
  vector<TH1D*> mcHistograms_matchJet_puWeightUp;
  vector<TH1D*> mcHistograms_matchJet_puWeightDown;
  vector<TH1D*> mcHistograms_matchJet_scaleUp;
  vector<TH1D*> mcHistograms_matchJet_scaleDown;
  vector<TH1D*> mcHistograms_matchJet_pdfUp;
  vector<TH1D*> mcHistograms_matchJet_pdfDown;
  vector<TH1D*> mcHistograms_matchJet_topPtUp;
  vector<TH1D*> mcHistograms_matchJet_topPtDown;
  vector<TH1D*> mcHistograms_matchJet_JECup;
  vector<TH1D*> mcHistograms_matchJet_JECdown;
  vector<TH1D*> mcHistograms_matchJet_leptonSFup;
  vector<TH1D*> mcHistograms_matchJet_leptonSFdown;
  vector<TH1D*> mcHistograms_matchJet_photonSFup;
  vector<TH1D*> mcHistograms_matchJet_photonSFdown;

  vector<TH1D*> mcQCDHistograms;

  vector<TH1D*> mcQCDHistograms_relIso_10;
  vector<TH1D*> mcQCDHistograms_relIso_m10;

  // h_xyz
  TH1D * h_gg;
  TH1D * h_qcd;
  TH1D * h_qcd_relIso_10;
  TH1D * h_qcd_relIso_m10;

  TTree * ggTree;
  TTree * qcdTree;
  
  
  Float_t value;

  TString variable;
  TString channel;
  int controlRegion;
  Int_t nBins;
  Float_t xlo, xhi;
  Float_t metCut;
    
  TFile * fLeptonSF;
  TH2D * sf_muon;
  TH2D * sf_electron;
  TH2D * sf_SingleElectronTrigger;

  TFile * fPhotonSF;
  TH2D * sf_photon_id;
  TH2D * sf_photon_veto;

  Int_t intLumi_int;

   // friend variables
  Float_t Ngamma, Nfake, Nphotons, met;

  Float_t leptonPt, leptonEta;
  Float_t leadGammaEt, leadGammaEta;
  Float_t trailGammaEt, trailGammaEta;

  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr;
  Float_t puWeightUp, puWeightDown, btagWeightUp, btagWeightDown;
  Float_t overlaps_whizard, overlaps_madgraph;
  Float_t topPtReweighting;

  Float_t leadMatchGamma, leadMatchElectron, leadMatchJet;

  Float_t relIso;

};

TemplateMaker::TemplateMaker(TString var, TString chan, int cRegion, Int_t n_Bins, Float_t x_lo, Float_t x_hi, Float_t cutOnMet) :
  variable(var),
  channel(chan),
  controlRegion(cRegion),
  nBins(n_Bins),
  xlo(x_lo),
  xhi(x_hi),
  metCut(cutOnMet)
{
  intLumi_int = 19712;

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

  mcHistograms_matchPhoton.clear();
  mcHistograms_matchPhoton_btagWeightUp.clear();
  mcHistograms_matchPhoton_btagWeightDown.clear();
  mcHistograms_matchPhoton_puWeightUp.clear();
  mcHistograms_matchPhoton_puWeightDown.clear();
  mcHistograms_matchPhoton_scaleUp.clear();
  mcHistograms_matchPhoton_scaleDown.clear();
  mcHistograms_matchPhoton_pdfUp.clear();
  mcHistograms_matchPhoton_pdfDown.clear();
  mcHistograms_matchPhoton_topPtUp.clear();
  mcHistograms_matchPhoton_topPtDown.clear();
  mcHistograms_matchPhoton_JECup.clear();
  mcHistograms_matchPhoton_JECdown.clear();
  mcHistograms_matchPhoton_leptonSFup.clear();
  mcHistograms_matchPhoton_leptonSFdown.clear();
  mcHistograms_matchPhoton_photonSFup.clear();
  mcHistograms_matchPhoton_photonSFdown.clear();

  mcHistograms_matchElectron.clear();
  mcHistograms_matchElectron_btagWeightUp.clear();
  mcHistograms_matchElectron_btagWeightDown.clear();
  mcHistograms_matchElectron_puWeightUp.clear();
  mcHistograms_matchElectron_puWeightDown.clear();
  mcHistograms_matchElectron_scaleUp.clear();
  mcHistograms_matchElectron_scaleDown.clear();
  mcHistograms_matchElectron_pdfUp.clear();
  mcHistograms_matchElectron_pdfDown.clear();
  mcHistograms_matchElectron_topPtUp.clear();
  mcHistograms_matchElectron_topPtDown.clear();
  mcHistograms_matchElectron_JECup.clear();
  mcHistograms_matchElectron_JECdown.clear();
  mcHistograms_matchElectron_leptonSFup.clear();
  mcHistograms_matchElectron_leptonSFdown.clear();
  mcHistograms_matchElectron_photonSFup.clear();
  mcHistograms_matchElectron_photonSFdown.clear();

  mcHistograms_matchJet.clear();
  mcHistograms_matchJet_btagWeightUp.clear();
  mcHistograms_matchJet_btagWeightDown.clear();
  mcHistograms_matchJet_puWeightUp.clear();
  mcHistograms_matchJet_puWeightDown.clear();
  mcHistograms_matchJet_scaleUp.clear();
  mcHistograms_matchJet_scaleDown.clear();
  mcHistograms_matchJet_pdfUp.clear();
  mcHistograms_matchJet_pdfDown.clear();
  mcHistograms_matchJet_topPtUp.clear();
  mcHistograms_matchJet_topPtDown.clear();
  mcHistograms_matchJet_JECup.clear();
  mcHistograms_matchJet_JECdown.clear();
  mcHistograms_matchJet_leptonSFup.clear();
  mcHistograms_matchJet_leptonSFdown.clear();
  mcHistograms_matchJet_photonSFup.clear();
  mcHistograms_matchJet_photonSFdown.clear();

  mcQCDHistograms.clear();

  mcQCDHistograms_relIso_10.clear();
  mcQCDHistograms_relIso_m10.clear();

}

TemplateMaker::~TemplateMaker() { 

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

  mcHistograms_matchPhoton.clear();
  mcHistograms_matchPhoton_btagWeightUp.clear();
  mcHistograms_matchPhoton_btagWeightDown.clear();
  mcHistograms_matchPhoton_puWeightUp.clear();
  mcHistograms_matchPhoton_puWeightDown.clear();
  mcHistograms_matchPhoton_scaleUp.clear();
  mcHistograms_matchPhoton_scaleDown.clear();
  mcHistograms_matchPhoton_pdfUp.clear();
  mcHistograms_matchPhoton_pdfDown.clear();
  mcHistograms_matchPhoton_topPtUp.clear();
  mcHistograms_matchPhoton_topPtDown.clear();
  mcHistograms_matchPhoton_JECup.clear();
  mcHistograms_matchPhoton_JECdown.clear();
  mcHistograms_matchPhoton_leptonSFup.clear();
  mcHistograms_matchPhoton_leptonSFdown.clear();
  mcHistograms_matchPhoton_photonSFup.clear();
  mcHistograms_matchPhoton_photonSFdown.clear();

  mcHistograms_matchElectron.clear();
  mcHistograms_matchElectron_btagWeightUp.clear();
  mcHistograms_matchElectron_btagWeightDown.clear();
  mcHistograms_matchElectron_puWeightUp.clear();
  mcHistograms_matchElectron_puWeightDown.clear();
  mcHistograms_matchElectron_scaleUp.clear();
  mcHistograms_matchElectron_scaleDown.clear();
  mcHistograms_matchElectron_pdfUp.clear();
  mcHistograms_matchElectron_pdfDown.clear();
  mcHistograms_matchElectron_topPtUp.clear();
  mcHistograms_matchElectron_topPtDown.clear();
  mcHistograms_matchElectron_JECup.clear();
  mcHistograms_matchElectron_JECdown.clear();
  mcHistograms_matchElectron_leptonSFup.clear();
  mcHistograms_matchElectron_leptonSFdown.clear();
  mcHistograms_matchElectron_photonSFup.clear();
  mcHistograms_matchElectron_photonSFdown.clear();

  mcHistograms_matchJet.clear();
  mcHistograms_matchJet_btagWeightUp.clear();
  mcHistograms_matchJet_btagWeightDown.clear();
  mcHistograms_matchJet_puWeightUp.clear();
  mcHistograms_matchJet_puWeightDown.clear();
  mcHistograms_matchJet_scaleUp.clear();
  mcHistograms_matchJet_scaleDown.clear();
  mcHistograms_matchJet_pdfUp.clear();
  mcHistograms_matchJet_pdfDown.clear();
  mcHistograms_matchJet_topPtUp.clear();
  mcHistograms_matchJet_topPtDown.clear();
  mcHistograms_matchJet_JECup.clear();
  mcHistograms_matchJet_JECdown.clear();
  mcHistograms_matchJet_leptonSFup.clear();
  mcHistograms_matchJet_leptonSFdown.clear();
  mcHistograms_matchJet_photonSFup.clear();
  mcHistograms_matchJet_photonSFdown.clear();

  mcQCDHistograms.clear();

  mcQCDHistograms_relIso_10.clear();
  mcQCDHistograms_relIso_m10.clear();

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

  removeWhizardOverlap.clear();
  removeMadgraphOverlap.clear();
  reweightTopPt.clear();

  fLeptonSF->Close();
  fPhotonSF->Close();

}

void TemplateMaker::SetTrees(TTree * gg, TTree * qcd) {

  ggTree = gg;
  qcdTree = qcd;
  
}

void TemplateMaker::SetAddresses() {

  ggTree->SetBranchAddress(variable, &value);
  ggTree->SetBranchAddress("Ngamma", &Ngamma);
  ggTree->SetBranchAddress("Nfake", &Nfake);
  ggTree->SetBranchAddress("Nphotons", &Nphotons);
  if(variable != "pfMET_t01") ggTree->SetBranchAddress("pfMET_t01", &met);
  if(channel.Contains("ele")) {
    ggTree->SetBranchAddress("ele_pt", &leptonPt);
    ggTree->SetBranchAddress("ele_eta", &leptonEta);
  }
  else {
    ggTree->SetBranchAddress("muon_pt", &leptonPt);
    ggTree->SetBranchAddress("muon_eta", &leptonEta);
  }
  ggTree->SetBranchAddress("leadPhotonEt", &leadGammaEt);
  ggTree->SetBranchAddress("leadPhotonEta", &leadGammaEta);
  ggTree->SetBranchAddress("trailPhotonEt", &trailGammaEt);
  ggTree->SetBranchAddress("trailPhotonEta", &trailGammaEta);

  qcdTree->SetBranchAddress(variable, &value);
  qcdTree->SetBranchAddress("Ngamma", &Ngamma);
  qcdTree->SetBranchAddress("Nfake", &Nfake);
  qcdTree->SetBranchAddress("Nphotons", &Nphotons);
  if(variable != "pfMET_t01") qcdTree->SetBranchAddress("pfMET_t01", &met);
  if(channel.Contains("ele")) {
    qcdTree->SetBranchAddress("ele_pt", &leptonPt);
    qcdTree->SetBranchAddress("ele_eta", &leptonEta);
    qcdTree->SetBranchAddress("ele_relIso", &relIso);
  }
  else {
    qcdTree->SetBranchAddress("muon_pt", &leptonPt);
    qcdTree->SetBranchAddress("muon_eta", &leptonEta);
    qcdTree->SetBranchAddress("muon_relIso", &relIso);
  }
  qcdTree->SetBranchAddress("leadPhotonEt", &leadGammaEt);
  qcdTree->SetBranchAddress("leadPhotonEta", &leadGammaEta);
  qcdTree->SetBranchAddress("trailPhotonEt", &trailGammaEt);
  qcdTree->SetBranchAddress("trailPhotonEta", &trailGammaEta);

  for(unsigned int i = 0; i < mcTrees.size(); i++) {
    mcTrees[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees[i]->SetBranchAddress("pileupWeightErr", &puWeightErr);
    mcTrees[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcTrees[i]->SetBranchAddress("btagWeightErr", &btagWeightErr);
    mcTrees[i]->SetBranchAddress("btagWeightUp", &btagWeightUp);
    mcTrees[i]->SetBranchAddress("btagWeightDown", &btagWeightDown);
    mcTrees[i]->SetBranchAddress("pileupWeightUp", &puWeightUp);
    mcTrees[i]->SetBranchAddress("pileupWeightDown", &puWeightDown);
    mcTrees[i]->SetBranchAddress(variable, &value);
    mcTrees[i]->SetBranchAddress("Ngamma", &Ngamma);
    mcTrees[i]->SetBranchAddress("Nfake", &Nfake);
    mcTrees[i]->SetBranchAddress("Nphotons", &Nphotons);
    if(variable != "pfMET_t01") mcTrees[i]->SetBranchAddress("pfMET_t01", &met);
    if(channel.Contains("ele")) {
      mcTrees[i]->SetBranchAddress("ele_pt", &leptonPt);
      mcTrees[i]->SetBranchAddress("ele_eta", &leptonEta);
    }
    else {
      mcTrees[i]->SetBranchAddress("muon_pt", &leptonPt);
      mcTrees[i]->SetBranchAddress("muon_eta", &leptonEta);
    }
    mcTrees[i]->SetBranchAddress("leadPhotonEt", &leadGammaEt);
    mcTrees[i]->SetBranchAddress("leadPhotonEta", &leadGammaEta);
    mcTrees[i]->SetBranchAddress("trailPhotonEt", &trailGammaEt);
    mcTrees[i]->SetBranchAddress("trailPhotonEta", &trailGammaEta);
    mcTrees[i]->SetBranchAddress("leadMatchGamma", &leadMatchGamma);
    mcTrees[i]->SetBranchAddress("leadMatchElectron", &leadMatchElectron);
    mcTrees[i]->SetBranchAddress("leadMatchJet", &leadMatchJet);

    mcTrees_JECup[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees_JECup[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcTrees_JECup[i]->SetBranchAddress(variable, &value);
    mcTrees_JECup[i]->SetBranchAddress("Ngamma", &Ngamma);
    mcTrees_JECup[i]->SetBranchAddress("Nfake", &Nfake);
    mcTrees_JECup[i]->SetBranchAddress("Nphotons", &Nphotons);
    if(variable != "pfMET_t01") mcTrees_JECup[i]->SetBranchAddress("pfMET_t01", &met);
    if(channel.Contains("ele")) {
      mcTrees_JECup[i]->SetBranchAddress("ele_pt", &leptonPt);
      mcTrees_JECup[i]->SetBranchAddress("ele_eta", &leptonEta);
    }
    else {
      mcTrees_JECup[i]->SetBranchAddress("muon_pt", &leptonPt);
      mcTrees_JECup[i]->SetBranchAddress("muon_eta", &leptonEta);
    }
    mcTrees_JECup[i]->SetBranchAddress("leadPhotonEt", &leadGammaEt);
    mcTrees_JECup[i]->SetBranchAddress("leadPhotonEta", &leadGammaEta);
    mcTrees_JECup[i]->SetBranchAddress("trailPhotonEt", &trailGammaEt);
    mcTrees_JECup[i]->SetBranchAddress("trailPhotonEta", &trailGammaEta);
    mcTrees_JECup[i]->SetBranchAddress("leadMatchGamma", &leadMatchGamma);
    mcTrees_JECup[i]->SetBranchAddress("leadMatchElectron", &leadMatchElectron);
    mcTrees_JECup[i]->SetBranchAddress("leadMatchJet", &leadMatchJet);

    mcTrees_JECdown[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcTrees_JECdown[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcTrees_JECdown[i]->SetBranchAddress(variable, &value);
    mcTrees_JECdown[i]->SetBranchAddress("Ngamma", &Ngamma);
    mcTrees_JECdown[i]->SetBranchAddress("Nfake", &Nfake);
    mcTrees_JECdown[i]->SetBranchAddress("Nphotons", &Nphotons);
    if(variable != "pfMET_t01") mcTrees_JECdown[i]->SetBranchAddress("pfMET_t01", &met);
    if(channel.Contains("ele")) {
      mcTrees_JECdown[i]->SetBranchAddress("ele_pt", &leptonPt);
      mcTrees_JECdown[i]->SetBranchAddress("ele_eta", &leptonEta);
    }
    else {
      mcTrees_JECdown[i]->SetBranchAddress("muon_pt", &leptonPt);
      mcTrees_JECdown[i]->SetBranchAddress("muon_eta", &leptonEta);
    }
    mcTrees_JECdown[i]->SetBranchAddress("leadPhotonEt", &leadGammaEt);
    mcTrees_JECdown[i]->SetBranchAddress("leadPhotonEta", &leadGammaEta);
    mcTrees_JECdown[i]->SetBranchAddress("trailPhotonEt", &trailGammaEt);
    mcTrees_JECdown[i]->SetBranchAddress("trailPhotonEta", &trailGammaEta);
    mcTrees_JECdown[i]->SetBranchAddress("leadMatchGamma", &leadMatchGamma);
    mcTrees_JECdown[i]->SetBranchAddress("leadMatchElectron", &leadMatchElectron);
    mcTrees_JECdown[i]->SetBranchAddress("leadMatchJet", &leadMatchJet);

    mcQCDTrees[i]->SetBranchAddress("pileupWeight", &puWeight);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightErr", &puWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeight", &btagWeight);
    mcQCDTrees[i]->SetBranchAddress("btagWeightErr", &btagWeightErr);
    mcQCDTrees[i]->SetBranchAddress("btagWeightUp", &btagWeightUp);
    mcQCDTrees[i]->SetBranchAddress("btagWeightDown", &btagWeightDown);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightUp", &puWeightUp);
    mcQCDTrees[i]->SetBranchAddress("pileupWeightDown", &puWeightDown);
    mcQCDTrees[i]->SetBranchAddress(variable, &value);
    mcQCDTrees[i]->SetBranchAddress("Ngamma", &Ngamma);
    mcQCDTrees[i]->SetBranchAddress("Nfake", &Nfake);
    mcQCDTrees[i]->SetBranchAddress("Nphotons", &Nphotons);
    if(variable != "pfMET_t01") mcQCDTrees[i]->SetBranchAddress("pfMET_t01", &met);
    if(channel.Contains("ele")) {
      mcQCDTrees[i]->SetBranchAddress("ele_pt", &leptonPt);
      mcQCDTrees[i]->SetBranchAddress("ele_eta", &leptonEta);
      mcQCDTrees[i]->SetBranchAddress("ele_relIso", &relIso);
    }
    else {
      mcQCDTrees[i]->SetBranchAddress("muon_pt", &leptonPt);
      mcQCDTrees[i]->SetBranchAddress("muon_eta", &leptonEta);
      mcQCDTrees[i]->SetBranchAddress("muon_relIso", &relIso);
    }
    mcQCDTrees[i]->SetBranchAddress("leadPhotonEt", &leadGammaEt);
    mcQCDTrees[i]->SetBranchAddress("leadPhotonEta", &leadGammaEta);
    mcQCDTrees[i]->SetBranchAddress("trailPhotonEt", &trailGammaEt);
    mcQCDTrees[i]->SetBranchAddress("trailPhotonEta", &trailGammaEta);

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

    mcTrees[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    mcTrees_JECup[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    mcTrees_JECdown[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
    mcQCDTrees[i]->SetBranchAddress("TopPtReweighting", &topPtReweighting);
  }

}

bool TemplateMaker::LoadMCBackground(TString fileName, TString scanName,
				     Double_t xsec, Double_t scaleErrorUp, Double_t scaleErrorDown, Double_t pdfErrorUp, Double_t pdfErrorDown,
				     bool remove_whizard, bool remove_madgraph, bool reweightTop) {

  mcFiles.push_back(new TFile(fileName, "READ"));
  if(!mcFiles.back()) {
    cout << "Could not load TFile " << fileName << endl;
    return false;
  }

  TString signalString, qcdString;
  if(controlRegion == kSR1 || controlRegion == kSR2 || controlRegion == kCR0 || controlRegion == kAny) {
    signalString = channel+"_signalTree";
    qcdString = (channel.Contains("ele")) ? channel+"_eQCDTree" : channel+"_muQCDTree";
  }
  else {
    signalString = channel+"_fakeTree";
    qcdString = (channel.Contains("ele")) ? channel+"_eQCDfakeTree" : channel+"_muQCDfakeTree";
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

  return true;
}

void TemplateMaker::BookTemplates() {
  
  TString suffix = channel;
  if(metCut > 0.) suffix += "_metCut_" + Form("%.1f", metCut);

  h_gg = new TH1D(variable+"_gg_"+suffix, variable, nBins, xlo, xhi);
  h_gg->Sumw2();

  h_qcd = new TH1D(variable+"_qcd_"+suffix, variable, nBins, xlo, xhi);
  h_qcd->Sumw2();

  h_qcd_relIso_10 = (TH1D*)h_qcd->Clone(variable+"_qcd_relIso_10_"+suffix);
  h_qcd_relIso_m10 = (TH1D*)h_qcd->Clone(variable+"_qcd_relIso_m10_"+suffix);

  for(unsigned int i = 0; i < mcTrees.size(); i++) {
  TH1D * h_bkg = new TH1D(variable+"_"+mcNames[i]+"_"+suffix, variable, nBins, xlo, xhi);
  h_bkg->Sumw2();
  mcHistograms.push_back(h_bkg);
  
  TH1D * h_bkg_btagWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_btagWeightUp");
  mcHistograms_btagWeightUp.push_back(h_bkg_btagWeightUp);
  
  TH1D * h_bkg_btagWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_btagWeightDown");
  mcHistograms_btagWeightDown.push_back(h_bkg_btagWeightDown);
  
  TH1D * h_bkg_puWeightUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_puWeightUp");
  mcHistograms_puWeightUp.push_back(h_bkg_puWeightUp);
  
  TH1D * h_bkg_puWeightDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_puWeightDown");
  mcHistograms_puWeightDown.push_back(h_bkg_puWeightDown);
  
  TH1D * h_bkg_scaleUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_scaleUp");
  mcHistograms_scaleUp.push_back(h_bkg_scaleUp);
  
  TH1D * h_bkg_scaleDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_scaleDown");
  mcHistograms_scaleDown.push_back(h_bkg_scaleDown);
  
  TH1D * h_bkg_pdfUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_pdfUp");
  mcHistograms_pdfUp.push_back(h_bkg_pdfUp);
  
  TH1D * h_bkg_pdfDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_pdfDown");
  mcHistograms_pdfDown.push_back(h_bkg_pdfDown);
  
  TH1D * h_bkg_topPtUp = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_topPtUp");
  mcHistograms_topPtUp.push_back(h_bkg_topPtUp);
  
  TH1D * h_bkg_topPtDown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_topPtDown");
  mcHistograms_topPtDown.push_back(h_bkg_topPtDown);
  
  TH1D * h_bkg_JECup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_JECUp");
  mcHistograms_JECup.push_back(h_bkg_JECup);
  
  TH1D * h_bkg_JECdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_JECDown");
  mcHistograms_JECdown.push_back(h_bkg_JECdown);
  
  TH1D * h_bkg_leptonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_leptonSFUp");
  mcHistograms_leptonSFup.push_back(h_bkg_leptonSFup);
  
  TH1D * h_bkg_leptonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_leptonSFDown");
  mcHistograms_leptonSFdown.push_back(h_bkg_leptonSFdown);
  
  TH1D * h_bkg_photonSFup = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_photonSFUp");
  mcHistograms_photonSFup.push_back(h_bkg_photonSFup);
  
  TH1D * h_bkg_photonSFdown = (TH1D*)h_bkg->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_photonSFDown");
  mcHistograms_photonSFdown.push_back(h_bkg_photonSFdown);

  TH1D * h_bkg_matchPhoton = new TH1D(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton", variable, nBins, xlo, xhi);
  h_bkg_matchPhoton->Sumw2();
  mcHistograms_matchPhoton.push_back(h_bkg_matchPhoton);
  
  TH1D * h_bkg_matchPhoton_btagWeightUp = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_btagWeightUp");
  mcHistograms_matchPhoton_btagWeightUp.push_back(h_bkg_matchPhoton_btagWeightUp);
  
  TH1D * h_bkg_matchPhoton_btagWeightDown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_btagWeightDown");
  mcHistograms_matchPhoton_btagWeightDown.push_back(h_bkg_matchPhoton_btagWeightDown);
  
  TH1D * h_bkg_matchPhoton_puWeightUp = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_puWeightUp");
  mcHistograms_matchPhoton_puWeightUp.push_back(h_bkg_matchPhoton_puWeightUp);
  
  TH1D * h_bkg_matchPhoton_puWeightDown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_puWeightDown");
  mcHistograms_matchPhoton_puWeightDown.push_back(h_bkg_matchPhoton_puWeightDown);
  
  TH1D * h_bkg_matchPhoton_scaleUp = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_scaleUp");
  mcHistograms_matchPhoton_scaleUp.push_back(h_bkg_matchPhoton_scaleUp);
  
  TH1D * h_bkg_matchPhoton_scaleDown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_scaleDown");
  mcHistograms_matchPhoton_scaleDown.push_back(h_bkg_matchPhoton_scaleDown);
  
  TH1D * h_bkg_matchPhoton_pdfUp = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_pdfUp");
  mcHistograms_matchPhoton_pdfUp.push_back(h_bkg_matchPhoton_pdfUp);
  
  TH1D * h_bkg_matchPhoton_pdfDown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_pdfDown");
  mcHistograms_matchPhoton_pdfDown.push_back(h_bkg_matchPhoton_pdfDown);
  
  TH1D * h_bkg_matchPhoton_topPtUp = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_topPtUp");
  mcHistograms_matchPhoton_topPtUp.push_back(h_bkg_matchPhoton_topPtUp);
  
  TH1D * h_bkg_matchPhoton_topPtDown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_topPtDown");
  mcHistograms_matchPhoton_topPtDown.push_back(h_bkg_matchPhoton_topPtDown);
  
  TH1D * h_bkg_matchPhoton_JECup = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_JECUp");
  mcHistograms_matchPhoton_JECup.push_back(h_bkg_matchPhoton_JECup);
  
  TH1D * h_bkg_matchPhoton_JECdown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_JECDown");
  mcHistograms_matchPhoton_JECdown.push_back(h_bkg_matchPhoton_JECdown);
  
  TH1D * h_bkg_matchPhoton_leptonSFup = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_leptonSFUp");
  mcHistograms_matchPhoton_leptonSFup.push_back(h_bkg_matchPhoton_leptonSFup);
  
  TH1D * h_bkg_matchPhoton_leptonSFdown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_leptonSFDown");
  mcHistograms_matchPhoton_leptonSFdown.push_back(h_bkg_matchPhoton_leptonSFdown);
  
  TH1D * h_bkg_matchPhoton_photonSFup = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_photonSFUp");
  mcHistograms_matchPhoton_photonSFup.push_back(h_bkg_matchPhoton_photonSFup);
  
  TH1D * h_bkg_matchPhoton_photonSFdown = (TH1D*)h_bkg_matchPhoton->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchPhoton"+"_photonSFDown");
  mcHistograms_matchPhoton_photonSFdown.push_back(h_bkg_matchPhoton_photonSFdown);

  TH1D * h_bkg_matchElectron = new TH1D(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron", variable, nBins, xlo, xhi);
  h_bkg_matchElectron->Sumw2();
  mcHistograms_matchElectron.push_back(h_bkg_matchElectron);
  
  TH1D * h_bkg_matchElectron_btagWeightUp = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_btagWeightUp");
  mcHistograms_matchElectron_btagWeightUp.push_back(h_bkg_matchElectron_btagWeightUp);
  
  TH1D * h_bkg_matchElectron_btagWeightDown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_btagWeightDown");
  mcHistograms_matchElectron_btagWeightDown.push_back(h_bkg_matchElectron_btagWeightDown);
  
  TH1D * h_bkg_matchElectron_puWeightUp = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_puWeightUp");
  mcHistograms_matchElectron_puWeightUp.push_back(h_bkg_matchElectron_puWeightUp);
  
  TH1D * h_bkg_matchElectron_puWeightDown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_puWeightDown");
  mcHistograms_matchElectron_puWeightDown.push_back(h_bkg_matchElectron_puWeightDown);
  
  TH1D * h_bkg_matchElectron_scaleUp = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_scaleUp");
  mcHistograms_matchElectron_scaleUp.push_back(h_bkg_matchElectron_scaleUp);
  
  TH1D * h_bkg_matchElectron_scaleDown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_scaleDown");
  mcHistograms_matchElectron_scaleDown.push_back(h_bkg_matchElectron_scaleDown);
  
  TH1D * h_bkg_matchElectron_pdfUp = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_pdfUp");
  mcHistograms_matchElectron_pdfUp.push_back(h_bkg_matchElectron_pdfUp);
  
  TH1D * h_bkg_matchElectron_pdfDown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_pdfDown");
  mcHistograms_matchElectron_pdfDown.push_back(h_bkg_matchElectron_pdfDown);
  
  TH1D * h_bkg_matchElectron_topPtUp = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_topPtUp");
  mcHistograms_matchElectron_topPtUp.push_back(h_bkg_matchElectron_topPtUp);
  
  TH1D * h_bkg_matchElectron_topPtDown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_topPtDown");
  mcHistograms_matchElectron_topPtDown.push_back(h_bkg_matchElectron_topPtDown);
  
  TH1D * h_bkg_matchElectron_JECup = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_JECUp");
  mcHistograms_matchElectron_JECup.push_back(h_bkg_matchElectron_JECup);
  
  TH1D * h_bkg_matchElectron_JECdown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_JECDown");
  mcHistograms_matchElectron_JECdown.push_back(h_bkg_matchElectron_JECdown);
  
  TH1D * h_bkg_matchElectron_leptonSFup = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_leptonSFUp");
  mcHistograms_matchElectron_leptonSFup.push_back(h_bkg_matchElectron_leptonSFup);
  
  TH1D * h_bkg_matchElectron_leptonSFdown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_leptonSFDown");
  mcHistograms_matchElectron_leptonSFdown.push_back(h_bkg_matchElectron_leptonSFdown);
  
  TH1D * h_bkg_matchElectron_photonSFup = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_photonSFUp");
  mcHistograms_matchElectron_photonSFup.push_back(h_bkg_matchElectron_photonSFup);
  
  TH1D * h_bkg_matchElectron_photonSFdown = (TH1D*)h_bkg_matchElectron->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchElectron"+"_photonSFDown");
  mcHistograms_matchElectron_photonSFdown.push_back(h_bkg_matchElectron_photonSFdown);

  TH1D * h_bkg_matchJet = new TH1D(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet", variable, nBins, xlo, xhi);
  h_bkg_matchJet->Sumw2();
  mcHistograms_matchJet.push_back(h_bkg_matchJet);
  
  TH1D * h_bkg_matchJet_btagWeightUp = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_btagWeightUp");
  mcHistograms_matchJet_btagWeightUp.push_back(h_bkg_matchJet_btagWeightUp);
  
  TH1D * h_bkg_matchJet_btagWeightDown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_btagWeightDown");
  mcHistograms_matchJet_btagWeightDown.push_back(h_bkg_matchJet_btagWeightDown);
  
  TH1D * h_bkg_matchJet_puWeightUp = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_puWeightUp");
  mcHistograms_matchJet_puWeightUp.push_back(h_bkg_matchJet_puWeightUp);
  
  TH1D * h_bkg_matchJet_puWeightDown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_puWeightDown");
  mcHistograms_matchJet_puWeightDown.push_back(h_bkg_matchJet_puWeightDown);
  
  TH1D * h_bkg_matchJet_scaleUp = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_scaleUp");
  mcHistograms_matchJet_scaleUp.push_back(h_bkg_matchJet_scaleUp);
  
  TH1D * h_bkg_matchJet_scaleDown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_scaleDown");
  mcHistograms_matchJet_scaleDown.push_back(h_bkg_matchJet_scaleDown);
  
  TH1D * h_bkg_matchJet_pdfUp = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_pdfUp");
  mcHistograms_matchJet_pdfUp.push_back(h_bkg_matchJet_pdfUp);
  
  TH1D * h_bkg_matchJet_pdfDown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_pdfDown");
  mcHistograms_matchJet_pdfDown.push_back(h_bkg_matchJet_pdfDown);
  
  TH1D * h_bkg_matchJet_topPtUp = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_topPtUp");
  mcHistograms_matchJet_topPtUp.push_back(h_bkg_matchJet_topPtUp);
  
  TH1D * h_bkg_matchJet_topPtDown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_topPtDown");
  mcHistograms_matchJet_topPtDown.push_back(h_bkg_matchJet_topPtDown);
  
  TH1D * h_bkg_matchJet_JECup = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_JECUp");
  mcHistograms_matchJet_JECup.push_back(h_bkg_matchJet_JECup);
  
  TH1D * h_bkg_matchJet_JECdown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_JECDown");
  mcHistograms_matchJet_JECdown.push_back(h_bkg_matchJet_JECdown);
  
  TH1D * h_bkg_matchJet_leptonSFup = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_leptonSFUp");
  mcHistograms_matchJet_leptonSFup.push_back(h_bkg_matchJet_leptonSFup);
  
  TH1D * h_bkg_matchJet_leptonSFdown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_leptonSFDown");
  mcHistograms_matchJet_leptonSFdown.push_back(h_bkg_matchJet_leptonSFdown);
  
  TH1D * h_bkg_matchJet_photonSFup = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_photonSFUp");
  mcHistograms_matchJet_photonSFup.push_back(h_bkg_matchJet_photonSFup);
  
  TH1D * h_bkg_matchJet_photonSFdown = (TH1D*)h_bkg_matchJet->Clone(variable+"_"+mcNames[i]+"_"+suffix+"_matchJet"+"_photonSFDown");
  mcHistograms_matchJet_photonSFdown.push_back(h_bkg_matchJet_photonSFdown);

  h_bkg = new TH1D(variable+"_qcd_"+mcNames[i]+"_"+suffix, variable, nBins, xlo, xhi);
  h_bkg->Sumw2();
  mcQCDHistograms.push_back(h_bkg);
  
  mcQCDHistograms_relIso_10.push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_10_"+mcNames[i]+"_"+suffix));
  mcQCDHistograms_relIso_m10.push_back((TH1D*)h_bkg->Clone(variable+"_qcd_relIso_m10_"+mcNames[i]+"_"+suffix));
  }

}

void TemplateMaker::FillData() {
  
  for(int i = 0; i < ggTree->GetEntries(); i++) {
    ggTree->GetEntry(i);

    if(!inControlRegion()) continue;

    if(variable == "pfMET_t01") met = value;
    if(metCut > 0. && met < metCut) continue;
    
    h_gg->Fill(value);

  }

  ggTree->ResetBranchAddresses();

}
  
void TemplateMaker::FillQCD() {

  for(int i = 0; i < qcdTree->GetEntries(); i++) {
    qcdTree->GetEntry(i);
    
    if(!inControlRegion()) continue;

    if(variable == "pfMET_t01") met = value;
    if(metCut > 0. && met < metCut) continue;

    h_qcd->Fill(value);

    if(relIso > 0.25 * 1.1) h_qcd_relIso_10->Fill(value);
    if(relIso > 0.25 * 0.9) h_qcd_relIso_m10->Fill(value);

  }

  qcdTree->ResetBranchAddresses();

}

void TemplateMaker::FillMCBackgrounds() {
 
  Float_t leptonSF, leptonSFup, leptonSFdown;
  Float_t photonSF, photonSFup, photonSFdown;

  for(unsigned int i = 0; i < mcTrees.size(); i++) {
    
    cout << " -- Starting on " << mcNames[i] << endl;

    for(int j = 0; j < mcTrees[i]->GetEntries(); j++) {
      mcTrees[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!inControlRegion()) continue;
      
      if(variable == "pfMET_t01") met = value;
      if(metCut > 0. && met < metCut) continue;

      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;

      Float_t addError2 = puWeightErr*puWeightErr/puWeight/puWeight + btagWeightErr*btagWeightErr/btagWeight/btagWeight;
      addError2 *= totalWeight*totalWeight;
      FillWithError(mcHistograms[i], value, totalWeight, addError2);
      if(leadMatchGamma == 1) FillWithError(mcHistograms_matchPhoton[i], value, totalWeight, addError2);
      if(leadMatchElectron == 1) FillWithError(mcHistograms_matchElectron[i], value, totalWeight, addError2);
      if(leadMatchJet == 1) FillWithError(mcHistograms_matchJet[i], value, totalWeight, addError2);

      totalWeight = puWeight * btagWeightUp * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_btagWeightUp[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_btagWeightUp[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_btagWeightUp[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_btagWeightUp[i]->Fill(value, totalWeight);

      totalWeight = puWeight * btagWeightDown * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_btagWeightDown[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_btagWeightDown[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_btagWeightDown[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_btagWeightDown[i]->Fill(value, totalWeight);

      totalWeight = puWeightUp * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_puWeightUp[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_puWeightUp[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_puWeightUp[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_puWeightUp[i]->Fill(value, totalWeight);

      totalWeight = puWeightDown * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_puWeightDown[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_puWeightDown[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_puWeightDown[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_puWeightDown[i]->Fill(value, totalWeight);

      totalWeight = puWeight * btagWeight * leptonSFup * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_leptonSFup[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_leptonSFup[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_leptonSFup[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_leptonSFup[i]->Fill(value, totalWeight);

      totalWeight = puWeight * btagWeight * leptonSFdown * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_leptonSFdown[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_leptonSFdown[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_leptonSFdown[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_leptonSFdown[i]->Fill(value, totalWeight);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFup;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_photonSFup[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_photonSFup[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_photonSFup[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_photonSFup[i]->Fill(value, totalWeight);

      totalWeight = puWeight * btagWeight * leptonSF * photonSFdown;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      mcHistograms_photonSFdown[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_photonSFdown[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_photonSFdown[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_photonSFdown[i]->Fill(value, totalWeight);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting * topPtReweighting;
      mcHistograms_topPtUp[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_topPtUp[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_topPtUp[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_topPtUp[i]->Fill(value, totalWeight);

      totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      mcHistograms_topPtDown[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_topPtDown[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_topPtDown[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_topPtDown[i]->Fill(value, totalWeight);

    }

  }

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    for(int k = 0; k < mcHistograms[i]->GetNbinsX(); k++) {
      Double_t content = mcHistograms[i]->GetBinContent(k+1);
      Double_t error = mcHistograms[i]->GetBinError(k+1);
      
      mcHistograms_scaleUp[i]->SetBinContent(k+1, content);
      mcHistograms_scaleDown[i]->SetBinContent(k+1, content);
      mcHistograms_pdfUp[i]->SetBinContent(k+1, content);
      mcHistograms_pdfDown[i]->SetBinContent(k+1, content);
      
      mcHistograms_scaleUp[i]->SetBinError(k+1, error);
      mcHistograms_scaleDown[i]->SetBinError(k+1, error);
      mcHistograms_pdfUp[i]->SetBinError(k+1, error);
      mcHistograms_pdfDown[i]->SetBinError(k+1, error);
    }

    mcHistograms[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_btagWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_btagWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_puWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_puWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_scaleUp[i]->Scale(intLumi_int * (crossSections[i] + scaleErrUp[i]) / mcNGen[i]);
    mcHistograms_scaleDown[i]->Scale(intLumi_int * (crossSections[i] - scaleErrDown[i]) / mcNGen[i]);
    mcHistograms_pdfUp[i]->Scale(intLumi_int * (crossSections[i] + pdfErrUp[i]) / mcNGen[i]);
    mcHistograms_pdfDown[i]->Scale(intLumi_int * (crossSections[i] - pdfErrDown[i]) / mcNGen[i]);
    mcHistograms_topPtUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_topPtDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_leptonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_leptonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_photonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_photonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
  }

  for(unsigned int i = 0; i < mcHistograms_matchPhoton.size(); i++) {
    for(int k = 0; k < mcHistograms_matchPhoton[i]->GetNbinsX(); k++) {
      Double_t content = mcHistograms_matchPhoton[i]->GetBinContent(k+1);
      Double_t error = mcHistograms_matchPhoton[i]->GetBinError(k+1);
      
      mcHistograms_matchPhoton_scaleUp[i]->SetBinContent(k+1, content);
      mcHistograms_matchPhoton_scaleDown[i]->SetBinContent(k+1, content);
      mcHistograms_matchPhoton_pdfUp[i]->SetBinContent(k+1, content);
      mcHistograms_matchPhoton_pdfDown[i]->SetBinContent(k+1, content);
      
      mcHistograms_matchPhoton_scaleUp[i]->SetBinError(k+1, error);
      mcHistograms_matchPhoton_scaleDown[i]->SetBinError(k+1, error);
      mcHistograms_matchPhoton_pdfUp[i]->SetBinError(k+1, error);
      mcHistograms_matchPhoton_pdfDown[i]->SetBinError(k+1, error);
    }

    mcHistograms_matchPhoton[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_btagWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_btagWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_puWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_puWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_scaleUp[i]->Scale(intLumi_int * (crossSections[i] + scaleErrUp[i]) / mcNGen[i]);
    mcHistograms_matchPhoton_scaleDown[i]->Scale(intLumi_int * (crossSections[i] - scaleErrDown[i]) / mcNGen[i]);
    mcHistograms_matchPhoton_pdfUp[i]->Scale(intLumi_int * (crossSections[i] + pdfErrUp[i]) / mcNGen[i]);
    mcHistograms_matchPhoton_pdfDown[i]->Scale(intLumi_int * (crossSections[i] - pdfErrDown[i]) / mcNGen[i]);
    mcHistograms_matchPhoton_topPtUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_topPtDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_leptonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_leptonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_photonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_photonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
  }

  for(unsigned int i = 0; i < mcHistograms_matchElectron.size(); i++) {
    for(int k = 0; k < mcHistograms_matchElectron[i]->GetNbinsX(); k++) {
      Double_t content = mcHistograms_matchElectron[i]->GetBinContent(k+1);
      Double_t error = mcHistograms_matchElectron[i]->GetBinError(k+1);
      
      mcHistograms_matchElectron_scaleUp[i]->SetBinContent(k+1, content);
      mcHistograms_matchElectron_scaleDown[i]->SetBinContent(k+1, content);
      mcHistograms_matchElectron_pdfUp[i]->SetBinContent(k+1, content);
      mcHistograms_matchElectron_pdfDown[i]->SetBinContent(k+1, content);
      
      mcHistograms_matchElectron_scaleUp[i]->SetBinError(k+1, error);
      mcHistograms_matchElectron_scaleDown[i]->SetBinError(k+1, error);
      mcHistograms_matchElectron_pdfUp[i]->SetBinError(k+1, error);
      mcHistograms_matchElectron_pdfDown[i]->SetBinError(k+1, error);
    }

    mcHistograms_matchElectron[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_btagWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_btagWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_puWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_puWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_scaleUp[i]->Scale(intLumi_int * (crossSections[i] + scaleErrUp[i]) / mcNGen[i]);
    mcHistograms_matchElectron_scaleDown[i]->Scale(intLumi_int * (crossSections[i] - scaleErrDown[i]) / mcNGen[i]);
    mcHistograms_matchElectron_pdfUp[i]->Scale(intLumi_int * (crossSections[i] + pdfErrUp[i]) / mcNGen[i]);
    mcHistograms_matchElectron_pdfDown[i]->Scale(intLumi_int * (crossSections[i] - pdfErrDown[i]) / mcNGen[i]);
    mcHistograms_matchElectron_topPtUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_topPtDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_leptonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_leptonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_photonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_photonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
  }

  for(unsigned int i = 0; i < mcHistograms_matchJet.size(); i++) {
    for(int k = 0; k < mcHistograms_matchJet[i]->GetNbinsX(); k++) {
      Double_t content = mcHistograms_matchJet[i]->GetBinContent(k+1);
      Double_t error = mcHistograms_matchJet[i]->GetBinError(k+1);
      
      mcHistograms_matchJet_scaleUp[i]->SetBinContent(k+1, content);
      mcHistograms_matchJet_scaleDown[i]->SetBinContent(k+1, content);
      mcHistograms_matchJet_pdfUp[i]->SetBinContent(k+1, content);
      mcHistograms_matchJet_pdfDown[i]->SetBinContent(k+1, content);
      
      mcHistograms_matchJet_scaleUp[i]->SetBinError(k+1, error);
      mcHistograms_matchJet_scaleDown[i]->SetBinError(k+1, error);
      mcHistograms_matchJet_pdfUp[i]->SetBinError(k+1, error);
      mcHistograms_matchJet_pdfDown[i]->SetBinError(k+1, error);
    }

    mcHistograms_matchJet[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_btagWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_btagWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_puWeightUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_puWeightDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_scaleUp[i]->Scale(intLumi_int * (crossSections[i] + scaleErrUp[i]) / mcNGen[i]);
    mcHistograms_matchJet_scaleDown[i]->Scale(intLumi_int * (crossSections[i] - scaleErrDown[i]) / mcNGen[i]);
    mcHistograms_matchJet_pdfUp[i]->Scale(intLumi_int * (crossSections[i] + pdfErrUp[i]) / mcNGen[i]);
    mcHistograms_matchJet_pdfDown[i]->Scale(intLumi_int * (crossSections[i] - pdfErrDown[i]) / mcNGen[i]);
    mcHistograms_matchJet_topPtUp[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_topPtDown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_leptonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_leptonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_photonSFup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_photonSFdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
  }
  
  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECup[i]->GetEntries(); j++) {
      mcTrees_JECup[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!inControlRegion()) continue;
      
      if(variable == "pfMET_t01") met = value;
      if(metCut > 0. && met < metCut) continue;

      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      
      mcHistograms_JECup[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_JECup[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_JECup[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_JECup[i]->Fill(value, totalWeight);
      
    }
    
    mcHistograms_JECup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_JECup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_JECup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_JECup[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    
  }
  
  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) {
    
    for(int j = 0; j < mcTrees_JECdown[i]->GetEntries(); j++) {
      mcTrees_JECdown[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!inControlRegion()) continue;
      
      if(variable == "pfMET_t01") met = value;
      if(metCut > 0. && met < metCut) continue;

      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;
      
      mcHistograms_JECdown[i]->Fill(value, totalWeight);
      if(leadMatchGamma == 1) mcHistograms_matchPhoton_JECdown[i]->Fill(value, totalWeight);
      if(leadMatchElectron == 1) mcHistograms_matchElectron_JECdown[i]->Fill(value, totalWeight);
      if(leadMatchJet == 1) mcHistograms_matchJet_JECdown[i]->Fill(value, totalWeight);
      
    }
    
    mcHistograms_JECdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchPhoton_JECdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchElectron_JECdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcHistograms_matchJet_JECdown[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
      
  }
  
  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) {

    for(int j = 0; j < mcQCDTrees[i]->GetEntries(); j++) {
      mcQCDTrees[i]->GetEntry(j);
      
      if(!checkBtagging()) continue;
      if(!inControlRegion()) continue;

      if(variable == "pfMET_t01") met = value;
      if(metCut > 0. && met < metCut) continue;

      if(removeWhizardOverlap[i] && overlaps_whizard > 0.001) continue;
      if(removeMadgraphOverlap[i] && overlaps_madgraph > 0.001) continue;
      if(reweightTopPt[i] && topPtReweighting < 0.) continue;

      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);
      
      GetLeptonSF(leptonSF, leptonSFup, leptonSFdown);
      GetPhotonSF(photonSF, photonSFup, photonSFdown);

      double totalWeight = puWeight * btagWeight * leptonSF * photonSF;
      if(reweightTopPt[i]) totalWeight *= topPtReweighting;

      Float_t addError2 = puWeightErr*puWeightErr/puWeight/puWeight + btagWeightErr*btagWeightErr/btagWeight/btagWeight;
      addError2 *= totalWeight*totalWeight;

      FillWithError(mcQCDHistograms[i], value, totalWeight, addError2);

      if(relIso > 0.2 * 1.1) mcQCDHistograms_relIso_10[i]->Fill(value, totalWeight);
      if(relIso > 0.2 * 0.9) mcQCDHistograms_relIso_m10[i]->Fill(value, totalWeight);
      
    }
    
    mcQCDHistograms[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);

    mcQCDHistograms_relIso_10[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    mcQCDHistograms_relIso_m10[i]->Scale(intLumi_int * crossSections[i] / mcNGen[i]);
    
  }

  for(unsigned int i = 0; i < mcTrees.size(); i++) mcTrees[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECup.size(); i++) mcTrees_JECup[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcTrees_JECdown.size(); i++) mcTrees_JECdown[i]->ResetBranchAddresses();
  for(unsigned int i = 0; i < mcQCDTrees.size(); i++) mcQCDTrees[i]->ResetBranchAddresses();

}

void TemplateMaker::SubtractMCFromQCD() {

  for(unsigned int i = 0; i < mcQCDHistograms.size(); i++) h_qcd->Add(mcQCDHistograms[i], -1.);

  for(unsigned int i = 0; i < mcQCDHistograms_relIso_10.size(); i++) h_qcd_relIso_10->Add(mcQCDHistograms_relIso_10[i], -1.);

  for(unsigned int i = 0; i < mcQCDHistograms_relIso_m10.size(); i++) h_qcd_relIso_m10->Add(mcQCDHistograms_relIso_m10[i], -1.);

  for(Int_t j = 0; j < h_qcd->GetNbinsX(); j++) {
    if(h_qcd->GetBinContent(j+1) < 0) {
      h_qcd->SetBinContent(j+1, 0.);
      h_qcd->SetBinError(j+1, 0.);
    }
  }

  for(Int_t j = 0; j < h_qcd_relIso_10->GetNbinsX(); j++) {
    if(h_qcd_relIso_10->GetBinContent(j+1) < 0) {
      h_qcd_relIso_10->SetBinContent(j+1, 0.);
      h_qcd_relIso_10->SetBinError(j+1, 0.);
    }
  }

  for(Int_t j = 0; j < h_qcd_relIso_m10->GetNbinsX(); j++) {
    if(h_qcd_relIso_m10->GetBinContent(j+1) < 0) {
      h_qcd_relIso_m10->SetBinContent(j+1, 0.);
      h_qcd_relIso_m10->SetBinError(j+1, 0.);
    }
  }

}

void TemplateMaker::SaveOutput() {

  TString outName = "fitTemplates.root";

  TFile * fOut = new TFile(outName, "UPDATE");

  h_gg->Write();
  h_qcd->Write();
  h_qcd_relIso_10->Write();
  h_qcd_relIso_m10->Write();

  for(unsigned int i = 0; i < mcHistograms.size(); i++) {
    mcHistograms[i]->Write();
    mcHistograms_btagWeightUp[i]->Write();
    mcHistograms_btagWeightDown[i]->Write();
    mcHistograms_puWeightUp[i]->Write();
    mcHistograms_puWeightDown[i]->Write();
    mcHistograms_scaleUp[i]->Write();
    mcHistograms_scaleDown[i]->Write();
    mcHistograms_pdfUp[i]->Write();
    mcHistograms_pdfDown[i]->Write();
    mcHistograms_topPtUp[i]->Write();
    mcHistograms_topPtDown[i]->Write();
    mcHistograms_JECup[i]->Write();
    mcHistograms_JECdown[i]->Write();
    mcHistograms_leptonSFup[i]->Write();
    mcHistograms_leptonSFdown[i]->Write();
    mcHistograms_photonSFup[i]->Write();
    mcHistograms_photonSFdown[i]->Write();

    mcHistograms_matchPhoton[i]->Write();
    mcHistograms_matchPhoton_btagWeightUp[i]->Write();
    mcHistograms_matchPhoton_btagWeightDown[i]->Write();
    mcHistograms_matchPhoton_puWeightUp[i]->Write();
    mcHistograms_matchPhoton_puWeightDown[i]->Write();
    mcHistograms_matchPhoton_scaleUp[i]->Write();
    mcHistograms_matchPhoton_scaleDown[i]->Write();
    mcHistograms_matchPhoton_pdfUp[i]->Write();
    mcHistograms_matchPhoton_pdfDown[i]->Write();
    mcHistograms_matchPhoton_topPtUp[i]->Write();
    mcHistograms_matchPhoton_topPtDown[i]->Write();
    mcHistograms_matchPhoton_JECup[i]->Write();
    mcHistograms_matchPhoton_JECdown[i]->Write();
    mcHistograms_matchPhoton_leptonSFup[i]->Write();
    mcHistograms_matchPhoton_leptonSFdown[i]->Write();
    mcHistograms_matchPhoton_photonSFup[i]->Write();
    mcHistograms_matchPhoton_photonSFdown[i]->Write();

    mcHistograms_matchElectron[i]->Write();
    mcHistograms_matchElectron_btagWeightUp[i]->Write();
    mcHistograms_matchElectron_btagWeightDown[i]->Write();
    mcHistograms_matchElectron_puWeightUp[i]->Write();
    mcHistograms_matchElectron_puWeightDown[i]->Write();
    mcHistograms_matchElectron_scaleUp[i]->Write();
    mcHistograms_matchElectron_scaleDown[i]->Write();
    mcHistograms_matchElectron_pdfUp[i]->Write();
    mcHistograms_matchElectron_pdfDown[i]->Write();
    mcHistograms_matchElectron_topPtUp[i]->Write();
    mcHistograms_matchElectron_topPtDown[i]->Write();
    mcHistograms_matchElectron_JECup[i]->Write();
    mcHistograms_matchElectron_JECdown[i]->Write();
    mcHistograms_matchElectron_leptonSFup[i]->Write();
    mcHistograms_matchElectron_leptonSFdown[i]->Write();
    mcHistograms_matchElectron_photonSFup[i]->Write();
    mcHistograms_matchElectron_photonSFdown[i]->Write();

    mcHistograms_matchJet[i]->Write();
    mcHistograms_matchJet_btagWeightUp[i]->Write();
    mcHistograms_matchJet_btagWeightDown[i]->Write();
    mcHistograms_matchJet_puWeightUp[i]->Write();
    mcHistograms_matchJet_puWeightDown[i]->Write();
    mcHistograms_matchJet_scaleUp[i]->Write();
    mcHistograms_matchJet_scaleDown[i]->Write();
    mcHistograms_matchJet_pdfUp[i]->Write();
    mcHistograms_matchJet_pdfDown[i]->Write();
    mcHistograms_matchJet_topPtUp[i]->Write();
    mcHistograms_matchJet_topPtDown[i]->Write();
    mcHistograms_matchJet_JECup[i]->Write();
    mcHistograms_matchJet_JECdown[i]->Write();
    mcHistograms_matchJet_leptonSFup[i]->Write();
    mcHistograms_matchJet_leptonSFdown[i]->Write();
    mcHistograms_matchJet_photonSFup[i]->Write();
    mcHistograms_matchJet_photonSFdown[i]->Write();
  }

  fOut->Close();

}

void TemplateMaker::GetLeptonSF(Float_t& central, Float_t& up, Float_t& down) {

  Float_t pt, eta, error;

  if(channel.Contains("ele")) {
    pt = min(leptonPt, (float)199.);
    pt = max(pt, (float)15.);
    eta = min(fabs(leptonEta), (double)2.39);

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
    pt = min(leptonPt, (float)499.);
    pt = max(pt, (float)10.);
    eta = min(fabs(leptonEta), (double)2.09);

    central = sf_muon->GetBinContent(sf_muon->FindBin(pt, eta));
    error = sf_muon->GetBinError(sf_muon->FindBin(pt, eta));

    up = error;
    down = 2. * central - error;
  }

  return;  
}

void TemplateMaker::GetPhotonSF(Float_t& central, Float_t& up, Float_t& down) {

  if(Nphotons == 0) {
    central = 1.;
    up = 1.;
    down = 1.;
    return;
  }

  Float_t et, eta;
  Float_t error = 0;

  if(Nphotons == 1) {
    et = min(leadGammaEt, (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(leadGammaEta), (double)1.4441);

    Float_t id_val = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    central = id_val * veto_val;
    error = central * sqrt(id_error*id_error/(id_val*id_val) + veto_error*veto_error/(veto_val*veto_val));
  }

  else if(Nphotons >= 2) {
    // lead photon
    et = min(leadGammaEt, (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(leadGammaEta), (double)1.4441);

    Float_t id_val_lead = sf_photon_id->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t id_error_lead = sf_photon_id->GetBinError(sf_photon_id->FindBin(et, eta));
    
    Float_t veto_val_lead = sf_photon_veto->GetBinContent(sf_photon_id->FindBin(et, eta));
    Float_t veto_error_lead = sf_photon_veto->GetBinError(sf_photon_id->FindBin(et, eta));

    // trail photon
    et = min(trailGammaEt, (float)999.);
    et = max(et, (float)15.);
    eta = min(fabs(trailGammaEta), (double)1.4441);

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
