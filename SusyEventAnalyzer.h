#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TArrayI.h>

#include <map>
#include <set>
#include <fstream>

#include "../src/SusyEvent.h"
#include "Muon.h"
#include "Electron.h"
#include "Photon.h"
#include "Jet.h"
#include "ScaleFactorInfo.h"
#include "BtagInfo.h"
#include "Metadata.h"

using namespace std;

template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}

bool CorrPtGreater(TLorentzVector p1, TLorentzVector p2) {
  return (p1.Pt() > p2.Pt());
}

class HLTInfo {
 public :
  HLTInfo() {
    paths.clear();
    types.clear();
  }
  HLTInfo(vector<TString> pathList, vector<int> eventType) {
    SetHLTPaths(pathList);
    SetEventTypes(eventType);
  }
  virtual ~HLTInfo() { ; }

  vector<TString> GetHLTNames() { return paths; }
  vector<int> GetEventTypes() { return types; }

  void SetHLTPaths(vector<TString> p) { 
    paths.clear();
    for(unsigned int i = 0; i < p.size(); i++) paths.push_back(p[i]);
  }
  void SetEventTypes(vector<int> t) { 
    types.clear();
    for(unsigned int i = 0; i < t.size(); i++) types.push_back(t[i]);
  }
  void AddEventTypes(int s) { types.push_back(s); }

 private :
  vector<TString> paths;
  vector<int> types;
};

class SusyEventAnalyzer {
 public:
  SusyEventAnalyzer(TTree&);
  virtual ~SusyEventAnalyzer();

  virtual void Data();
  virtual void Acceptance();
  virtual void GeneratorInfo();
  virtual void LeptonInfo();
  virtual void CalculateBtagEfficiency();
  virtual void PileupWeights(TString puFile, TString puFile_up, TString puFile_down);
  virtual void ZGammaData();
  virtual void ZGammaAcceptance();

  // utility functions
  float deltaR(TLorentzVector& p1, TLorentzVector& p2);
  float deltaR(TLorentzVector& p1, TVector3& p2);
  float deltaR(TVector3& p1, TVector3& p2);
  TVector3 FindJetVertex(susy::PFJet jet, vector<susy::Track> tracks);
  double d0correction(TVector3& beamSpot, susy::Track& track) const;
  double dZcorrection(TVector3& beamSpot, susy::Track& track) const;
  void AddHlt(vector<TString> v, vector<int> eventType) {
    hltInfos.push_back(HLTInfo(v, eventType));
  }

  // metadata and btag stuff
  void SetUseSyncFile(bool v) { useSyncFile = v; }
  void IncludeSyncFile(char* file);
  void SetCheckSingleEvent(bool v) { singleEvent = v; }
  void AddCheckSingleEvent(Int_t run, Int_t lumi, ULong_t eventnum) {
    single_run = run;
    single_lumi = lumi;
    single_event = eventnum;
  }
  void SetScanName(TString v) { scan = v; }
  void SetIsMC(bool v)        { isMC = v; }
  void SetIsFastSim(bool v) {
    if(v) {
      isMC = v;
      isFastSim = v;
    }
    else isFastSim = v;
  }
  void SetRejectFakeElectrons(bool v)  { rejectFakeElectrons = v;}
  void SetDoBtagScaling(bool v)        { doBtagScaling = v; }
  void SetBtagTechnicalStop(TString t) { btagTechnicalStop = t; }
  void SetDoPileupReweighting(bool v)  { doPileupReweighting = v; }
  void SetBtagger(TString t) { btagger = t; }
  void AddValidTagger(TString t) { validTags.push_back(t); }
  pair<unsigned int, Float_t> tagEncoder(TString tagger) {
    if(tagger == "CSVL") return make_pair(5, 0.244);
    if(tagger == "CSVM") return make_pair(5, 0.679);
    if(tagger == "CSVT") return make_pair(5, 0.898);
    return make_pair(0, 1.e6);
  }
  void SetUseTrigger(bool v) { useTrigger = v; }
  void SetUseJson(bool v)    { useJson = v; }

  void SetUseDPhiCut(bool v) { useDPhiCut = v; }

  void SetUseDeltaRCutsOnPhotons(bool v) { useDeltaRCutsOnPhotons = v; }

  void IncludeAJson(TString const&);
  void SetOutput(TString const& v) { outputName = v; }
  void SetPrintInterval(int v) { printInterval = v; }
  void SetPrintLevel(int v) { printLevel = v; }
  void SetProcessNEvents(int v) { processNEvents = v; }
  void CopyEvents(bool v) { copyEvents = v; }

  // major analysis logic
  void findMuons(susy::Event& ev, 
		 vector<susy::Muon*>& tightMuons, vector<susy::Muon*>& looseMuons, 
		 float& HT, int mode);
  void findElectrons(susy::Event& ev, 
		     vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons, 
		     vector<susy::Electron*>& tightEles, vector<susy::Electron*>& looseEles, 
		     float& HT, int mode);
  void findPhotons(susy::Event& ev, 
		   vector<susy::Photon*>& photons,
		   vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons,
		   vector<susy::Electron*> tightEles, vector<susy::Electron*> looseEles,
		   float& HT, 
		   int photonMode);

  // in data
  void findJets(susy::Event& ev, 
		vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons,
		vector<susy::Electron*> tightEles, vector<susy::Electron*> looseEles,
		vector<susy::Photon*> photons,
		vector<susy::PFJet*>& pfJets, vector<susy::PFJet*>& btags, 
		ScaleFactorInfo sf,
		vector<BtagInfo>& tagInfos, vector<float>& csvValues, 
		vector<TLorentzVector>& pfJets_corrP4, vector<TLorentzVector>& btags_corrP4, 
		float& HT, TLorentzVector& hadronicSystem);

  void findJets_inMC(susy::Event& ev, 
		     vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons,
		     vector<susy::Electron*> tightEles, vector<susy::Electron*> looseEles,
		     vector<susy::Photon*> photons,
		     vector<susy::PFJet*>& pfJets, vector<susy::PFJet*>& btags, 
		     ScaleFactorInfo sf,
		     vector<BtagInfo>& tagInfos, vector<float>& csvValues, 
		     vector<TLorentzVector>& pfJets_corrP4, vector<TLorentzVector>& btags_corrP4, 
		     float& HT, TLorentzVector& hadronicSystem,
		     int systematic);
     
  bool GetDiJetPt(susy::Event& ev, vector<susy::Photon*> candidates, float& diJetPt, float& leadpt, float& trailpt);
  bool PhotonMatchesElectron(susy::Event& ev, vector<susy::Photon*> candidates, int& bothMatchCounter);
  int FigureTTbarDecayMode();

  double TopPtReweighting(susy::Event& ev);
  double TopPtReweighting_ttHbb(susy::Event& ev);

  void fill_whizard_phaseSpace(TH2D*& h);
  void fill_madgraph_phaseSpace(TH2D*& h);
  bool overlaps_whizard();
  bool overlaps_madgraph();
  
  void SetTreeValues(map<TString, float>& treeMap,
		     vector<susy::Muon*> tightMuons, vector<susy::Electron*> tightEles, 
		     vector<susy::PFJet*> pfJets, vector<susy::PFJet*> btags,
		     vector<susy::Photon*> photons,
		     vector<TLorentzVector> pfJets_corrP4, vector<TLorentzVector> btags_corrP4,
		     vector<float> csvValues,
		     TLorentzVector hadronicSystem,
		     float HT, float HT_jets,
		     int nPVertex,
		     float eventWeight, float eventWeightErr, float eventWeight_up, float eventWeight_down,
		     Long64_t jentry);

  void SetTreeValues_ZGamma(map<TString, float>& treeMap,
			    vector<susy::Muon*> tightMuons, vector<susy::Electron*> tightEles, 
			    vector<susy::PFJet*> pfJets, vector<susy::PFJet*> btags,
			    vector<susy::Photon*> photons,
			    vector<TLorentzVector> pfJets_corrP4, vector<TLorentzVector> btags_corrP4,
			    vector<float> csvValues,
			    TLorentzVector hadronicSystem,
			    float HT, float HT_jets,
			    int nPVertex,
			    float eventWeight, float eventWeightErr, float eventWeight_up, float eventWeight_down,
			    Long64_t jentry);

  // lazy junk
  void FillMetFilter2D(susy::Event& ev, TH2F*& h);

 protected:
  bool IsGoodLumi(UInt_t, UInt_t) const;
  bool PassTrigger(TString path) const;
  bool PassTriggers(int eventType);

  susy::Event event;
  TTree *fTree;
  TString outputName;
  int printLevel;
  unsigned printInterval;
  int processNEvents;
  bool copyEvents;
  map<unsigned, set<unsigned> > goodLumiList;
  mutable pair<unsigned, unsigned> currentLumi;
  mutable bool currentLumiIsGood;

 private:
  vector<HLTInfo> hltInfos;

  vector<UInt_t> syncRuns;
  vector<UInt_t> syncLumi;
  vector<ULong_t> syncEvents;
  UInt_t single_run;
  UInt_t single_lumi;
  ULong_t single_event;
  bool useSyncFile;
  bool singleEvent;
  
  vector<TString> validTags;
  TString btagger;

  bool isMC;
  bool isFastSim;
  bool rejectFakeElectrons;
  bool doBtagScaling;
  TString btagTechnicalStop;
  bool doPileupReweighting;
  bool useTrigger;
  bool useJson;

  bool useDPhiCut;

  bool useDeltaRCutsOnPhotons;

  TString scan;

};

SusyEventAnalyzer::SusyEventAnalyzer(TTree& tree) :
  event(),
  fTree(&tree),
  outputName("analysis"),
  printLevel(0),
  printInterval(1000),
  processNEvents(-1),
  copyEvents(false),
  goodLumiList(),
  currentLumi(0, 0),
  currentLumiIsGood(true),
  rejectFakeElectrons(false),
  useDeltaRCutsOnPhotons(true)  
{
  event.setInput(tree);
}

SusyEventAnalyzer::~SusyEventAnalyzer() {}

void SusyEventAnalyzer::IncludeAJson(TString const& _fileName) {
  if(_fileName == "") return;

  ifstream inputFile(_fileName);
  if(!inputFile.is_open()){
    cerr << "Cannot open JSON file " << _fileName << endl;
    return;
  }

  string line;
  TString jsonText;
  while(true){
    getline(inputFile, line);
    if(!inputFile.good()) break;
    jsonText += line;
  }
  inputFile.close();

  TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
  TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

  TArrayI positions(2);
  positions[1] = 0;
  while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
    TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
    TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

    unsigned run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
    set<unsigned>& lumis(goodLumiList[run]);

    TArrayI lumiPos(2);
    lumiPos[1] = 0;
    while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
      TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
      int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
      int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
      for(int lumi(begin); lumi <= end; ++lumi)
        lumis.insert(lumi);
    }
  }
}

bool SusyEventAnalyzer::IsGoodLumi(UInt_t run, UInt_t lumi) const {
  if(goodLumiList.size() == 0) return true;
  if(run == currentLumi.first && lumi == currentLumi.second) return currentLumiIsGood;
  currentLumi.first = run;
  currentLumi.second = lumi;
  currentLumiIsGood = false;

  map<unsigned, set<unsigned> >::const_iterator rItr(goodLumiList.find(run));
  if(rItr != goodLumiList.end()){
    set<unsigned>::const_iterator lItr(rItr->second.find(lumi));
    if(lItr != rItr->second.end()) currentLumiIsGood = true;
  }

  return currentLumiIsGood;
}

bool SusyEventAnalyzer::PassTrigger(TString path) const {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event.hltMap.begin(); it != event.hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second))) {
      pass = true;
      break;
    }
  }
  return pass;
}

bool SusyEventAnalyzer::PassTriggers(int eventType) {

  // Types: gg = 1
  //	    eg = 2
  // 	    ee = 3
  // 	    ff = 4
  //        gf = 5

  bool pass = false;

  for(unsigned int i = 0; i < hltInfos.size(); i++) {

    bool correctType = false;
    for(unsigned int j = 0; j < ((hltInfos[i]).GetEventTypes()).size(); j++) {
      if(eventType == ((hltInfos[i]).GetEventTypes())[j]) {
	correctType = true;
	break;
      }
    }
    if(!correctType) continue;

    for(unsigned int j = 0; j < ((hltInfos[i]).GetHLTNames()).size(); j++) {
      if(PassTrigger( ((hltInfos[i]).GetHLTNames())[j]) ) {
	pass = true;
	break;
      }
    }

    if(pass == true) break;
  }

  return pass;
}

float SusyEventAnalyzer::deltaR(TLorentzVector& p1, TLorentzVector& p2) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

float SusyEventAnalyzer::deltaR(TLorentzVector& p1, TVector3& p2) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

float SusyEventAnalyzer::deltaR(TVector3& p1, TVector3& p2) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

TVector3 SusyEventAnalyzer::FindJetVertex(susy::PFJet jet, vector<susy::Track> tracks) {

  vector<TVector3> vertices;
  vector<float> sumPt;
  vector<int> nTracks;

  float epsilon = 1.e-6;

  for(unsigned int i = 0; i < jet.tracklist.size(); i++) {
    TVector3 vtx = tracks[jet.tracklist[i]].vertex;
    bool newVertex = true;
    unsigned int oldVertexIndex = -1;
    for(unsigned int j = 0; j < vertices.size(); j++) {
      if(abs(vtx.X() - vertices[j].X()) < epsilon &&
	 abs(vtx.Y() - vertices[j].Y()) < epsilon &&
	 abs(vtx.Z() - vertices[j].Z()) < epsilon) {
	newVertex = false;
	oldVertexIndex = j;
	break;
      }
    }

    if(newVertex) {
      vertices.push_back(vtx);
      sumPt.push_back(tracks[jet.tracklist[i]].momentum.Pt());
      nTracks.push_back(1);
    }
    else {
      sumPt[oldVertexIndex] += tracks[jet.tracklist[i]].momentum.Pt();
      nTracks[oldVertexIndex]++;
    }

  }

  int maxNtracks = -1;
  float maxSumPt = -1.;
  unsigned int matchByN = 0;
  unsigned int matchBySumPt = 0;

  for(unsigned int i = 0; i < vertices.size(); i++) {
    if(nTracks[i] > maxNtracks) {
      matchByN = i;
      maxNtracks = nTracks[i];
    }
    if(sumPt[i] > maxSumPt) {
      matchBySumPt = i;
      maxSumPt = sumPt[i];
    }
  }

  if(matchByN != matchBySumPt) cout << "Disagreeing nTracks and sumPt vertex-finding!" << endl;

  if(vertices.size() > 0) return vertices[matchByN];

  return TVector3(0.);

}

double SusyEventAnalyzer::d0correction(TVector3& beamspot, susy::Track& track) const {
  double d0 = track.d0() - beamspot.X()*sin(track.phi()) + beamspot.Y()*cos(track.phi());
  return d0;
}

double SusyEventAnalyzer::dZcorrection(TVector3& beamSpot, susy::Track& track) const {

  if(&(beamSpot) == 0x0 || &(track) == 0x0) {
    cout << endl << endl << "Something majorly broken in dZcorrection!!!" <<	endl <<	endl;
    return 1.e6;
  }

  if(track.momentum.Pt() == 0.) return 1.e6;

  double dz = (track.vertex.Z() - beamSpot.Z()) - ((track.vertex.X() - beamSpot.X())*track.momentum.Px() + (track.vertex.Y() - beamSpot.Y())*track.momentum.Py()) / track.momentum.Pt() * (track.momentum.Pz() / track.momentum.Pt());
  return dz;
}

void SusyEventAnalyzer::findPhotons(susy::Event& ev, 
				    vector<susy::Photon*>& photons,
				    vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons,
				    vector<susy::Electron*> tightEles, vector<susy::Electron*> looseEles,
				    float& HT, 
				    int photonMode) {
  
  map<TString, vector<susy::Photon> >::iterator phoMap = ev.photons.find("photons");
  if(phoMap != event.photons.end()) {
    for(vector<susy::Photon>::iterator it = phoMap->second.begin();
	it != phoMap->second.end(); it++) {
      
      // kSignalPhotons only looks for passing photons
      // kFakePhotons considers both passing and fake photons
      if((photonMode == kSignalPhotons && is_loosePhoton(*it, event.rho25)) ||
	 (photonMode == kFakePhotons && (is_loosePhoton(*it, event.rho25) || is_fakePhoton(*it, event.rho25)))) {

	bool overlap = false;

	float this_dR;

	if(useDeltaRCutsOnPhotons) { // default true
	  for(unsigned int i = 0; i < tightMuons.size(); i++) {
	    this_dR = deltaR(tightMuons[i]->momentum, it->momentum);
	    if(this_dR < 0.7) overlap = true;
	  }
	  
	  for(unsigned int i = 0; i < looseMuons.size(); i++) {
	    if(deltaR(looseMuons[i]->momentum, it->momentum) <= 0.7) overlap = true;
	  }
	  
	  for(unsigned int i = 0; i < tightEles.size(); i++) {
	    this_dR = deltaR(tightEles[i]->momentum, it->momentum);
	    if(this_dR < 0.7) overlap = true;
	  }
	  
	  for(unsigned int i = 0; i < looseEles.size(); i++) {
	    if(deltaR(looseEles[i]->momentum, it->momentum) <= 0.7) overlap = true;
	  }
	}
	  
	for(unsigned int i = 0; i < photons.size(); i++) {
	  this_dR = deltaR(photons[i]->momentum, it->momentum);
	  if(this_dR < 0.1) overlap = true;
	}
	
	if(overlap) continue;

	photons.push_back(&*it);
	HT += it->momentum.Pt();
	
      }
      
    } // for photon
  } // if

  sort(photons.begin(), photons.end(), EtGreater<susy::Photon>);
    
  return;

}

void SusyEventAnalyzer::findJets(susy::Event& ev,
				 vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons,
				 vector<susy::Electron*> tightEles, vector<susy::Electron*> looseEles,
				 vector<susy::Photon*> photons,
				 vector<susy::PFJet*>& pfJets, vector<susy::PFJet*>& btags, 
				 ScaleFactorInfo sf,
				 vector<BtagInfo>& tagInfos, vector<float>& csvValues, 
				 vector<TLorentzVector>& pfJets_corrP4, vector<TLorentzVector>& btags_corrP4, 
				 float& HT, TLorentzVector& hadronicSystem) {

  map<TString, susy::PFJetCollection>::iterator pfJets_it = ev.pfJets.find("ak5");
  if(pfJets_it != ev.pfJets.end()) {
    susy::PFJetCollection& jetColl = pfJets_it->second;
      
    for(vector<susy::PFJet>::iterator it = jetColl.begin();
	it != jetColl.end(); it++) {
	
      map<TString, Float_t>::iterator s_it = it->jecScaleFactors.find("L1FastL2L3");
      float scale = s_it->second;
	  
      TLorentzVector corrP4 = scale * it->momentum;

      if(!isGoodJet(*it, corrP4)) continue;
      if(JetOverlapsElectron(corrP4, tightEles, looseEles)) continue;
      if(JetOverlapsMuon(corrP4, tightMuons, looseMuons)) continue;
      if(useDeltaRCutsOnPhotons && JetsOverlapsPhoton(corrP4, photons)) continue;
      if(!useDeltaRCutsOnPhotons && btagger == "CSVM" && it->bTagDiscriminators[susy::kCSV] <= 0.679 && JetsOverlapsPhoton(corrP4, photons)) continue;

      pfJets.push_back(&*it);
      csvValues.push_back(it->bTagDiscriminators[susy::kCSV]);
      HT += corrP4.Pt();
      hadronicSystem += corrP4;
      pfJets_corrP4.push_back(corrP4);
	
      if(fabs(corrP4.Eta()) < 2.4) {
	if(isMC) {
	  BtagInfo info((*it), corrP4, btagger, 1., isMC, isFastSim, sf);
	  tagInfos.push_back(info);
	}

	if((btagger == "CSVL" && it->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && it->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && it->bTagDiscriminators[susy::kCSV] > 0.898)) {
	  btags.push_back(&*it);
	  btags_corrP4.push_back(corrP4);
	}
	 
      }
   
    } // loop over jet coll
  } // if the jet coll exists
  sort(pfJets.begin(), pfJets.end(), EtGreater<susy::PFJet>);
  sort(btags.begin(), btags.end(), EtGreater<susy::PFJet>);
  sort(csvValues.begin(), csvValues.end(), greater<float>());
  sort(pfJets_corrP4.begin(), pfJets_corrP4.end(), CorrPtGreater);
  sort(btags_corrP4.begin(), btags_corrP4.end(), CorrPtGreater);

}

void SusyEventAnalyzer::findJets_inMC(susy::Event& ev, 
				      vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons,
				      vector<susy::Electron*> tightEles, vector<susy::Electron*> looseEles,
				      vector<susy::Photon*> photons,
				      vector<susy::PFJet*>& pfJets, vector<susy::PFJet*>& btags, 
				      ScaleFactorInfo sf,
				      vector<BtagInfo>& tagInfos, vector<float>& csvValues, 
				      vector<TLorentzVector>& pfJets_corrP4, vector<TLorentzVector>& btags_corrP4, 
				      float& HT, TLorentzVector& hadronicSystem,
				      int systematic) {
  
  map<TString, susy::PFJetCollection>::iterator pfJets_it = ev.pfJets.find("ak5");
  if(pfJets_it != ev.pfJets.end()) {
    susy::PFJetCollection& jetColl = pfJets_it->second;
      
    for(vector<susy::PFJet>::iterator it = jetColl.begin();
	it != jetColl.end(); it++) {
	
      map<TString, Float_t>::iterator s_it = it->jecScaleFactors.find("L1FastL2L3");
      float scale = s_it->second;
      float jecScaleUncertainty = it->jecUncertainty;

      if(systematic == kJECup) scale *= (1.0 + jecScaleUncertainty);
      else if(systematic == kJECdown) scale *= (1.0 - jecScaleUncertainty);

      TLorentzVector corrP4 = scale * it->momentum;

      if(!isGoodJet(*it, corrP4)) continue;
      if(JetOverlapsElectron(corrP4, tightEles, looseEles)) continue;
      if(JetOverlapsMuon(corrP4, tightMuons, looseMuons)) continue;
      if(useDeltaRCutsOnPhotons && JetsOverlapsPhoton(corrP4, photons)) continue;
      if(!useDeltaRCutsOnPhotons && btagger == "CSVM" && it->bTagDiscriminators[susy::kCSV] <= 0.679 && JetsOverlapsPhoton(corrP4, photons)) continue;

      pfJets.push_back(&*it);
      csvValues.push_back(it->bTagDiscriminators[susy::kCSV]);
      HT += corrP4.Pt();
      hadronicSystem += corrP4;
      pfJets_corrP4.push_back(corrP4);
	
      if(fabs(corrP4.Eta()) < 2.4) {
	if(isMC) {
	  BtagInfo info((*it), corrP4, btagger, 1., isMC, isFastSim, sf);
	  tagInfos.push_back(info);
	}

	if((btagger == "CSVL" && it->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && it->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && it->bTagDiscriminators[susy::kCSV] > 0.898)) {
	  btags.push_back(&*it);
	  btags_corrP4.push_back(corrP4);
	}
	 
      }
   
    } // loop over jet coll
  } // if the jet coll exists
  sort(pfJets.begin(), pfJets.end(), EtGreater<susy::PFJet>);
  sort(btags.begin(), btags.end(), EtGreater<susy::PFJet>);
  sort(csvValues.begin(), csvValues.end(), greater<float>());
  sort(pfJets_corrP4.begin(), pfJets_corrP4.end(), CorrPtGreater);
  sort(btags_corrP4.begin(), btags_corrP4.end(), CorrPtGreater);

}


void SusyEventAnalyzer::findMuons(susy::Event& ev, vector<susy::Muon*>& tightMuons, vector<susy::Muon*>& looseMuons, float& HT, int mode) {

  map<TString, vector<susy::Muon> >::iterator muMap = ev.muons.find("muons");
  if(muMap != ev.muons.end()) {
    for(vector<susy::Muon>::iterator mu_it = muMap->second.begin(); mu_it != muMap->second.end(); mu_it++) {

      if((int)mu_it->bestTrackIndex() >= (int)(event.tracks).size() || (int)mu_it->bestTrackIndex() < 0) continue;

      // don't consider a muon that's already been picked. this shouldn't occur in the reconstruction anyway
      bool overlapsMuon = false;
      for(unsigned int i = 0; i < tightMuons.size(); i++) {
	if(deltaR(tightMuons[i]->momentum, mu_it->momentum) <= 0.1) {
	  overlapsMuon = true;
	  break;
	}
      }
      for(unsigned int i = 0; i < looseMuons.size(); i++) {
	if(deltaR(looseMuons[i]->momentum, mu_it->momentum) <= 0.1) {
	  overlapsMuon = true;
	  break;
	}
      }
      if(overlapsMuon) continue;

      bool passesTight = isTightMuon(*mu_it, 
				     event.tracks, 
				     d0correction(event.vertices[0].position, event.tracks[mu_it->bestTrackIndex()]), 
				     dZcorrection(event.vertices[0].position, event.tracks[mu_it->bestTrackIndex()]));

      if(mode != kMuonQCD) {
	if(passesTight && isIsolatedMuon(*mu_it)) {
	  tightMuons.push_back(&*mu_it);
	  HT += mu_it->momentum.Pt();
	}
	else if(isVetoMuon(*mu_it)) {
	  looseMuons.push_back(&*mu_it);
	  HT += mu_it->momentum.Pt();
	}
      }

      else {
	if(passesTight && isAntiIsolatedMuon(*mu_it)) {
	  tightMuons.push_back(&*mu_it);
	  HT += mu_it->momentum.Pt();
	}
	else if(isVetoMuon(*mu_it)) {
	  looseMuons.push_back(&*mu_it);
	  HT += mu_it->momentum.Pt();
	}
      }

    }

  }

}

void SusyEventAnalyzer::findElectrons(susy::Event& ev, vector<susy::Muon*> tightMuons, vector<susy::Muon*> looseMuons, vector<susy::Electron*>& tightEles, vector<susy::Electron*>& looseEles, float& HT, int mode) {

  map<TString, vector<susy::Electron> >::iterator eleMap = ev.electrons.find("gsfElectrons");
  if(eleMap != ev.electrons.end()) {
    for(vector<susy::Electron>::iterator ele_it = eleMap->second.begin(); ele_it != eleMap->second.end(); ele_it++) {

      if((int)ele_it->gsfTrackIndex >= (int)(event.tracks).size() || (int)ele_it->gsfTrackIndex < 0) continue;

      // this should never really happen...
      bool overlapsMuon = false;
      for(unsigned int i = 0; i < tightMuons.size(); i++) {
	if(deltaR(tightMuons[i]->momentum, ele_it->momentum) <= 0.1) {
	  overlapsMuon = true;
	  break;
	}
      }
      for(unsigned int i = 0; i < looseMuons.size(); i++) {
	if(deltaR(looseMuons[i]->momentum, ele_it->momentum) <= 0.1) {
	  overlapsMuon = true;
	  break;
	}
      }
      if(overlapsMuon) continue;

      // don't consider electrons already picked. this shouldn't happen in the reconstruction anyway
      bool overlapsElectron = false;
      for(unsigned int i = 0; i < tightEles.size(); i++) {
	if(deltaR(tightEles[i]->momentum, ele_it->momentum) <= 0.1) {
	  overlapsElectron = true;
	  break;
	}
      }
      for(unsigned int i = 0; i < looseEles.size(); i++) {
	if(deltaR(looseEles[i]->momentum, ele_it->momentum) <= 0.1) {
	  overlapsElectron = true;
	  break;
	}
      }
      if(overlapsElectron) continue;

      bool passesTight = isTightElectron(*ele_it, 
					 event.superClusters, 
					 event.rho25, 
					 d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]), 
					 dZcorrection(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));

      bool passesLoose = isLooseElectron(*ele_it,
					 event.superClusters, 
					 event.rho25, 
					 d0correction(event.vertices[0].position, event.tracks[ele_it->gsfTrackIndex]));

      if(mode != kElectronQCD) {
	if(passesTight && isIsolatedElectron(*ele_it, event.superClusters, event.rho25) && ele_it->mvaTrig > 0.5) {
	  tightEles.push_back(&*ele_it);
	  HT += ele_it->momentum.Pt();
	}
	else if(passesLoose) {
	  looseEles.push_back(&*ele_it);
	  HT += ele_it->momentum.Pt();
	}
      }

      else {
	if(passesTight && isAntiIsolatedElectron(*ele_it, event.superClusters, event.rho25) && ele_it->mvaTrig < -0.1 && ele_it->mvaTrig > -1.0) {
	  tightEles.push_back(&*ele_it);
	  HT += ele_it->momentum.Pt();
	}
	else if(passesLoose) {
	  looseEles.push_back(&*ele_it);
	  HT += ele_it->momentum.Pt();
	}
      }

    }
  }

}

void SusyEventAnalyzer::FillMetFilter2D(susy::Event& ev, TH2F*& h) {

  for(int i = 0; i < susy::nMetFilters; i++) {
    for(int j = 0; j < susy::nMetFilters; j++) {

      if(!ev.passMetFilter(i) && !ev.passMetFilter(j)) h->Fill(i, j);

    }
  }

}

bool SusyEventAnalyzer::GetDiJetPt(susy::Event& ev, vector<susy::Photon*> candidates, float& diJetPt, float& leadpt, float& trailpt) {

  vector<TLorentzVector> jetP4s;

  map<TString, susy::PFJetCollection>::iterator iJets = ev.pfJets.find("ak5");
  if(iJets == ev.pfJets.end()) return false;
  susy::PFJetCollection& jetColl = iJets->second;

  bool worked = true;

  for(vector<susy::Photon*>::iterator it = candidates.begin(); it != candidates.end(); it++) {
      
    bool matched = false;
      
    for(vector<susy::PFJet>::iterator iJet = jetColl.begin(); iJet != jetColl.end(); iJet++) {
      float theJES = 1.0;
      map<TString, Float_t>::const_iterator iCorr1 = iJet->jecScaleFactors.find("L1FastL2L3");
      map<TString, Float_t>::const_iterator iCorr2 = iJet->jecScaleFactors.find("L2L3");
      TLorentzVector corrP4 = iJet->momentum;
      if(iCorr1 != iJet->jecScaleFactors.end()) {
	if(iCorr2 != iJet->jecScaleFactors.end()) {
	  theJES = iCorr1->second / iCorr2->second;
	  corrP4 = theJES * iJet->momentum;
	}
	else {
	  theJES = iCorr1->second;
	  corrP4 = theJES * iJet->momentum;
	}
      }

      float dEta = corrP4.Eta() - (*it)->momentum.Eta();
      float dPhi = TVector2::Phi_mpi_pi(corrP4.Phi() - (*it)->momentum.Phi());
      float dR = sqrt(dEta*dEta + dPhi*dPhi);

      if(corrP4.Et() > 20. &&
	 fabs(corrP4.Eta()) < 2.6 &&
	 dR < 0.3) {
	  
	jetP4s.push_back(corrP4);
	matched = true;
	break;
      }

    }
      
    if(!matched) {
      worked = false;
      jetP4s.clear();
      jetP4s.push_back(candidates[0]->momentum);
      jetP4s.push_back(candidates[1]->momentum);
      break;
    }
      
  }
  
  sort(jetP4s.begin(), jetP4s.end(), CorrPtGreater);

  diJetPt = (jetP4s[0] + jetP4s[1]).Pt();
  leadpt = jetP4s[0].Pt();
  trailpt = jetP4s[1].Pt();

  return worked;
}

bool SusyEventAnalyzer::PhotonMatchesElectron(susy::Event& ev, vector<susy::Photon*> candidates, int& bothMatchCounter) {

  if(ev.isRealData) return false;

  bool matchesLead = false;
  bool matchesTrail = false;

  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

    if(it->status == 3 && fabs(it->pdgId) == 11) {
      if(deltaR(it->momentum, candidates[0]->momentum) < 0.1) matchesLead = true;
      if(deltaR(it->momentum, candidates[1]->momentum) < 0.1) matchesTrail = true;
    }

  }

  if(matchesLead && matchesTrail) bothMatchCounter++;

  return (matchesLead || matchesTrail);
}

int SusyEventAnalyzer::FigureTTbarDecayMode() {

  int decayMode = -1;

  int nElectronicWs = 0;
  int nMuonicWs = 0;
  int nTauonicWs = 0;
  
  int firstWindex = -1;
  
  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {
    
    //if(it->status != 3) continue;
    //if(it->momentum.Pt() < 20.) continue;
    
    bool isFromWfromStopOrTop = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
      (
       fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 6 || 
       fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000006
       );
    
    if(!isFromWfromStopOrTop) continue;
    
    if(firstWindex == it->motherIndex) continue;
    
    if(fabs(it->pdgId) == 11) nElectronicWs++;
    else if(fabs(it->pdgId) == 13) nMuonicWs++;
    else if(fabs(it->pdgId) == 15) nTauonicWs++;
    
    if(firstWindex < 0) firstWindex = it->motherIndex;
    
  }
  
  /* decayMode:
     0 hadronic
     1 semi-ele
     2 semi-mu
     3 semi-tau
     4 di-ele
     5 di-mu
     6 di-tau
     7 ele-mu
     8 ele-tau
     9 mu-tau
  */
  
  if((nElectronicWs + nMuonicWs + nTauonicWs) == 0) decayMode = 0;
  else if(nElectronicWs == 1 && (nMuonicWs + nTauonicWs) == 0) decayMode = 1;
  else if(nMuonicWs == 1 && (nElectronicWs + nTauonicWs) == 0) decayMode = 2;
  else if(nTauonicWs == 1 && (nElectronicWs + nMuonicWs) == 0) decayMode = 3;
  else if(nElectronicWs == 2 && (nMuonicWs + nTauonicWs) == 0) decayMode = 4;
  else if(nMuonicWs == 2 && (nElectronicWs + nTauonicWs) == 0) decayMode = 5;
  else if(nTauonicWs == 2 && (nElectronicWs + nMuonicWs) == 0) decayMode = 6;
  else if(nElectronicWs == 1 && nMuonicWs == 1 && nTauonicWs == 0) decayMode = 7;
  else if(nElectronicWs == 1 && nMuonicWs == 0 && nTauonicWs == 1) decayMode = 8;
  else if(nElectronicWs == 0 && nMuonicWs == 1 && nTauonicWs == 1) decayMode = 9;

  return decayMode;
}

double SusyEventAnalyzer::TopPtReweighting(susy::Event& ev) {

  susy::Particle * top     = 0;
  susy::Particle * antitop = 0;

  for(vector<susy::Particle>::iterator it = ev.genParticles.begin(); it != ev.genParticles.end(); it++) {
    if(abs(it->pdgId) != 6 || it->status != 3) continue;
    if(it->pdgId == 6 && !top) top = &*it;
    if(it->pdgId == -6 && !antitop) antitop = &*it;
    if(top && antitop) break;
  }

  if(!top || !antitop) return -1;

  double weight;

  if(scan == "ttJetsHadronic") {
    weight = 0.156 - 0.00137*top->momentum.Pt();
    weight += 0.156 - 0.00137*antitop->momentum.Pt();
    weight = exp(weight / 2.);
  }

  else if(scan == "ttJetsSemiLep") {
    weight = 0.159 - 0.00141*top->momentum.Pt();
    weight += 0.159 - 0.00141*antitop->momentum.Pt();
    weight = exp(weight / 2.);
  }

  else if(scan == "ttJetsFullLep") {
    weight = 0.148 - 0.00129*top->momentum.Pt();
    weight += 0.148 - 0.00129*antitop->momentum.Pt();
    weight = exp(weight / 2.);
  }

  else weight = -1;

  return weight;
}

double SusyEventAnalyzer::TopPtReweighting_ttHbb(susy::Event& ev) {

  susy::Particle * top     = 0;

  for(vector<susy::Particle>::iterator it = ev.genParticles.begin(); it != ev.genParticles.end(); it++) {
    if(it->pdgId == 6 && it->status == 3 && !top) {
      top = &*it;
      break;
    }
  }

  if(!top) return -1;

  double p0 = 1.18246e+00;
  double p1 = 4.63312e+02;
  double p2 = 2.10061e-06;

  double topPt = top->momentum.Pt();

  if(topPt > p1) topPt = p1;

  double weight = p0 + p2 * topPt * (topPt - 2 * p1);

  return weight;
}

void SusyEventAnalyzer::fill_whizard_phaseSpace(TH2D*& h) {

  vector<susy::Particle*> photons;
  vector<susy::Particle*> legs;

  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

    if(it->momentum.Pt() < 2.0) continue;

    int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

    bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
    bool fromWboson = (mom_id == 24);

    if(abs(it->pdgId) == 22 && (fromQCD || fromWboson)) photons.push_back(&*it);
    if(abs(it->pdgId) == 5 && mom_id == 6) legs.push_back(&*it);

  }

  for(unsigned int i = 0; i < photons.size(); i++) {

    double minDR = 999.0;
    double dr;

    for(unsigned int j = 0; j < legs.size(); j++) {
      dr = deltaR(photons[i]->momentum, legs[j]->momentum);
      if(dr < minDR) minDR = dr;
    }

    h->Fill(photons[i]->momentum.Pt(), minDR);
    
  }

}

void SusyEventAnalyzer::fill_madgraph_phaseSpace(TH2D*& h) {

  double eta_cut = 3.0;
  
  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

    if(it->momentum.Pt() < 2.0) continue;

    int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

    bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
    bool fromWboson = (mom_id == 24);
    bool fromLepton = (mom_id == 11 || mom_id == 13 || mom_id == 15);

    if(abs(it->pdgId) == 22 &&
       fabs(it->momentum.Eta()) < eta_cut &&
       (fromQCD || fromWboson || fromLepton)) {

      double minDR = 999.0;
      double dR;

      for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	if(oit->momentum.Pt() < 2.0) continue;
	if(oit == it) continue;
	if(oit->momentum.M() > 10) continue;
	if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	dR = deltaR(it->momentum, oit->momentum);

	if(dR < minDR) minDR = dR;

      } // other gen particle loop
	
      h->Fill(it->momentum.Pt(), minDR);

    } // if candidate photon

  } // photon loop

}

bool SusyEventAnalyzer::overlaps_whizard() {

  double et_cut = 20.0;
  double dR_cut = 0.1;

  bool haveOverlap = false;

  vector<susy::Particle*> photons;
  vector<susy::Particle*> legs;

  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

    if(it->momentum.Pt() < 2.0) continue;

    int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

    bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
    bool fromWboson = (mom_id == 24);

    if(abs(it->pdgId) == 22 && it->momentum.Et() > et_cut && (fromQCD || fromWboson)) photons.push_back(&*it);
    if(abs(it->pdgId) == 5 && mom_id == 6) legs.push_back(&*it);

  }

  for(unsigned int i = 0; i < photons.size(); i++) {

    bool haveLeg = false;
    bool closeToLeg = false;

    for(unsigned int j = 0; j < legs.size(); j++) {
      haveLeg = true;
      if(deltaR(photons[i]->momentum, legs[j]->momentum) < dR_cut) closeToLeg = true;
    }
    
    if(haveLeg && !closeToLeg) haveOverlap = true;
    
  }

  return haveOverlap;
}

bool SusyEventAnalyzer::overlaps_madgraph() {

  double et_cut = 13.0;
  double eta_cut = 3.0;
  double dR_cut = 0.2;

  bool haveOverlap = false;
  
  for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

    if(it->momentum.Pt() < 2.0) continue;

    int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

    bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
    bool fromWboson = (mom_id == 24);
    bool fromLepton = (mom_id == 11 || mom_id == 13 || mom_id == 15);

    if(abs(it->pdgId) == 22 && 
       it->momentum.Et() > et_cut &&
       fabs(it->momentum.Eta()) < eta_cut &&
       (fromQCD || fromWboson || fromLepton)) {

      double minDR = 999.0;
      double dR;

      for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	if(oit->momentum.Pt() < 2.0) continue;
	if(oit == it) continue;
	if(oit->momentum.M() > 10) continue;
	if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	dR = deltaR(it->momentum, oit->momentum);

	if(dR < minDR) minDR = dR;

      } // other gen particle loop
	
      if(minDR > dR_cut) haveOverlap = true;

    } // if candidate photon

  } // photon loop

  return haveOverlap;

}
  
void SusyEventAnalyzer::SetTreeValues(map<TString, float>& treeMap,
				      vector<susy::Muon*> tightMuons, vector<susy::Electron*> tightEles, 
				      vector<susy::PFJet*> pfJets, vector<susy::PFJet*> btags,
				      vector<susy::Photon*> photons,
				      vector<TLorentzVector> pfJets_corrP4, vector<TLorentzVector> btags_corrP4,
				      vector<float> csvValues,
				      TLorentzVector hadronicSystem,
				      float HT, float HT_jets,
				      int nPVertex,
				      float eventWeight, float eventWeightErr, float eventWeight_up, float eventWeight_down,
				      Long64_t jentry) {

  treeMap["Nmuons"] = tightMuons.size();
  treeMap["Nelectrons"] = tightEles.size();
  treeMap["nPV"] = nPVertex;
  treeMap["metFilterBit"] = event.metFilterBit;
  if(isMC && scan == "stop-bino") treeMap["ttbarDecayMode"] = FigureTTbarDecayMode();
  if(isMC) {
    treeMap["overlaps_whizard"] = overlaps_whizard();
    treeMap["overlaps_madgraph"] = overlaps_madgraph();
    treeMap["TopPtReweighting"] = TopPtReweighting(event);
    treeMap["TopPtReweighting_ttHbb"] = TopPtReweighting_ttHbb(event);
  }
  treeMap["Nphotons"] = photons.size();

  int nGamma = 0;
  int nFake = 0;
  for(unsigned int i = 0; i < photons.size(); i++) {
    if(is_loosePhoton(*photons[i], event.rho25)) nGamma++;
    if(is_fakePhoton(*photons[i], event.rho25)) nFake++;
  }
  treeMap["Ngamma"] = nGamma;
  treeMap["Nfake"] = nFake;

  treeMap["Njets"] = pfJets.size();
  treeMap["Nbtags"] = btags.size();
  treeMap["HT_jets"] = HT_jets;
  treeMap["HT"] = HT + HT_jets;
  treeMap["hadronic_pt"] = hadronicSystem.Pt();
  
  treeMap["max_csv"] = (csvValues.size() >= 1) ? csvValues[0] : -1.;
  treeMap["submax_csv"] = (csvValues.size() >= 2) ? csvValues[1] : -1.;
  treeMap["min_csv"] = (csvValues.size() >= 1) ? csvValues.back() : -1.;
  
  treeMap["jet1_pt"] = (pfJets_corrP4.size() >= 1) ? pfJets_corrP4[0].Pt() : -1.;
  treeMap["jet2_pt"] = (pfJets_corrP4.size() >= 2) ? pfJets_corrP4[1].Pt() : -1.;
  treeMap["jet3_pt"] = (pfJets_corrP4.size() >= 3) ? pfJets_corrP4[2].Pt() : -1.;
  treeMap["jet4_pt"] = (pfJets_corrP4.size() >= 4) ? pfJets_corrP4[3].Pt() : -1.;
  
  treeMap["btag1_pt"] = (btags_corrP4.size() >= 1) ? btags_corrP4[0].Pt() : -1.;
  treeMap["btag2_pt"] = (btags_corrP4.size() >= 2) ? btags_corrP4[1].Pt() : -1.;

  if(isMC) treeMap["pileupWeight"] = eventWeight;
  if(isMC) treeMap["pileupWeightErr"] = eventWeightErr;
  if(isMC) treeMap["pileupWeightUp"] = eventWeight_up;
  if(isMC) treeMap["pileupWeightDown"] = eventWeight_down;
        
  susy::MET* pfMet         = &(event.metMap.find("pfMet")->second);
  susy::MET* pfMetType1    = &(event.metMap.find("pfType1CorrectedMet")->second);
  susy::MET* pfMetType1p2  = &(event.metMap.find("pfType1p2CorrectedMet")->second);
  susy::MET* pfMetType01   = &(event.metMap.find("pfType01CorrectedMet")->second);
  susy::MET* pfMetType01p2 = &(event.metMap.find("pfType01p2CorrectedMet")->second);
  susy::MET* pfNoPileUpMet = &(event.metMap.find("pfNoPileUpMet")->second);
  susy::MET* pfMVAMet      = &(event.metMap.find("pfMVAMet")->second);
  susy::MET* genMet        = &(event.metMap.find("genMetTrue")->second);
  
  treeMap["pfMET"]     = pfMet->met();
  treeMap["pfMET_x"]   = pfMet->metX();
  treeMap["pfMET_y"]   = pfMet->metY();
  treeMap["pfMET_phi"] = pfMet->mEt.Phi();
  
  TVector2 sysShiftCorr(1.62861e-01 - 2.38517e-02*nPVertex, 3.60860e-01 - 1.30335e-01*nPVertex);
  if(!isMC) sysShiftCorr = TVector2(4.83642e-02 + 2.48870e-01*nPVertex, -1.50135e-01 - 8.27917e-02*nPVertex);
  treeMap["pfMET_sysShift"]     = (pfMet->mEt - sysShiftCorr).Mod();
  treeMap["pfMET_sysShift_phi"] = (pfMet->mEt - sysShiftCorr).Phi();
  
  treeMap["pfMET_t1"]    = pfMetType1->met();
  treeMap["pfMET_t1p2"]  = pfMetType1p2->met();
  treeMap["pfMET_t01"]   = pfMetType01->met();
  treeMap["pfMET_t01p2"] = pfMetType01p2->met();
  treeMap["pfNoPUMET"]   = pfNoPileUpMet->met();
  treeMap["pfMVAMET"]    = pfMVAMet->met();
  if(isMC) treeMap["genMET"]      = genMet->met();

  float mLepGammaLead = -1.;
  float mLepGammaTrail = -1.;
  float mLepGammaGamma = -1.;
  if(photons.size() > 0) {
    if(tightEles.size() == 1) mLepGammaLead = (tightEles[0]->momentum + photons[0]->momentum).M();
    else if(tightMuons.size() == 1) mLepGammaLead = (tightMuons[0]->momentum + photons[0]->momentum).M();
  }
  if(photons.size() > 1) {
    if(tightEles.size() == 1) {
      mLepGammaTrail = (tightEles[0]->momentum + photons[1]->momentum).M();
      mLepGammaGamma = (tightEles[0]->momentum + photons[0]->momentum + photons[1]->momentum).M();
    }
    else if(tightMuons.size() == 1) {
      mLepGammaTrail = (tightMuons[0]->momentum + photons[1]->momentum).M();
      mLepGammaGamma = (tightMuons[0]->momentum + photons[0]->momentum + photons[1]->momentum).M();
    }
  }
  treeMap["mLepGammaLead"] = (photons.size() > 0) ? mLepGammaLead : -1;
  treeMap["mLepGammaTrail"] = (photons.size() > 1) ? mLepGammaTrail : -1;
  treeMap["mLepGammaGamma"] = (photons.size() > 1) ? mLepGammaGamma : -1;
  
  // Photon angles
  if(photons.size() > 0) {
    if(tightEles.size() == 1) {
      treeMap["dR_leadPhoton_l"] = deltaR(tightEles[0]->momentum, photons[0]->momentum);
      treeMap["dEta_leadPhoton_l"] = fabs(tightEles[0]->momentum.Eta() - photons[0]->momentum.Eta());
      treeMap["dPhi_leadPhoton_l"] = TVector2::Phi_mpi_pi(tightEles[0]->momentum.Phi() - photons[0]->momentum.Phi());
      treeMap["cosTheta_leadPhoton_l"] = TMath::Cos(tightEles[0]->momentum.Angle(photons[0]->momentum.Vect()));
    }
    else if(tightMuons.size() == 1) {
      treeMap["dR_leadPhoton_l"] = deltaR(tightMuons[0]->momentum, photons[0]->momentum);
      treeMap["dEta_leadPhoton_l"] = fabs(tightMuons[0]->momentum.Eta() - photons[0]->momentum.Eta());
      treeMap["dPhi_leadPhoton_l"] = TVector2::Phi_mpi_pi(tightMuons[0]->momentum.Phi() - photons[0]->momentum.Phi());
      treeMap["cosTheta_leadPhoton_l"] = TMath::Cos(tightMuons[0]->momentum.Angle(photons[0]->momentum.Vect()));
    }

    if(btags_corrP4.size() > 0) {
      double min_b_angle = 999.;
      int min_b_index = -1;
      for(unsigned int i = 0; i < btags_corrP4.size(); i++) {
	min_b_angle = min(min_b_angle, (double)deltaR(btags_corrP4[i], photons[0]->momentum));
	min_b_index = i;
      }
      treeMap["dR_leadPhoton_b_min"] = min_b_angle;
      treeMap["dEta_leadPhoton_b_min"] = fabs(btags_corrP4[min_b_index].Eta() - photons[0]->momentum.Eta());
      treeMap["dPhi_leadPhoton_b_min"] = TVector2::Phi_mpi_pi(btags_corrP4[min_b_index].Phi() - photons[0]->momentum.Phi());
      treeMap["cosTheta_leadPhoton_b_min"] = TMath::Cos(btags_corrP4[min_b_index].Angle(photons[0]->momentum.Vect()));
    }
    else {
      treeMap["dR_leadPhoton_b_min"] = -2.;
      treeMap["dEta_leadPhoton_b_min"] = -2.;
      treeMap["dPhi_leadPhoton_b_min"] = -200.;
      treeMap["cosTheta_leadPhoton_b_min"] = -200.;
    }

  }
  else {
    treeMap["dR_leadPhoton_l"] = -1.;
    treeMap["dEta_leadPhoton_l"] = -1.;
    treeMap["dPhi_leadPhoton_l"] = -100.;
    treeMap["cosTheta_leadPhoton_l"] = -100.;
    treeMap["dR_leadPhoton_b_min"] = -1.;
    treeMap["dEta_leadPhoton_b_min"] = -1.;
    treeMap["dPhi_leadPhoton_b_min"] = -100.;
    treeMap["cosTheta_leadPhoton_b_min"] = -100.;
  }

  if(photons.size() > 1) {
    if(tightEles.size() == 1) {
      treeMap["dR_trailPhoton_l"] = deltaR(tightEles[0]->momentum, photons[1]->momentum);
      treeMap["dEta_trailPhoton_l"] = fabs(tightEles[0]->momentum.Eta() - photons[1]->momentum.Eta());
      treeMap["dPhi_trailPhoton_l"] = TVector2::Phi_mpi_pi(tightEles[0]->momentum.Phi() - photons[1]->momentum.Phi());
      treeMap["cosTheta_trailPhoton_l"] = TMath::Cos(tightEles[0]->momentum.Angle(photons[1]->momentum.Vect()));
    }
    else if(tightMuons.size() == 1) {
      treeMap["dR_trailPhoton_l"] = deltaR(tightMuons[0]->momentum, photons[1]->momentum);
      treeMap["dEta_trailPhoton_l"] = fabs(tightMuons[0]->momentum.Eta() - photons[1]->momentum.Eta());
      treeMap["dPhi_trailPhoton_l"] = TVector2::Phi_mpi_pi(tightMuons[0]->momentum.Phi() - photons[1]->momentum.Phi());
      treeMap["cosTheta_trailPhoton_l"] = TMath::Cos(tightMuons[0]->momentum.Angle(photons[1]->momentum.Vect()));
    }

    if(btags_corrP4.size() > 0) {
      double min_b_angle = 999.;
      int min_b_index = -1;
      for(unsigned int i = 0; i < btags_corrP4.size(); i++) {
	min_b_angle = min(min_b_angle, (double)deltaR(btags_corrP4[i], photons[1]->momentum));
	min_b_index = i;
      }
      treeMap["dR_trailPhoton_b_min"] = min_b_angle;
      treeMap["dEta_trailPhoton_b_min"] = fabs(btags_corrP4[min_b_index].Eta() - photons[1]->momentum.Eta());
      treeMap["dPhi_trailPhoton_b_min"] = TVector2::Phi_mpi_pi(btags_corrP4[min_b_index].Phi() - photons[1]->momentum.Phi());
      treeMap["cosTheta_trailPhoton_b_min"] = TMath::Cos(btags_corrP4[min_b_index].Angle(photons[1]->momentum.Vect()));
    }
    else {
      treeMap["dR_trailPhoton_b_min"] = -2.;
      treeMap["dEta_trailPhoton_b_min"] = -2.;
      treeMap["dPhi_trailPhoton_b_min"] = -200.;
      treeMap["cosTheta_trailPhoton_b_min"] = -200.;
    }

  }
  else {
    treeMap["dR_trailPhoton_l"] = -1.;
    treeMap["dEta_trailPhoton_l"] = -1.;
    treeMap["dPhi_trailPhoton_l"] = -100.;
    treeMap["cosTheta_trailPhoton_l"] = -100.;
    treeMap["dR_trailPhoton_b_min"] = -1.;
    treeMap["dEta_trailPhoton_b_min"] = -1.;
    treeMap["dPhi_trailPhoton_b_min"] = -100.;
    treeMap["cosTheta_trailPhoton_b_min"] = -100.;
  }

  // photon gen match
  if(isMC) {
    if(photons.size() > 0) {

      bool foundMatch = false;

      for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

	if(it->momentum.Pt() < 2.0) continue;

	bool et_eta_match = deltaR(it->momentum, photons[0]->momentum) < 0.2 &&
	  (fabs(photons[0]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 1.0;

	if(!et_eta_match) continue;

	int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

	bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
	bool fromBoson = (mom_id == 23 || mom_id == 24 || mom_id == 25);
	bool fromLepton = (mom_id == 11 || mom_id == 13 || mom_id == 15);

	if(et_eta_match && abs(it->pdgId) == 22) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {
	    
	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if((fromQCD || fromBoson || fromLepton) &&
	     (fabs(photons[0]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[0]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[0]->momentum) < 0.01) {
	    treeMap["leadMatchGamma"] = 1;
	    treeMap["leadMatchElectron"] = 0;
	    treeMap["leadMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["leadMatchGamma"] = 0;
	    treeMap["leadMatchElectron"] = 0;
	    treeMap["leadMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to a photon

	else if(et_eta_match && abs(it->pdgId) == 11) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if(fromBoson &&
	     (fabs(photons[0]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[0]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[0]->momentum) < 0.04) {
	    treeMap["leadMatchGamma"] = 0;
	    treeMap["leadMatchElectron"] = 1;
	    treeMap["leadMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["leadMatchGamma"] = 0;
	    treeMap["leadMatchElectron"] = 0;
	    treeMap["leadMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to an electron

      } // gen particles loop

      if(!foundMatch) {
	treeMap["leadMatchGamma"] = 0;
	treeMap["leadMatchElectron"] = 0;
	treeMap["leadMatchJet"] = 1;
      }
	
    } // nPhotons > 0

    if(photons.size() > 1) {
      
      bool foundMatch = false;

      for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

	if(it->momentum.Pt() < 2.0) continue;

	bool et_eta_match = deltaR(it->momentum, photons[1]->momentum) < 0.2 &&
	  (fabs(photons[1]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 1.0;

	if(!et_eta_match) continue;

	int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

	bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
	bool fromBoson = (mom_id == 23 || mom_id == 24 || mom_id == 25);
	bool fromLepton = (mom_id == 11 || mom_id == 13 || mom_id == 15);

	if(et_eta_match && abs(it->pdgId) == 22) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if((fromQCD || fromBoson || fromLepton) &&
	     (fabs(photons[1]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[1]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[1]->momentum) < 0.01) {
	    treeMap["trailMatchGamma"] = 1;
	    treeMap["trailMatchElectron"] = 0;
	    treeMap["trailMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["trailMatchGamma"] = 0;
	    treeMap["trailMatchElectron"] = 0;
	    treeMap["trailMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to a photon

	else if(et_eta_match && abs(it->pdgId) == 11) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if(fromBoson &&
	     (fabs(photons[1]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[1]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[1]->momentum) < 0.04) {
	    treeMap["trailMatchGamma"] = 0;
	    treeMap["trailMatchElectron"] = 1;
	    treeMap["trailMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["trailMatchGamma"] = 0;
	    treeMap["trailMatchElectron"] = 0;
	    treeMap["trailMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to an electron

	if(foundMatch) break;

      } // gen particles loop

      if(!foundMatch) {
	treeMap["trailMatchGamma"] = 0;
	treeMap["trailMatchElectron"] = 0;
	treeMap["trailMatchJet"] = 1;
      }

    } // if nphotons > 1

  } // ifMC


  // Transverse W mass

  float metphi = (pfMet->mEt - sysShiftCorr).Phi();
  float metphi_t1 = (pfMetType1->mEt - sysShiftCorr).Phi();
  float metphi_t1p2 = (pfMetType1p2->mEt - sysShiftCorr).Phi();
  float metphi_t01 = (pfMetType01->mEt - sysShiftCorr).Phi();
  float metphi_t01p2 = (pfMetType01p2->mEt - sysShiftCorr).Phi();
  float metphi_nopumet = (pfNoPileUpMet->mEt - sysShiftCorr).Phi();
  float metphi_mvamet = (pfMVAMet->mEt - sysShiftCorr).Phi();
  float metphi_genmet = (isMC) ? (genMet->mEt - sysShiftCorr).Phi() : -100.;

  float leptonphi = -100.;
  float leptonpt = -1.;
  if(tightEles.size() == 1) {
    leptonphi = tightEles[0]->momentum.Phi();
    leptonpt = tightEles[0]->momentum.Pt();
  }
  if(tightMuons.size() == 1) {
    leptonphi = tightMuons[0]->momentum.Phi();
    leptonpt = tightMuons[0]->momentum.Pt();
  }

  if(tightMuons.size() + tightEles.size() == 1) {

    treeMap["dPhi_met_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi);
    treeMap["dPhi_met_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi);

    treeMap["dPhi_met_t1_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi_t1);
    treeMap["dPhi_met_t1_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi_t1);

    treeMap["dPhi_met_t1p2_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi_t1p2);
    treeMap["dPhi_met_t1p2_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi_t1p2);

    treeMap["dPhi_met_t01_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi_t01);
    treeMap["dPhi_met_t01_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi_t01);

    treeMap["dPhi_met_t01p2_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi_t01p2);
    treeMap["dPhi_met_t01p2_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi_t01p2);

    treeMap["dPhi_nopumet_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi_nopumet);
    treeMap["dPhi_nopumet_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi_nopumet);
    
    treeMap["dPhi_mvamet_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi_mvamet);
    treeMap["dPhi_mvamet_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi_mvamet);

    if(isMC) {
      treeMap["dPhi_genmet_l"] = TVector2::Phi_mpi_pi(leptonphi - metphi_genmet);
      treeMap["dPhi_genmet_ht"] = TVector2::Phi_mpi_pi(hadronicSystem.Phi() - metphi_genmet);
    }

    float w_mT = 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi));
    w_mT *= 2. * leptonpt * treeMap["pfMET_sysShift"];

    float w_mT_t1 = 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi_t1));
    w_mT_t1 *= 2. * leptonpt * treeMap["pfMET_t1"];

    float w_mT_t1p2 = 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi_t1p2));
    w_mT_t1p2 *= 2. * leptonpt * treeMap["pfMET_t1p2"];

    float w_mT_t01 = 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi_t01));
    w_mT_t01 *= 2. * leptonpt * treeMap["pfMET_t01"];

    float w_mT_t01p2 = 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi_t01p2));
    w_mT_t01p2 *= 2. * leptonpt * treeMap["pfMET_t01p2"];

    float w_mT_nopumet = 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi_nopumet));
    w_mT_nopumet *= 2. * leptonpt * treeMap["pfMET_nopumet"];

    float w_mT_mvamet = 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi_mvamet));
    w_mT_mvamet *= 2. * leptonpt * treeMap["pfMET_mvamet"];

    float w_mT_genmet = (isMC) ? 1. - TMath::Cos(TVector2::Phi_mpi_pi(leptonphi - metphi_genmet)) : -1.;
    if(isMC) w_mT_genmet *= 2. * leptonpt * treeMap["pfMET_genmet"];

    treeMap["w_mT"] = sqrt(w_mT);
    treeMap["w_mT_t1"] = sqrt(w_mT_t1);
    treeMap["w_mT_t1p2"] = sqrt(w_mT_t1p2);
    treeMap["w_mT_t01"] = sqrt(w_mT_t01);
    treeMap["w_mT_t01p2"] = sqrt(w_mT_t01p2);
    treeMap["w_mT_nopumet"] = sqrt(w_mT_nopumet);
    treeMap["w_mT_mvamet"] = sqrt(w_mT_mvamet);
    treeMap["w_mT_genmet"] = (isMC) ? sqrt(w_mT_genmet) : -1.;
  }
  else {
    treeMap["w_mT"] = -1.;
    treeMap["w_mT_t1"] = -1.;
    treeMap["w_mT_t1p2"] = -1.;
    treeMap["w_mT_t01"] = -1.;
    treeMap["w_mT_t01p2"] = -1.;
    treeMap["w_mT_nopumet"] = -1.;
    treeMap["w_mT_mvamet"] = -1.;
    treeMap["w_mT_genmet"] = -1.;
  }
  
  // M3 calculation
  if(pfJets_corrP4.size() > 2) {
    int maxSumPt = 0.;
    int max_pos = 0;
    int comb = 1 << pfJets_corrP4.size();
    
    for(int i = 0; i < comb; i++) {

      int npicked = 0; 
      for(unsigned int j = 0; j < pfJets_corrP4.size(); j++) { 
	if(((i >> j) & 0x1) == 1) npicked++;
      } 

      if(npicked == 3) {
	TLorentzVector thisCombination(0., 0., 0., 0.);
	for(unsigned int j = 0; j < pfJets_corrP4.size(); j++) {
	  if(((i >> j ) & 0x1) == 1) thisCombination += pfJets_corrP4[j];
	}
	if(thisCombination.Pt() > maxSumPt) {
	  maxSumPt = thisCombination.Pt();
	  max_pos = i;
	}
      }

    }
    
    TLorentzVector maxCombination(0., 0., 0., 0.);
    for(unsigned int i = 0; i < pfJets_corrP4.size(); i++) {
      if(((max_pos >> i) & 0x1) == 1) maxCombination += pfJets_corrP4[i];
    }
    
    treeMap["m3"] = maxCombination.M();
  }
  else treeMap["m3"] = -1.;

  // M3 calculation -- uncorrected jets
  if(pfJets.size() > 2) {
    int maxSumPt = 0.;
    int max_pos = 0;
    int comb = 1 << pfJets.size();
    
    for(int i = 0; i < comb; i++) {

      int npicked = 0; 
      for(unsigned int j = 0; j < pfJets.size(); j++) { 
	if(((i >> j) & 0x1) == 1) npicked++; 
      } 

      if(npicked == 3) {
	TLorentzVector thisCombination(0., 0., 0., 0.);
	for(unsigned int j = 0; j < pfJets.size(); j++) {
	  if(((i >> j ) & 0x1) == 1) thisCombination += pfJets[j]->momentum;
	}
	if(thisCombination.Pt() > maxSumPt) {
	  maxSumPt = thisCombination.Pt();
	  max_pos = i;
	}
      }

    }
    
    TLorentzVector maxCombination(0., 0., 0., 0.);
    for(unsigned int i = 0; i < pfJets.size(); i++) {
      if(((max_pos >> i) & 0x1) == 1) maxCombination += pfJets[i]->momentum;
    }
    
    treeMap["m3_uncorr"] = maxCombination.M();
  }
  else treeMap["m3_uncorr"] = -1.;

  treeMap["ele_pt"] = (tightEles.size() > 0) ? tightEles[0]->momentum.Pt() : -1.;
  treeMap["ele_phi"] = (tightEles.size() > 0) ? tightEles[0]->momentum.Phi() : -100.;
  treeMap["ele_eta"] = (tightEles.size() > 0) ? tightEles[0]->momentum.Eta() : -100.;
  treeMap["ele_mvaTrigV0"] = (tightEles.size() > 0) ? tightEles[0]->mvaTrig : -100.;

  if(tightEles.size() > 0) {
    float ele_eta = fabs(event.superClusters[tightEles[0]->superClusterIndex].position.Eta());
    float ea;
    if(ele_eta < 1.0) ea = 0.13;
    else if(ele_eta < 1.479) ea = 0.14;
    else if(ele_eta < 2.0) ea = 0.07;
    else if(ele_eta < 2.2) ea = 0.09;
    else if(ele_eta < 2.3) ea = 0.11;
    else if(ele_eta < 2.4) ea = 0.11;
    else ea = 0.14;
      
    float ele_iso = max(0., (double)(tightEles[0]->photonIso + tightEles[0]->neutralHadronIso - event.rho25*ea));
    ele_iso += tightEles[0]->chargedHadronIso;
    
    treeMap["ele_relIso"] = ele_iso / tightEles[0]->momentum.Pt();
  }
  else treeMap["ele_relIso"] = -100.;

  treeMap["muon_pt"] = (tightMuons.size() > 0) ? tightMuons[0]->momentum.Pt() : -1.;
  treeMap["muon_phi"] = (tightMuons.size() > 0) ? tightMuons[0]->momentum.Phi() : -100.;
  treeMap["muon_eta"] = (tightMuons.size() > 0) ? tightMuons[0]->momentum.Eta() : -100.;

  if(tightMuons.size() > 0) {
    float mu_iso = max(0., (double)(tightMuons[0]->sumNeutralHadronEt04 + tightMuons[0]->sumPhotonEt04 - 0.5*(tightMuons[0]->sumPUPt04)));
    mu_iso += tightMuons[0]->sumChargedHadronPt04;
    float mu_pt = tightMuons[0]->momentum.Pt();
    treeMap["muon_relIso"] = mu_iso / mu_pt;
  }
  else treeMap["muon_relIso"] = -10.;
  
  treeMap["leadPhotonEt"] = (photons.size() > 0) ? photons[0]->momentum.Et() : -1.;
  treeMap["leadPhotonEta"] = (photons.size() > 0) ? photons[0]->momentum.Eta() : -100.;
  treeMap["leadPhotonPhi"] = (photons.size() > 0) ? photons[0]->momentum.Phi() : -100.;
  treeMap["leadChargedHadronIso"] = (photons.size() > 0) ? chargedHadronIso_corrected(*photons[0], event.rho25) : -100.;
  treeMap["leadNeutralHadronIso"] = (photons.size() > 0) ? neutralHadronIso_corrected(*photons[0], event.rho25) - 0.04*(photons[0]->momentum.Pt()) : -100.;
  treeMap["leadPhotonIso"] = (photons.size() > 0) ? photonIso_corrected(*photons[0], event.rho25) - 0.005*(photons[0]->momentum.Pt()) : -100.;
  treeMap["leadSigmaIetaIeta"] = (photons.size() > 0) ? photons[0]->sigmaIetaIeta : -100.;
  treeMap["lead_nPixelSeeds"] = (photons.size() > 0) ? photons[0]->nPixelSeeds : -10.;
  treeMap["leadMVAregEnergy"] = (photons.size() > 0) ? photons[0]->MVAregEnergy : -10.;
  treeMap["leadMVAregErr"] = (photons.size() > 0) ? photons[0]->MVAregErr : -10.;
  
  treeMap["trailPhotonEt"] = (photons.size() > 1) ? photons[1]->momentum.Et() : -1.;
  treeMap["trailPhotonEta"] = (photons.size() > 1) ? photons[1]->momentum.Eta() : -100.;
  treeMap["trailPhotonPhi"] = (photons.size() > 1) ? photons[1]->momentum.Phi() : -100.;
  treeMap["trail_nPixelSeeds"] = (photons.size() > 1) ? photons[1]->nPixelSeeds : -100.;
  treeMap["trailChargedHadronIso"] = (photons.size() > 1) ? chargedHadronIso_corrected(*photons[1], event.rho25) : -100.;
  treeMap["trailNeutralHadronIso"] = (photons.size() > 1) ? neutralHadronIso_corrected(*photons[1], event.rho25) - 0.04*(photons[1]->momentum.Pt()) : -100.;
  treeMap["trailPhotonIso"] = (photons.size() > 1) ? photonIso_corrected(*photons[1], event.rho25) - 0.005*(photons[1]->momentum.Pt()) : -100.;
  treeMap["trailSigmaIetaIeta"] = (photons.size() > 1) ? photons[1]->sigmaIetaIeta : -100.;
  treeMap["trailMVAregEnergy"] = (photons.size() > 1) ? photons[1]->MVAregEnergy : -10.;
  treeMap["trailMVAregErr"] = (photons.size() > 1) ? photons[1]->MVAregErr : -10.;
  
  treeMap["diEMpT"] = (photons.size() > 1) ? (photons[0]->momentum + photons[1]->momentum).Pt() : -1.;
  treeMap["photon_invmass"] = (photons.size() > 1) ? (photons[0]->momentum + photons[1]->momentum).M() : -10.;
  
  float diJetPt, lead_matched_jetpt, trail_matched_jetpt;
  bool matchingWorked = (photons.size() > 1 && GetDiJetPt(event, photons, diJetPt, lead_matched_jetpt, trail_matched_jetpt));
  treeMap["diJetPt"] = matchingWorked ? diJetPt : -1.;
  
  float dEta = (photons.size() > 1) ? photons[0]->momentum.Eta() - photons[1]->momentum.Eta() : -100;
  float photon_dPhi = (photons.size() > 1) ? TVector2::Phi_mpi_pi(photons[0]->momentum.Phi() - photons[1]->momentum.Phi()) : -100;
  float photon_dR = (photons.size() > 1) ? sqrt(dEta*dEta + photon_dPhi*photon_dPhi) : -100;
  treeMap["photon_dPhi"] = photon_dPhi;
  treeMap["photon_dR"] = photon_dR;
  
  if(!isMC) {
    treeMap["runNumber"] = event.runNumber;
    treeMap["eventNumber"] = event.eventNumber;
    treeMap["lumiBlock"] = event.luminosityBlockNumber;
    treeMap["jentry"] = jentry;
  }

}

void SusyEventAnalyzer::IncludeSyncFile(char* file) {
  FILE* syncList = fopen(file, "r");
  if(!syncList) return;

  char line[256];
  Int_t _run, _lumi;
  ULong_t _event;

  while(fgets(line, 255, syncList)) {
    sscanf(line, "%d %d %lu", &_run, &_lumi, &_event);
    syncRuns.push_back(_run);
    syncLumi.push_back(_lumi);
    syncEvents.push_back(_event);
  }
    
  fclose(syncList);
}

#endif

void SusyEventAnalyzer::SetTreeValues_ZGamma(map<TString, float>& treeMap,
					     vector<susy::Muon*> tightMuons, vector<susy::Electron*> tightEles, 
					     vector<susy::PFJet*> pfJets, vector<susy::PFJet*> btags,
					     vector<susy::Photon*> photons,
					     vector<TLorentzVector> pfJets_corrP4, vector<TLorentzVector> btags_corrP4,
					     vector<float> csvValues,
					     TLorentzVector hadronicSystem,
					     float HT, float HT_jets,
					     int nPVertex,
					     float eventWeight, float eventWeightErr, float eventWeight_up, float eventWeight_down,
					     Long64_t jentry) {

  treeMap["nPV"] = nPVertex;
  treeMap["metFilterBit"] = event.metFilterBit;
  if(isMC && scan == "stop-bino") treeMap["ttbarDecayMode"] = FigureTTbarDecayMode();
  if(isMC) {
    treeMap["overlaps_whizard"] = overlaps_whizard();
    treeMap["overlaps_madgraph"] = overlaps_madgraph();
    treeMap["TopPtReweighting"] = TopPtReweighting(event);
    treeMap["TopPtReweighting_ttHbb"] = TopPtReweighting_ttHbb(event);
  }
  treeMap["Nphotons"] = photons.size();

  int nGamma = 0;
  int nFake = 0;
  for(unsigned int i = 0; i < photons.size(); i++) {
    if(is_loosePhoton(*photons[i], event.rho25)) nGamma++;
    if(is_fakePhoton(*photons[i], event.rho25)) nFake++;
  }
  treeMap["Ngamma"] = nGamma;
  treeMap["Nfake"] = nFake;

  treeMap["Njets"] = pfJets.size();
  treeMap["Nbtags"] = btags.size();
  treeMap["HT_jets"] = HT_jets;
  treeMap["HT"] = HT + HT_jets;
  treeMap["hadronic_pt"] = hadronicSystem.Pt();
  
  treeMap["max_csv"] = (csvValues.size() >= 1) ? csvValues[0] : -1.;
  treeMap["submax_csv"] = (csvValues.size() >= 2) ? csvValues[1] : -1.;
  treeMap["min_csv"] = (csvValues.size() >= 1) ? csvValues.back() : -1.;
  
  treeMap["jet1_pt"] = (pfJets_corrP4.size() >= 1) ? pfJets_corrP4[0].Pt() : -1.;
  treeMap["jet2_pt"] = (pfJets_corrP4.size() >= 2) ? pfJets_corrP4[1].Pt() : -1.;
  treeMap["jet3_pt"] = (pfJets_corrP4.size() >= 3) ? pfJets_corrP4[2].Pt() : -1.;
  treeMap["jet4_pt"] = (pfJets_corrP4.size() >= 4) ? pfJets_corrP4[3].Pt() : -1.;
  
  treeMap["btag1_pt"] = (btags_corrP4.size() >= 1) ? btags_corrP4[0].Pt() : -1.;
  treeMap["btag2_pt"] = (btags_corrP4.size() >= 2) ? btags_corrP4[1].Pt() : -1.;

  if(isMC) treeMap["pileupWeight"] = eventWeight;
  if(isMC) treeMap["pileupWeightErr"] = eventWeightErr;
  if(isMC) treeMap["pileupWeightUp"] = eventWeight_up;
  if(isMC) treeMap["pileupWeightDown"] = eventWeight_down;
        
  susy::MET* pfMet         = &(event.metMap.find("pfMet")->second);
  susy::MET* pfMetType1    = &(event.metMap.find("pfType1CorrectedMet")->second);
  susy::MET* pfMetType1p2  = &(event.metMap.find("pfType1p2CorrectedMet")->second);
  susy::MET* pfMetType01   = &(event.metMap.find("pfType01CorrectedMet")->second);
  susy::MET* pfMetType01p2 = &(event.metMap.find("pfType01p2CorrectedMet")->second);
  susy::MET* pfNoPileUpMet = &(event.metMap.find("pfNoPileUpMet")->second);
  susy::MET* pfMVAMet      = &(event.metMap.find("pfMVAMet")->second);
  susy::MET* genMet        = &(event.metMap.find("genMetTrue")->second);
  
  treeMap["pfMET"]     = pfMet->met();
  treeMap["pfMET_x"]   = pfMet->metX();
  treeMap["pfMET_y"]   = pfMet->metY();
  treeMap["pfMET_phi"] = pfMet->mEt.Phi();
  
  TVector2 sysShiftCorr(1.62861e-01 - 2.38517e-02*nPVertex, 3.60860e-01 - 1.30335e-01*nPVertex);
  if(!isMC) sysShiftCorr = TVector2(4.83642e-02 + 2.48870e-01*nPVertex, -1.50135e-01 - 8.27917e-02*nPVertex);
  treeMap["pfMET_sysShift"]     = (pfMet->mEt - sysShiftCorr).Mod();
  treeMap["pfMET_sysShift_phi"] = (pfMet->mEt - sysShiftCorr).Phi();
  
  treeMap["pfMET_t1"]    = pfMetType1->met();
  treeMap["pfMET_t1p2"]  = pfMetType1p2->met();
  treeMap["pfMET_t01"]   = pfMetType01->met();
  treeMap["pfMET_t01p2"] = pfMetType01p2->met();
  treeMap["pfNoPUMET"]   = pfNoPileUpMet->met();
  treeMap["pfMVAMET"]    = pfMVAMet->met();
  if(isMC) treeMap["genMET"]      = genMet->met();

  // photon gen match
  if(isMC) {
    if(photons.size() > 0) {

      bool foundMatch = false;

      for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

	if(it->momentum.Pt() < 2.0) continue;

	bool et_eta_match = deltaR(it->momentum, photons[0]->momentum) < 0.2 &&
	  (fabs(photons[0]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 1.0;

	if(!et_eta_match) continue;

	int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

	bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
	bool fromBoson = (mom_id == 23 || mom_id == 24 || mom_id == 25);
	bool fromLepton = (mom_id == 11 || mom_id == 13 || mom_id == 15);

	if(et_eta_match && abs(it->pdgId) == 22) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {
	    
	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if((fromQCD || fromBoson || fromLepton) &&
	     (fabs(photons[0]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[0]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[0]->momentum) < 0.01) {
	    treeMap["leadMatchGamma"] = 1;
	    treeMap["leadMatchElectron"] = 0;
	    treeMap["leadMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["leadMatchGamma"] = 0;
	    treeMap["leadMatchElectron"] = 0;
	    treeMap["leadMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to a photon

	else if(et_eta_match && abs(it->pdgId) == 11) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if(fromBoson &&
	     (fabs(photons[0]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[0]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[0]->momentum) < 0.04) {
	    treeMap["leadMatchGamma"] = 0;
	    treeMap["leadMatchElectron"] = 1;
	    treeMap["leadMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["leadMatchGamma"] = 0;
	    treeMap["leadMatchElectron"] = 0;
	    treeMap["leadMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to an electron

      } // gen particles loop

      if(!foundMatch) {
	treeMap["leadMatchGamma"] = 0;
	treeMap["leadMatchElectron"] = 0;
	treeMap["leadMatchJet"] = 1;
      }
	
    } // nPhotons > 0

    if(photons.size() > 1) {
      
      bool foundMatch = false;

      for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

	if(it->momentum.Pt() < 2.0) continue;

	bool et_eta_match = deltaR(it->momentum, photons[1]->momentum) < 0.2 &&
	  (fabs(photons[1]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 1.0;

	if(!et_eta_match) continue;

	int mom_id = abs(event.genParticles[it->motherIndex].pdgId);

	bool fromQCD = ((mom_id >= 1 && mom_id <= 6) || mom_id == 2212 || mom_id == 21);
	bool fromBoson = (mom_id == 23 || mom_id == 24 || mom_id == 25);
	bool fromLepton = (mom_id == 11 || mom_id == 13 || mom_id == 15);

	if(et_eta_match && abs(it->pdgId) == 22) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if((fromQCD || fromBoson || fromLepton) &&
	     (fabs(photons[1]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[1]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[1]->momentum) < 0.01) {
	    treeMap["trailMatchGamma"] = 1;
	    treeMap["trailMatchElectron"] = 0;
	    treeMap["trailMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["trailMatchGamma"] = 0;
	    treeMap["trailMatchElectron"] = 0;
	    treeMap["trailMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to a photon

	else if(et_eta_match && abs(it->pdgId) == 11) {

	  double minDR = 999.0;
	  double dR;

	  for(vector<susy::Particle>::iterator oit = event.genParticles.begin(); oit != event.genParticles.end(); oit++) {

	    if(oit->momentum.Pt() < 2.0) continue;
	    if(oit == it) continue;
	    if(oit->momentum.M() > 10) continue;
	    if(abs(it->pdgId) == 12 || abs(it->pdgId) == 14 || abs(it->pdgId) == 16) continue; // skip neutrinos

	    dR = deltaR(it->momentum, oit->momentum);

	    if(dR < minDR) minDR = dR;

	  } // other gen particle loop

	  if(fromBoson &&
	     (fabs(photons[1]->momentum.Pt() - it->momentum.Pt()) / it->momentum.Pt()) < 0.1 &&
	     minDR > 0.2 &&
	     fabs(photons[1]->momentum.Eta() - it->momentum.Eta()) < 0.005 &&
	     deltaR(it->momentum, photons[1]->momentum) < 0.04) {
	    treeMap["trailMatchGamma"] = 0;
	    treeMap["trailMatchElectron"] = 1;
	    treeMap["trailMatchJet"] = 0;
	    foundMatch = true;
	    break;
	  }
	  else {
	    treeMap["trailMatchGamma"] = 0;
	    treeMap["trailMatchElectron"] = 0;
	    treeMap["trailMatchJet"] = 1;
	    foundMatch = true;
	    break;
	  }

	} // if matched to an electron

	if(foundMatch) break;

      } // gen particles loop

      if(!foundMatch) {
	treeMap["trailMatchGamma"] = 0;
	treeMap["trailMatchElectron"] = 0;
	treeMap["trailMatchJet"] = 1;
      }

    } // if nphotons > 1

  } // ifMC

  treeMap["ele1_pt"] = (tightEles.size() > 0) ? tightEles[0]->momentum.Pt() : -1.;
  treeMap["ele1_phi"] = (tightEles.size() > 0) ? tightEles[0]->momentum.Phi() : -100.;
  treeMap["ele1_eta"] = (tightEles.size() > 0) ? tightEles[0]->momentum.Eta() : -100.;
  treeMap["ele1_mvaTrigV0"] = (tightEles.size() > 0) ? tightEles[0]->mvaTrig : -100.;

  treeMap["ele2_pt"] = (tightEles.size() > 1) ? tightEles[1]->momentum.Pt() : -1.;
  treeMap["ele2_phi"] = (tightEles.size() > 1) ? tightEles[1]->momentum.Phi() : -100.;
  treeMap["ele2_eta"] = (tightEles.size() > 1) ? tightEles[1]->momentum.Eta() : -100.;
  treeMap["ele2_mvaTrigV0"] = (tightEles.size() > 1) ? tightEles[1]->mvaTrig : -100.;

  if(tightEles.size() > 0) {
    float ele_eta = fabs(event.superClusters[tightEles[0]->superClusterIndex].position.Eta());
    float ea;
    if(ele_eta < 1.0) ea = 0.13;
    else if(ele_eta < 1.479) ea = 0.14;
    else if(ele_eta < 2.0) ea = 0.07;
    else if(ele_eta < 2.2) ea = 0.09;
    else if(ele_eta < 2.3) ea = 0.11;
    else if(ele_eta < 2.4) ea = 0.11;
    else ea = 0.14;
      
    float ele_iso = max(0., (double)(tightEles[0]->photonIso + tightEles[0]->neutralHadronIso - event.rho25*ea));
    ele_iso += tightEles[0]->chargedHadronIso;
    
    treeMap["ele1_relIso"] = ele_iso / tightEles[0]->momentum.Pt();
  }
  else treeMap["ele1_relIso"] = -100.;

  if(tightEles.size() > 1) {
    float ele_eta = fabs(event.superClusters[tightEles[1]->superClusterIndex].position.Eta());
    float ea;
    if(ele_eta < 1.0) ea = 0.13;
    else if(ele_eta < 1.479) ea = 0.14;
    else if(ele_eta < 2.0) ea = 0.07;
    else if(ele_eta < 2.2) ea = 0.09;
    else if(ele_eta < 2.3) ea = 0.11;
    else if(ele_eta < 2.4) ea = 0.11;
    else ea = 0.14;
      
    float ele_iso = max(0., (double)(tightEles[1]->photonIso + tightEles[1]->neutralHadronIso - event.rho25*ea));
    ele_iso += tightEles[1]->chargedHadronIso;
    
    treeMap["ele1_relIso"] = ele_iso / tightEles[1]->momentum.Pt();
  }
  else treeMap["ele2_relIso"] = -100.;

  treeMap["muon1_pt"] = (tightMuons.size() > 0) ? tightMuons[0]->momentum.Pt() : -1.;
  treeMap["muon1_phi"] = (tightMuons.size() > 0) ? tightMuons[0]->momentum.Phi() : -100.;
  treeMap["muon1_eta"] = (tightMuons.size() > 0) ? tightMuons[0]->momentum.Eta() : -100.;

  treeMap["muon2_pt"] = (tightMuons.size() > 1) ? tightMuons[1]->momentum.Pt() : -1.;
  treeMap["muon2_phi"] = (tightMuons.size() > 1) ? tightMuons[1]->momentum.Phi() : -100.;
  treeMap["muon2_eta"] = (tightMuons.size() > 1) ? tightMuons[1]->momentum.Eta() : -100.;

  if(tightMuons.size() > 0) {
    float mu_iso = max(0., (double)(tightMuons[0]->sumNeutralHadronEt04 + tightMuons[0]->sumPhotonEt04 - 0.5*(tightMuons[0]->sumPUPt04)));
    mu_iso += tightMuons[0]->sumChargedHadronPt04;
    float mu_pt = tightMuons[0]->momentum.Pt();
    treeMap["muon_relIso"] = mu_iso / mu_pt;
  }
  else treeMap["muon2_relIso"] = -10.;

  if(tightMuons.size() > 1) {
    float mu_iso = max(0., (double)(tightMuons[1]->sumNeutralHadronEt04 + tightMuons[1]->sumPhotonEt04 - 0.5*(tightMuons[1]->sumPUPt04)));
    mu_iso += tightMuons[1]->sumChargedHadronPt04;
    float mu_pt = tightMuons[1]->momentum.Pt();
    treeMap["muon_relIso"] = mu_iso / mu_pt;
  }
  else treeMap["muon2_relIso"] = -10.;
  
  TLorentzVector z_vector = (tightEles.size() == 2) ? tightEles[0]->momentum + tightEles[1]->momentum : tightMuons[0]->momentum + tightMuons[1]->momentum;

  treeMap["z_mass"] = z_vector.M();
  treeMap["z_pt"] = z_vector.Pt();
  treeMap["z_eta"] = z_vector.Eta();
  treeMap["z_phi"] = z_vector.Phi();

  treeMap["leadPhotonEt"] = (photons.size() > 0) ? photons[0]->momentum.Et() : -1.;
  treeMap["leadPhotonEta"] = (photons.size() > 0) ? photons[0]->momentum.Eta() : -100.;
  treeMap["leadPhotonPhi"] = (photons.size() > 0) ? photons[0]->momentum.Phi() : -100.;
  treeMap["leadChargedHadronIso"] = (photons.size() > 0) ? chargedHadronIso_corrected(*photons[0], event.rho25) : -100.;
  treeMap["leadNeutralHadronIso"] = (photons.size() > 0) ? neutralHadronIso_corrected(*photons[0], event.rho25) - 0.04*(photons[0]->momentum.Pt()) : -100.;
  treeMap["leadPhotonIso"] = (photons.size() > 0) ? photonIso_corrected(*photons[0], event.rho25) - 0.005*(photons[0]->momentum.Pt()) : -100.;
  treeMap["leadSigmaIetaIeta"] = (photons.size() > 0) ? photons[0]->sigmaIetaIeta : -100.;
  treeMap["lead_nPixelSeeds"] = (photons.size() > 0) ? photons[0]->nPixelSeeds : -10.;
  treeMap["leadMVAregEnergy"] = (photons.size() > 0) ? photons[0]->MVAregEnergy : -10.;
  treeMap["leadMVAregErr"] = (photons.size() > 0) ? photons[0]->MVAregErr : -10.;
  
  treeMap["trailPhotonEt"] = (photons.size() > 1) ? photons[1]->momentum.Et() : -1.;
  treeMap["trailPhotonEta"] = (photons.size() > 1) ? photons[1]->momentum.Eta() : -100.;
  treeMap["trailPhotonPhi"] = (photons.size() > 1) ? photons[1]->momentum.Phi() : -100.;
  treeMap["trail_nPixelSeeds"] = (photons.size() > 1) ? photons[1]->nPixelSeeds : -100.;
  treeMap["trailChargedHadronIso"] = (photons.size() > 1) ? chargedHadronIso_corrected(*photons[1], event.rho25) : -100.;
  treeMap["trailNeutralHadronIso"] = (photons.size() > 1) ? neutralHadronIso_corrected(*photons[1], event.rho25) - 0.04*(photons[1]->momentum.Pt()) : -100.;
  treeMap["trailPhotonIso"] = (photons.size() > 1) ? photonIso_corrected(*photons[1], event.rho25) - 0.005*(photons[1]->momentum.Pt()) : -100.;
  treeMap["trailSigmaIetaIeta"] = (photons.size() > 1) ? photons[1]->sigmaIetaIeta : -100.;
  treeMap["trailMVAregEnergy"] = (photons.size() > 1) ? photons[1]->MVAregEnergy : -10.;
  treeMap["trailMVAregErr"] = (photons.size() > 1) ? photons[1]->MVAregErr : -10.;
  
  treeMap["diEMpT"] = (photons.size() > 1) ? (photons[0]->momentum + photons[1]->momentum).Pt() : -1.;
  treeMap["photon_invmass"] = (photons.size() > 1) ? (photons[0]->momentum + photons[1]->momentum).M() : -10.;
  
  float diJetPt, lead_matched_jetpt, trail_matched_jetpt;
  bool matchingWorked = (photons.size() > 1 && GetDiJetPt(event, photons, diJetPt, lead_matched_jetpt, trail_matched_jetpt));
  treeMap["diJetPt"] = matchingWorked ? diJetPt : -1.;
  
  float dEta = (photons.size() > 1) ? photons[0]->momentum.Eta() - photons[1]->momentum.Eta() : -100;
  float photon_dPhi = (photons.size() > 1) ? TVector2::Phi_mpi_pi(photons[0]->momentum.Phi() - photons[1]->momentum.Phi()) : -100;
  float photon_dR = (photons.size() > 1) ? sqrt(dEta*dEta + photon_dPhi*photon_dPhi) : -100;
  treeMap["photon_dPhi"] = photon_dPhi;
  treeMap["photon_dR"] = photon_dR;
  
  if(!isMC) {
    treeMap["runNumber"] = event.runNumber;
    treeMap["eventNumber"] = event.eventNumber;
    treeMap["lumiBlock"] = event.luminosityBlockNumber;
    treeMap["jentry"] = jentry;
  }

}
