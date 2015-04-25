#include "math.h"
#include "TMath.h"

#include <vector>

#include "../src/SusyEvent.h"

using namespace std;

bool isIsolatedMuon(susy::Muon mu) {
  float mu_iso = max(0., (mu.sumNeutralHadronEt04 + mu.sumPhotonEt04 - 0.5*(mu.sumPUPt04)));
  mu_iso += mu.sumChargedHadronPt04;

  return (mu_iso / mu.momentum.Pt() < 0.12);
}

bool isAntiIsolatedMuon(susy::Muon mu) {
  float mu_iso = max(0., (mu.sumNeutralHadronEt04 + mu.sumPhotonEt04 - 0.5*(mu.sumPUPt04)));
  mu_iso += mu.sumChargedHadronPt04;

  return (mu_iso / mu.momentum.Pt() > 0.25 * 0.9 && mu_iso / mu.momentum.Pt() < 1.0); // hard-cut on -10%, and require 0.25 in CreateHistograms.h
}

// This doesn't check for relIso! Muons passing this are either signal or QCD muon candidates
bool isTightMuon(susy::Muon mu, vector<susy::Track> tracks, double d0, double dz) {

  float mu_iso = max(0., (mu.sumNeutralHadronEt04 + mu.sumPhotonEt04 - 0.5*(mu.sumPUPt04)));
  mu_iso += mu.sumChargedHadronPt04;

  bool hasTracks = (int)mu.trackIndex < (int)tracks.size() && 
    (int)mu.standAloneTrackIndex < (int)tracks.size() && 
    (int)mu.combinedTrackIndex < (int)tracks.size() && 
    (int)mu.bestTrackIndex() < (int)tracks.size() && 
    (int)mu.bestTrackIndex() >= 0;

  if(!hasTracks) return false;
  
  bool passes = mu.isGlobalMuon() && 
    mu.isPFMuon() && 
    tracks[mu.combinedTrackIndex].normChi2() < 10. && 
    mu.nValidMuonHits > 0 && 
    mu.nMatchedStations > 1 &&
    fabs(d0) < 0.2 &&
    fabs(dz) < 0.5 &&
    tracks[mu.trackIndex].numberOfValidPixelHits > 0 && 
    (mu.nPixelLayersWithMeasurement + mu.nStripLayersWithMeasurement) > 5 && 
    mu.momentum.Pt() > 30. &&
    fabs(mu.momentum.Eta()) < 2.1;
  
  return passes;
  
}

bool isVetoMuon(susy::Muon mu) {

  float mu_iso = max(0., (mu.sumNeutralHadronEt04 + mu.sumPhotonEt04 - 0.5*(mu.sumPUPt04)));
  mu_iso += mu.sumChargedHadronPt04;

  bool passes = (mu.isPFMuon() &&
		 (mu.isGlobalMuon() || mu.isTrackerMuon()) &&
		 mu.momentum.Pt() > 10. && // (ttH(bb) for now)
		 fabs(mu.momentum.Eta()) < 2.5 &&
		 mu_iso / mu.momentum.Pt() < 0.2);

  return passes;
}

