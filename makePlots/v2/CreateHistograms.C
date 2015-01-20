#include "CreateHistograms.h"

void CreateHistograms(int nPhotons_req = 1) {

  // negative cut means invert cut
  double sigmaIetaIetaCut = 0.012;

  // negative cut means no cut
  double metCut = -1.;

  // negative cut means no cut
  int nPhotons = 1;

  HistogramMaker * hMaker = new HistogramMaker("ele_bjj", "noSigmaIetaIeta", nPhotons, metCut, sigmaIetaIetaCut);
  hMaker->LoadLeptonSFs("../../data/lepton_SF_8TeV_53x_baseline.root");
  hMaker->LoadPhotonSFs("../../data/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root");

  const int nMetBins_0g = 17;
  Double_t xbins_met_0g[nMetBins_0g+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 300, 650};
  const int nKinematicBins_0g = 28;
  Double_t xbins_kinematic_0g[nKinematicBins_0g+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000};

  const int nMetBins_1g = 10;
  Double_t xbins_met_1g[nMetBins_1g+1] = {0, 10, 20, 30, 40, 50, 75, 100, 150, 300, 650};
  const int nKinematicBins_1g = 20;
  Double_t xbins_kinematic_1g[nKinematicBins_1g+1] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000};

  const int nMetBins_2g = 6;
  Double_t xbins_met_2g[nMetBins_2g+1] = {0, 20, 50, 75, 100, 150, 1000};
  const int nKinematicBins_2g = 10;
  Double_t xbins_kinematic_2g[nKinematicBins_2g+1] = {0, 25, 50, 100, 150, 200, 400, 600, 1000, 1500, 2000};

  if(nPhotons_req < 1) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("pfMET", nMetBins_0g, xbins_met_0g);
    hMaker->BookHistogram("HT", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("Njets", 20, 0., 20.);
    hMaker->BookHistogram("Nbtags", 20, 0., 20.);
    hMaker->BookHistogram("max_csv", 20, 0., 1.);
    hMaker->BookHistogram("submax_csv", 20, 0., 1.);
    hMaker->BookHistogram("HT_jets", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("hadronic_pt", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("jet1_pt", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("jet2_pt", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("jet3_pt", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("btag1_pt", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("w_mT", nKinematicBins_0g, xbins_kinematic_0g);    // 13
    hMaker->BookHistogram("m3", nKinematicBins_0g, xbins_kinematic_0g);
    hMaker->BookHistogram("ele_pt", nKinematicBins_0g, xbins_kinematic_0g);  // 15
    hMaker->BookHistogram("ele_eta", 60, -2.5, 2.5);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_0g, xbins_kinematic_0g); // 17
    hMaker->BookHistogram("muon_eta", 60, -2.5, 2.5);                  // 18

    hMaker->BookHistogram("nPV", 70, 0., 70.);

    hMaker->BookHistogram2D("Nphotons", "pfMET", 3, 0., 3., 20, 0., 350., 1.e-2, 1.e5);

    /*
      hMaker->BookHistogram2D("Njets", "Nbtags", 15, 0., 15., 7, 0., 7.);
      hMaker->BookHistogram2D("HT", "pfMET", 20, 0., 1200., 20, 0., 350.);
      hMaker->BookHistogram2D("w_mT", "Njets", 60, 0., 600., 15, 0., 15.);
      hMaker->BookHistogram2D("w_mT", "Nbtags", 60, 0., 600., 7, 0., 7.);
      hMaker->BookHistogram2D("w_mT", "pfMET", 60, 0., 600., 20, 0., 350.);
      hMaker->BookHistogram2D("w_mT", "HT", 60, 0., 600., 20, 0., 1200.);
    */
  }

  if(nPhotons_req == 1) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("pfMET", nMetBins_1g, xbins_met_1g);
    hMaker->BookHistogram("HT", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("Njets", 20, 0., 20.);
    hMaker->BookHistogram("Nbtags", 20, 0., 20.);
    hMaker->BookHistogram("max_csv", 20, 0., 1.);
    hMaker->BookHistogram("submax_csv", 20, 0., 1.);
    hMaker->BookHistogram("HT_jets", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("hadronic_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet2_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("jet3_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("btag1_pt", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("w_mT", nKinematicBins_1g, xbins_kinematic_1g);    // 13
    hMaker->BookHistogram("m3", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("ele_pt", nKinematicBins_1g, xbins_kinematic_1g);  // 15
    hMaker->BookHistogram("ele_eta", 60, 0., 2.5, true);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_1g, xbins_kinematic_1g); // 17
    hMaker->BookHistogram("muon_eta", 60, 0., 2.5, true);
    hMaker->BookHistogram("leadPhotonEt", nKinematicBins_1g, xbins_kinematic_1g); // 19
    hMaker->BookHistogram("leadPhotonEta", 40, 0., 1.5, true);                  // 20
    hMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159); // 21
    hMaker->BookHistogram("leadSigmaIetaIeta", 200, 0., 0.1);      // 22
    hMaker->BookHistogram("leadChargedHadronIso", 35, 0, 15.0);
    hMaker->BookHistogram("mLepGammaLead", nKinematicBins_1g, xbins_kinematic_1g);
    hMaker->BookHistogram("cosTheta_leadPhoton_l", 50, -1., 1.);
    hMaker->BookHistogram("cosTheta_leadPhoton_b_min", 50, -1., 1.);
    hMaker->BookHistogram("dR_leadPhoton_l", 120, 0., 6.);
    hMaker->BookHistogram("dR_leadPhoton_b_min", 120, 0., 6.);
    hMaker->BookHistogram("dEta_leadPhoton_l", 120, 0., 6.);
    hMaker->BookHistogram("dEta_leadPhoton_b_min", 120, 0., 6.);
    hMaker->BookHistogram("dPhi_leadPhoton_l", 120, 0., 6.2);
    hMaker->BookHistogram("dPhi_leadPhoton_b_min", 120, 0., 6.2);
    
    hMaker->BookHistogram2D("leadSigmaIetaIeta", "pfMET", 200, 0., 0.1, 20, 0., 350.);
    hMaker->BookHistogram2D("leadChargedHadronIso", "pfMET", 70, 0., 15., 20, 0., 350.);
  }
  
  if(nPhotons_req == 2) {
    hMaker->BookHistogram("Nphotons", 4, 0., 4.);
    hMaker->BookHistogram("pfMET", nMetBins_2g, xbins_met_2g);
    hMaker->BookHistogram("HT", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("Njets", 20, 0., 20.);
    hMaker->BookHistogram("Nbtags", 20, 0., 20.);
    hMaker->BookHistogram("max_csv", 20, 0., 1.);
    hMaker->BookHistogram("submax_csv", 20, 0., 1.);
    hMaker->BookHistogram("HT_jets", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("hadronic_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("jet1_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("jet2_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("jet3_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("btag1_pt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("w_mT", nKinematicBins_2g, xbins_kinematic_2g);    // 13
    hMaker->BookHistogram("m3", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("ele_pt", nKinematicBins_2g, xbins_kinematic_2g);  // 15
    hMaker->BookHistogram("ele_eta", 60, 0., 2.5, true);                   // 16
    hMaker->BookHistogram("muon_pt", nKinematicBins_2g, xbins_kinematic_2g); // 17
    hMaker->BookHistogram("muon_eta", 60, 0., 2.5, true);
    hMaker->BookHistogram("leadPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 19
    hMaker->BookHistogram("leadPhotonEta", 40, 0., 1.5, true);                  // 20
    hMaker->BookHistogram("leadPhotonPhi", 63, -3.14159, 3.14159);
    hMaker->BookHistogram("leadSigmaIetaIeta", 200, 0., 0.1);
    hMaker->BookHistogram("leadChargedHadronIso", 35, 0, 15.0);
    hMaker->BookHistogram("mLepGammaLead", nKinematicBins_2g, xbins_kinematic_2g);
        
    hMaker->BookHistogram2D("leadSigmaIetaIeta", "pfMET", 200, 0., 0.1, 20, 0., 350.);
    hMaker->BookHistogram2D("leadChargedHadronIso", "pfMET", 70, 0., 15., 20, 0., 350.);

    hMaker->BookHistogram("trailPhotonEt", nKinematicBins_2g, xbins_kinematic_2g); // 24
    hMaker->BookHistogram("trailPhotonPhi", 63, -3.14159, 3.14159);          // 25
    hMaker->BookHistogram("trailPhotonEta", 40, 0., 1.5, true);                  // 26
    hMaker->BookHistogram("trailSigmaIetaIeta", 200, 0, 0.1);                // 27
    hMaker->BookHistogram("trailChargedHadronIso", 35, 0, 15.0);
    hMaker->BookHistogram("diEMpT", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("diJetPt", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("photon_invmass", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("photon_dR", 50, 0., 5.);
    hMaker->BookHistogram("photon_dPhi", 35, 0., 3.14159);
    hMaker->BookHistogram("mLepGammaTrail", nKinematicBins_2g, xbins_kinematic_2g);
    hMaker->BookHistogram("mLepGammaGamma", nKinematicBins_2g, xbins_kinematic_2g);

    hMaker->BookHistogram("cosTheta_leadPhoton_l", 50, -1., 1.);
    hMaker->BookHistogram("cosTheta_leadPhoton_b_min", 50, -1., 1.);
    hMaker->BookHistogram("dR_leadPhoton_l", 120, 0., 6.);
    hMaker->BookHistogram("dR_leadPhoton_b_min", 120, 0., 6.);
    hMaker->BookHistogram("dEta_leadPhoton_l", 120, 0., 6.);
    hMaker->BookHistogram("dEta_leadPhoton_b_min", 120, 0., 6.);
    hMaker->BookHistogram("dPhi_leadPhoton_l", 120, 0., 6.2);
    hMaker->BookHistogram("dPhi_leadPhoton_b_min", 120, 0., 6.2);

    hMaker->BookHistogram("cosTheta_trailPhoton_l", 50, -1., 1.);
    hMaker->BookHistogram("cosTheta_trailPhoton_b_min", 50, -1., 1.);
    hMaker->BookHistogram("dR_trailPhoton_l", 120, 0., 6.);
    hMaker->BookHistogram("dR_trailPhoton_b_min", 120, 0., 6.);
    hMaker->BookHistogram("dEta_trailPhoton_l", 120, 0., 6.);
    hMaker->BookHistogram("dEta_trailPhoton_b_min", 120, 0., 6.);
    hMaker->BookHistogram("dPhi_trailPhoton_l", 120, 0., 6.2);
    hMaker->BookHistogram("dPhi_trailPhoton_b_min", 120, 0., 6.2);
  }



  hMaker->HistogramData("../inputs/SingleElectron.root");
  hMaker->HistogramQCD("../inputs/SingleElectron.root");
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttA_2to5.root", "ttA_2to5",
				0.9081 * 2, false, true);
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep",
				245.8 * 0.105, true, true);
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep",
				245.8 * 0.438, true, true);
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic",
				245.8 * 0.0457, true, true);

  hMaker->SetSigmaIetaIetaCut(-1. * sigmaIetaIetaCut);
  hMaker->HistogramData("../inputs/SingleElectron.root");
  hMaker->HistogramQCD("../inputs/SingleElectron.root");
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttA_2to5.root", "ttA_2to5",
				0.9081 * 2, false, true);
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttJetsFullLep.root", "ttJetsFullLep",
				245.8 * 0.105, true, true);
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttJetsSemiLep.root", "ttJetsSemiLep",
				245.8 * 0.438, true, true);
  hMaker->HistogramMCBackground("../inputs/signal_contamination_ttJetsHadronic.root", "ttJetsHadronic",
				245.8 * 0.0457, true, true);

  delete hMaker;

}
