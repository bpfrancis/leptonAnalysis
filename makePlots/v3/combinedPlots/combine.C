#include <vector>

using namespace std;

const int numSpecialEleNames = 10;
TString specialEleNames[numSpecialEleNames] = {"eleSFUp", "eleSFDown",
					       "ele_qcdDefUp", "ele_qcdDefDown",
					       "userSystA_eleUp", "userSystA_eleDown",
					       "userSystB_eleUp", "userSystB_eleDown",
					       "userSystC_eleUp", "userSystC_eleDown"};

const int numSpecialMuonNames = 10;
TString specialMuonNames[numSpecialMuonNames] = {"muonSFUp", "muonSFDown",
						 "muon_qcdDefUp", "muon_qcdDefDown",
						 "userSystA_muonUp", "userSystA_muonDown",
						 "userSystB_muonUp", "userSystB_muonDown",
						 "userSystC_muonUp", "userSystC_muonDown"};

TH1D * DivideByBinWidth(TH1D * h) {

  for(Int_t i = 0; i < h->GetNbinsX(); i++) {
    Double_t val = h->GetBinContent(i+1);
    Double_t err = h->GetBinError(i+1);
    Double_t width = h->GetBinWidth(i+1);

    h->SetBinContent(i+1, val / width);
    h->SetBinError(i+1, err / width);
  }

  return h;
}


void addChannels() {

  TString channelNames[4] = {"SR1", "SR2", "CR1", "CR2"};
  
  TFile * input_bjj = new TFile("../save_paper/limitInputs_bjj.root", "READ");

  TFile * fCombined = new TFile("limitInputs_combined_bjj.root", "RECREATE");
  
  for(int channel = 0; channel < 4; channel++) {

    TDirectoryFile * eleDir = (TDirectoryFile*)input_bjj->Get("ele_" + channelNames[channel]);
    TDirectoryFile * muonDir = (TDirectoryFile*)input_bjj->Get("muon_" + channelNames[channel]);
    fCombined->mkdir("comb_" + channelNames[channel]);
    
    TObject * obj;
    TKey * key;
    
    TH1D * hEle;
    TH1D * hMuon;

    TIter next(eleDir->GetListOfKeys());
    
    while((key = (TKey*)next())) {

      TString name = key->GetName();
      
      if(!name.Contains("userSyst") && !name.Contains("signal_") && key->GetCycle() != 2) continue;
      
      TString eleName = name;
      TString muonName = name;

      obj = eleDir->Get(name);
      if(!(obj->InheritsFrom("TH1D"))) {
	cout << "Skipping (" << obj->ClassName() << ") " << name << endl;
	continue;
      }

      if(name.Contains("_stat_")) continue;

      if(name.Contains("signal_")) {
	if(!(name.Contains("signal_mst_460_m1_175")) && !(name.Contains("signal_mst_560_m1_325"))) continue;
      }

      hEle = (TH1D*)eleDir->Get(eleName);
      hMuon = (TH1D*)muonDir->Get(muonName);

      // if a systematic only exists for one channel, then get the central value for the other channel
      for(int i = 0; i < numSpecialEleNames; i++) {
	if(name.Contains(specialEleNames[i])) {
	  muonName = muonName.ReplaceAll("_" + specialEleNames[i], "");
	  hMuon = (TH1D*)muonDir->Get(muonName);
	}
      }

      if(!hEle) {
	cout << "Did not find ele version of: " << name << endl;
	continue;
      }
      if(!hMuon) {
	cout << "Did not find muon version of: " << name << endl;
	continue;
      }

      if(!muonName.Contains("qcd")) hEle->Add(hMuon);
      fCombined->cd("comb_" + channelNames[channel]);
      hEle->Write(eleName);
      
    }

    TIter next_muon(muonDir->GetListOfKeys());
    
    while((key = (TKey*)next_muon())) {

      TString name = key->GetName();
      
      if(!name.Contains("userSyst") && !name.Contains("signal_") && key->GetCycle() != 2) continue;
      
      TString eleName = name;
      TString muonName = name;
      
      bool hasSpecialMuonName = false;
      for(int i = 0; i < numSpecialMuonNames; i++) {
	if(name.Contains(specialMuonNames[i])) {
	  hasSpecialMuonName = true;
	  break;
	}
      }

      if(!hasSpecialMuonName) continue;

      obj = muonDir->Get(name);
      if(!(obj->InheritsFrom("TH1D"))) {
	cout << "Skipping (" << obj->ClassName() << ") " << name << endl;
	continue;
      }

      if(name.Contains("_stat_")) continue;

      if(name.Contains("signal_")) {
	if(!(name.Contains("signal_mst_460_m1_175")) && !(name.Contains("signal_mst_560_m1_325"))) continue;
      }

      eleName = eleName.ReplaceAll("_" + specialMuonNames[i], "");
      
      TString newEleName = name.ReplaceAll("_" + specialMuonNames[i], "");
      hEle = (TH1D*)eleDir->Get(eleName);
      hMuon = (TH1D*)muonDir->Get(muonName);

      if(!hEle) {
	cout << "Did not find ele version of: " << name << endl;
	continue;
      }
      if(!hMuon) {
	cout << "Did not find muon version of: " << name << endl;
	continue;
      }

      hMuon->Add(hEle);
      fCombined->cd("comb_" + channelNames[channel]);
      hMuon->Write(muonName);
      
    }

    fCombined->Write();

  }
    
  fCombined->Close();
  
  input_bjj->Close();
  //input_jjj->Close();
  
}

void calculateErrors() {

  Double_t xbins_met_2g[5] = {0, 20, 50, 100, 800};
  
  TString channelNames[4] = {"SR1", "SR2", "CR1", "CR2"};

  bool needsQCD[4] = {false, false, true, false};

  const int numBackgrounds = 10;
  TString backgrounds[numBackgrounds] = {"qcd", "ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};
  unsigned int layers[numBackgrounds] = {0, 1, 2, 3, 4, 5, 5, 6, 6, 7};

  const int numSaveNames = 8;
  TString saveNames[numSaveNames] = {"qcd", "ttjets", "wjets", "zjets", "singleTop", "vv", "ttV", "ttgamma"};

  const int numSystematics = 21;
  TString systematicNames[numSystematics] = {"btagWeight", "puWeight", "topPt",
					     "scale_tt", "scale_V", "scale_VV",
					     "pdf_gg", "pdf_qq", "pdf_qg",
					     "JEC", "eleSF", "muonSF", "photonSF",
					     "ele_qcdDef", "muon_qcdDef",
					     "userSystA_ele", "userSystA_muon",
					     "userSystB_ele", "userSystB_muon",
					     "userSystC_ele", "userSystC_muon"};

  TFile * input = new TFile("limitInputs_combined_bjj.root", "READ");

  TFile * output = new TFile("combined.root", "RECREATE");

  TH1D * diff_CR1;
  
  for(int channel = 0; channel < 4; channel++) {

    TDirectoryFile * dir = (TDirectoryFile*)input->Get("comb_" + channelNames[channel]);

    TH1D * data = (TH1D*)dir->Get("data_obs");
    TH1D * siga = (TH1D*)dir->Get("signal_mst_460_m1_175");
    TH1D * sigb = (TH1D*)dir->Get("signal_mst_560_m1_325");
    
    TH1D * central;
    TH1D * up;
    TH1D * down;

    vector<TH1D*> layerHistograms;    // key [iBkg]
    vector<TH1D*> systHistogramsUp;   // key [iSyst]
    vector<TH1D*> systHistogramsDown; // key [iSyst]

    vector<TH1D*> layerHistograms_systUp;   // key [iBkg] out of 10! not 8
    vector<TH1D*> layerHistograms_systDown; // key [iBkg] out of 10! not 8
    
    TH1D * central;
    TH1D * h_syst_up;
    TH1D * h_syst_down;

    TH1D * h_error_syst;

    output->mkdir(channelNames[channel]);
    
    for(int iBkg = 0; iBkg < numBackgrounds; iBkg++) {

      if(!needsQCD[channel] && iBkg == 0) continue;
      
      central = (TH1D*)dir->Get(backgrounds[iBkg]);

      layerHistograms_systUp.push_back((TH1D*)central->Clone("central_"+backgrounds[iBkg]+"_systUp"));
      layerHistograms_systUp.back()->Reset();
      
      layerHistograms_systDown.push_back((TH1D*)central->Clone("central_"+backgrounds[iBkg]+"_systDown"));
      layerHistograms_systDown.back()->Reset();
      
      // if first layer (qcd ttjets)
      if((needsQCD[channel] && iBkg == 0) ||
	 (!needsQCD[channel] && iBkg == 1)) {
	
	layerHistograms.push_back((TH1D*)central->Clone("central_"+backgrounds[iBkg]));

	// clone an empty h_error_syst for later
	h_error_syst = (TH1D*)central->Clone("error_syst");
	h_error_syst->Reset();
      }

      if((needsQCD[channel] && iBkg > 0) ||
	 (!needsQCD[channel] && iBkg > 1)) {

	for(unsigned int i = 0; i < layerHistograms.size(); i++) layerHistograms[i]->Add(central);

	// if this is a new layer, store it
	// if not, then we've already added the second hist into its layer
	if(layers[iBkg] != layers[iBkg-1]) layerHistograms.push_back((TH1D*)central->Clone("central_"+backgrounds[iBkg]));

      }

    } // for iBkg

    // now set up systematics histograms
    // here we just need to claim stuff for each "_systName" rather than all of "bkg_systName"
    for(int iSyst = 0; iSyst < numSystematics; iSyst++) {
      systHistogramsUp.push_back((TH1D*)layerHistograms[0]->Clone("systUp_" + systematicNames[iSyst]));
      systHistogramsDown.push_back((TH1D*)layerHistograms[0]->Clone("systUp_" + systematicNames[iSyst]));

      systHistogramsUp.back()->Reset();
      systHistogramsDown.back()->Reset();
    }

    // now add up the systematics
    for(int iBkg = 0; iBkg < numBackgrounds; iBkg++) {

      if(!needsQCD[channel] && iBkg == 0) continue;
      
      for(int iSyst = 0; iSyst < numSystematics; iSyst++) {

	if(!needsQCD[channel] && systematicNames[iSyst].Contains("qcd")) continue;
	
	h_syst_up = (TH1D*)dir->Get(backgrounds[iBkg] + "_" + systematicNames[iSyst] + "Up");
	h_syst_down = (TH1D*)dir->Get(backgrounds[iBkg] + "_" + systematicNames[iSyst] + "Down");

	if(!h_syst_up || !h_syst_down) continue;

	systHistogramsUp[iSyst]->Add(h_syst_up);
	systHistogramsDown[iSyst]->Add(h_syst_down);
	
      } // for syst
    } // for bkg

    // add up systematics keyed for each bkg (ie ttjets_systUp)
    for(int iBkg = 0; iBkg < layerHistograms_systUp.size(); iBkg++) {
      for(int iSyst = 0; iSyst < numSystematics; iSyst++) {
	TH1D * h_temp_up = (TH1D*)dir->Get(backgrounds[iBkg] + "_" + systematicNames[iSyst] + "Up");
	TH1D * h_temp_down = (TH1D*)dir->Get(backgrounds[iBkg] + "_" + systematicNames[iSyst] + "Down");

	if(!h_temp_up || !h_temp_down) continue;

	//cout << "Durp adding for " << backgrounds[iBkg] << " " << systematicNames[iSyst] << endl;
	layerHistograms_systUp[iBkg]->Add(h_temp_up);
	layerHistograms_systDown[iBkg]->Add(h_temp_down);
      }
    }

    /*
    // print out some info to check
    cout << endl << "Channel: " << channelNames[channel] << endl;
    cout << "Data: " << data->Integral() << endl;
    cout << "============" << endl;
    for(unsigned int i = 0; i < layerHistograms.size(); i++) {
      for(int jBin = 0; jBin < layerHistograms[i]->GetNbinsX(); jBin++) {
	cout << saveNames[i] << " bin " << jBin << " -- " << layerHistograms[i]->GetBinContent(jBin+1);
	cout << " +/- " << layerHistograms[i]->GetBinError(jBin+1) << endl;
      }
    }
    cout << "systematics: " << endl;
    for(unsigned int i = 0; i < layerHistograms_systUp.size(); i++) {
      cout << backgrounds[i] << " -- ";
      for(int jBin = 0; jBin < layerHistograms_systUp[i]->GetNbinsX(); jBin++) {
	cout << ", " << fabs(layerHistograms_systUp[i]->GetBinContent(jBin+1) - layerHistograms_systDown[i]->GetBinContent(jBin+1)) / 2.;
      }
    }
    cout << endl << "============" << endl;
    cout << "sig a: " << siga->Integral() << endl;
    cout << "sig b: " << sigb->Integral() << endl;
    cout << endl << endl;
    */
    
    // divide by bin widths

    data = (TH1D*)DivideByBinWidth(data);
    siga = (TH1D*)DivideByBinWidth(siga);
    sigb = (TH1D*)DivideByBinWidth(sigb);

    for(unsigned int i = 0; i < layerHistograms.size(); i++) layerHistograms[i] = (TH1D*)DivideByBinWidth(layerHistograms[i]);
    
    for(unsigned int i = 0; i < numSystematics; i++) {
      systHistogramsUp[i] = (TH1D*)DivideByBinWidth(systHistogramsUp[i]);
      systHistogramsDown[i] = (TH1D*)DivideByBinWidth(systHistogramsDown[i]);
    }

    // durp! calculate diff_CR1
    if(channelNames[channel] == "CR1") {
      diff_CR1 = (TH1D*)data->Clone("diff_CR1");

      for(int iBin = 0; iBin < diff_CR1->GetNbinsX(); iBin++) {
	double durp_diff = data->GetBinContent(iBin+1);
	double durp_bkg = layerHistograms[0]->GetBinContent(iBin+1);
	diff_CR1->SetBinContent(iBin+1, fabs(durp_diff - durp_bkg) / durp_bkg);
      }
    }
    
    // now calculate systematic errors based off of layerHistograms[0], which is the sum total background now
    for(int iSyst = 0; iSyst < numSystematics; iSyst++) {
      for(int iBin = 0; iBin < layerHistograms[0]->GetNbinsX(); iBin++) {
	
	double value_up = systHistogramsUp[iSyst]->GetBinContent(iBin+1);
	double value_down = systHistogramsDown[iSyst]->GetBinContent(iBin+1);
	
	double current = h_error_syst->GetBinError(iBin+1);
	double delta = fabs(value_up - value_down) / 2.;

	h_error_syst->SetBinContent(iBin+1, layerHistograms[0]->GetBinContent(iBin+1));
	h_error_syst->SetBinError(iBin+1, sqrt(current*current + delta*delta));
      }
    }

    if(channelNames[channel] == "CR2") {
      for(int iBin = 0; iBin < h_error_syst->GetNbinsX(); iBin++) {
	double old_ = h_error_syst->GetBinError(iBin+1);
	double new_ = diff_CR1->GetBinContent(iBin+1);
	h_error_syst->SetBinError(iBin+1, sqrt(old_*old_ + new_*new_));
      }
    }

    // rebin CR2
    if(channelNames[channel] == "CR2") {
      data = (TH1D*)data->Rebin(4, "data_reb", xbins_met_2g);
      siga = (TH1D*)siga->Rebin(4, "siga_reb", xbins_met_2g);
      sigb = (TH1D*)sigb->Rebin(4, "sigb_reb", xbins_met_2g);
      for(unsigned int i = 0; i < layerHistograms.size(); i++) layerHistograms[i] = (TH1D*)layerHistograms[i]->Rebin(4, "bkg_reb", xbins_met_2g);
      for(unsigned int i = 0; i < layerHistograms_systUp.size(); i++) {
	//cout << "Durp rebinning " << layerHistograms_systUp[i]->GetTitle() << endl;
	layerHistograms_systUp[i] = (TH1D*)layerHistograms_systUp[i]->Rebin(4, "bkg_systUp_reb", xbins_met_2g);
	layerHistograms_systDown[i] = (TH1D*)layerHistograms_systDown[i]->Rebin(4, "bkg_systDown_reb", xbins_met_2g);
      }
      for(unsigned int i = 0; i < systHistogramsUp.size(); i++) {
	systHistogramsUp[i] = (TH1D*)systHistogramsUp[i]->Rebin(4, "bkg_rebUp", xbins_met_2g);
	systHistogramsDown[i] = (TH1D*)systHistogramsDown[i]->Rebin(4, "bkg_rebDown", xbins_met_2g);
      }
      h_error_syst = (TH1D*)h_error_syst->Rebin(4, "error_syst_reb", xbins_met_2g);

    }

    // now work out the ratio plots
    TH1D * ratio = (TH1D*)data->Clone("ratio");
    ratio->Reset();
    ratio->SetTitle("Data / Background");
    for(int iBin = 0; iBin < ratio->GetNbinsX(); iBin++) {
      if(layerHistograms[0]->GetBinContent(iBin+1) == 0.) continue;
      ratio->SetBinContent(iBin+1, data->GetBinContent(iBin+1) / layerHistograms[0]->GetBinContent(iBin+1));
      ratio->SetBinError(iBin+1, data->GetBinError(iBin+1) / layerHistograms[0]->GetBinContent(iBin+1));
    }

    TH1D * ratio_stat = (TH1D*)layerHistograms[0]->Clone("ratio_stat");
    ratio_stat->Reset();
    for(int iBin = 0; iBin < ratio_stat->GetNbinsX(); iBin++) {
      ratio_stat->SetBinContent(iBin+1, 1.);
      if(layerHistograms[0]->GetBinContent(iBin+1) == 0.) ratio_stat->SetBinError(iBin+1, 0.);
      else ratio_stat->SetBinError(iBin+1, layerHistograms[0]->GetBinError(iBin+1) / layerHistograms[0]->GetBinContent(iBin+1));
    }

    TH1D * ratio_sys = (TH1D*)layerHistograms[0]->Clone("ratio_sys");
    ratio_sys->Reset();
    for(int iBin = 0; iBin < ratio_sys->GetNbinsX(); iBin++) {
      ratio_sys->SetBinContent(iBin+1, 1.);
      if(layerHistograms[0]->GetBinContent(iBin+1) == 0.) ratio_sys->SetBinError(iBin+1, 0.);
      else {
	ratio_sys->SetBinError(iBin+1, h_error_syst->GetBinError(iBin+1) / layerHistograms[0]->GetBinContent(iBin+1));
	if(channelNames[channel] == "CR2") ratio_sys->SetBinError(iBin+1, ratio_sys->GetBinError(iBin+1) * diff_CR1->GetBinContent(iBin+1));
      }
    }

    // check we got everything
    
    if((needsQCD[channel] && layerHistograms.size() != numSaveNames) ||
       (!needsQCD[channel] && layerHistograms.size() != numSaveNames - 1)) {
      cout << endl << "Didn't calculate something right, expected " << numSaveNames << " layers and got " << layerHistograms.size() << endl;
      return;
    }
    
    // now save the output
    output->cd(channelNames[channel]);

    data->Write("data");
    siga->Write("signal_a");
    sigb->Write("signal_b");
    h_error_syst->Write("error_syst");
    ratio->Write("ratio");
    ratio_stat->Write("ratio_stat");
    ratio_sys->Write("ratio_sys");
    
    for(unsigned int i = 0; i < layerHistograms.size(); i++) {
      if(needsQCD[channel]) layerHistograms[i]->Write(saveNames[i]);
      else layerHistograms[i]->Write(saveNames[i+1]);
    }
    
    layerHistograms[0]->Write("error_stat");

    
  } // for channels
  
  output->Close();
  input->Close();
}
		   
void combine() {
  //addChannels();
  calculateErrors();
}
