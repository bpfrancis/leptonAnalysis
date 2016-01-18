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
#include <iomanip>
#include <fstream>

#include <map>
#include <vector>
#include <stdio.h>

using namespace std;

void remakeANTables() {

  TString backgrounds[10] = {"qcd", "ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};
  TString bkgLabels[10] = {"QCD", "t#bar{t} + Jets", "W + Jets", "Z/#gamma* + Jets", "Single top", 
			   "VV, V#gamma", "do not draw", "t#bar{t} + V", "do not draw", "t#bar{t} + #gamma"};

  TString systematics[15] = {"btagWeight", "puWeight", "JEC", "eleSF", "muonSF", "photonSF", "topPt",
                             "scale_tt", "scale_V", "scale_VV",
                             "pdf_gg", "pdf_qq", "pdf_gq",
                             "ele_qcdDef", "muon_qcdDef"};
  
  TString latexLabels[14] = {"QCD", "$t\\bar{t}$ + jets", "$W$ + jets", "$Z$ + jets", "Single $t$",
			     "Diboson", "V$\\gamma$", "$t\\bar{t} + W$", "$t\\bar{t} + Z$", "$t\\bar{t} + \\gamma$",
			     "Total Background", "GMSB (460\\_175)", "GMSB (560\\_325)", "Data"};

  int bkgColors[10] = {kSpring-6, kGray, kOrange-3, kYellow, kRed, kViolet-2, kAzure-2, kCyan, kOrange-5, 8};

  TString channels[10] = {"ele_Any", "muon_Any",
			  "ele_SR1", "muon_SR1",
			  "ele_SR2", "muon_SR2",
			  "ele_CR1", "muon_CR1",
			  "ele_CR2", "muon_CR2"};

  TString latexChannels[10] = {"Pre-selection", "Pre-selection",
			       "SR1", "SR1",
			       "SR2", "SR2",
			       "CR1", "CR1",
			       "CR2", "CR2"};

  TFile * input = new TFile("/uscms_data/d2/bfrancis/stopBino/limits/CMSSW_7_1_5/src/LimitSetting/condor_submit/limitInputs_bjj.root", "READ");

  double value, error;

  vector< vector<double> > tableValues, tableErrors;
  tableValues.resize(10); // tableValues[channel][row]
  tableErrors.resize(10);
  for(int i = 0; i < 10; i++) {
    tableValues[i].resize(14, 0.); // 10 backgrounds, total, siga, sigb, data
    tableErrors[i].resize(14, 0.);
  }

  for(int i = 0; i < 10; i++) {

    TH1D * h = (TH1D*)input->Get(channels[i]+"/data_obs");
    value = h->IntegralAndError(0, -1, error);
    tableValues[i][13] = value;
    tableErrors[i][13] = error;

    for(int j = 0; j < 10; j++) {

      h = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]);
      value = h->IntegralAndError(0, -1, error);
      tableValues[i][j] = value;

      for(int k = 0; k < 15; k++) {
	TH1D * h_up = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]+"_"+systematics[k]+"Up");
	TH1D * h_down = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]+"_"+systematics[k]+"Down");

	if(!h_up || !h_down) continue;

	double value_up = fabs(h_up->Integral() - value);
	double value_down = fabs(h_down->Integral() - value);
	double value_avg = (value_up + value_down)/2.;
	error = sqrt(error*error + value_avg*value_avg);
      }

      tableErrors[i][j] = error;
      
      if(i > 1 && j == 0) continue; // don't add QCD to total for SR1/2
      
      if(value > 0) {
	tableValues[i][10] += value;
	tableErrors[i][10] = sqrt(tableErrors[i][10]*tableErrors[i][10] + error*error);
      }
    }

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175");
    value = h->IntegralAndError(0, -1, error);
    tableValues[i][11] = value;
    
    for(int k = 0; k < 15; k++) {
	TH1D * h_up = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175"+"_"+systematics[k]+"Up");
	TH1D * h_down = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175"+"_"+systematics[k]+"Down");

	if(!h_up || !h_down) continue;

	double value_up = fabs(h_up->Integral() - value);
	double value_down = fabs(h_down->Integral() - value);
	double value_avg = (value_up + value_down)/2.;
	error = sqrt(error*error + value_avg*value_avg);
    }

    tableErrors[i][11] = error;

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325");
    value = h->IntegralAndError(0, -1, error);
    tableValues[i][12] = value;

    for(int k = 0; k < 15; k++) {
	TH1D * h_up = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325"+"_"+systematics[k]+"Up");
	TH1D * h_down = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325"+"_"+systematics[k]+"Down");

	if(!h_up || !h_down) continue;

	double value_up = fabs(h_up->Integral() - value);
	double value_down = fabs(h_down->Integral() - value);
	double value_avg = (value_up + value_down)/2.;
	error = sqrt(error*error + value_avg*value_avg);
    }

    tableErrors[i][12] = error;
    
  }

  int channelsInTables[4][3] = {{0, 2, 4}, // ele any/sr1/sr2
				{1, 3, 5}, // mu  any/sr1/sr2
				{6, 8, -1}, // ele cr1/cr2
				{7, 9, -1}};

  TString tableHeaders[4] = {"Channel & Pre-selection & SR1 & SR2",
			     "Channel & Pre-selection & SR1 & SR2",
			     "Channel & CR1 & CR2",
			     "Channel & CR1 & CR2"};

  TString hline = "\\hline";
  TString double_hline = "\\hline\\hline";
  TString rowBreak = " \\\\";

  for(int iTable = 0; iTable < 4; iTable++) { // make 4 tables

    cout << endl << endl;

    cout << hline << endl << tableHeaders[iTable] << rowBreak << endl;
    cout << double_hline << endl;

    for(int iRow = 0; iRow < 14; iRow++) {

      cout << latexLabels[iRow];
      
      for(int jTable = 0; jTable < 3; jTable++) { // with up to 3 columns each

	if(channelsInTables[iTable][jTable] >= 0) {
	  int chan = channelsInTables[iTable][jTable];
	  double x = tableValues[chan][iRow];
	  double y = tableErrors[chan][iRow];

	  if(iRow == 13) { // data row
	    cout << " & " << setprecision(0) << fixed << x;
	    continue;
	  }
	  
	  if(x < 1.e-4) cout << " & --";
	  else {
	    if(x >= 20.0) cout << " & " << setprecision(0) << fixed << x << " $\\pm$ " << sqrt(x) << " $\\pm$ " << y;
	    else if(x >= 2.0) cout << " & " << setprecision(1) << fixed << x << " $\\pm$ " << sqrt(x) << " $\\pm$ " << y;
	    else cout << " & " << setprecision(2) << fixed << x << " $\\pm$ " << sqrt(x) << " $\\pm$ " << y;
	  }
	}

      } // for columns in a row
      cout << rowBreak << endl;
      if(iRow == 9 || iRow == 10 || iRow == 12) cout << double_hline << endl;
      else cout << hline << endl;
      
    } // for rows in the table

  } // for tables

  input->Close();

}

void remakePasTables() {

  const int nBkgs = 9;
  const int nLines = 11;
  const int nSystematics = 15;
  const int nChannels = 4;
  
  TString backgrounds[nBkgs] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};
  int bkgGroups[nBkgs] = {0, 1, 1, 2, 3, 4, 5, 5, 6};
  
  TString groupedLabels[nLines] = {"$t\\bar{t}$",
				   "$W/Z$ + jets",
				   "Single $t$",
				   "Diboson", 
				   "$V\\gamma$",
				   "$t\\bar{t} + W/Z$",
				   "$t\\bar{t} + \\gamma$",
				   "Total Background",
				   "GMSB (460\\_175)",
				   "GMSB (560\\_325)",
				   "Data"};

  const int index_total = 7;
  const int index_siga = 8;
  const int index_sigb = 9;
  const int index_data = 10;

  TString systematics[nSystematics] = {"btagWeight", "puWeight", "JEC", "eleSF", "muonSF", "photonSF", "topPt",
				       "scale_tt", "scale_V", "scale_VV",
				       "pdf_gg", "pdf_qq", "pdf_gq",
				       "ele_qcdDef", "muon_qcdDef"};
  
  TString channels[nChannels] = {"ele_SR1", "muon_SR1",
				 "ele_SR2", "muon_SR2"}

  TFile * input = new TFile("../save_paper/limitInputs_bjj.root", "READ");

  double value, error;

  vector< vector<double> > tableValues, tableErrors;
  tableValues.resize(nChannels); // tableValues[channel][row]
  tableErrors.resize(nChannels);
  for(int i = 0; i < nChannels; i++) {
    tableValues[i].resize(nLines, 0.); // 7 backgrounds, total, siga, sigb, data
    tableErrors[i].resize(nLines, 0.);
  }

  for(int i = 0; i < nChannels; i++) {

    TH1D * h = (TH1D*)input->Get(channels[i]+"/data_obs");
    value = h->IntegralAndError(0, -1, error);
    tableValues[i][index_data] = value;
    tableErrors[i][index_data] = error;

    int index_line = 0; // j iterates backgrounds, but only iterate index_line if bkgGroups[j] != bkgGroups[j-1]
    
    for(int j = 0; j < nBkgs; j++) {

      if(j > 0 && bkgGroups[j] != bkgGroups[j-1]) index_line++;
      
      h = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]);
      value = h->IntegralAndError(0, -1, error);
      tableValues[i][index_line] += value;

      for(int k = 0; k < nSystematics; k++) {
	TH1D * h_up = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]+"_"+systematics[k]+"Up");
	TH1D * h_down = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]+"_"+systematics[k]+"Down");

	if(!h_up || !h_down) continue;

	double value_up = fabs(h_up->Integral() - value);
	double value_down = fabs(h_down->Integral() - value);
	double value_avg = (value_up + value_down)/2.;
	error = sqrt(error*error + value_avg*value_avg);
      }

      tableErrors[i][index_line] = sqrt(tableErrors[i][index_line]*tableErrors[i][index_line] + error*error);
      
      if(value > 0) {
	tableValues[i][index_total] += value;
	tableErrors[i][index_total] = sqrt(tableErrors[i][index_total]*tableErrors[i][index_total] + error*error);
      }
      
    }

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175");
    value = h->IntegralAndError(0, -1, error);
    tableValues[i][index_siga] = value;
    
    for(int k = 0; k < 15; k++) {
	TH1D * h_up = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175"+"_"+systematics[k]+"Up");
	TH1D * h_down = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175"+"_"+systematics[k]+"Down");

	if(!h_up || !h_down) continue;

	double value_up = fabs(h_up->Integral() - value);
	double value_down = fabs(h_down->Integral() - value);
	double value_avg = (value_up + value_down)/2.;
	error = sqrt(error*error + value_avg*value_avg);
    }

    tableErrors[i][index_siga] = error;

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325");
    value = h->IntegralAndError(0, -1, error);
    tableValues[i][index_sigb] = value;

    for(int k = 0; k < 15; k++) {
	TH1D * h_up = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325"+"_"+systematics[k]+"Up");
	TH1D * h_down = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325"+"_"+systematics[k]+"Down");

	if(!h_up || !h_down) continue;

	double value_up = fabs(h_up->Integral() - value);
	double value_down = fabs(h_down->Integral() - value);
	double value_avg = (value_up + value_down)/2.;
	error = sqrt(error*error + value_avg*value_avg);
    }

    tableErrors[i][index_sigb] = error;
    
  }

  TString tableHeader = "\\multirow{2}{*}{Channel} & \\multicolumn{2}{c|}{SR1} & \\multicolumn{2}{c|}{SR2}";
  TString channelHeader = " & e & $\\mu$ & e & $\\mu$";
  
  TString hline = "\\hline";
  TString double_hline = "\\hline\\hline";
  TString rowBreak = " \\\\";

  cout << endl << endl;

  cout << double_hline << endl;
  cout << tableHeader << rowBreak << endl;
  cout << "\\cline{2-5}" << endl;
  cout << channelHeader << rowBreak << endl;
  cout << hline << endl;

  for(int iRow = 0; iRow < nLines; iRow++) {

    cout << groupedLabels[iRow];

    for(int jColumn = 0; jColumn < nChannels; jColumn++) {
    
      double x = tableValues[jColumn][iRow];
      double y = tableErrors[jColumn][iRow];

      if(iRow == index_data) {
	cout << " & " << setprecision(0) << fixed << x;
	continue;
      }

      if(x < 1.e-4) cout << " & --";
      else {
	if(x >= 20.0) cout << " & " << setprecision(0) << fixed << x << " $\\pm$ " << sqrt(x) << " $\\pm$ " << y;
	else if(x >= 2.0) cout << " & " << setprecision(1) << fixed << x << " $\\pm$ " << sqrt(x) << " $\\pm$ " << y;
	else cout << " & " << setprecision(2) << fixed << x << " $\\pm$ " << sqrt(x) << " $\\pm$ " << y;
      }

    } // for columns in a row

    cout << rowBreak << endl;
    if(iRow == index_total-1 || iRow == index_sigb-1 || iRow == index_data-1) cout << hline << endl;
    if(iRow == index_data) cout << double_hline << endl;
    
  } // for rows in the table

  input->Close();

}

void remakeTables() {

  remakeANTables();
  remakePasTables();

}
