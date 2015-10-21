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

void metTable() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TString backgrounds[10] = {"qcd", "ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};
  TString bkgLabels[10] = {"QCD", "t#bar{t} + Jets", "W + Jets", "Z/#gamma* + Jets", "Single top", 
			   "VV, V#gamma", "do not draw", "t#bar{t} + V", "do not draw", "t#bar{t} + #gamma"};

  TString latexLabels[14] = {"QCD", "$t\\bar{t}$ + jets", "$W$ + jets", "$Z$ + jets", "Single $t$",
			     "Diboson", "V$\\gamma$", "$t\\bar{t} + W$", "$t\\bar{t} + Z$", "$t\\bar{t} + \\gamma$",
			     "Total Background", "GMSB (460\\_175)", "GMSB (560\\_325)", "Data"};

  const int nChannels = 6;
  TString channels[nChannels] = {"ele_Any", "muon_Any",
				 "ele_SR1", "muon_SR1",
				 "ele_SR2", "muon_SR2"};

  const int nMetBins = 5;
  double xbins_met[nMetBins+1] = {0, 25, 50, 100, 300, -1};
  double metBins[nMetBins+1];
  
  TFile * input = new TFile("../limitInputs_bjj.root", "READ");

  double value, error;

  vector< vector< vector<double> > > tableValues;
  vector< vector< vector<double> > > tableErrors;
  // tableValues[ele_SR1][ttgamma][met 0-25]
  
  tableValues.resize((unsigned int)nChannels);
  tableErrors.resize((unsigned int)nChannels);
  
  for(int i = 0; i < nChannels; i++) {
    tableValues[i].resize(14);
    tableErrors[i].resize(14);
  }
  
  for(int i = 0; i < nChannels; i++) {
    for(int j = 0; j < 14; j++) {
      tableValues[i][j].resize(nMetBins, 0.);
      tableErrors[i][j].resize(nMetBins, 0.);
    }
  }
  
  for(int i = 0; i < nChannels; i++) {

    // data
    TH1D * h = (TH1D*)input->Get(channels[i]+"/data_obs");

    for(int k = 0; k < nMetBins+1; k++) metBins[k] = h->GetXaxis()->FindBin(xbins_met[k]);
    
    for(int k = 0; k < nMetBins; k++) {
      value = h->IntegralAndError(metBins[k], metBins[k+1]-1, error);
      tableValues[i][13][k] = value;
      tableErrors[i][13][k] = error;
    }

    // backgrounds

    vector<double> total(nMetBins, 0.);
    vector<double> totalError(nMetBins, 0.);
    
    for(int j = 0; j < 10; j++) {

      h = (TH1D*)input->Get(channels[i]+"/"+backgrounds[j]);
      for(int k = 0; k < nMetBins; k++) {
	value = h->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	tableValues[i][j][k] = value;
	tableErrors[i][j][k] = error;

	if(i > 1 && j == 0) continue; // don't add QCD to total for SR1/2
	
	if(value > 0) {
	  total[k] += value;
	  totalError[k] = sqrt(totalError[k]*totalError[k] + error*error);
	  
	  tableValues[i][10][k] += value;
	  tableErrors[i][10][k] = sqrt(tableErrors[i][10][k]*tableErrors[i][10][k] + error*error);
	}

      }
	
    }

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175");
    for(int k = 0; k < nMetBins; k++) {
	value = h->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	tableValues[i][11][k] = value;
	tableErrors[i][11][k] = error;
    }
	
    h = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325");
    for(int k = 0; k < nMetBins; k++) {
	value = h->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	tableValues[i][12][k] = value;
	tableErrors[i][12][k] = error;
    }

  }

  // now set up the latex

  TString preamble = "\\begin{table}\\footnotesize\n\\begin{center}\n\\begin{tabular}{|l|c|c|c|c|c|}\n";
  TString postamble = "\\end{tabular}\n\\caption{durp}\n\\end{center}\n\\end{table}\n";

  TString tableHeader = "Channel & \\multicolumn{5}{|c|}{MET (GeV)}";
  TString binHeader = " & 0--25 & 25--50 & 50--100 & 100--300 & $\\geq$ 300";
  TString hline = "\\hline";
  TString double_hline = "\\hline\\hline";
  TString rowBreak = " \\\\";

  for(int iTable = 0; iTable < nChannels; iTable++) {

    cout << endl << endl;

    cout << channels[iTable] << endl << endl;

    cout << preamble;
    cout << hline << endl << tableHeader << rowBreak << endl;
    cout << hline << endl << binHeader << rowBreak << endl;
    cout << double_hline << endl;

    for(int iRow = 0; iRow < 14; iRow++) {

      cout << latexLabels[iRow];

      for(int jColumn = 0; jColumn < nMetBins; jColumn++) {
      
	double x = tableValues[iTable][iRow][jColumn];
	double y = tableErrors[iTable][iRow][jColumn];

	if(iRow == 13) { // data row
	  cout << " & " << setprecision(0) << fixed << x;
	  continue;
	}
	
	if(x < 1.e-4) cout << " & --";
	else {
	  if(x >= 20.0) cout << " & " << setprecision(0) << fixed << x << " $\\pm$ " << y;
	  else if(x >= 2.0) cout << " & " << setprecision(1) << fixed << x << " $\\pm$ " << y;
	  else cout << " & " << setprecision(2) << fixed << x << " $\\pm$ " << y;
	}

      } // for columns in a row
      cout << rowBreak << endl;
      if(iRow == 9 || iRow == 10 || iRow == 12) cout << double_hline << endl;
      else cout << hline << endl;
      
    } // for rows in the table

    cout << postamble;
    
  } // for tables


  input->Close();

}
