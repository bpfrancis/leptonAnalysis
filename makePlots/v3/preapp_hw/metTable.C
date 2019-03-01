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

  const int nBackgrounds = 10;
  TString backgrounds[nBackgrounds] = {"qcd", "ttjets", "wjets", "zjets", "singleTop", "diboson", "vgamma", "ttW", "ttZ", "ttgamma"};
  TString bkgLabels[nBackgrounds] = {"QCD", "t#bar{t} + Jets", "W + Jets", "Z/#gamma* + Jets", "Single top", 
				     "VV, V#gamma", "do not draw", "t#bar{t} + V", "do not draw", "t#bar{t} + #gamma"};

  TString latexLabels[14] = {"QCD", "$t\\bar{t}$ + jets", "$W$ + jets", "$Z$ + jets", "Single $t$",
			     "Diboson", "V$\\gamma$", "$t\\bar{t} + W$", "$t\\bar{t} + Z$", "$t\\bar{t} + \\gamma$",
			     "Total Background", "GMSB (460\\_175)", "GMSB (560\\_325)", "Data"};

  const int nChannels = 10;
  TString channels[nChannels] = {"ele_Any", "muon_Any",
				 "ele_SR1", "muon_SR1",
				 "ele_SR2", "muon_SR2",
				 "ele_CR1", "muon_CR1",
				 "ele_CR2", "muon_CR2"};

  const int numSystematics = 21;
  TString systematicNames[numSystematics] = {"btagWeight", "puWeight", "topPt",
					     "scale_tt", "scale_V", "scale_VV",
					     "pdf_gg", "pdf_qq", "pdf_qg",
					     "JEC", "eleSF", "muonSF", "photonSF",
					     "ele_qcdDef", "muon_qcdDef",
					     "userSystA_ele", "userSystA_muon",
					     "userSystB_ele", "userSystB_muon",
					     "userSystC_ele", "userSystC_muon"};
  
  const int nMetBins = 4;
  double xbins_met[nMetBins+1] = {0, 25, 50, 100, -1};
  double metBins[nMetBins+1];
  
  TFile * input = new TFile("../save_paper/limitInputs_bjj.root", "READ");

  double value, error;

  vector< vector< vector<double> > > tableValues;
  vector< vector< vector<double> > > tableErrors;
  vector< vector< vector<double> > > tableSystematics;
  // tableValues[ele_SR1][ttgamma][met 0-25]
  
  tableValues.resize((unsigned int)nChannels);
  tableErrors.resize((unsigned int)nChannels);
  tableSystematics.resize((unsigned int)nChannels);
  
  for(int i = 0; i < nChannels; i++) {
    tableValues[i].resize(14);
    tableErrors[i].resize(14);
    tableSystematics[i].resize(14);
  }
  
  for(int i = 0; i < nChannels; i++) {
    for(int j = 0; j < 14; j++) {
      tableValues[i][j].resize(nMetBins, 0.);
      tableErrors[i][j].resize(nMetBins, 0.);
      tableSystematics[i][j].resize(nMetBins, 0.);
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
      tableSystematics[i][13][k] = 0.;
    }

    // backgrounds (total is jBkg ~ 10)
    for(int jBkg = 0; jBkg < nBackgrounds; jBkg++) {

      h = (TH1D*)input->Get(channels[i]+"/"+backgrounds[jBkg]);
      for(int k = 0; k < nMetBins; k++) {
	value = h->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	tableValues[i][jBkg][k] = value;
	tableErrors[i][jBkg][k] = error;
	tableSystematics[i][jBkg][k] = 0.;

	if(channels[i] != "ele_CR1" && jBkg == 0) continue; // don't add QCD to total except for ele_CR1
	
	if(value > 0) {
	  tableValues[i][10][k] += value;
	  tableErrors[i][10][k] = sqrt(tableErrors[i][10][k]*tableErrors[i][10][k] + error*error);
	
	  for(int iSyst = 0; iSyst < numSystematics; iSyst++) {
	    TH1D * h_up = (TH1D*)input->Get(channels[i]+"/"+backgrounds[jBkg]+"_"+systematicNames[iSyst]+"Up");
	    TH1D * h_down = (TH1D*)input->Get(channels[i]+"/"+backgrounds[jBkg]+"_"+systematicNames[iSyst]+"Down");
	    
	    if(!h_up || !h_down) continue;
	    
	    double val_up = h_up->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	    double val_down = h_down->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	    double systNew = fabs(val_up - val_down) / 2.;
	    double systOld = tableSystematics[i][jBkg][k];
	    
	    tableSystematics[i][jBkg][k] = sqrt(systOld*systOld + systNew*systNew);
	    tableSystematics[i][10][k] = sqrt(tableSystematics[i][10][k]*tableSystematics[i][10][k] + systNew*systNew);
	  }

	} // if value > 0, track the total
	  
      }
	
    }

    h = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175");
    for(int k = 0; k < nMetBins; k++) {
	value = h->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	tableValues[i][11][k] = value;
	tableErrors[i][11][k] = error;

	for(int iSyst = 0; iSyst < numSystematics; iSyst++) {
	  TH1D * h_up = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175_"+systematicNames[iSyst]+"Up");
	  TH1D * h_down = (TH1D*)input->Get(channels[i]+"/signal_mst_460_m1_175_"+systematicNames[iSyst]+"Down");

	  if(!h_up || !h_down) continue;

	  double val_up = h_up->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	  double val_down = h_down->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	  double systNew = fabs(val_up - val_down) / 2.;
	  double systOld = tableSystematics[i][11][k];

	  tableSystematics[i][11][k] = sqrt(systOld*systOld + systNew*systNew);
	}
    }
	
    h = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325");
    for(int k = 0; k < nMetBins; k++) {
	value = h->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	tableValues[i][12][k] = value;
	tableErrors[i][12][k] = error;

	for(int iSyst = 0; iSyst < numSystematics; iSyst++) {
	  TH1D * h_up = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325_"+systematicNames[iSyst]+"Up");
	  TH1D * h_down = (TH1D*)input->Get(channels[i]+"/signal_mst_560_m1_325_"+systematicNames[iSyst]+"Down");

	  if(!h_up || !h_down) continue;

	  double val_up = h_up->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	  double val_down = h_down->IntegralAndError(metBins[k], metBins[k+1]-1, error);
	  double systNew = fabs(val_up - val_down) / 2.;
	  double systOld = tableSystematics[i][12][k];

	  tableSystematics[i][12][k] = sqrt(systOld*systOld + systNew*systNew);
	}
    }

  }

  // now set up the latex

  TString preamble = "\\begin{table}\\footnotesize\n\\begin{center}\n\\begin{tabular}{|l|c|c|c|c|}\n";
  TString postamble = "\\end{tabular}\n\\caption{durp}\n\\end{center}\n\\end{table}\n";

  TString tableHeader = "Channel & \\multicolumn{4}{|c|}{MET (GeV)}";
  TString binHeader = " & 0--25 & 25--50 & 50--100 & $\\geq$ 100";
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
	double z = tableSystematics[iTable][iRow][jColumn];

	if(iRow == 13) { // data row
	  cout << " & " << setprecision(0) << fixed << x;
	  continue;
	}

	if(x < 1.e-3) cout << " & --";
	else {
	  int prec = 0;

	  int multiplier = 1;
	  while(true) {
	    if((z == 0. && (int)(x*multiplier) > 0 && (int)(y*multiplier) > 0) ||
	       ((int)(x*multiplier) > 0 && (int)(y*multiplier) > 0 && (int)(z*multiplier) > 0)) break;
	    else {
	      multiplier *= 10;
	      prec++;
	    }
	  }

	  if(prec > 2) prec = 2;

	  cout << " & " << setprecision(prec) << fixed << x << " $\\pm$ " << y << " $\\pm$ " << z;
	}

      } // for columns in a row
      cout << rowBreak << endl;
      if(iRow == 9 || iRow == 10 || iRow == 12) cout << double_hline << endl;
      else cout << hline << endl;
      
    } // for rows in the table

    cout << postamble;
    
  } // for tables

  // combination tables
  for(int iTable = 0; iTable < nChannels - 1; iTable += 2) {

    cout << endl << endl;

    cout << "combined: " << channels[iTable] << " + " << channels[iTable+1] << endl;

    cout << preamble;
    cout << hline << endl << tableHeader << rowBreak << endl;
    cout << hline << endl << binHeader << rowBreak << endl;
    cout << double_hline << endl;

    for(int iRow = 0; iRow < 14; iRow++) {

      cout << latexLabels[iRow];

      for(int jColumn = 0; jColumn < nMetBins; jColumn++) {
      
	double x = tableValues[iTable][iRow][jColumn] + tableValues[iTable+1][iRow][jColumn];

	double ya = tableErrors[iTable][iRow][jColumn];
	double yb = tableErrors[iTable+1][iRow][jColumn];
	double y = sqrt(ya*ya + yb*yb);
	
	double za = tableSystematics[iTable][iRow][jColumn];
	double zb = tableSystematics[iTable+1][iRow][jColumn];
	double z = sqrt(za*za + zb*zb);

	if(iRow == 13) { // data row
	  cout << " & " << setprecision(0) << fixed << x;
	  continue;
	}

	if(x < 1.e-3) cout << " & --";
	else {
	  int prec = 0;

	  int multiplier = 1;
	  while(true) {
	    if((z == 0. && (int)(x*multiplier) > 0 && (int)(y*multiplier) > 0) ||
	       ((int)(x*multiplier) > 0 && (int)(y*multiplier) > 0 && (int)(z*multiplier) > 0)) break;
	    else {
	      multiplier *= 10;
	      prec++;
	    }
	  }

	  if(prec > 2) prec = 2;
	  
	  cout << " & " << setprecision(prec) << fixed << x << " $\\pm$ " << y << " $\\pm$ " << z;
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
