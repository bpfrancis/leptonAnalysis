#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

void dibosonSplit() {

  TString channel = "muon_SR1/";
  
  TFile * input = new TFile("../save_paper/limitInputs_bjj.root");

  const int nHistos = 5;
  TString names[nHistos] = {"ww", "wz", "zz", "wg", "zg"};
  TString latexNames[nHistos] = {"$\\quad WW$", "$\\quad WZ$", "$\\quad ZZ$", "$\\quad W\\gamma$", "$\\quad Z\\gamma$"};

  const int nMetBins = 5;
  double xbins_met[nMetBins+1] = {0, 20, 50, 100, 300, -1};
  double metBins[nMetBins+1];
  double realbins[nMetBins+1] = {0, 20, 50, 100, 300, 800};

  double fullWidths[nMetBins] = {20., 30., 50., 200., 500.};

  const int numSystematics = 21;
  TString systematicNames[numSystematics] = {"btagWeight", "puWeight", "topPt",
                                             "scale_tt", "scale_V", "scale_VV",
                                             "pdf_gg", "pdf_qq", "pdf_qg",
                                             "JEC", "eleSF", "muonSF", "photonSF",
                                             "ele_qcdDef", "muon_qcdDef",
                                             "userSystA_ele", "userSystA_muon",
                                             "userSystB_ele", "userSystB_muon",
                                             "userSystC_ele", "userSystC_muon"};

  double value, error;
  double sysError;

  TFile * output = new TFile("durp.root", "RECREATE");
  
  for(int i = 0; i < nHistos; i++) {

    cout << endl << endl;
    cout << names[i] << endl;
    cout << "------" << endl;

    TH1D * h = (TH1D*)input->Get(channel+names[i]);

    TH1D * h_up;
    TH1D * h_down;

    TH1D * h_out = (TH1D*)h->Clone(names[i]);
    h_out = (TH1D*)h_out->Rebin(nMetBins, "h_out_reb", realbins);
    
    for(int jBin = 0; jBin < nMetBins; jBin++) metBins[jBin] = h->GetXaxis()->FindBin(xbins_met[jBin]);
    
    for(int jBin = 0; jBin < nMetBins; jBin++) {

      value = h->IntegralAndError(metBins[jBin], metBins[jBin+1]-1, error);
      if(jBin == 0) cout << latexNames[i] << " & ";
      cout << setprecision(2) << value << " $\\pm$ " << error << " $\\pm$ ";
      //cout << value << " +/- " << error << " +/- ";

      sysError = 0.;
      
      for(int jSyst = 0; jSyst < numSystematics; jSyst++) {
	if(channel.Contains("ele") && systematicNames[jSyst].Contains("muon")) continue;
	if(channel.Contains("muon") && systematicNames[jSyst].Contains("ele")) continue;
	
	h_up = (TH1D*)input->Get(channel+names[i]+"_"+systematicNames[jSyst]+"Up");
	h_down = (TH1D*)input->Get(channel+names[i]+"_"+systematicNames[jSyst]+"Down");
	if(!h_up || !h_down) continue;

	double newErrorUp = h_up->Integral(metBins[jBin], metBins[jBin+1]-1);
	double newErrorDown = h_down->Integral(metBins[jBin], metBins[jBin+1]-1);
	double newError = fabs(newErrorUp - newErrorDown) / 2.;
	
	sysError = sqrt(sysError*sysError + newError*newError);
      }

      double oldError = h_out->GetBinError(jBin+1);
      h_out->SetBinError(jBin+1, sqrt(oldError*oldError + sysError*sysError));

      cout << setprecision(2) << sysError;
      if(jBin == nMetBins-1) cout << " \\\\" << endl;
      else cout << " & ";
      //cout << sysError << endl;

    } // for bins

    output->cd();
    h_out->Write(names[i]);

  } // for backgrounds

  output->Close();
  
  input->Close();
  
}

