#include <vector>

using namespace std;

void plot(vector<TString> names, TString channel, TString variable, TString photon, bool logy) {

  if(names.size() == 0) return;

  vector<TFile*> files_0g, files_1g, files_2g;
  
  TH1D * h_0g;
  TH1D * h_1g;
  TH1D * h_2g;

  for(unsigned int i = 0; i < names.size(); i++) {
    TFile * f = new TFile("0g/"+names[i]+"_"+channel+"_"+photon+".root", "READ");
    if(i == 0) {
      TH1D * h = (TH1D*)f->Get(variable);
      h_0g = (TH1D*)h->Clone();
    }
    else h_0g->Add((TH1D*)f->Get(variable));

    f = new TFile("1g/"+names[i]+"_"+channel+"_"+photon+".root", "READ");
    if(i == 0) {
      TH1D * h = (TH1D*)f->Get(variable);
      h_1g = (TH1D*)h->Clone();
    }
    else h_1g->Add((TH1D*)f->Get(variable));
		   
    f = new TFile("2g/"+names[i]+"_"+channel+"_"+photon+".root", "READ");
    if(i == 0) {
      TH1D * h = (TH1D*)f->Get(variable);
      h_2g = (TH1D*)h->Clone();
    }
    else h_2g->Add((TH1D*)f->Get(variable));

  }

  h_0g->Scale(1./h_0g->Integral());
  h_1g->Scale(1./h_1g->Integral());
  h_2g->Scale(1./h_2g->Integral());

  h_1g->SetLineColor(kBlue);
  h_2g->SetLineColor(kRed);

  TCanvas * can = new TCanvas("canvas_"+variable, "Plot", 10, 10, 800, 800);
  can->SetLogy(logy);

  h_0g->Draw();
  h_1g->Draw("same");
  h_2g->Draw("same");

  can->SaveAs("durp.png");

}


void compareShapes() {

  TString channel = "ele_bjj";
  TString variable = "pfMET";
  TString photon = "gamma";
  bool logy = true;

  vector<TString> names;
  names.push_back("ttJetsSemiLep");
  names.push_back("ttJetsFullLep");
  names.push_back("ttJetsHadronic");

  plot(names, channel, variable, photon, logy);
}
  
