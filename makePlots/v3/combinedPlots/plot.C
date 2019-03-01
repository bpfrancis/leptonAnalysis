#include <vector>

using namespace std;

void plot() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  const int numChannels = 4;
  TString channels[numChannels] = {"SR1", "SR2", "CR1", "CR2"};
  
  const int numBackgrounds = 8;
  TString backgroundNames[numBackgrounds] = {"qcd", "ttjets", "wjets", "zjets", "singleTop", "vv", "ttV", "ttgamma"};
  TString legendNames[numBackgrounds] = {"QCD", "t#bar{t} + Jets", "W + Jets", "Z/#gamma* + Jets", "Single top", "VV, V#gamma", "t#bar{t} + V", "t#bar{t} + #gamma"};
  int backgroundColors[numBackgrounds] = {kSpring-6, kGray, kOrange-3, kYellow, kRed, kViolet-2, kCyan, 8};
  
  TFile * input = new TFile("combined.root", "READ");

  TCanvas * pasCan = new TCanvas("pasCan", "Plot", 50, 50, 800, 600);
  pasPadhi = new TPad("pasPadhi", "pasPadhi", 0, 0.3, 1, 1);
  pasPadlo = new TPad("pasPadlo", "pasPadlo", 0, 0, 1, 0.3);
  pasPadhi->SetLogy(true);
  pasPadhi->SetTickx(true);
  pasPadhi->SetTicky(true);
  
  pasPadhi->SetLeftMargin(0.12);
  pasPadhi->SetRightMargin(0.04);
  pasPadhi->SetTopMargin(0.08/(1.0 - 0.3));
  pasPadhi->SetBottomMargin(0);
  
  pasPadlo->SetLeftMargin(0.12);
  pasPadlo->SetRightMargin(0.04);
  pasPadlo->SetTopMargin(0);
  pasPadlo->SetBottomMargin(0.12/0.3);
  
  pasPadhi->Draw();
  pasPadlo->Draw();

  TLegend * pasLegDrawSignal;
  pasLegDrawSignal = new TLegend(0.55, 0.53, 0.9, 0.83, NULL, "brNDC");
  pasLegDrawSignal->SetNColumns(2);
  pasLegDrawSignal->SetFillColor(0);
  pasLegDrawSignal->SetBorderSize(0);
  pasLegDrawSignal->SetTextSize(0.028 / 0.7);
  
  for(int channel = 0; channel < numChannels; channel++) {
  
    TH1D * data = (TH1D*)input->Get(channels[channel] + "/data");
    TH1D * error_syst = (TH1D*)input->Get(channels[channel] + "/error_syst");

    vector<TH1D*> backgrounds;

    bool drawnFirst = false;

    pasPadhi->cd();

    // durp
    pasLegDrawSignal->Clear();
    pasLegDrawSignal->AddEntry(data, "Data", "LP");
    pasLegDrawSignal->AddEntry((TObject*)0, "", "");
    pasLegDrawSignal->AddEntry(error_syst, "Stat. #oplus Syst. Errors", "F");
    pasLegDrawSignal->AddEntry((TObject*)0, "", "");
    //pasLegDrawSignal->AddEntry(siga, "GGM (460_175)", "L");
    //pasLegDrawSignal->AddEntry(sigb, "GGM (560_325)", "L");
    
    for(int iBkg = 0; iBkg < numBackgrounds; iBkg++) {
      if(channels[channel] != "CR1" && iBkg == 0) continue;

      backgrounds.push_back((TH1D*)input->Get(channels[channel] + "/" + backgroundNames[iBkg]));
      backgrounds.back()->SetFillColor(backgroundColors[iBkg]);
      backgrounds.back()->SetMarkerSize(0);
      backgrounds.back()->SetLineColor(1);

      pasLegDrawSignal->AddEntry(backgrounds.back(), legendNames[iBkg], "F");
      
      if(!drawnFirst) {
	backgrounds.back()->SetTitle("pfMet_t01");
	backgrounds.back()->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");
	backgrounds.back()->GetYaxis()->SetTitle("Number of Events / GeV");
	if(channels[channel] == "CR1") backgrounds.back()->GetYaxis()->SetRangeUser(2.e-4, 1.e2);
	if(channels[channel] == "CR2") backgrounds.back()->GetYaxis()->SetRangeUser(3.e-4, 40.);
	if(channels[channel] == "SR1") backgrounds.back()->GetYaxis()->SetRangeUser(7.e-4, 1.e2);
	if(channels[channel] == "SR2") backgrounds.back()->GetYaxis()->SetRangeUser(3.e-4, 2.0);
	
	backgrounds.back()->Draw("hist");
	drawnFirst = true;
      }
      else backgrounds.back()->Draw("hist same");
    }

    if(channels[channel] != "CR1") pasLegDrawSignal->AddEntry((TObject*)0, "", "");

    data->SetLineColor(kBlack);
  
    error_syst->SetFillColor(kOrange+10);
    error_syst->SetFillStyle(3154);
    error_syst->SetMarkerSize(0);
    
    error_syst->Draw("same e2");
    data->Draw("same e1");
    backgrounds.front()->Draw("same axis");
    // draw siga/sigb
    pasLegDrawSignal->Draw("same");

    TMathText pasLumiLatex;
    pasLumiLatex.SetNDC();
    pasLumiLatex.SetTextAngle(0);
    pasLumiLatex.SetTextColor(kBlack);
    pasLumiLatex.SetTextFont(42);
    pasLumiLatex.SetTextAlign(31);
    pasLumiLatex.SetTextSize(0.048 / 0.7);
    pasLumiLatex.DrawMathText(0.96, 0.9, "19.7 fb^{-1} (8 TeV) #ell#gamma+bjj");

    TLatex pasCMSLatex;
    pasCMSLatex.SetNDC();
    pasCMSLatex.SetTextAngle(0);
    pasCMSLatex.SetTextColor(kBlack);
    pasCMSLatex.SetTextFont(61);
    pasCMSLatex.SetTextAlign(11);
    pasCMSLatex.SetTextSize(0.06 / 0.7);
    pasCMSLatex.DrawLatex(0.12, 0.9, "CMS");

    TLatex pasPrelimLatex;
    pasPrelimLatex.SetNDC();
    pasPrelimLatex.SetTextAngle(0);
    pasPrelimLatex.SetTextColor(kBlack);
    pasPrelimLatex.SetTextFont(52);
    pasPrelimLatex.SetTextAlign(11);
    pasPrelimLatex.SetTextSize(0.0456 / 0.7);
    pasPrelimLatex.DrawLatex(0.2178, 0.9, "Preliminary");

    pasPadlo->cd();

    TH1D * ratio = (TH1D*)input->Get(channels[channel] + "/ratio");
    TH1D * ratio_stat = (TH1D*)input->Get(channels[channel] + "/ratio_stat");
    TH1D * ratio_sys = (TH1D*)input->Get(channels[channel] + "/ratio_sys");

    ratio->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetTitle("Data / Background");
    ratio->GetYaxis()->SetLabelSize(0.06);
    ratio->GetYaxis()->SetTitleSize(0.06);
    ratio->GetYaxis()->SetTitleOffset(0.5);
    ratio->GetYaxis()->SetNdivisions(508);
    ratio->SetLineColor(kBlack);
    
    ratio_stat->SetFillStyle(1001);
    ratio_stat->SetFillColor(kGray+1);
    ratio_stat->SetLineColor(kGray+1);
    ratio_stat->SetMarkerColor(kGray+1);

    ratio_sys->SetFillStyle(1001);
    ratio_sys->SetFillColor(kGray);
    ratio_sys->SetLineColor(kGray);
    ratio_sys->SetMarkerColor(kGray);
    
    ratio->Draw("e1");
    ratio_sys->Draw("e2 same");
    //ratio_stat->Draw("e2 same");
    ratio->Draw("e1 same");
    ratio->Draw("axis same");
    // not in pas style?
    //ratioLeg->Draw("same");

    TLine * oneLine = new TLine(0, 1, 800, 1);
    oneLine->SetLineStyle(2);
    oneLine->Draw();
    
    pasCan->SaveAs(channels[channel] + ".pdf");

    backgrounds.clear();


  } // for channels

  input->Close();
  
}
      

  
