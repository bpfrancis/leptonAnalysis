void overlay_limits() {

  gROOT->ForceStyle();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  //For the temperature plots
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleOffset(1.4, "xz");
  gStyle->SetTitleOffset(1.9, "y");

  gStyle->SetNdivisions(505);
  gStyle->SetTitleFont(43, "xyz");
  gStyle->SetTitleSize(32, "xyz");
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetLabelSize(0.04, "xyz");

  TCanvas * can = new TCanvas("can", "can", 10, 10, 900, 800);

  int xMin = 222.5;
  int xMax = 960;
  int yMin = 137.5;
  int yMax = 775;

  TString xLabel = "m_{Stop} [GeV]";
  TString yLabel = "m_{Bino} [GeV]";

  TString scan = "stop-bino";

  const int nX = 16;
  const int nY = 16;

  Double_t mst[nX] = {235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910};
  Double_t mBino[nY] = {150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725};

  Double_t xBins[nX+1];
  xBins[0] = 222.5;
  for(int i = 1; i < nX; i++) xBins[i] = (mst[i] + mst[i-1])/2.;
  xBins[nX] = 960;

  Double_t yBins[nY+1];
  yBins[0] = 137.5;
  for(int i = 1; i < nY; i++) yBins[i] = (mBino[i] + mBino[i-1])/2.;
  yBins[nY] = 775;

  double leg_xmin = 0.58;
  double leg_xmax = 0.9;
  double leg_ymin = 0.64;
  double leg_ymax = 0.87;

  leg_xmin -= 0.35;
  leg_xmax -= 0.35;

  TFile * file1 = new TFile("limits_sr1.root");
  TFile * file2 = new TFile("limits_sr2.root");
  TFile * fileBoth = new TFile("limits_both.root");

  TH2D * h_back = new TH2D("h_back_"+scan,";"+xLabel+";"+yLabel,100,xMin,xMax,100,yMin,yMax);

  TGraph * sr1 = (TGraph*)file1->Get("contour_obs");
  TGraph * sr2 = (TGraph*)file2->Get("contour_obs");
  TGraph * both = (TGraph*)fileBoth->Get("contour_obs");

  TGraph * exp_band_sr1 = (TGraph*)file1->Get("contour_exp_1s_band");
  TGraph * exp_band_sr2 = (TGraph*)file2->Get("contour_exp_1s_band");
  TGraph * exp_band_both = (TGraph*)fileBoth->Get("contour_exp_1s_band");

  sr1->SetLineColor(kBlue);
  sr1->SetLineWidth(3);

  sr2->SetLineColor(kRed);
  sr2->SetLineWidth(3);

  both->SetLineColor(kBlack);
  both->SetLineWidth(3);

  exp_band_sr1->SetFillColor(kBlue-9);

  exp_band_sr2->SetFillStyle(3154);
  exp_band_sr2->SetFillColor(kRed-3);

  exp_band_both->SetFillColor(kOrange-3);

  TGraph * upperDiagonalRegion = new TGraph(4);
  upperDiagonalRegion->SetPoint(0, yMin + 172.5, yMin);
  upperDiagonalRegion->SetPoint(1, xMin, yMin);
  upperDiagonalRegion->SetPoint(2, xMin, yBins[nY]);
  upperDiagonalRegion->SetPoint(3, yBins[nX] + 172.5, yBins[nY]);
  upperDiagonalRegion->SetFillColor(16);

  TLine * nlspLine = new TLine(222.5, 222.5, 775, 775);
  nlspLine->SetLineStyle(2);
  nlspLine->SetLineWidth(2);

  TLine * virtualLine = new TLine(310, 137.5, 947.5, 775);
  virtualLine->SetLineStyle(2);
  virtualLine->SetLineWidth(2);

  TLatex * nlspComment = new TLatex(250, 275, "mStop < mBino");
  nlspComment->SetTextAngle(45);
  nlspComment->SetTextSize(0.02);
  
  TLatex * virtualComment = new TLatex(300, 150, "mStop - mBino < m_{t}");
  virtualComment->SetTextAngle(45);
  virtualComment->SetTextSize(0.02);

  TLegend* leg2 = new TLegend(leg_xmin, leg_ymin, leg_xmax - 0.05, leg_ymax, NULL, "brNDC");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.03);

  leg2->AddEntry(sr1, "Obs SR1", "L");
  leg2->AddEntry(exp_band_sr1, "Exp SR1 #pm1#sigma", "F");
  leg2->AddEntry(sr2, "Obs SR2", "L");
  leg2->AddEntry(exp_band_sr2, "Exp SR2 #pm1#sigma", "F");
  leg2->AddEntry(both, "Obs Both", "L");
  leg2->AddEntry(exp_band_both, "Exp Both #pm1#sigma", "F");
		 
  h_back->Draw("axis");
  exp_band_both->Draw("same F");
  exp_band_sr1->Draw("same F");
  exp_band_sr2->Draw("same F");
  sr1->Draw("same L");
  sr2->Draw("same L");
  both->Draw("same L");

  upperDiagonalRegion->Draw("same f");
  virtualLine->Draw("same");
  nlspLine->Draw("same");
  virtualComment->Draw("same");
  nlspComment->Draw("same");
  leg2->Draw("same");

  can->SaveAs("overlay.pdf");

}
