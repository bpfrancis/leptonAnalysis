void fillPotHoles(TH2D *h) {

  // fill cells which have empty value

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  double epsilon = 1e-10;

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {
      double val = h->GetBinContent(ix,iy);
      if(val != val) h->SetBinContent(ix,iy,0); // checking for NAN
    }
  }

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {

      if(h->GetXaxis()->GetBinCenter(ix) < h->GetYaxis()->GetBinCenter(iy)) continue;

      double val = h->GetBinContent(ix,iy);
      if(fabs(val) > epsilon) continue;
      int ncnt = 0;
      double sum = 0;
      double sumErr = 0;
      double up    = h->GetBinContent(ix,iy+1);
      if(fabs(up) > epsilon && iy < nbinsY){
	sum += up;
	sumErr += h->GetBinError(ix,iy+1)*h->GetBinError(ix,iy+1);
	ncnt++;
      }
      double down  = h->GetBinContent(ix,iy-1);
      if(fabs(down) > epsilon && iy > 1){
	sum += down;
	sumErr += h->GetBinError(ix,iy-1)*h->GetBinError(ix,iy-1);
	ncnt++;
      }
      double left  = h->GetBinContent(ix-1,iy);
      if(fabs(left) > epsilon && ix > 1){
	sum += left;
	sumErr += h->GetBinError(ix-1,iy)*h->GetBinError(ix-1,iy);
	ncnt++;
      }
      double right = h->GetBinContent(ix+1,iy);
      if(fabs(right) > epsilon && ix < nbinsX){
	sum += right;
	sumErr += h->GetBinError(ix+1,iy)*h->GetBinError(ix+1,iy);
	ncnt++;
      }
      if(ncnt > 0) {
	h->SetBinContent(ix,iy,sum/ncnt);
	h->SetBinError(ix,iy,sqrt(sumErr)/ncnt);
      }
    } // for iy
  } // for ix

}

void go(TString channel, TString title) {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  Double_t mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 
		      385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 
		      1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  Double_t mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 
			275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 
			725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};
  
  Double_t xbins[31];
  xbins[0] = 0;
  xbins[1] = 55;
  for(int i = 1; i < 29; i++) xbins[i+1] = (mst[i] + mst[i-1])/2.;
  xbins[30] = 6510;

  Double_t ybins[33];
  ybins[0] = 0;
  ybins[1] = 12.5;
  for(int i = 1; i < 31; i++) ybins[i+1] = (mBino[i] + mBino[i-1])/2.;
  ybins[32] = 2175;

  char code[100];
  int index1, index2;

  TH2D * h_acc = new TH2D("acc_"+channel, "acc_"+channel, 30, xbins, 32, ybins);

  TFile * input = new TFile("../limitInputs_bjj.root");

  TFile * f_xsec = new TFile("/eos/uscms/store/user/bfrancis/data/stop-bino_xsecs.root", "READ");
  TH2D * h_xsec = (TH2D*)f_xsec->Get("real_xsec");

  TH1D * h;

  for(int i = 0; i < 899; i++) {

    index1 = mst[int(i)/31];
    index2 = mBino[int(i)%31];

    sprintf(code, "_mst_%d_m1_%d", index1, index2);
    TString code_t = code;

    h = (TH1D*)input->Get(channel+"/signal"+code_t);
    if(!h) continue;

    double xsec = h_xsec->GetBinContent(h_xsec->FindBin(index1, index2)) * 19712.;
    double branchingRatio = 0.438 / 3.;

    double value = h->Integral() / xsec / branchingRatio;
    
    h_acc->SetBinContent(h_acc->FindBin(index1, index2), value);

  }

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);
  fillPotHoles(h_acc);

  h_acc->GetXaxis()->SetTitle("Stop mass (GeV/c^{2})");
  h_acc->GetXaxis()->SetRangeUser(0, 1600);
  h_acc->GetXaxis()->SetLabelSize(0.025);

  h_acc->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_acc->GetYaxis()->SetTitleOffset(1.2);
  h_acc->GetYaxis()->SetLabelSize(0.025);
  h_acc->GetYaxis()->SetRangeUser(0, 1600);

  h_acc->GetZaxis()->SetLabelSize(0.02);
  if(channel.Contains("_Any")) h_acc->GetZaxis()->SetRangeUser(0, 0.52);
  else if(channel.Contains("_SR1")) h_acc->GetZaxis()->SetRangeUser(0, 0.28);
  else if(channel.Contains("_SR2")) h_acc->GetZaxis()->SetRangeUser(0, 0.19);

  h_acc->Draw("colz");

  TPaveText * labelText = new TPaveText(0.2, 0.6, 0.5, 0.75, "NDC");
  labelText->SetFillColor(0);
  labelText->SetFillStyle(0);
  labelText->SetLineColor(0);
  labelText->SetBorderSize(0);
  labelText->AddText("A #times #scale[1.5]{#varepsilon} / BR(t#bar{t} #rightarrow l#nubbjj)");
  labelText->Draw("same");

  TPaveText * reqText = new TPaveText(0.2, 0.55, 0.5, 0.6, "NDC");
  reqText->SetFillColor(0);
  reqText->SetFillStyle(0);
  reqText->SetLineColor(0);
  reqText->SetBorderSize(0);
  reqText->AddText(title);
  reqText->Draw("same");

  TLine * nlspLine = new TLine(0, 0, 1600, 1600);
  nlspLine->SetLineStyle(2);
  nlspLine->SetLineWidth(2);

  TLine * virtualLine = new TLine(172.5, 0, 1600, 1600 - 172.5);
  virtualLine->SetLineStyle(2);
  virtualLine->SetLineWidth(2);

  TLatex * nlspComment = new TLatex(1100 + 172.5, 1125 + 172.5, "mStop < mBino");
  nlspComment->SetTextAngle(45);
  nlspComment->SetTextSize(0.02);
  
  TLatex * virtualComment = new TLatex(172.5 + 1150, 1100, "mStop - mBino < m_{t}");
  virtualComment->SetTextAngle(45);
  virtualComment->SetTextSize(0.02);

  virtualLine->Draw("same");
  nlspLine->Draw("same");
  virtualComment->Draw("same");
  nlspComment->Draw("same");

  can->SaveAs("acceptance_"+channel+".pdf");

  delete can;
  input->Close();

}

void plotAcceptance() {

  go("ele_Any", "ele Pre-selection");
  go("ele_SR1", "ele SR1");
  go("ele_SR2", "ele SR2");

  go("muon_Any", "muon Pre-selection");
  go("muon_SR1", "muon SR1");
  go("muon_SR2", "muon SR2");

}
