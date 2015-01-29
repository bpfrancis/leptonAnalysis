#include "TH1.h"
#include "TTree.h"
#include "TString.h"

#include <vector>
#include <stdarg.h>

using namespace std;

bool LargerHistogram(const TH1D* h1, const TH1D* h2) {
  return (h1->GetMaximum() > h2->GetMaximum());
}

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

void fillPotHoles(TH2D *h) {

  // fill cells which have empty value

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  double epsilon = 1e-10;

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {
      double val = h->GetBinContent(ix,iy);
      if(isnan(val)) h->SetBinContent(ix,iy,0); // checking for NAN
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

void fixBadCells(TH2D* h) {

  // fix bad cells which have wrong sign compared to 4 surrounding cells.
  // then assign average value from 4 cells

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  double epsilon = 0;

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {
      double val = h->GetBinContent(ix,iy);
      if(val < epsilon) {
        int ncnt = 0;
        double up    = h->GetBinContent(ix,iy+1);
        if(up > epsilon) ncnt++;
        double down  = h->GetBinContent(ix,iy-1);
        if(down > epsilon) ncnt++;
        double left  = h->GetBinContent(ix-1,iy);
        if(left > epsilon) ncnt++;
        double right = h->GetBinContent(ix+1,iy);
        if(right > epsilon) ncnt++;
        if(ncnt == 4){
          val = (up+down+left+right)/ncnt;
          h->SetBinContent(ix,iy,val);
          up    = h->GetBinError(ix,iy+1);
          down  = h->GetBinError(ix,iy-1);
          left  = h->GetBinError(ix-1,iy);
          right = h->GetBinError(ix+1,iy);
          val = sqrt(up*up + down*down + left*left + right*right)/ncnt;
          h->SetBinError(ix,iy,val);
        }
      }
      else {
        int ncnt = 0;
        double up    = h->GetBinContent(ix,iy+1);
        if(up < epsilon) ncnt++;
        double down  = h->GetBinContent(ix,iy-1);
        if(down < epsilon) ncnt++;
        double left  = h->GetBinContent(ix-1,iy);
        if(left < epsilon) ncnt++;
        double right = h->GetBinContent(ix+1,iy);
        if(right < epsilon) ncnt++;
        if(ncnt == 4){
          val = (up+down+left+right)/ncnt;
          h->SetBinContent(ix,iy,val);
          up    = h->GetBinError(ix,iy+1);
          down  = h->GetBinError(ix,iy-1);
          left  = h->GetBinError(ix-1,iy);
          right = h->GetBinError(ix+1,iy);
          val = sqrt(up*up + down*down + left*left + right*right)/ncnt;
          h->SetBinError(ix,iy,val);
        }
      }
    } // for iy
  } // for ix
}
