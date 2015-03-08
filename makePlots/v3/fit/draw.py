import ROOT
from ROOT import gROOT
from ROOT import gStyle

def drawPlots(data, signal, signalSF, signalName, background, backgroundSF, backgroundName, xlo, xhi, name):

    signal.Scale(signalSF)
    background.Scale(backgroundSF)

    signal.SetLineColor(2)
    signal.SetLineWidth(3)

    background.SetLineColor(3)
    background.SetLineWidth(3)
    
    sumHist = signal.Clone('sumHist')
    sumHist.Add(background)
    sumHist.SetLineColor(1)
    sumHist.SetLineWidth(3)
    sumHist.GetXaxis().SetRangeUser(xlo, xhi)

    sumErrors = sumHist.Clone('sumErrors')
    sumErrors.SetFillColor(ROOT.kOrange+10)
    sumErrors.SetFillStyle(3154)
    sumErrors.SetMarkerSize(0)

    ratio = data.Clone('ratio')
    ratio.Reset()
    ratio.SetTitle('Data / Background')
    for ibin in range(ratio.GetNbinsX()):
        if(sumHist.GetBinContent(ibin+1) == 0):
            continue
        ratio.SetBinContent(ibin+1, data.GetBinContent(ibin+1) / sumHist.GetBinContent(ibin+1))
        ratio.SetBinError(ibin+1, data.GetBinError(ibin+1) / sumHist.GetBinContent(ibin+1))

    ratio.GetXaxis().SetRangeUser(xlo, xhi)
    ratio.GetYaxis().SetRangeUser(0.0, 2.0)

    ratioStat = sumHist.Clone('ratioStat')
    for ibin in range(ratioStat.GetNbinsX()):
        ratioStat.SetBinContent(ibin+1, 1.0)
        if(sumHist.GetBinContent(ibin+1) == 0.0):
            ratioStat.SetBinError(ibin+1, 0.0)
        else:
            ratioStat.SetBinError(ibin+1, sumHist.GetBinError(ibin+1) / sumHist.GetBinContent(ibin+1))


    ratioStat.SetFillStyle(1001)
    ratioStat.SetFillColor(ROOT.kGray+1)
    ratioStat.SetLineColor(ROOT.kGray+1)
    ratioStat.SetMarkerColor(ROOT.kGray+1)
            
    leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.95, '', 'brNDC')
    leg.AddEntry(data, 'Data', 'ELP')
    leg.AddEntry(signal, signalName, 'L')
    leg.AddEntry(background, backgroundName, 'L')
    
    can = ROOT.TCanvas()
    padhi = ROOT.TPad('padhi', 'padhi', 0, 0.3, 1, 1)
    padlo = ROOT.TPad('padlo', 'padlo', 0, 0, 1, 0.3)
    padhi.SetLogy(False)
    padhi.SetTickx(True)
    padhi.SetTicky(True)
    padhi.SetBottomMargin(0)
    padlo.SetTopMargin(0)
    padlo.SetBottomMargin(0.2)
    padhi.Draw()
    padlo.Draw()

    padhi.cd()
    sumHist.Draw('hist')
    signal.Draw('hist same')
    background.Draw('hist same')
    sumErrors.Draw('e2 same')
    data.Draw('e1 same')
    leg.Draw()
    
    padlo.cd()
    ratio.Draw('e1')
    ratioStat.Draw('e2 same')
    ratio.Draw('e1 same')
    ratio.Draw('axis same')
    
    can.SaveAs(name+'.png')
    #can.SaveAs(name+'.C')
    
def drawPlots_threeParameters(a, aSF, b, bSF, c, cSF, data, name):

    a.Scale(aSF)
    a.SetLineColor(2)
    a.SetLineWidth(3)

    b.Scale(bSF)
    b.SetLineColor(3)
    b.SetLineWidth(3)

    c.Scale(cSF)
    c.SetLineColor(4)
    c.SetLineWidth(3)

    sumHist = a.Clone('sumHist')
    sumHist.Add(b)
    sumHist.Add(c)
    sumHist.SetLineColor(1)
    sumHist.SetLineWidth(3)

    leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.95, '', 'brNDC')
    leg.AddEntry(data, 'Data', 'ELP')
    leg.AddEntry(a, a.GetName(), 'L')
    leg.AddEntry(b, b.GetName(), 'L')
    leg.AddEntry(c, c.GetName(), 'L')

    can = ROOT.TCanvas()
    can.SetLogy(False)

    sumHist.Draw('hist')
    a.Draw('hist same')
    b.Draw('hist same')
    c.Draw('hist same')
    data.Draw('e1 same')
    leg.Draw()
    can.SaveAs(name+'.png')

def drawShapes(a, b, c, data, name):

    a.Scale(1. / a.Integral())
    a.SetLineColor(2)
    a.SetLineWidth(3)

    b.Scale(1. / b.Integral())
    b.SetLineColor(3)
    b.SetLineWidth(3)

    c.Scale(1. / c.Integral())
    c.SetLineColor(4)
    c.SetLineWidth(3)

    data.Scale(1. / data.Integral())

    leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.95, '', 'brNDC')
    leg.AddEntry(data, 'Data', 'ELP')
    leg.AddEntry(a, a.GetName(), 'L')
    leg.AddEntry(b, b.GetName(), 'L')
    leg.AddEntry(c, c.GetName(), 'L')

    can = ROOT.TCanvas()
    can.SetLogy(False)

    a.Draw('hist')
    b.Draw('hist same')
    c.Draw('hist same')
    data.Draw('e1 same')
    leg.Draw()
    can.SaveAs(name+'.png')
