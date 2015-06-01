import ROOT
from ROOT import gROOT
from ROOT import gStyle
import math

openfiles = {}

def get1DHist(filename, histname):
    if filename not in openfiles:
        openfiles[filename] = ROOT.TFile(filename,'READ')
    thisFile = openfiles[filename]

    hist = ROOT.TH1D()
    hist = thisFile.Get(histname)
    hist.SetDirectory(0)
    hist.SetFillColor(0)
    return hist

def integrateError(h, xlo, xhi):
    lowbin = h.FindBin(xlo)
    highbin = h.FindBin(xhi) - 1
    
    err = ROOT.Double(0.0)
    integral = h.IntegralAndError(lowbin, highbin, err)
    return integral, err

def ScaleWithError(h, sf, err):

    if sf == 1.0:
        return

    if sf == 0.0:
        h.Reset()
        return

    for ibin in range(h.GetNbinsX()):
        oldcontent = h.GetBinContent(ibin+1)
        olderror = h.GetBinError(ibin+1)

        if oldcontent == 0:
            continue

        newcontent = sf * oldcontent
        newerror = newcontent * math.sqrt((err*err/sf/sf) + (olderror*olderror/oldcontent/oldcontent))

        h.SetBinContent(ibin+1, newcontent)
        h.SetBinError(ibin+1, newerror)

def drawPlots(data, signal, signalSF, signalName, background, backgroundSF, backgroundName, xlo, xhi, name, xaxisLabel):

    gROOT.SetStyle('Plain')
    gStyle.SetOptStat(0000)
    gStyle.SetOptTitle(0);

    signal.Scale(signalSF)
    background.Scale(backgroundSF)

    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.5)
    data.SetFillColor(1)
    data.SetLineColor(1)

    signal.SetLineColor(2)
    signal.SetLineWidth(3)

    background.SetLineColor(3)
    background.SetLineWidth(3)
    
    sumHist = signal.Clone('sumHist')
    sumHist.Add(background)
    sumHist.SetLineColor(1)
    sumHist.SetLineWidth(3)
    sumHist.GetXaxis().SetRangeUser(xlo, xhi)

    axisMin = signal.GetBinContent(signal.GetMinimumBin()) * .5
    if data.GetBinContent(data.GetMinimumBin()) * .5 < axisMin:
        axisMin = data.GetBinContent(data.GetMinimumBin()) * .5
    if background.GetBinContent(background.GetMinimumBin()) * .5 < axisMin:
        axisMin = background.GetBinContent(background.GetMinimumBin()) * .5        

    axisMax = sumHist.GetBinContent(sumHist.GetMaximumBin()) * 5
    sumHist.GetYaxis().SetRangeUser(axisMin, axisMax)

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

    ratio.GetYaxis().SetRangeUser(0.0, 1.9)
    ratio.GetYaxis().SetTitle('Data / Background')
    ratio.GetYaxis().SetLabelFont(63)
    ratio.GetYaxis().SetLabelSize(48)
    ratio.GetYaxis().SetNdivisions(508)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleOffset(0.6)

    ratio.GetXaxis().SetRangeUser(xlo, xhi)
    ratio.GetXaxis().SetLabelFont(63)
    ratio.GetXaxis().SetLabelSize(48)
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetTitleOffset(0.6)
    ratio.GetXaxis().SetTitle(xaxisLabel)

    ratio.SetLineWidth(2)

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
    leg.AddEntry(data, 'Data', 'LP')
    leg.AddEntry(signal, signalName, 'L')
    leg.AddEntry(background, backgroundName, 'L')
    
    can = ROOT.TCanvas('can', 'plot', 10, 10, 2000, 2000)
    padhi = ROOT.TPad('padhi', 'padhi', 0, 0.3, 1, 1)
    padlo = ROOT.TPad('padlo', 'padlo', 0, 0, 1, 0.3)
    padhi.SetLogy(True)
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
    #sumErrors.Draw('e2 same')
    data.Draw('e1 same')
    leg.Draw()
    
    padlo.cd()
    ratio.Draw('e1')
    ratioStat.Draw('e2 same')
    ratio.Draw('e1 same')
    ratio.Draw('axis same')
    
    can.SaveAs('plots/'+name+'.pdf')
    #can.SaveAs('canvases/'+name+'.C')

def makeFit(varname, varmin, varmax, signalHist, backgroundHist, dataHist):

    # RooFit variables
    var = ROOT.RooRealVar(varname, varname, varmin, varmax)
    argList = ROOT.RooArgList()
    argList.add(var)
    argSet = ROOT.RooArgSet()
    argSet.add(var)

    # create PDFs
    signalDataHist = ROOT.RooDataHist('signalDataHist', 'signal RooDataHist', argList, signalHist)
    signalPdf = ROOT.RooHistPdf('signalPdf', varname+' of signal', argSet, signalDataHist)

    backgroundDataHist = ROOT.RooDataHist('backgroundDataHist', 'background RooDataHist', argList, backgroundHist)
    backgroundPdf = ROOT.RooHistPdf('backgroundPdf', varname+' of background', argSet, backgroundDataHist)

    # data
    dataDataHist = ROOT.RooDataHist('data '+varname, varname+' in Data', argList, dataHist)

    # signal fraction parameter
    signalFractionVar = ROOT.RooRealVar('signal fraction', 'signal fraction', 0.5, 0.0, 1.0)
    sumPdf = ROOT.RooAddPdf('totalPdf', 'signal and background', signalPdf, backgroundPdf, signalFractionVar)

    # fit
    sumPdf.fitTo(dataDataHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1))
    #sumPdf.fitTo(dataDataHist, ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.PrintLevel(-1))

    print 'fit returned value ', signalFractionVar.getVal(), ' +/- ', signalFractionVar.getError()
    return (signalFractionVar.getVal(), signalFractionVar.getError())


