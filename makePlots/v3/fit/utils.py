import ROOT
from ROOT import gROOT
from ROOT import gStyle
import math
import CMS_lumi, tdrstyle

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

def integrateErrorSum(h1, h2):
    err1 = ROOT.Double(0.0)
    integral1 = h1.IntegralAndError(0, -1, err1)

    err2 = ROOT.Double(0.0)
    integral2 = h2.IntegralAndError(0, -1, err2)

    integral = integral1 + integral2
    err = ( err1**2 + err2**2 )**0.5
    
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

def drawPlots(data, signal, signalSF, signalName, background, backgroundSF, backgroundName, xlo, xhi, name, xaxisLabel, axisMin, axisMax, doLogy):

    gROOT.SetStyle('Plain')
    gStyle.SetOptStat(0000)
    gStyle.SetOptTitle(0);

    signal.Scale(signalSF)
    background.Scale(backgroundSF)

    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.5)
    data.SetFillColor(1)
    data.SetLineColor(1)

    background.SetLineColor(1)
    background.SetFillColor(8)
    
    sumHist = signal.Clone('sumHist')
    sumHist.Add(background)
    sumHist.SetLineColor(1)
    sumHist.SetFillColor(ROOT.kGray)
    sumHist.GetXaxis().SetRangeUser(xlo, xhi)

    sumHist.GetYaxis().SetRangeUser(axisMin, axisMax)

    ratio = data.Clone('ratio')
    ratio.Reset()
    for ibin in range(ratio.GetNbinsX()):
        if(sumHist.GetBinContent(ibin+1) == 0):
            continue
        ratio.SetBinContent(ibin+1, data.GetBinContent(ibin+1) / sumHist.GetBinContent(ibin+1))
        ratio.SetBinError(ibin+1, data.GetBinError(ibin+1) / sumHist.GetBinContent(ibin+1))

    ratio.GetYaxis().SetRangeUser(0.0, 1.9)
    ratio.GetYaxis().SetTitle('Data / Background')
    ratio.GetYaxis().SetLabelFont(63)
    ratio.GetYaxis().SetLabelSize(48)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleOffset(0.6)

    ratio.GetXaxis().SetRangeUser(xlo, xhi)
    ratio.GetXaxis().SetLabelFont(63)
    ratio.GetXaxis().SetLabelSize(48)
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetTitleOffset(0.6)
    ratio.GetXaxis().SetTitle(xaxisLabel)

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

    legtitle = ''
    if 'ele_bjj' in name or 'ele_jjj' in name:
        legtitle = 'ele'
    else:
        legtitle = 'muon'

    if 'jjj' in name:
        legtitle = legtitle + ' (no b-tag)'

    leg = ROOT.TLegend(0.7, 0.6, 0.85, 0.85, legtitle, 'brNDC')
    leg.AddEntry(data, 'Data', 'LP')
    leg.AddEntry(sumHist, signalName, 'F')
    leg.AddEntry(background, backgroundName, 'F')
    leg.SetFillColor(0)
    leg.SetTextSize(0.028)

    lumiHeader = ROOT.TPaveText(0.1, 0.901, 0.9, 0.94, 'NDC')
    lumiHeader.SetFillColor(0)
    lumiHeader.SetFillStyle(0)
    lumiHeader.SetLineColor(0)
    lumiHeader.SetBorderSize(0)
    lumiHeader.AddText('CMS Preliminary 2015     #sqrt{s} = 8 TeV     #intL = 19.7 fb^{-1}')

    can = ROOT.TCanvas('can', 'plot', 10, 10, 2000, 2000)
    padhi = ROOT.TPad('padhi', 'padhi', 0, 0.3, 1, 1)
    padlo = ROOT.TPad('padlo', 'padlo', 0, 0, 1, 0.3)
    padhi.SetLogy(doLogy)
    padhi.SetTickx(True)
    padhi.SetTicky(True)
    padhi.SetBottomMargin(0)
    padlo.SetTopMargin(0)
    padlo.SetBottomMargin(0.2)
    padhi.Draw()
    padlo.Draw()

    padhi.cd()
    sumHist.Draw('hist')
    background.Draw('hist same')
    data.Draw('e1 same')
    data.Draw('axis same')
    leg.Draw()
    lumiHeader.Draw()
    
    padlo.cd()
    ratio.Draw('e1')
    ratioStat.Draw('e2 same')
    ratio.Draw('e1 same')
    ratio.Draw('axis same')
    
    can.SaveAs('plots/'+name+'.pdf')
    #can.SaveAs('canvases/'+name+'.C')

def drawPlots_tdrstyle(data, signal, signalSF, signalName, background, backgroundSF, backgroundName, xlo, xhi, name, xaxisLabel, axisMin, axisMax, doLogy):

    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_8TeV = '19.7 fb^{-1}'
    CMS_lumi.extraText = 'Preliminary'
    CMS_lumi.writeExtraText = 1
    CMS_lumi.lumi_sqrtS = '8 TeV'

    iPos = 11
    H = 600
    W = 800
    # iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
    iPeriod = 2

    # references for T, B, L, R
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    
    can = ROOT.TCanvas('can', 'plot', 50, 50, W, H)

    padFraction = 0.3
    
    padhi = ROOT.TPad('padhi', 'padhi', 0, padFraction, 1, 1)
    padhi.SetLeftMargin( L/W )
    padhi.SetRightMargin( R/W )
    padhi.SetTopMargin(T/H/(1.0 - padFraction))
    padhi.SetBottomMargin(0)
    padhi.SetTickx(0)
    padhi.SetTicky(0)    
    padhi.SetLogy(doLogy)
    
    padlo = ROOT.TPad('padlo', 'padlo', 0, 0, 1, padFraction)
    padlo.SetLeftMargin( L/W )
    padlo.SetRightMargin( R/W )
    padlo.SetTopMargin(0)
    padlo.SetBottomMargin(B/H/padFraction)
    padlo.SetTickx(0)
    padlo.SetTicky(0)

    padhi.Draw()
    padlo.Draw()
        
    signal.Scale(signalSF)
    background.Scale(backgroundSF)

    background.SetFillColor(8)
    
    sumHist = signal.Clone('sumHist')
    sumHist.Add(background)
    sumHist.SetFillColor(ROOT.kGray)
    sumHist.GetXaxis().SetRangeUser(xlo, xhi)
    sumHist.GetYaxis().SetRangeUser(axisMin, axisMax)
    sumHist.GetYaxis().SetTitle('Events / 10 GeV')

    ratio = data.Clone('ratio')
    ratio.Reset()
    for ibin in range(ratio.GetNbinsX()):
        if(sumHist.GetBinContent(ibin+1) == 0):
            continue
        ratio.SetBinContent(ibin+1, data.GetBinContent(ibin+1) / sumHist.GetBinContent(ibin+1))
        ratio.SetBinError(ibin+1, data.GetBinError(ibin+1) / sumHist.GetBinContent(ibin+1))

    ratio.GetYaxis().SetRangeUser(0.45, 1.55)
    ratio.GetYaxis().SetTitle('Data / Background')
    ratio.GetYaxis().SetLabelSize(0.06)
    ratio.GetYaxis().SetTitleSize(0.06)
    ratio.GetYaxis().SetTitleOffset(0.5)
    #ratio.GetYaxis().SetNdivisions(512)

    ratio.GetXaxis().SetRangeUser(xlo, xhi)
    ratio.GetXaxis().SetTitle(xaxisLabel)
    ratio.GetXaxis().SetLabelSize(0.10)
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetTitleOffset(1.2)

    oneLine = ROOT.TLine(xlo, 1.0, xhi + 10.0, 1.0)
    oneLine.SetLineStyle(2)

    reqtitle = ''
    if 'ele_bjj' in name or 'ele_jjj' in name:
        if 'z_mass' in name:
            reqtitle = 'ee'
        else:
            reqtitle = 'e#gamma'
    else:
        if 'z_mass' in name:
            reqtitle = '#mu#mu'
        else:
            reqtitle = '#mu#gamma'

    if 'bjj' in name:
        reqtitle = reqtitle + '+bjj'
    else:
        reqtitle = reqtitle + '+jjj'

    #leg = ROOT.TLegend(0.7, 0.6, 0.85, 0.85, '', 'brNDC')
    leg = ROOT.TLegend(0.72, 0.6, 0.89, 0.85, '', 'brNDC')
    leg.SetFillColor(0)
    #leg.SetTextSize(0.028)
    leg.SetTextSize(0.028 / 0.6)
    leg.SetBorderSize(0)
    leg.AddEntry(data, 'Data', 'LP')
    leg.AddEntry(sumHist, signalName, 'F')
    leg.AddEntry(background, backgroundName, 'F')

    padhi.cd()
    sumHist.Draw('hist')
    background.Draw('hist same')
    data.Draw('e1 same')
    sumHist.Draw('axis same')
    leg.Draw()

    padlo.cd()
    ratio.Draw('e1')
    oneLine.Draw()
    ratio.Draw('axis same')

    padhi.cd()
    CMS_lumi.CMS_lumi(padhi, iPeriod, iPos, reqtitle)

    can.cd()
    
    padFix = ROOT.TPad('pad4', 'pad4', .09, .295, .118, .33)
    if not doLogy:
        padFix.Draw()

    can.SaveAs('plots/'+name+'.pdf')
    can.SaveAs('canvases/'+name+'.C')
    
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


