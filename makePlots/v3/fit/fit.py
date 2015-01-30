import ROOT
from ROOT import gROOT
from ROOT import gStyle

def drawPlots(data, signal, signalSF, signalName, background, backgroundSF, backgroundName, name):

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

    leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.95, '', 'brNDC')
    leg.AddEntry(data, 'Data', 'ELP')
    leg.AddEntry(signal, signalName, 'L')
    leg.AddEntry(background, backgroundName, 'L')

    can = ROOT.TCanvas()
    can.SetLogy(True)

    sumHist.Draw('hist')
    signal.Draw('hist same')
    background.Draw('hist same')
    data.Draw('e1 same')
    leg.Draw()
    can.SaveAs(name+'.pdf')

    can.SetLogy(False)
    data.Divide(sumHist)
    data.Draw('e1')
    can.SaveAs(name+'_ratio.pdf')

def makeFit(varname, varmin, varmax, signalHist, backgroundHist, dataHist, plotName):

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

    print 'fit returned value ', signalFractionVar.getVal(), ' +/- ', signalFractionVar.getError()
    return (signalFractionVar.getVal(), signalFractionVar.getError())

openfiles = {}

def get1DHist(filename, histname):
    if filename not in openfiles:
        openfiles[filename] = ROOT.TFile(filename,'READ')
    file = openfiles[filename]

    hist = ROOT.TH1D()
    hist = file.Get(histname)
    hist.SetDirectory(0)
    hist.SetFillColor(0)
    return hist

def doM3Fit(channel, controlRegion, xlo, xhi):

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    dataHist = get1DHist(input, 'm3_gg_'+channel)

    topHist = get1DHist(input, 'm3_ttJetsHadronic_'+channel)
    topHist.Add(get1DHist(input, 'm3_ttJetsFullLep_'+channel))
    topHist.Add(get1DHist(input, 'm3_ttJetsSemiLep_'+channel))
    topHist.Add(get1DHist(input, 'm3_ttA_2to5_'+channel))

    wjetsHist = get1DHist(qcdFile, 'm3_W1JetsToLNu_'+channel)
    wjetsHist.Add(get1DHist(qcdFile, 'm3_W2JetsToLNu_'+channel))
    wjetsHist.Add(get1DHist(qcdFile, 'm3_W3JetsToLNu_'+channel))
    wjetsHist.Add(get1DHist(qcdFile, 'm3_W4JetsToLNu_'+channel))

    bkgHist = get1DHist(input, 'm3_qcd_'+channel)
    bkgHist.Add(get1DHist(input, 'm3_TBar_s_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_TBar_t_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_TBar_tW_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_T_s_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_T_t_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_T_tW_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_dy1JetsToLL_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_dy2JetsToLL_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_dy3JetsToLL_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_dy4JetsToLL_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_WW_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_WZ_'+channel))
    bkgHist.Add(get1DHist(input, 'm3_ZZ_'+channel))
    bkgHist.Add(get1DHist(input, 'ttWJets'+channel))
    bkgHist.Add(get1DHist(input, 'm3_ttZJets_'+channel))
    
    dataHist.Add(bkgHist, -1.0)

    (fitFrac, fitFracErr) = makeFit('m3', xlo, xhi, topHist, wjetsHist, dataHist, 'm3_fit.png')

    dataInt = dataHist.Integral()
    topInt = topHist.Integral()
    wjetsInt = wjetsHist.Integral()
    topSF = fitFrac * dataInt / topInt
    topSFerror = fitFracErr * dataInt / topInt
    print '#'*80
    print 'Correction to the ttJets scale factor: ', topSF, ' +-', topSFerror, '(fit error only)'
    ttgammaSF = (1.0-fitFrac) * dataInt / wjetsInt
    ttgammaSFerror = fitFracErr * dataInt / wjetsInt
    print 'Correction to ttgamma scale factor: ', wjetsSF, ' +-', wjetsSFerror,'(fit error only)'
    print '#'*80

    drawPlots(dataHist, topHist, topSF, 'ttbar', wjetsHist, wjetsSF, 'wjets', 'm3')
