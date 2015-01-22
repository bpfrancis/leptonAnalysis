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

def doFit(variable, channel, xlo, xhi):

    dataHist = get1DHist('data_'+channel+'_fake.root', variable)

    ttJetsHist = get1DHist('ttJetsSemiLep_'+channel+'_fake.root', variable)
    ttJetsHist.Add(get1DHist('ttJetsFullLep_'+channel+'_fake.root', variable))
    ttJetsHist.Add(get1DHist('ttJetsHadronic_'+channel+'_fake.root', variable))

    ttgammaHist = get1DHist('ttA_2to5_'+channel+'_fake.root', variable)

    (fitFrac, fitFracErr) = makeFit(variable, xlo, xhi, ttJetsHist, ttgammaHist, dataHist, variable+'_fit.png')

    dataInt = dataHist.Integral()
    topInt = ttJetsHist.Integral()
    ttgammaInt = ttgammaHist.Integral()
    ttJetsSF = fitFrac * dataInt / topInt
    ttJetsSFerror = fitFracErr * dataInt / topInt
    print '#'*80
    print 'Correction to the ttJets scale factor: ', ttJetsSF, ' +-', ttJetsSFerror, '(fit error only)'
    ttgammaSF = (1.0-fitFrac) * dataInt / ttgammaInt
    ttgammaSFerror = fitFracErr * dataInt / ttgammaInt
    print 'Correction to ttgamma scale factor: ', ttgammaSF, ' +-', ttgammaSFerror,'(fit error only)'
    print '#'*80

    drawPlots(dataHist, ttJetsHist, ttJetsSF, 'ttJets', ttgammaHist, ttgammaSF, 'ttgamma', variable)
