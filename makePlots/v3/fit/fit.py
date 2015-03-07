from draw import *
import math

def ScaleWithError(h, sf, err):

    if sf == 1.0:
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

###########################################################
def doQCDFit(channel, controlRegion, systematic, xlo, xhi):

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    dataHist = get1DHist(input, 'pfMET_gg_'+channel)
    qcdHist = get1DHist(input, 'pfMET_qcd_'+channel)

    MCHist = get1DHist(input, 'pfMET_ttJetsHadronic_'+channel+systematic)
    MCHist.Add(get1DHist(input, 'pfMET_ttJetsSemiLep_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_ttJetsFullLep_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_W1JetsToLNu_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_W2JetsToLNu_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_W3JetsToLNu_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_W4JetsToLNu_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_TBar_s_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_TBar_t_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_TBar_tW_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_T_s_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_T_t_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_T_tW_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_dy1JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_dy2JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_dy3JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_dy4JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_ttA_2to5_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_WW_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_WZ_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_ZZ_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_TTWJets_'+channel+systematic))
    MCHist.Add(get1DHist(input, 'pfMET_TTZJets_'+channel+systematic))

    (qcdFrac, qcdFracErr) = makeFit('pfMET', xlo, xhi, qcdHist, MCHist, dataHist)

    #lowbin = dataHist.FindBin(xlo)
    #highbin = dataHist.FindBin(xhi) - 1

    lowbin = 1
    highbin = dataHist.GetNbinsX()+1

    dataInt = dataHist.Integral(lowbin, highbin)
    qcdInt = qcdHist.Integral(lowbin, highbin)
    mcInt = MCHist.Integral(lowbin, highbin)

    QCDSF = qcdFrac*dataInt/qcdInt
    QCDSFerror = qcdFracErr*dataInt/qcdInt
    MCSF = (1-qcdFrac)*dataInt/mcInt
    MCSFerror = qcdFracErr*dataInt/mcInt

    drawPlots(dataHist, qcdHist, QCDSF, 'QCD', MCHist, MCSF, 'MC', xlo, xhi, 'pfMET_'+channel+systematic)

    return (QCDSF, QCDSFerror, MCSF, MCSFerror)

def doM3Fit(channel, controlRegion, systematic, output, xlo, xhi):

    (QCDSF, QCDSFerror, MCSF, MCSFerror) = doQCDFit(channel, controlRegion, systematic, 0.0, 300.0)

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    dataHist = get1DHist(input, 'm3_gg_'+channel)

    topHist = get1DHist(input, 'm3_ttJetsHadronic_'+channel+systematic)
    topHist.Add(get1DHist(input, 'm3_ttJetsFullLep_'+channel+systematic))
    topHist.Add(get1DHist(input, 'm3_ttJetsSemiLep_'+channel+systematic))
    ScaleWithError(topHist, MCSF, MCSFerror)

    wjetsHist = get1DHist(input, 'm3_W1JetsToLNu_'+channel+systematic)
    wjetsHist.Add(get1DHist(input, 'm3_W2JetsToLNu_'+channel+systematic))
    wjetsHist.Add(get1DHist(input, 'm3_W3JetsToLNu_'+channel+systematic))
    wjetsHist.Add(get1DHist(input, 'm3_W4JetsToLNu_'+channel+systematic))
    ScaleWithError(wjetsHist, MCSF, MCSFerror)

    qcdHist = get1DHist(input, 'm3_qcd_'+channel)
    ScaleWithError(qcdHist, QCDSF, QCDSFerror)

    bkgHist = get1DHist(input, 'm3_TBar_s_'+channel+systematic)
    bkgHist.Add(get1DHist(input, 'm3_TBar_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_TBar_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_T_s_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_T_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_T_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_dy1JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_dy2JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_dy3JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_dy4JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_WW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_WZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_ZZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_TTWJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_TTZJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_ttA_2to5_'+channel+systematic))
    ScaleWithError(bkgHist, MCSF, MCSFerror)

    dataHist.Add(qcdHist, -1.0)
    dataHist.Add(bkgHist, -1.0)

    (fitFrac, fitFracErr) = makeFit('m3', xlo, xhi, topHist, wjetsHist, dataHist)

    lowbin = dataHist.FindBin(xlo)
    highbin = dataHist.FindBin(xhi) - 1

    dataInt = dataHist.Integral(lowbin, highbin)
    topInt = topHist.Integral(lowbin, highbin)
    wjetsInt = wjetsHist.Integral(lowbin, highbin)
    
    topSF = fitFrac * dataInt / topInt
    topSFerror = fitFracErr * dataInt / topInt
    
    wjetsSF = (1.0-fitFrac) * dataInt / wjetsInt
    wjetsSFerror = fitFracErr * dataInt / wjetsInt

    output.write(systematic+'\t'+
            str(topSF)+'\t'+
            str(topSFerror)+'\t'+
            str(wjetsSF)+'\t'+
            str(wjetsSFerror)+'\t'+
            str(QCDSF)+'\t'+
            str(QCDSFerror)+'\t'+
            str(MCSF)+'\t'+
            str(MCSFerror)+'\n')

    drawPlots(dataHist, topHist, topSF, 'ttbar', wjetsHist, wjetsSF, 'wjets', xlo, xhi, 'm3_'+channel+systematic)

def doSigmaFit(channel, controlRegion, systematic, output, xlo, xhi):

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    dataHist = get1DHist(input, 'leadSigmaIetaIeta_gg_'+channel)

    topHist = get1DHist(input, 'leadSigmaIetaIeta_ttJetsHadronic_'+channel+systematic)
    topHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ttJetsFullLep_'+channel+systematic))
    topHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ttJetsSemiLep_'+channel+systematic))

    ttgammaHist = get1DHist(input, 'leadSigmaIetaIeta_ttA_2to5_'+channel+systematic)

    bkgHist = get1DHist(input, 'leadSigmaIetaIeta_TBar_s_'+channel+systematic)
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_TBar_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_TBar_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_T_s_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_T_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_T_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_dy1JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_dy2JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_dy3JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_dy4JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_W1JetsToLNu_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_W2JetsToLNu_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_W3JetsToLNu_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_W4JetsToLNu_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_WW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_WZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ZZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_TTWJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_TTZJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ttA_2to5_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_qcd_'+channel)

    dataHist.Add(bkgHist, -1.0)

    (fitFrac, fitFracErr) = makeFit('leadSigmaIetaIeta', xlo, xhi, topHist, ttgammaHist, dataHist)

    lowbin = dataHist.FindBin(xlo)
    highbin = dataHist.FindBin(xhi) - 1

    dataInt = dataHist.Integral(lowbin, highbin)
    topInt = topHist.Integral(lowbin, highbin)
    ttgammaInt = ttgammaHist.Integral(lowbin, highbin)
    
    topSF = fitFrac * dataInt / topInt
    topSFerror = fitFracErr * dataInt / topInt
    
    ttgammaSF = (1.0-fitFrac) * dataInt / ttgammaInt
    ttgammaSFerror = fitFracErr * dataInt / ttgammaInt

    output.write(systematic+'\t'+
            str(topSF)+'\t'+
            str(topSFerror)+'\t'+
            str(ttgammaSF)+'\t'+
            str(ttgammaSFerror)+'\t'+
            str(QCDSF)+'\t'+
            str(QCDSFerror)+'\t'+
            str(MCSF)+'\t'+
            str(MCSFerror)+'\n')

    drawPlots(dataHist, topHist, topSF, 'ttbar', ttgammaHist, ttgammaSF, 'ttgamma', xlo, xhi, 'leadSigmaIetaIeta_'+channel+systematic)
