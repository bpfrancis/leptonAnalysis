from utils import *

def normalizeQCD(dataHist, qcdHist, MCHist):

    lowbin = 1
    highbin = dataHist.FindBin(20.0) - 1

    ndata = dataHist.Integral(lowbin, highbin)
    nqcd = qcdHist.Integral(lowbin, highbin)
    nmc = MCHist.Integral(lowbin, highbin)

    sigma_data = 0
    sigma_qcd = 0
    sigma_mc = 0

    for ibin in range(highbin):
        sigma_data += pow(dataHist.GetBinError(ibin+1), 2)
        sigma_qcd += pow(qcdHist.GetBinError(ibin+1), 2)
        sigma_mc += pow(MCHist.GetBinError(ibin+1), 2)

    sigma_data = math.sqrt(sigma_data)
    sigma_qcd = math.sqrt(sigma_qcd)
    sigma_mc = math.sqrt(sigma_mc)

    qcdScale = (ndata - nmc) / nqcd;
    if qcdScale < 0:
        return (0.0, 0.0)

    qcdScaleError = (pow(sigma_data, 2) + pow(sigma_mc, 2)) / pow(ndata - nmc, 2)
    qcdScaleError += pow(sigma_qcd, 2) / pow(nqcd, 2)
    qcdScaleError = qcdScale * math.sqrt(qcdScaleError)
    
    return (qcdScale, qcdScaleError)

def doQCDFit(channel, controlRegion, systematic, output, xlo, xhi):

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    varName = 'pfMET_t01'

    dataHist = get1DHist(input, varName+'_gg_'+channel)
    qcdHist = get1DHist(input, varName+'_qcd_'+channel)

    MCHist = get1DHist(input, varName+'_ttJetsHadronic_'+channel+systematic)
    MCHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_W3JetsToLNu_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_TBar_s_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_T_s_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_T_t_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_dy1JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_TTGamma_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_WW_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_WZ_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systematic))
    MCHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systematic))

    (qcdFrac, qcdFracErr) = makeFit(varName+'', xlo, xhi, qcdHist, MCHist, dataHist)

    lowbin = dataHist.FindBin(xlo)
    highbin = dataHist.FindBin(xhi) - 1

    dataInt = dataHist.Integral(lowbin, highbin)
    qcdInt = qcdHist.Integral(lowbin, highbin)

    QCDSF = qcdFrac*dataInt/qcdInt
    QCDSFerror = qcdFracErr*dataInt/qcdInt

    (QCDSF, QCDSFerror) = normalizeQCD(channel, controlRegion, systematic)

    output.write(systematic+'\t'+
                 str(QCDSF)+'\t'+
                 str(QCDSFerror)+'\n')

    drawPlots(dataHist, qcdHist, QCDSF, 'QCD', MCHist, 1.0, 'MC', xlo, xhi, varName+'_'+channel+systematic)

    return (QCDSF, QCDSFerror)

def doM3Fit(channel, controlRegion, systematic, output_wjets, output_ttbar, xlo, xhi, qcdSFPair):

    (QCDSF, QCDSFerror) = qcdSFPair

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    dataHist = get1DHist(input, 'm3_gg_'+channel)

    topHist = get1DHist(input, 'm3_ttJetsHadronic_'+channel+systematic)
    topHist.Add(get1DHist(input, 'm3_ttJetsFullLep_'+channel+systematic))
    topHist.Add(get1DHist(input, 'm3_ttJetsSemiLep_'+channel+systematic))

    wjetsHist = get1DHist(input, 'm3_W3JetsToLNu_'+channel+systematic)
    wjetsHist.Add(get1DHist(input, 'm3_W4JetsToLNu_'+channel+systematic))

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
    bkgHist.Add(get1DHist(input, 'm3_TTGamma_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_ZGToLLG_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'm3_WGToLNuG_'+channel+systematic))

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

    output_wjets.write(systematic+'\t'+
                       str(wjetsSF)+'\t'+
                       str(wjetsSFerror)+'\n')

    output_ttbar.write(systematic+'\t'+
                       str(topSF)+'\t'+
                       str(topSFerror)+'\n')

    drawPlots(dataHist, topHist, topSF, 'ttbar', wjetsHist, wjetsSF, 'wjets', xlo, xhi, 'm3_'+channel+systematic)

    return (topSF, topSFerror, wjetsSF, wjetsSFerror)

def doSigmaFit(channel, controlRegion, systematic, output_ttbar, output_ttgamma, xlo, xhi, wjetsResults, topM3Results):

    (wjetsSF, wjetsSFerror) = wjetsResults
    (topM3sf, topM3sfError) = topM3Results

    input = '../histograms_'+channel+'_SigmaPlot.root'

    dataHist = get1DHist(input, 'leadSigmaIetaIeta_gg_'+channel)

    topHist = get1DHist(input, 'leadSigmaIetaIeta_ttJetsHadronic_'+channel+systematic)
    topHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ttJetsFullLep_'+channel+systematic))
    topHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ttJetsSemiLep_'+channel+systematic))
    ScaleWithError(topHist, topM3sf, topM3sfError)

    ttgammaHist = get1DHist(input, 'leadSigmaIetaIeta_TTGamma_'+channel+systematic)

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
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_WW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_WZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ZZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_TTWJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_TTZJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_ZGToLLG_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadSigmaIetaIeta_WGToLNuG_'+channel+systematic))

    wjetsHist = get1DHist(input, 'leadSigmaIetaIeta_W3JetsToLNu_'+channel+systematic)
    wjetsHist.Add(get1DHist(input, 'leadSigmaIetaIeta_W4JetsToLNu_'+channel+systematic))
    ScaleWithError(wjetsHist, wjetsSF, wjetsSFerror)

    qcdHist = get1DHist(input, 'leadSigmaIetaIeta_qcd_'+channel)
    (qcdSF, qcdSFerror) = normalizeQCD(dataHist, qcdHist, MCHist)
    ScaleWithError(qcdHist, qcdSF, qcdSFerror)
    
    dataHist.Add(bkgHist, -1.0)
    dataHist.Add(qcdHist, -1.0)
    dataHist.Add(wjetsHist, -1.0)

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

    output_ttbar.write(systematic+'\t'+
                       str(topSF)+'\t'+
                       str(topSFerror)+'\n')

    output_ttgamma.write(systematic+'\t'+
                         str(ttgammaSF)+'\t'+
                         str(ttgammaSFerror)+'\n')

    drawPlots(dataHist, topHist, topSF, 'ttbar', ttgammaHist, ttgammaSF, 'ttgamma', xlo, xhi, 'leadSigmaIetaIeta_'+channel+systematic)

    return (topSF, topSFerror, ttgammaSF, ttgammaSFerror)

def doChHadIsoFit(channel, controlRegion, systematic, output, xlo, xhi, version):

    input = '../histograms_'+channel+'_'+version+'.root'

    dataHist = get1DHist(input, 'leadChargedHadronIso_gg_'+channel)

    topHist = get1DHist(input, 'leadChargedHadronIso_ttJetsHadronic_'+channel+systematic)
    topHist.Add(get1DHist(input, 'leadChargedHadronIso_ttJetsFullLep_'+channel+systematic))
    topHist.Add(get1DHist(input, 'leadChargedHadronIso_ttJetsSemiLep_'+channel+systematic))

    ttgammaHist = get1DHist(input, 'leadChargedHadronIso_TTGamma_'+channel+systematic)

    bkgHist = get1DHist(input, 'leadChargedHadronIso_TBar_s_'+channel+systematic)
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_TBar_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_TBar_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_T_s_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_T_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_T_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_WW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_WZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_ZZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_TTWJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, 'leadChargedHadronIso_TTZJets_'+channel+systematic))

    wjetsHist = get1DHist(input, 'leadChargedHadronIso_W1JetsToLNu_'+channel+systematic)
    wjetsHist.Add(get1DHist(input, 'leadChargedHadronIso_W2JetsToLNu_'+channel+systematic))
    wjetsHist.Add(get1DHist(input, 'leadChargedHadronIso_W3JetsToLNu_'+channel+systematic))
    wjetsHist.Add(get1DHist(input, 'leadChargedHadronIso_W4JetsToLNu_'+channel+systematic))

    zjetsHist = get1DHist(input, 'leadChargedHadronIso_dy1JetsToLL_'+channel+systematic)
    zjetsHist.Add(get1DHist(input, 'leadChargedHadronIso_dy2JetsToLL_'+channel+systematic))
    zjetsHist.Add(get1DHist(input, 'leadChargedHadronIso_dy3JetsToLL_'+channel+systematic))
    zjetsHist.Add(get1DHist(input, 'leadChargedHadronIso_dy4JetsToLL_'+channel+systematic))

    qcdHist = get1DHist(input, 'leadChargedHadronIso_qcd_'+channel)
    if channel == 'ele_bjj':
        ScaleWithError(qcdHist, 0.57381, 0.079112)
        ScaleWithError(wjetsHist, 1.72687920543, 0.0272087469619)
    if channel == 'muon_bjj':
        ScaleWithError(qcdHist, 0.0144609, 0.0192772)
        ScaleWithError(wjetsHist, 1.50657303747, 0.0263282690732)

    dataHist.Add(bkgHist, -1.0)
    dataHist.Add(qcdHist, -1.0)
    dataHist.Add(wjetsHist, -1.0)

    ttgammaHist.Add(zjetsHist)

    (fitFrac, fitFracErr) = makeFit('leadChargedHadronIso', xlo, xhi, topHist, ttgammaHist, dataHist)

    lowbin = dataHist.FindBin(xlo)
    highbin = dataHist.FindBin(xhi) - 1

    dataInt = dataHist.Integral()
    topInt = topHist.Integral()
    ttgammaInt = ttgammaHist.Integral()

    topSF = fitFrac * dataInt / topInt
    topSFerror = fitFracErr * dataInt / topInt

    ttgammaSF = (1.0-fitFrac) * dataInt / ttgammaInt
    ttgammaSFerror = fitFracErr * dataInt / ttgammaInt

    if systematic == '':
        topInt = topInt	* topSF
        ttgammaInt = ttgammaInt	* ttgammaSF
        #print '\n\nN(top) = '+repr(topInt)+'\nN(ttgamma) = '+repr(ttgammaInt)+'\n\n'

    output.write(systematic+'\t'+
            str(topSF)+'\t'+
            str(topSFerror)+'\t'+
            str(ttgammaSF)+'\t'+
            str(ttgammaSFerror)+'\n')

    drawPlots(dataHist, topHist, topSF, 'ttbar', ttgammaHist, ttgammaSF, 'ttgamma', xlo, xhi, 'leadChargedHadronIso_'+channel+systematic)

    return (fitFrac, topSF, ttgammaSF)
