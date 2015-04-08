from utils import *

def normalizeQCD(input, channel, systematic):

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

def doSigmaFit(varName, channel, controlRegion, systematic, output_ttbar, output_ttgamma, xlo, xhi, wjetsResults, topM3Results):

    (wjetsSF, wjetsSFerror) = wjetsResults
    (topM3sf, topM3sfError) = topM3Results

    input = '../histograms_'+channel+'_SigmaPlot.root'

    dataHist = get1DHist(input, varName+'_gg_'+channel)

    topHist = get1DHist(input, varName+'_ttJetsHadronic_'+channel+systematic)
    topHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systematic))
    topHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systematic))
    ScaleWithError(topHist, topM3sf, topM3sfError)

    ttgammaHist = get1DHist(input, varName+'_TTGamma_'+channel+systematic)

    bkgHist = get1DHist(input, varName+'_TBar_s_'+channel+systematic)
    bkgHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_s_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_dy1JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_WW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_WZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systematic))

    wjetsHist = get1DHist(input, varName+'_W3JetsToLNu_'+channel+systematic)
    wjetsHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systematic))
    ScaleWithError(wjetsHist, wjetsSF, wjetsSFerror)

    MCHist = bkgHist.Clone('bkgClone')
    MCHist.Add(topHist)
    MCHist.Add(ttgammaHist)
    MCHist.Add(wjetsHist)

    qcdHist = get1DHist(input, varName+'_qcd_'+channel)
    (qcdSF, qcdSFerror) = normalizeQCD(input, channel, systematic)
    ScaleWithError(qcdHist, qcdSF, qcdSFerror)
    
    dataHist.Add(bkgHist, -1.0)
    dataHist.Add(qcdHist, -1.0)
    dataHist.Add(wjetsHist, -1.0)

    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, topHist, ttgammaHist, dataHist)

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

    drawPlots(dataHist, topHist, topSF, 'ttbar', ttgammaHist, ttgammaSF, 'ttgamma', xlo, xhi, varName+'_'+channel+systematic)

    return (topSF, topSFerror, ttgammaSF, ttgammaSFerror)

def doElectronFit(channel, controlRegion, systematic, output_z, xlo, xhi):

    varName = 'mLepGammaLead'

    input = '../histograms_'+channel+'_SigmaPlot.root'

    dataHist = get1DHist(input, varName+'_gg_'+channel)

    zHist = get1DHist(input, varName+'_dy1JetsToLL_'+channel+systematic)
    zHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systematic))

    bkgHist = get1DHist(input, varName+'_TBar_s_'+channel+systematic)
    bkgHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_s_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_WW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_WZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsHadronic_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTGamma_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_W3JetsToLNu_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systematic))

    MCHist = bkgHist.Clone('bkgClone')
    MCHist.Add(zHist)

    qcdHist = get1DHist(input, varName+'_qcd_'+channel)
    (qcdSF, qcdSFerror) = normalizeQCD(input, channel, systematic)
    ScaleWithError(qcdHist, qcdSF, qcdSFerror)

    bkgHist.Add(qcdHist)
    
    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, zHist, bkgHist, dataHist)

    lowbin = dataHist.FindBin(xlo)
    highbin = dataHist.FindBin(xhi) - 1
    
    dataInt = dataHist.Integral(lowbin, highbin)
    zInt = zHist.Integral(lowbin, highbin)
    bkgInt = bkgHist.Integral(lowbin, highbin)

    zSF = fitFrac * dataInt / zInt
    zSFerror = fitFracErr * dataInt / zInt

    bkgSF = (1.0-fitFrac) * dataInt / bkgInt
    bkgSFerror = fitFracErr * dataInt / bkgInt

    output_z.write(systematic+'\t'+
                   str(zSF)+'\t'+
                   str(zSFerror)+'\n')

    drawPlots(dataHist, zHist, zSF, 'ttbar', bkgHist, bkgSF, 'bkg', xlo, xhi, varName+'_'+channel+systematic)

    return (zSF, zSFerror)

