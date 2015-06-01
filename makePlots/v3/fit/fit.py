from utils import *

def normalizeQCD(input, channel, systematic, wjetsResults, topM3Results, dilepResults, eleFakeRateResults):

    (wjetsSF, wjetsSFerror) = wjetsResults
    (topM3sf, topM3sfError) = topM3Results
    (dilepSF, dilepSFerror) = dilepResults
    (eleFakeRateSF, eleFakeRateSFerror) = eleFakeRateResults

    varName = 'pfMET_t01'

    qcdName = '_qcd_'
    systName = systematic

    if systematic == '_qcdDefUp':
        qcdName = '_qcd_relIso_10_'
        systName = ''
    elif systematic == '_qcdDefDown':
        qcdName = '_qcd_relIso_m10_'
        systName = ''

    dataHist = get1DHist(input, varName+'_gg_'+channel)
    qcdHist = get1DHist(input, varName+qcdName+channel)

    topHist = get1DHist(input, varName+'_ttJetsHadronic_'+channel+systName)
    topHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systName))
    topHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systName))
    ScaleWithError(topHist, topM3sf, topM3sfError)

    wjetsHist = get1DHist(input, varName+'_W3JetsToLNu_'+channel+systName)
    wjetsHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systName))
    ScaleWithError(wjetsHist, wjetsSF, wjetsSFerror)

    zHist = get1DHist(input, varName+'_dy1JetsToLL_'+channel+systName)
    zHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systName))
    ScaleWithError(zHist, dilepSF, dilepSFerror)
    ScaleWithError(zHist, eleFakeRateSF, eleFakeRateSFerror)

    MCHist = get1DHist(input, varName+'_TBar_s_'+channel+systName)
    MCHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_T_s_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_T_t_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TTGamma_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_WW_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_WZ_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systName))

    MCHist.Add(topHist)
    MCHist.Add(wjetsHist)
    MCHist.Add(zHist)

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

def doQCDFit(channel, controlRegion, systematic, output, output_allMC_qcd, xlo, xhi, dilepResults):

    (dilepSF, dilepSFerror) = dilepResults

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    varName = 'pfMET_t01'

    qcdName = '_qcd_'
    systName = systematic

    if systematic == '_qcdDefUp':
        qcdName = '_qcd_relIso_10_'
        systName = ''
    elif systematic == '_qcdDefDown':
        qcdName = '_qcd_relIso_m10_'
        systName = ''

    dataHist = get1DHist(input, varName+'_gg_'+channel)
    qcdHist = get1DHist(input, varName+qcdName+channel)

    MCHist = get1DHist(input, varName+'_dy1JetsToLL_'+channel+systName)
    MCHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systName))
    ScaleWithError(MCHist, dilepSF, dilepSFerror)

    MCHist.Add(get1DHist(input, varName+'_ttJetsHadronic_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_W3JetsToLNu_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TBar_s_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_T_s_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_T_t_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TTGamma_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_WW_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_WZ_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systName))
    MCHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systName))

    (dataInt, dataIntError) = integrateError(dataHist, xlo, xhi)
    (qcdInt, qcdIntError) = integrateError(qcdHist, xlo, xhi)
    (bkgInt, bkgIntError) = integrateError(MCHist, xlo, xhi)

    (qcdFrac, qcdFracErr) = makeFit(varName+'', xlo, xhi, qcdHist, MCHist, dataHist)

    QCDSF = qcdFrac * dataInt / qcdInt
    QCDSFerror = QCDSF * ( (qcdFracErr/qcdFrac)**2 + (dataIntError/dataInt)**2 + (qcdIntError/qcdInt)**2 )**0.5

    bkgSF = (1.0-qcdFrac) * dataInt / bkgInt
    bkgSFerror = bkgSF * ( (qcdFracErr/(1.0-qcdFrac))**2 + (dataIntError/dataInt)**2 + (bkgIntError/bkgInt)**2 )**0.5

    if systematic == '':
        output.write('central\t'+
                     str(QCDSF)+'\t'+
                     str(QCDSFerror)+'\n')
        output_allMC_qcd.write('central\t'+
                               str(bkgSF)+'\t'+
                               str(bkgSFerror)+'\n')

        #drawPlots(dataHist, qcdHist, QCDSF, 'QCD', MCHist, 1.0, 'MC', xlo, xhi, varName+'_'+channel+systematic, '#slash{E}_{T} (GeV)')
        drawPlots(dataHist, qcdHist, QCDSF, 'QCD', MCHist, bkgSF, 'MC', xlo, xhi, varName+'_'+channel+systematic, '#slash{E}_{T} (GeV)')
    else:
        output.write(systematic+'\t'+
                     str(QCDSF)+'\t'+
                     str(QCDSFerror)+'\n')
        output_allMC_qcd.write(systematic+'\t'+
                               str(bkgSF)+'\t'+
                               str(bkgSFerror)+'\n')

    return (QCDSF, QCDSFerror)

def doM3Fit(channel, controlRegion, systematic, output_wjets, output_ttbar, xlo, xhi, qcdSFPair, dilepPair):

    (QCDSF, QCDSFerror) = qcdSFPair
    (dilepSF, dilepSFerror) = dilepPair

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    qcdName = '_qcd_'
    systName = systematic

    if systematic == '_qcdDefUp':
        qcdName = '_qcd_relIso_10_'
        systName = ''
    elif systematic == '_qcdDefDown':
        qcdName = '_qcd_relIso_m10_'
        systName = ''

    dataHist = get1DHist(input, 'm3_gg_'+channel)

    topHist = get1DHist(input, 'm3_ttJetsHadronic_'+channel+systName)
    topHist.Add(get1DHist(input, 'm3_ttJetsFullLep_'+channel+systName))
    topHist.Add(get1DHist(input, 'm3_ttJetsSemiLep_'+channel+systName))

    wjetsHist = get1DHist(input, 'm3_W3JetsToLNu_'+channel+systName)
    wjetsHist.Add(get1DHist(input, 'm3_W4JetsToLNu_'+channel+systName))

    qcdHist = get1DHist(input, 'm3'+qcdName+channel)
    ScaleWithError(qcdHist, QCDSF, QCDSFerror)

    bkgHist = get1DHist(input, 'm3_dy1JetsToLL_'+channel+systName)
    bkgHist.Add(get1DHist(input, 'm3_dy2JetsToLL_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_dy3JetsToLL_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_dy4JetsToLL_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_ZGToLLG_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_WGToLNuG_'+channel+systName))
    ScaleWithError(bkgHist, dilepSF, dilepSFerror)

    bkgHist.Add(get1DHist(input, 'm3_TBar_s_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_TBar_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_TBar_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_T_s_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_T_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_T_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_WW_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_WZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_ZZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_TTWJets_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_TTZJets_'+channel+systName))
    bkgHist.Add(get1DHist(input, 'm3_TTGamma_'+channel+systName))
   
    dataHist.Add(qcdHist, -1.0)
    dataHist.Add(bkgHist, -1.0)

    (dataInt, dataIntError) = integrateError(dataHist, xlo, xhi)
    (topInt, topIntError) = integrateError(topHist, xlo, xhi)
    (wjetsInt, wjetsIntError) = integrateError(wjetsHist, xlo, xhi)

    (fitFrac, fitFracErr) = makeFit('m3', xlo, xhi, topHist, wjetsHist, dataHist)

    topSF = fitFrac * dataInt / topInt
    topSFerror = topSF * ( (fitFracErr/fitFrac)**2 + (dataIntError/dataInt)**2 + (topIntError/topInt)**2 )**0.5

    wjetsSF = (1.0-fitFrac) * dataInt / wjetsInt
    wjetsSFerror = wjetsSF * ( (fitFracErr/(1.0-fitFrac))**2 + (dataIntError/dataInt)**2 + (wjetsIntError/wjetsInt)**2 )**0.5

    if systematic == '':
        output_wjets.write('central\t'+
                           str(wjetsSF)+'\t'+
                           str(wjetsSFerror)+'\n')

        output_ttbar.write('central\t'+
                           str(topSF)+'\t'+
                           str(topSFerror)+'\n')

        drawPlots(dataHist, topHist, topSF, 't#bar{t} + Jets', wjetsHist, wjetsSF, 'W + Jets', xlo, xhi, 'm3_'+channel+systematic, 'M3 (GeV/c^2)')

    else:
        output_wjets.write(systematic+'\t'+
                           str(wjetsSF)+'\t'+
                           str(wjetsSFerror)+'\n')

        output_ttbar.write(systematic+'\t'+
                           str(topSF)+'\t'+
                           str(topSFerror)+'\n')

    return (topSF, topSFerror, wjetsSF, wjetsSFerror)

def doSigmaFit(varName, channel, controlRegion, systematic, output_ttbar, output_ttgamma, xlo, xhi, wjetsResults, topM3Results, dilepResults, eleFakeRateResults):

    (wjetsSF, wjetsSFerror) = wjetsResults
    (topM3sf, topM3sfError) = topM3Results
    (dilepSF, dilepSFerror) = dilepResults
    (eleFakeRateSF, eleFakeRateSFerror) = eleFakeRateResults

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    qcdName = '_qcd_'
    systName = systematic

    if systematic == '_qcdDefUp':
        qcdName = '_qcd_relIso_10_'
        systName = ''
    elif systematic == '_qcdDefDown':
        qcdName = '_qcd_relIso_m10_'
        systName = ''

    dataHist = get1DHist(input, varName+'_gg_'+channel)

    topHist = get1DHist(input, varName+'_ttJetsHadronic_'+channel+systName)
    topHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systName))
    topHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systName))
    ScaleWithError(topHist, topM3sf, topM3sfError)

    ttgammaHist = get1DHist(input, varName+'_TTGamma_'+channel+systName)

    bkgHist = get1DHist(input, varName+'_TBar_s_'+channel+systName)
    bkgHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_s_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_WW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_WZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systName))
    
    wjetsHist = get1DHist(input, varName+'_W3JetsToLNu_'+channel+systName)
    wjetsHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systName))
    ScaleWithError(wjetsHist, wjetsSF, wjetsSFerror)

    zHist = get1DHist(input, varName+'_dy1JetsToLL_'+channel+systName)
    zHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systName))
    ScaleWithError(zHist, dilepSF, dilepSFerror)
    ScaleWithError(zHist, eleFakeRateSF, eleFakeRateSFerror)

    qcdHist = get1DHist(input, varName+qcdName+channel)
    (qcdSF, qcdSFerror) = normalizeQCD(input, channel, systematic, wjetsResults, topM3Results, dilepResults, eleFakeRateResults)
    ScaleWithError(qcdHist, qcdSF, qcdSFerror)
    
    dataHist.Add(bkgHist, -1.0)
    dataHist.Add(qcdHist, -1.0)
    dataHist.Add(wjetsHist, -1.0)
    dataHist.Add(zHist, -1.0)

    (dataInt, dataIntError) = integrateError(dataHist, xlo, xhi)
    (topInt, topIntError) = integrateError(topHist, xlo, xhi)
    (ttgammaInt, ttgammaIntError) = integrateError(ttgammaHist, xlo, xhi)

    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, topHist, ttgammaHist, dataHist)

    topSF = fitFrac * dataInt / topInt
    topSFerror = topSF * ( (fitFracErr/fitFrac)**2 + (dataIntError/dataInt)**2 + (topIntError/topInt)**2 )**0.5

    ttgammaSF = (1.0-fitFrac) * dataInt / ttgammaInt
    ttgammaSFerror = ttgammaSF * ( (fitFracErr/(1.0-fitFrac))**2 + (dataIntError/dataInt)**2 + (ttgammaIntError/ttgammaInt)**2 )**0.5

    if systematic == '':
        output_ttbar.write('central\t'+
                           str(topSF)+'\t'+
                           str(topSFerror)+'\n')

        output_ttgamma.write('central\t'+
                             str(ttgammaSF)+'\t'+
                             str(ttgammaSFerror)+'\n')

        xaxisLabel = '#sigma_{i#eta i#eta}'
        if varName == 'leadChargedHadronIso':
            xaxisLabel = 'Ch. Hadron Iso. (GeV)'

        drawPlots(dataHist, topHist, topSF, 't#bar{t} + Jets', ttgammaHist, ttgammaSF, 't#bar{t} + #gamma', xlo, xhi, varName+'_'+channel+systematic, xaxisLabel)
    else:
        output_ttbar.write(systematic+'\t'+
                           str(topSF)+'\t'+
                           str(topSFerror)+'\n')

        output_ttgamma.write(systematic+'\t'+
                             str(ttgammaSF)+'\t'+
                             str(ttgammaSFerror)+'\n')

    return (topSF, topSFerror, ttgammaSF, ttgammaSFerror)

def doElectronFit(channel, controlRegion, systematic, output_z, xlo, xhi, dilepResults):

    (dilepSF, dilepSFerror) = dilepResults

    varName = 'mLepGammaLead'

    inputMatched = '../fitTemplates.root'
    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    qcdName = '_qcd_'
    systName = systematic

    if systematic == '_qcdDefUp':
        qcdName = '_qcd_relIso_10_'
        systName = ''
    elif systematic == '_qcdDefDown':
        qcdName = '_qcd_relIso_m10_'
        systName = ''

    dataHist = get1DHist(input, varName+'_gg_'+channel)

    zHist = get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_matchElectron'+systName)
    zHist.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_matchElectron'+systName))
    ScaleWithError(zHist, dilepSF, dilepSFerror)

    bkgHist = get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_matchPhoton'+systName)
    bkgHist.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_matchPhoton'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_matchPhoton'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_matchPhoton'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_matchPhoton'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_matchPhoton'+systName))

    bkgHist.Add(get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_matchJet'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_matchJet'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_matchJet'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_matchJet'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_matchJet'+systName))
    bkgHist.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_matchJet'+systName))

    ScaleWithError(bkgHist, dilepSF, dilepSFerror)
    
    bkgHist.Add(get1DHist(input, varName+'_TBar_s_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_s_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_WW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_WZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsHadronic_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTGamma_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_W3JetsToLNu_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systName))

    qcdHist = get1DHist(input, varName+qcdName+channel)

    dummyResults = (1.0, 0.0)

    (qcdSF, qcdSFerror) = normalizeQCD(input, channel, systematic, dummyResults, dummyResults, dilepResults, dummyResults)
    ScaleWithError(qcdHist, qcdSF, qcdSFerror)

    bkgHist.Add(qcdHist)

    (dataInt, dataIntError) = integrateError(dataHist, xlo, xhi)
    (zInt, zIntError) = integrateError(zHist, xlo, xhi)
    (bkgInt, bkgIntError) = integrateError(bkgHist, xlo, xhi)
    
    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, zHist, bkgHist, dataHist)

    zSF = fitFrac * dataInt / zInt
    zSFerror = zSF * ( (fitFracErr/fitFrac)**2 + (dataIntError/dataInt)**2 + (zIntError/zInt)**2 )**0.5

    bkgSF = (1.0-fitFrac) * dataInt / bkgInt
    bkgSFerror = bkgSF * ( (fitFracErr/(1.0-fitFrac))**2 + (dataIntError/dataInt)**2 + (bkgIntError/bkgInt)**2 )**0.5

    if systematic == '':
        output_z.write('central\t'+
                       str(zSF)+'\t'+
                       str(zSFerror)+'\n')

        drawPlots(dataHist, zHist, zSF, 'Z + Jets', bkgHist, bkgSF, 'Bkg', xlo, xhi, varName+'_'+channel+systematic, 'm(l, #gamma) (GeV/c^2)')
    else:
        output_z.write(systematic+'\t'+
                       str(zSF)+'\t'+
                       str(zSFerror)+'\n')

    return (zSF, zSFerror)

def doDileptonFit(channel, controlRegion, systematic, output, xlo, xhi):

    input = '../../zgamma/histograms_'+channel+'_'+controlRegion+'.root'

    varName = 'z_mass'

    dataHist = get1DHist(input, varName+'_gg_'+channel)

    systName = systematic
    if systematic == '_qcdDefUp' or systematic == '_qcdDefDown':
        systName = ''

    zHist = get1DHist(input, varName+'_dy1JetsToLL_'+channel+systName)
    zHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systName))
    zHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systName))

    bkgHist = get1DHist(input, varName+'_TBar_s_'+channel+systName)
    bkgHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_s_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_t_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_WW_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_WZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsHadronic_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_TTGamma_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_W3JetsToLNu_'+channel+systName))
    bkgHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systName))

    (dataInt, dataIntError) = integrateError(dataHist, xlo, xhi)
    (zInt, zIntError) = integrateError(zHist, xlo, xhi)
    (bkgInt, bkgIntError) = integrateError(bkgHist, xlo, xhi)

    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, zHist, bkgHist, dataHist)

    zSF = fitFrac * dataInt / zInt
    zSFerror = zSF * ( (fitFracErr/fitFrac)**2 + (dataIntError/dataInt)**2 + (zIntError/zInt)**2 )**0.5

    bkgSF = (1.0-fitFrac) * dataInt / bkgInt
    bkgSFerror = bkgSF * ( (fitFracErr/(1.0-fitFrac))**2 + (dataIntError/dataInt)**2 + (bkgIntError/bkgInt)**2 )**0.5

    if systematic == '':
        output.write('central\t'+
                     str(zSF)+'\t'+
                     str(zSFerror)+'\n')

        xaxisLabel = 'm(ee) (GeV/c^2)'
        if channel == 'muon_bjj':
            xaxisLabel = 'm(#mu#mu) (GeV/c^2)'

        drawPlots(dataHist, zHist, zSF, 'Z + Jets', bkgHist, bkgSF, 'Bkg', xlo, xhi, varName+'_'+channel+systName, xaxisLabel)

    else:
        output.write(systematic+'\t'+
                     str(zSF)+'\t'+
                     str(zSFerror)+'\n')

    return (zSF, zSFerror)
