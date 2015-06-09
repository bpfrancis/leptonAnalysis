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

def doQCDFit(channel, controlRegion, systematic, output, output_allMC_qcd, xlo, xhi, dilepResults, axisMin, axisMax, doLogy):

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

        #drawPlots(dataHist, qcdHist, QCDSF, 'QCD', MCHist, 1.0, 'MC', xlo, xhi, varName+'_'+channel+systematic, '#slash{E}_{T} (GeV)', axisMin, axisMax, doLogy)
        drawPlots(dataHist, qcdHist, QCDSF, 'QCD', MCHist, bkgSF, 'MC', xlo, xhi, varName+'_'+channel+systematic, '#slash{E}_{T} (GeV)', axisMin, axisMax, doLogy)
    else:
        output.write(systematic+'\t'+
                     str(QCDSF)+'\t'+
                     str(QCDSFerror)+'\n')
        output_allMC_qcd.write(systematic+'\t'+
                               str(bkgSF)+'\t'+
                               str(bkgSFerror)+'\n')

    return (QCDSF, QCDSFerror)

def doM3Fit(channel, controlRegion, systematic, output_wjets, output_ttbar, xlo, xhi, qcdSFPair, dilepPair, axisMin, axisMax, doLogy):

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

        drawPlots(dataHist, topHist, topSF, 't#bar{t} + Jets', wjetsHist, wjetsSF, 'W + Jets', xlo, xhi, 'm3_'+channel+systematic, 'M3 (GeV/c^2)', axisMin, axisMax, doLogy)

    else:
        output_wjets.write(systematic+'\t'+
                           str(wjetsSF)+'\t'+
                           str(wjetsSFerror)+'\n')

        output_ttbar.write(systematic+'\t'+
                           str(topSF)+'\t'+
                           str(topSFerror)+'\n')

    return (topSF, topSFerror, wjetsSF, wjetsSFerror)

def doSigmaFit(varName, channel, controlRegion, systematic, output_ttbar, output_ttgamma, output_purity, xlo, xhi, wjetsResults, topM3Results, dilepResults, eleFakeRateResults, axisMin, axisMax, doLogy):

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
    #dataHist.Add(qcdHist, -1.0)
    dataHist.Add(wjetsHist, -1.0)
    dataHist.Add(zHist, -1.0)

    (dataInt, dataIntError) = integrateError(dataHist, xlo, xhi)
    (topInt, topIntError) = integrateError(topHist, xlo, xhi)
    (ttgammaInt, ttgammaIntError) = integrateError(ttgammaHist, xlo, xhi)

    regularCut = 0.012
    if varName == 'leadChargedHadronIso':
        regularCut = 2.6
    (dataInt_passCut, dataIntError_passCut) = integrateError(dataHist, 0, regularCut)
    (topInt_passCut, topIntError_passCut) = integrateError(topHist, 0, regularCut)
    (ttgammaInt_passCut, ttgammaIntError_passCut) = integrateError(ttgammaHist, 0, regularCut)

    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, ttgammaHist, topHist, dataHist)

    ttgammaSF = fitFrac * dataInt / ttgammaInt
    ttgammaSFerror = ttgammaSF * ( (fitFracErr/fitFrac)**2 + (dataIntError/dataInt)**2 + (ttgammaIntError/ttgammaInt)**2 )**0.5

    topSF = (1.0-fitFrac) * dataInt / topInt
    topSFerror = topSF * ( (fitFracErr/(1.0-fitFrac))**2 + (dataIntError/dataInt)**2 + (topIntError/topInt)**2 )**0.5

    purity = ttgammaSF * ttgammaInt_passCut / (ttgammaSF * ttgammaInt_passCut + topSF * topInt_passCut)
    purityError = ttgammaSFerror*ttgammaSFerror/ttgammaSF/ttgammaSF + topSFerror*topSFerror/topSF/topSF + ttgammaIntError_passCut*ttgammaIntError_passCut/ttgammaInt_passCut/ttgammaInt_passCut + topIntError_passCut*topIntError_passCut/topInt_passCut/topInt_passCut
    purityError = math.sqrt(purityError) * topSF * topInt_passCut
    purityError = purityError / (ttgammaSF * ttgammaInt_passCut + topSF * topInt_passCut)

    if systematic == '':
        output_ttbar.write('central\t'+
                           str(topSF)+'\t'+
                           str(topSFerror)+'\n')

        output_ttgamma.write('central\t'+
                             str(ttgammaSF)+'\t'+
                             str(ttgammaSFerror)+'\n')

        output_purity.write('central\t'+
                            str(purity)+'\t'+
                            str(purityError)+'\n')

        xaxisLabel = '#sigma_{i#eta i#eta}'
        if varName == 'leadChargedHadronIso':
            xaxisLabel = 'Ch. Hadron Iso. (GeV)'

        drawPlots(dataHist, ttgammaHist, ttgammaSF, 't#bar{t} + #gamma', topHist, topSF, 't#bar{t} + Jets', xlo, xhi, varName+'_'+channel+systematic, xaxisLabel, axisMin, axisMax, doLogy)
    else:
        output_ttbar.write(systematic+'\t'+
                           str(topSF)+'\t'+
                           str(topSFerror)+'\n')

        output_ttgamma.write(systematic+'\t'+
                             str(ttgammaSF)+'\t'+
                             str(ttgammaSFerror)+'\n')

        output_purity.write(systematic+'\t'+
                            str(1.0-fitFrac)+'\t'+
                            str(fitFracErr)+'\n')

    return (topSF, topSFerror, ttgammaSF, ttgammaSFerror)

def doElectronFit(channel, controlRegion, systematic, output_z, xlo, xhi, dilepResults, axisMin, axisMax, doLogy):

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

        drawPlots(dataHist, zHist, zSF, 'Z(#gamma) + Jets', bkgHist, bkgSF, 'Bkg', xlo, xhi, varName+'_'+channel+systematic, 'm(l, #gamma) (GeV/c^2)', axisMin, axisMax, doLogy)
    else:
        output_z.write(systematic+'\t'+
                       str(zSF)+'\t'+
                       str(zSFerror)+'\n')

    return (zSF, zSFerror)

def doDileptonFit(channel, controlRegion, systematic, output, xlo, xhi, axisMin, axisMax, doLogy):

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

        drawPlots(dataHist, zHist, zSF, 'Z(#gamma) + Jets', bkgHist, bkgSF, 'Bkg', xlo, xhi, varName+'_'+channel+systName, xaxisLabel, axisMin, axisMax, doLogy)

    else:
        output.write(systematic+'\t'+
                     str(zSF)+'\t'+
                     str(zSFerror)+'\n')

    return (zSF, zSFerror)

def doSigmaFitWithMatching(varName, channel, controlRegion, systematic, output_jet, output_photon, output_purity, xlo, xhi, wjetsResults, topM3Results, dilepResults, eleFakeRateResults, axisMin, axisMax, doLogy):

    (wjetsSF, wjetsSFerror) = wjetsResults
    (topM3sf, topM3sfError) = topM3Results
    (dilepSF, dilepSFerror) = dilepResults
    (eleFakeRateSF, eleFakeRateSFerror) = eleFakeRateResults

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

    dataHist = get1DHist(inputMatched, varName+'_gg_'+channel+'_metCut_50')

    photonHist = get1DHist(inputMatched, varName+'_ttJetsHadronic_'+channel+'_metCut_50_matchPhoton'+systName)
    photonHist.Add(get1DHist(inputMatched, varName+'_ttJetsFullLep_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_ttJetsSemiLep_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_ttJetsHadronic_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_ttJetsFullLep_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_ttJetsSemiLep_'+channel+'_metCut_50_matchElectron'+systName))
    ScaleWithError(photonHist, topM3sf, topM3sfError)

    photonHist.Add(get1DHist(inputMatched, varName+'_TTGamma_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TBar_s_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TBar_t_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TBar_tW_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_T_s_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_T_t_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_T_tW_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_WW_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_WZ_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_ZZ_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TTWJets_'+channel+'_metCut_50_matchPhoton'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TTZJets_'+channel+'_metCut_50_matchPhoton'+systName))

    photonHist.Add(get1DHist(inputMatched, varName+'_TTGamma_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TBar_s_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TBar_t_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TBar_tW_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_T_s_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_T_t_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_T_tW_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_WW_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_WZ_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_ZZ_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TTWJets_'+channel+'_metCut_50_matchElectron'+systName))
    photonHist.Add(get1DHist(inputMatched, varName+'_TTZJets_'+channel+'_metCut_50_matchElectron'+systName))

    wjetsHist = get1DHist(inputMatched, varName+'_W3JetsToLNu_'+channel+'_metCut_50_matchPhoton'+systName)
    wjetsHist.Add(get1DHist(inputMatched, varName+'_W4JetsToLNu_'+channel+'_metCut_50_matchPhoton'+systName))
    wjetsHist.Add(get1DHist(inputMatched, varName+'_W3JetsToLNu_'+channel+'_metCut_50_matchElectron'+systName))
    wjetsHist.Add(get1DHist(inputMatched, varName+'_W4JetsToLNu_'+channel+'_metCut_50_matchElectron'+systName))
    ScaleWithError(wjetsHist, wjetsSF, wjetsSFerror)
    photonHist.Add(wjetsHist)

    zHist = get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_metCut_50_matchPhoton'+systName)
    zHist.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_metCut_50_matchPhoton'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_metCut_50_matchPhoton'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_metCut_50_matchPhoton'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_metCut_50_matchPhoton'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_metCut_50_matchPhoton'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_metCut_50_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_metCut_50_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_metCut_50_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_metCut_50_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_metCut_50_matchElectron'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_metCut_50_matchElectron'+systName))
    ScaleWithError(zHist, dilepSF, dilepSFerror)
    ScaleWithError(zHist, eleFakeRateSF, eleFakeRateSFerror)
    photonHist.Add(zHIst)

    jetHist = get1DHist(inputMatched, varName+'_ttJetsHadronic_'+channel+'_metCut_50_matchJet'+systName)
    jetHist.Add(get1DHist(inputMatched, varName+'_ttJetsFullLep_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_ttJetsSemiLep_'+channel+'_metCut_50_matchJet'+systName))
    ScaleWithError(jetHist, topM3sf, topM3sfError)

    jetHist.Add(get1DHist(inputMatched, varName+'_TTGamma_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_TBar_s_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_TBar_t_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_TBar_tW_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_T_s_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_T_t_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_T_tW_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_WW_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_WZ_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_ZZ_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_TTWJets_'+channel+'_metCut_50_matchJet'+systName))
    jetHist.Add(get1DHist(inputMatched, varName+'_TTZJets_'+channel+'_metCut_50_matchJet'+systName))

    wjetsHist = get1DHist(inputMatched, varName+'_W3JetsToLNu_'+channel+'_metCut_50_matchJet'+systName)
    wjetsHist.Add(get1DHist(inputMatched, varName+'_W4JetsToLNu_'+channel+'_metCut_50_matchJet'+systName))
    ScaleWithError(wjetsHist, wjetsSF, wjetsSFerror)
    jetHist.Add(wjetsHist)

    zHist = get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_metCut_50_matchJet'+systName)
    zHist.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_metCut_50_matchJet'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_metCut_50_matchJet'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_metCut_50_matchJet'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_metCut_50_matchJet'+systName))
    zHist.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_metCut_50_matchJet'+systName))
    ScaleWithError(zHist, dilepSF, dilepSFerror)
    ScaleWithError(zHist, eleFakeRateSF, eleFakeRateSFerror)
    jetHist.Add(zHist)

    (dataInt, dataIntError) = integrateError(dataHist, xlo, xhi)
    (jetInt, jetIntError) = integrateError(jetHist, xlo, xhi)
    (photonInt, photonIntError) = integrateError(photonHist, xlo, xhi)

    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, photonHist, jetHist, dataHist)

    photonSF = fitFrac * dataInt / photonInt
    photonSFerror = photonSF * ( (fitFracErr/fitFrac)**2 + (dataIntError/dataInt)**2 + (photonIntError/photonInt)**2 )**0.5

    jetSF = (1.0-fitFrac) * dataInt / jetInt
    jetSFerror = jetSF * ( (fitFracErr/(1.0-fitFrac))**2 + (dataIntError/dataInt)**2 + (jetIntError/jetInt)**2 )**0.5

    # purity measurement

    regularCut = 0.012
    if varName == 'leadChargedHadronIso':
        regularCut = 2.6

    jetHistScaled = jetHist.Clone('jetHistScaled_'+channel+systName)
    ScaleWithError(jetHistScaled, jetSF, jetSFerror)

    photonHistScaled = photonHist.Clone('photonHistScaled_'+channel+systName)
    ScaleWithError(photonHistScaled, photonSF, photonSFerror)

    (nJetMC, nJetMCErr) = integrateError(jetHist, xlo, regularCut)
    (nJetData, nJetDataErr) = integrateError(jetHistScaled, xlo, regularCut)

    (nPhotonMC, nPhotonMCErr) = integrateError(photonHist, xlo, regularCut)
    (nPhotonData, nPhotonDataErr) = integrateError(photonHistScaled, xlo, regularCut)

    purityMC = nPhotonMC / (nPhotonMC + nJetMC)
    purityMCError = math.sqrt(nPhotonMCErr*nPhotonMCErr/nPhotonMC/nPhotonMC + nJetMCErr*nJetMCErr/nJetMC/nJetMC)
    purityMCError = nPhotonMC*nJetMC*purityMCError/(nPhotonMC + nJetMC)/(nPhotonMC + nJetMC)

    purityData = nPhotonData / (nPhotonData + nJetData)
    purityDataError = math.sqrt(nPhotonDataErr*nPhotonDataErr/nPhotonData/nPhotonData + nJetDataErr*nJetDataErr/nJetData/nJetData)
    purityDataError = nPhotonData*nJetData*purityDataError/(nPhotonData + nJetData)/(nPhotonData + nJetData)

    if systematic == '':
        output_jet.write('central\t'+
                           str(jetSF)+'\t'+
                           str(jetSFerror)+'\n')

        output_photon.write('central\t'+
                             str(photonSF)+'\t'+
                             str(photonSFerror)+'\n')

        output_purity.write('central\t'+
                            str(purityMC)+'\t'+
                            str(purityMCError)+'\t'+
                            str(purityData)+'\t'+
                            str(purityDataError)+'\n')

        xaxisLabel = '#sigma_{i#eta i#eta}'
        if varName == 'leadChargedHadronIso':
            xaxisLabel = 'Ch. Hadron Iso. (GeV)'

        if varName == 'leadSigmaIetaIeta':
            dataHist.Rebin(2)
            photonHist.Rebin(2)
            jetHist.Rebin(2)
        if varName == 'leadChargedHadronIso':
            dataHist.Rebin(4)
            photonHist.Rebin(4)
            jetHist.Rebin(4)

        drawPlots(dataHist, photonHist, photonSF, 'Prompt #gamma', jetHist, jetSF, 'Background', xlo, xhi, varName+'_'+channel+systematic, xaxisLabel, axisMin, axisMax, doLogy)
    else:
        output_jet.write(systematic+'\t'+
                           str(jetSF)+'\t'+
                           str(jetSFerror)+'\n')

        output_photon.write(systematic+'\t'+
                             str(photonSF)+'\t'+
                             str(photonSFerror)+'\n')

        output_purity.write(systematic+'\t'+
                            str(purityMC)+'\t'+
                            str(purityMCError)+'\t'+
                            str(purityData)+'\t'+
                            str(purityDataError)+'\n')


    return (jetSF, jetSFerror, photonSF, photonSFerror)

def wigglePurity(varName, outName, channel, systematic, wjetsResults, topM3Results, dilepResults, eleFakeRateResults, promptResults, nonpromptResults):

    (wjetsSF, wjetsSFerror) = wjetsResults
    (topM3sf, topM3sfError) = topM3Results
    (dilepSF, dilepSFerror) = dilepResults
    (eleFakeRateSF, eleFakeRateSFerror) = eleFakeRateResults

    (promptSF, promptSFerror) = promptResults
    (jetSF, jetSFerror) = nonpromptResults

    inputMatched = '../fitTemplates.root'

    if systematic == '_qcdDefUp':
        systName = ''
    elif systematic == '_qcdDefDown':
        systName = ''

    # ttjets

    photons = get1DHist(inputMatched, varName+'_ttJetsHadronic_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_ttJetsFullLep_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_ttJetsSemiLep_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_ttJetsHadronic_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_ttJetsFullLep_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_ttJetsSemiLep_'+channel+'_matchElectron'+systName))
    ScaleWithError(photons, topM3sf, topM3sfError)

    jets = get1DHist(inputMatched, varName+'_ttJetsHadronic_'+channel+'_matchJet'+systName)
    jets.Add(get1DHist(inputMatched, varName+'_ttJetsFullLep_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_ttJetsSemiLep_'+channel+'_matchJet'+systName))
    ScaleWithError(jets, topM3sf, topM3sfError)

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_ttjets = after / before
    scaleError_ttjets = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # wjets
    
    photons = get1DHist(inputMatched, varName+'_W3JetsToLNu_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_W4JetsToLNu_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_W3JetsToLNu_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_W4JetsToLNu_'+channel+'_matchElectron'+systName))
    ScaleWithError(photons, wjetsSF, wjetsSFerror)
    
    jets = get1DHist(inputMatched, varName+'_W3JetsToLNu_'+channel+'_matchJet'+systName)
    jets.Add(get1DHist(inputMatched, varName+'_W4JetsToLNu_'+channel+'_matchJet'+systName))
    ScaleWithError(jets, wjetsSF, wjetsSFerror)

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_wjets = after / before
    scaleError_wjets = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # zjets
    
    photons = get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_matchElectron'+systName))
    ScaleWithError(photons, dilepSF, dilepSFerror)
    ScaleWithError(photons, eleFakeRateSF, eleFakeRateSFerror)
    
    jets = get1DHist(inputMatched, varName+'_dy1JetsToLL_'+channel+'_matchJet'+systName)
    jets.Add(get1DHist(inputMatched, varName+'_dy2JetsToLL_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_dy3JetsToLL_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_dy4JetsToLL_'+channel+'_matchJet'+systName))
    ScaleWithError(jets, dilepSF, dilepSFerror)
    ScaleWithError(jets, eleFakeRateSF, eleFakeRateSFerror)

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_zjets = after / before
    scaleError_zjets = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # singleTop

    photons = get1DHist(inputMatched, varName+'_TBar_s_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_TBar_t_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_TBar_tW_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_T_s_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_T_t_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_T_tW_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_TBar_s_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_TBar_t_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_TBar_tW_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_T_s_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_T_t_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_T_tW_'+channel+'_matchElectron'+systName))

    jets = get1DHist(inputMatched, varName+'_TBar_s_'+channel+'_matchJet'+systName)
    jets.Add(get1DHist(inputMatched, varName+'_TBar_t_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_TBar_tW_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_T_s_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_T_t_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_T_tW_'+channel+'_matchJet'+systName))

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_singleTop = after / before
    scaleError_singleTop = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # diboson

    photons = get1DHist(inputMatched, varName+'_WW_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_WZ_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_ZZ_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_WW_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_WZ_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_ZZ_'+channel+'_matchElectron'+systName))

    jets = get1DHist(inputMatched, varName+'_WW_'+channel+'_matchJet'+systName)
    jets.Add(get1DHist(inputMatched, varName+'_WZ_'+channel+'_matchJet'+systName))
    jets.Add(get1DHist(inputMatched, varName+'_ZZ_'+channel+'_matchJet'+systName))

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_diboson = after / before
    scaleError_diboson = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # ttW

    photons = get1DHist(inputMatched, varName+'_TTWJets_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_TTWJets_'+channel+'_matchElectron'+systName))

    jets = get1DHist(inputMatched, varName+'_TTWJets_'+channel+'_matchJet'+systName)

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_ttW = after / before
    scaleError_ttW = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # ttZ

    photons = get1DHist(inputMatched, varName+'_TTZJets_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_TTZJets_'+channel+'_matchElectron'+systName))

    jets = get1DHist(inputMatched, varName+'_TTZJets_'+channel+'_matchJet'+systName)

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_ttZ = after / before
    scaleError_ttZ = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # ttgamma

    photons = get1DHist(inputMatched, varName+'_TTGamma_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_TTGamma_'+channel+'_matchElectron'+systName))

    jets = get1DHist(inputMatched, varName+'_TTGamma_'+channel+'_matchJet'+systName)

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_ttgamma = after / before
    scaleError_ttgamma = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # vgamma

    photons = get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_matchPhoton'+systName)
    photons.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_matchPhoton'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_matchElectron'+systName))
    photons.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_matchElectron'+systName))
    ScaleWithError(photons, dilepSF, dilepSFerror)
    ScaleWithError(photons, eleFakeRateSF, eleFakeRateSFerror)

    jets = get1DHist(inputMatched, varName+'_ZGToLLG_'+channel+'_matchJet'+systName)
    jets.Add(get1DHist(inputMatched, varName+'_WGToLNuG_'+channel+'_matchJet'+systName))
    ScaleWithError(jets, dilepSF, dilepSFerror)
    ScaleWithError(jets, eleFakeRateSF, eleFakeRateSFerror)

    (before, beforeErr) = integrateErrorAll(photons, jets)
    ScaleWithError(photons, promptSF, promptSFerror)
    ScaleWithError(jets, jetSF, jetSFerror)
    (after, afterErr) = integrateErrorAll(photons, jets)

    scale_vgamma = after / before
    scaleError_vgamma = scale * ( (afterErr/after)**2 + (beforeErr/before)**2 )**0.5

    # write output

    output_ttjets = open('puritySF_ttjets_'+outName+'_'+channel+'.txt', 'a')
    output_wjets = open('puritySF_wjets_'+outName+'_'+channel+'.txt', 'a')
    output_zjets = open('puritySF_zjets_'+outName+'_'+channel+'.txt', 'a')
    output_singleTop = open('puritySF_singleTop_'+outName+'_'+channel+'.txt', 'a')
    output_diboson = open('puritySF_diboson_'+outName+'_'+channel+'.txt', 'a')
    output_ttW = open('puritySF_ttW_'+outName+'_'+channel+'.txt', 'a')
    output_ttZ = open('puritySF_ttZ_'+outName+'_'+channel+'.txt', 'a')
    output_ttgamma = open('puritySF_ttgamma_'+outName+'_'+channel+'.txt', 'a')
    output_vgamma = open('puritySF_vgamma_'+outName+'_'+channel+'.txt', 'a')

    firstLine = systematic
    if systematic == '':
        firstLine = 'central'
        output_ttjets.write('systematic\tSF\tError\n')
        output_wjets.write('systematic\tSF\tError\n')
        output_zjets.write('systematic\tSF\tError\n')
        output_singleTop.write('systematic\tSF\tError\n')
        output_diboson.write('systematic\tSF\tError\n')
        output_ttW.write('systematic\tSF\tError\n')
        output_ttZ.write('systematic\tSF\tError\n')
        output_ttgamma.write('systematic\tSF\tError\n')
        output_vgamma.write('systematic\tSF\tError\n')

    output_ttjets.write(firstLine+'\t'+
                        str(scale_ttjets)+'\t'+
                        str(scaleError_ttjets)+'\n')

    output_wjets.write(firstLine+'\t'+
                       str(scale_wjets)+'\t'+
                       str(scaleError_wjets)+'\n')

    output_zjets.write(firstLine+'\t'+
                       str(scale_zjets)+'\t'+
                       str(scaleError_zjets)+'\n')

    output_singleTop.write(firstLine+'\t'+
                           str(scale_singleTop)+'\t'+
                           str(scaleError_singleTop)+'\n')

    output_diboson.write(firstLine+'\t'+
                         str(scale_diboson)+'\t'+
                         str(scaleError_diboson)+'\n')

    output_ttW.write(firstLine+'\t'+
                     str(scale_ttW)+'\t'+
                     str(scaleError_ttW)+'\n')

    output_ttZ.write(firstLine+'\t'+
                     str(scale_ttZ)+'\t'+
                     str(scaleError_ttZ)+'\n')

    output_ttgamma.write(firstLine+'\t'+
                         str(scale_ttgamma)+'\t'+
                         str(scaleError_ttgamma)+'\n')

    output_vgamma.write(firstLine+'\t'+
                        str(scale_vgamma)+'\t'+
                        str(scaleError_vgamma)+'\n')
