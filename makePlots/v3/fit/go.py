from fit import *

channels = ['ele_bjj', 'muon_bjj']
channels_noTag = ['ele_jjj', 'muon_jjj']

systematics = ['',
               '_btagWeightUp', '_btagWeightDown',
               '_puWeightUp', '_puWeightDown',
               '_scaleUp', '_scaleDown',
               '_pdfUp', '_pdfDown',
               '_topPtUp', '_topPtDown',
               '_JECUp', '_JECDown',
               '_leptonSFUp', '_leptonSFDown',
               '_photonSFUp', '_photonSFDown',
               '_qcdDefUp', '_qcdDefDown']

dilepRegion = 'Any'
zFitRegion = 'SR1' # maybe SigmaPlot
qcdFitRegion = 'Any'
m3FitRegion = 'Any'
purityFitRegion = 'SigmaPlot'

def goPhotonPurity(version, versionShort, metCutName, channel, xlo, xhi, wjetsResults, topM3Results, dilepResults, eleFakeRateResults):

    promptResults = []
    nonpromptResults = []

    output_prompt = open('photonSigSF_'+versionShort+metCutName+'_'+channel+'.txt', 'w')
    output_nonprompt = open('jetBkgSF_'+versionShort+metCutName+'_'+channel+'.txt', 'w')
    output_purity = open('photonPurity_'+versionShort+metCutName+'_'+channel+'.txt', 'w')

    output_prompt.write('systematic\tSF\tError\n')
    output_nonprompt.write('systematic\tSF\tError\n')
    output_purity.write('systematic\tpurityMC\tError\tpurityData\tError\n')
    
    jsyst = 0
    for systematic in systematics:
        (jetSF, jetSFerror, photonSF, photonSFerror) = doSigmaFitWithMatching(version, channel, purityFitRegion, systematic, metCutName,
                                                                              output_nonprompt, output_prompt, output_purity,
                                                                              xlo, xhi,
                                                                              wjetsResults[jsyst], topM3Results[jsyst], dilepResults[jsyst], eleFakeRateResults[jsyst],
                                                                              2, 2e3, True)
        nonpromptResults.append((jetSF, jetSFerror))
        promptResults.append((photonSF, photonSFerror))
        jsyst += 1

    output_prompt.close()
    output_nonprompt.close()
    output_purity.close()

    jsyst = 0
    for systematic in systematics:
        wigglePurity('pfMET_t01', versionShort, channel, 'SR1', systematic, metCutName,
                     wjetsResults[jsyst], topM3Results[jsyst], dilepResults[jsyst], eleFakeRateResults[jsyst],
                     promptResults[jsyst], nonpromptResults[jsyst])
        fixLimitInputs(channel, 'SR1', systematic, versionShort, metCutName,
                       wjetsResults, topM3Results, dilepResults, eleFakeRateResults, promptResults, nonpromptResults)

        wigglePurity('pfMET_t01', versionShort, channel, 'SR2', systematic, metCutName,
                     wjetsResults[jsyst], topM3Results[jsyst], dilepResults[jsyst], eleFakeRateResults[jsyst],
                     promptResults[jsyst], nonpromptResults[jsyst])
        fixLimitInputs(channel, 'SR2', systematic, versionShort, metCutName,
                       wjetsResults, topM3Results, dilepResults, eleFakeRateResults, promptResults, nonpromptResults)
        
        jsyst += 1

ichan = 0
for channel in channels:

    # Get the di-lepton zmass fit results
    dilepResults = []
    dilepResults_noTag = []

    output_dilepton = open('dilepSF_'+channel+'.txt', 'w')
    output_dilepton.write('systematic\tSF\tError\n')
    
    output_dilepton_noTag = open('dilepSF_'+channels_noTag[ichan]+'.txt', 'w')
    output_dilepton_noTag.write('systematic\tSF\tError\n')

    for systematic in systematics:
        dilepResults.append(doDileptonFit(channel, dilepRegion, systematic, output_dilepton, 0.0, 180.0, 0.9, 2e4, True))
        dilepResults_noTag.append(doDileptonFit(channels_noTag[ichan], dilepRegion, systematic, output_dilepton_noTag, 0.0, 180.0, 0.9, 2e5, True))

    output_dilepton.close()
    output_dilepton_noTag.close()

    # Calculate e --> gamma fake rate sf; do this in jjj

    eleFakeRateResults = []

    output_eleFakeRate = open('eleFakeRateSF_'+channels_noTag[ichan]+'.txt', 'w')
    output_eleFakeRate.write('systematic\tSF\tError\n')

    isyst = 0
    for systematic in systematics:
        if channel.find('ele') == 0:
            eleFakeRateResults.append(doElectronFit(channels_noTag[ichan], zFitRegion, systematic, output_eleFakeRate, 20.0, 180.0, dilepResults_noTag[isyst], 0, 3500, False))
        else:
            eleFakeRateResults.append((1.0, 0.0))
        isyst += 1

    output_eleFakeRate.close()

    # Now calculate QCD sf in kAny for the WJets fit

    qcdResults_kAny = []
    
    output_qcd = open('qcdSF_kAny_'+channel+'.txt', 'w')
    output_qcd.write('systematic\tSF\tError\n')

    output_allMC_qcd = open('qcd_allMCSF_kAny_'+channel+'.txt', 'w')
    output_allMC_qcd.write('systematic\tSF\tError\n')

    isyst = 0
    for systematic in systematics:
        qcdResults_kAny.append(doQCDFit(channel, qcdFitRegion, systematic, output_qcd, output_allMC_qcd, 0.0, 300.0, dilepResults[isyst], 50, 2e5, True))
        isyst += 1

    output_qcd.close()
    output_allMC_qcd.close()

    # Now calculate WJets sf in kAny

    wjetsResults = []
    ttbarResults_M3 = []

    output_wjets = open('wjetsSF_'+channel+'.txt', 'w')
    output_wjets.write('systematic\tSF\tError\n')

    output_ttbar = open('ttbarSF_M3_'+channel+'.txt', 'w')
    output_ttbar.write('systematic\tSF\tError\n')

    isyst = 0
    for systematic in systematics:
        (topSF, topSFerror, wjetsSF, wjetsSFerror) = doM3Fit(channel, m3FitRegion, systematic, output_wjets, output_ttbar, 40.0, 660.0, qcdResults_kAny[isyst], dilepResults[isyst], 1, 8e4, True)
        wjetsResults.append((wjetsSF, wjetsSFerror))
        ttbarResults_M3.append((topSF, topSFerror))
        isyst += 1
                
    output_wjets.close()
    output_ttbar.close()
 
    # Now calculate TTGamma sf in kSigmaPlot for lead sigma ieta ieta
    # using the M3 results above, and just normalizing QCD in MET < 20
    goPhotonPurity('leadChargedHadronIso', 'chHadIso', '', channel, 0.0, 20.0, wjetsResults, ttbarResults_M3, dilepResults, eleFakeRateResults)
    goPhotonPurity('leadChargedHadronIso', 'chHadIso', '_metCut_50', channel, 0.0, 20.0, wjetsResults, ttbarResults_M3, dilepResults, eleFakeRateResults)

    goPhotonPurity('leadSigmaIetaIeta', 'sigma', '', channel, 0.006, 0.02, wjetsResults, ttbarResults_M3, dilepResults, eleFakeRateResults)
    goPhotonPurity('leadSigmaIetaIeta', 'sigma', '_metCut_50', channel, 0.006, 0.02, wjetsResults, ttbarResults_M3, dilepResults, eleFakeRateResults)
    
    ichan += 1


