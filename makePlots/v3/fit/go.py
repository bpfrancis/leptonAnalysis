from fit import *

channels = ['ele_bjj', 'muon_bjj']
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

zFitRegion = 'SR1' # maybe SigmaPlot
qcdFitRegion = 'Any'
m3FitRegion = 'Any'
sigmaFitRegion = 'SigmaPlot' # MET < 50, g/f
chHadIsoFitRegion = 'SigmaPlot'

for channel in channels:

    # First calculate e --> gamma fake rate sf

    zResults = []

    channel_noTag = channel
    if channel.find('ele') == 0:
        channel_noTag = 'ele_jjj'
    else:
        channel_noTag = 'muon_jjj'

    output_z = open('zSF_'+channel_noTag+'.txt', 'w')
    output_z.write('systematic\tSF\tError\n')

    for systematic in systematics:
        if channel.find('ele') == 0:
            zResults.append(doElectronFit(channel_noTag, zFitRegion, systematic, output_z, 20.0, 180.0))
        else:
            zResults.append((1.0, 0.0))

    output_z.close()

    # Now calculate QCD sf in kAny for the WJets fit

    qcdResults_kAny = []
    
    output_qcd = open('qcdSF_kAny_'+channel+'.txt', 'w')
    output_qcd.write('systematic\tSF\tError\n')

    for systematic in systematics:
        qcdResults_kAny.append(doQCDFit(channel, qcdFitRegion, systematic, output_qcd, 0.0, 300.0))

    output_qcd.close()

    # Now calculate WJets sf in kAny

    wjetsResults = []
    ttbarResults_M3 = []

    output_wjets = open('wjetsSF_'+channel+'.txt', 'w')
    output_wjets.write('systematic\tSF\tError\n')

    output_ttbar = open('ttbarSF_M3_'+channel+'.txt', 'w')
    output_ttbar.write('systematic\tSF\tError\n')

    isyst = 0
    for systematic in systematics:
        (topSF, topSFerror, wjetsSF, wjetsSFerror) = doM3Fit(channel, m3FitRegion, systematic, output_wjets, output_ttbar, 40.0, 660.0, qcdResults_kAny[isyst])
        wjetsResults.append((wjetsSF, wjetsSFerror))
        ttbarResults_M3.append((topSF, topSFerror))
        isyst += 1
                
    output_wjets.close()
    output_ttbar.close()
 
    # Now calculate TTGamma sf in kSigmaPlot for lead sigma ieta ieta
    # using the M3 results above, and just normalizing QCD in MET < 20

    ttbarResults_sigma = []
    ttgammaResults_sigma = []

    output_ttjets = open('ttbarSF_sigma_'+channel+'.txt', 'w')
    output_ttjets.write('systematic\tSF\tError\n')

    output_ttgamma = open('ttgammaSF_sigma_'+channel+'.txt', 'w')
    output_ttgamma.write('systematic\tSF\tError\n')

    isyst = 0
    for systematic in systematics:
        (topSF, topSFerror, ttgammaSF, ttgammaSFerror) = doSigmaFit('leadSigmaIetaIeta', channel, sigmaFitRegion, systematic, output_ttjets, output_ttgamma, 0.006, 0.02, wjetsResults[isyst], ttbarResults_M3[isyst], zResults[isyst])
        ttbarResults_sigma.append((topSF, topSFerror))
        ttgammaResults_sigma.append((ttgammaSF, ttgammaSFerror))
        isyst += 1

    output_ttjets.close()
    output_ttgamma.close()

    # Now calculate TTGamma sf in kSigmaPlot for charged hadron isolation
    # using the M3 results above, and just normalizing QCD in MET < 20

    ttbarResults_chHadIso = []
    ttgammaResults_chHadIso = []

    output_ttjets = open('ttbarSF_chHadIso_'+channel+'.txt', 'w')
    output_ttjets.write('systematic\tSF\tError\n')

    output_ttgamma = open('ttgammaSF_chHadIso_'+channel+'.txt', 'w')
    output_ttgamma.write('systematic\tSF\tError\n')

    isyst = 0
    for systematic in systematics:
        (topSF, topSFerror, ttgammaSF, ttgammaSFerror) = doSigmaFit('leadChargedHadronIso', channel, chHadIsoFitRegion, systematic, output_ttjets, output_ttgamma, 0.0, 20.0, wjetsResults[isyst], ttbarResults_M3[isyst], zResults[isyst])
        ttbarResults_chHadIso.append((topSF, topSFerror))
        ttgammaResults_chHadIso.append((ttgammaSF, ttgammaSFerror))
        isyst += 1

    output_ttjets.close()
    output_ttgamma.close()
