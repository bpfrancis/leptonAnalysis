from fit import *
import math

controlRegion = 'SigmaPlot'

#sigmaVersion = 'noSigma'
#sigmaVersion = 'noChHadIso_fromSuperFake'
sigmaVersion = 'noNeither_fromSuperFake'
#sigmaVersion = 'noSigma_fromSuperFake'
#sigmaVersion = 'superFake'

#chHadVersion = 'noChHadIso_fromSuperFake'
#chHadVersion = 'noNeither_fromSuperFake'
chHadVersion = 'superFake'

channels = ['ele_bjj', 'muon_bjj']
systematics = ['',
               '_btagWeightUp', '_btagWeightDown',
               '_puWeightUp', '_puWeightDown',
               '_scaleUp', '_scaleDown',
               '_pdfUp', '_pdfDown',
               '_topPtUp', '_topPtDown',
               '_JECUp', '_JECDown',
               '_leptonSFUp', '_leptonSFDown',
               '_photonSFUp', '_photonSFDown']

systematics = ['']

sigmaFitFrac = 0
sigmaTTJetsSF = 0
sigmaTTGammaSF = 0

chHadFitFrac = 0
chHadTTJetsSF = 0
chHadTTGammaSF = 0

for channel in channels:
    sigmaOutput = open('sigmaFitResults_'+channel+'.txt', 'w')
    sigmaOutput.write('systematic\ttopSF\ttopSFerror\tttgammaSF\tttgammaSFerror\n')

    for systematic in systematics:
        if systematic == '':
		(sigmaFitFrac, sigmaTTJetsSF, sigmaTTGammaSF) = doSigmaFit(channel, controlRegion, systematic, sigmaOutput, 0.006, 0.02, sigmaVersion)
	else :
		(fitFrac, ttjetsSF, ttgammaSF) = doSigmaFit(channel, controlRegion, systematic, sigmaOutput, 0.006, 0.02, sigmaVersion)

    chHadOutput = open('chHadIsoFitResults_'+channel+'.txt', 'w')
    chHadOutput.write('systematic\ttopSF\ttopSFerror\tttgammaSF\tttgammaSFerror\n')

    for systematic in systematics:
	if systematic == '':
		(chHadFitFrac, chHadTTJetsSF, chHadTTGammaSF) = doChHadIsoFit(channel, controlRegion, systematic, chHadOutput, 0.0, 20.0, chHadVersion)
	else :
		(fitFrac, ttjetsSF, ttgammaSF) = doChHadIsoFit(channel, controlRegion, systematic, chHadOutput, 0.0, 20.0, chHadVersion)

    sigmaOutput.close()
    chHadOutput.close()

    print '='*80
    print 'Systematics for ', channel, ':'
    if sigmaTTJetsSF < chHadTTJetsSF :
        print 'ttjets  SF = ', sigmaTTJetsSF, ' +(', chHadTTJetsSF, ') -(', (2.*sigmaTTJetsSF - chHadTTJetsSF), ')'
        print 'ttgamma SF = ', sigmaTTGammaSF, ' +(', chHadTTGammaSF, ') -(', (2.*sigmaTTGammaSF - chHadTTGammaSF), ')'
    else :
        print 'ttjets  SF = ', sigmaTTJetsSF, ' +(', (2.*sigmaTTJetsSF - chHadTTJetsSF), ') -(', chHadTTJetsSF, ')'
        print 'ttgamma SF = ', sigmaTTGammaSF, ' +(', (2.*sigmaTTGammaSF - chHadTTGammaSF), ') -(', chHadTTGammaSF, ')'
    print '='*80
    print '\n\n'
