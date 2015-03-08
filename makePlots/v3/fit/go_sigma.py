from fit import *

controlRegion = 'SigmaPlot'

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

for channel in channels:
    output = open('sigmaFitResults_'+channel+'.txt', 'w')
    output.write('systematic\ttopSF\ttopSFerror\tttgammaSF\tttgammaSFerror\n')

    for systematic in systematics:
        doSigmaFit(channel, controlRegion, systematic, output, 0.005, 0.02)

    output.close()

