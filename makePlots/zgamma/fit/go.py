from fit import *

channels = ['ele_bjj', 'muon_bjj', 'ele_jjj', 'muon_jjj']
systematics = ['',
               '_btagWeightUp', '_btagWeightDown',
               '_puWeightUp', '_puWeightDown',
               '_scaleUp', '_scaleDown',
               '_pdfUp', '_pdfDown',
               '_topPtUp', '_topPtDown',
               '_JECUp', '_JECDown',
               '_leptonSFUp', '_leptonSFDown',
               '_photonSFUp', '_photonSFDown']

region = 'Any'

for channel in channels:

    results = []

    output = open('zMassSF_'+channel+'.txt', 'w')
    output.write('systematic\tSF\tError\n')

    for systematic in systematics:
        results.append(doFit(channel, region, systematic, output, 20.0, 180.0))

    output.close()
