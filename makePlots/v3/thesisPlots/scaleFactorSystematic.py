import math

def go(path):

    syst = 0.0
    values = []

    f = open(path, 'r')
    
    for line in f:
        if 'systematic' in line:
            continue
        if 'scale' in line or 'pdf' in line:
            continue
        name, val, err = line.split()
        values.append(float(val))
        
    for i in range(1, len(values), 2):
        up = abs(values[i] - values[0])
        down = abs(values[i+1] - values[0])
        avg = (up+down)/2
            
        syst = math.sqrt(syst*syst + avg*avg)
            
    print path, ' -- ', syst

    f.close()

def purityError(path):

    systMC = 0.0
    systData = 0.0
    
    mc = []
    mcStat = []
    data = []
    dataStat = []

    f = open(path, 'r')
    
    for line in f:
        if 'systematic' in line:
            continue
        name, x, xerr, y, yerr = line.split()
        mc.append(float(x))
        mcStat.append(float(xerr))
        data.append(float(y))
        dataStat.append(float(yerr))
        
    for i in range(1, len(mc), 2):
        up = abs(mc[i] - mc[0])
        down = abs(mc[i+1] - mc[0])
        avg = (up+down)/2
            
        systMC = math.sqrt(systMC*systMC + avg*avg)

    for i in range(1, len(data), 2):
        up = abs(data[i] - data[0])
        down = abs(data[i+1] - data[0])
        avg = (up+down)/2
            
        systData = math.sqrt(systData*systData + avg*avg)
            
    channel = 'e'
    if 'muon' in path:
        channel = '$\\mu$'

    plusMinus = ' $\\pm$ '
    columnBreak = ' & '

    x = round(mc[0], 3)
    y = round(mcStat[0], 3)
    z = round(systMC, 3)
    firstBit = str(x) + plusMinus + str(y) + plusMinus + str(z)

    x = round(data[0], 3)
    y = round(dataStat[0], 3)
    z = round(systData, 3)
    secondBit = str(x) + plusMinus + str(y) + plusMinus + str(z)

    print channel, columnBreak, firstBit, columnBreak, secondBit, ' \\\\'

    f.close()

def purityTable(purityPath, photonPath, jetPath):

    value = [0.0, 0.0, 0.0, 0.0]
    statistic = [0.0, 0.0, 0.0, 0.0]
    systematic = [0.0, 0.0, 0.0, 0.0]

    mcPurities = []
    dataPurities = []

    photonSFs = []
    jetSFs = []
    
    purityFile = open(purityPath, 'r')
    photonFile = open(photonPath, 'r')
    jetFile = open(jetPath, 'r')

    for line in purityFile:
        if 'systematic' in line:
            continue
        if 'scale' in line or 'pdf' in line:
            continue
        
        name, x, xerr, y, yerr = line.split()

        if value[0] == 0.0:
            value[0] = float(x)
            statistic[0] = float(xerr)
            value[1] = float(y)
            statistic[1] = float(yerr)

        mcPurities.append(float(x))
        dataPurities.append(float(y))
        
    for i in range(1, len(mcPurities), 2):
        up = abs(mcPurities[i] - mcPurities[0])
        down = abs(mcPurities[i+1] - mcPurities[0])
        avg = (up+down)/2
            
        systematic[0] = math.sqrt(systematic[0]*systematic[0] + avg*avg)

    for i in range(1, len(dataPurities), 2):
        up = abs(dataPurities[i] - dataPurities[0])
        down = abs(dataPurities[i+1] - dataPurities[0])
        avg = (up+down)/2
            
        systematic[1] = math.sqrt(systematic[1]*systematic[1] + avg*avg)

    for line in photonFile:
        if 'systematic' in line:
            continue
        if 'scale' in line or 'pdf' in line:
            continue
        
        name, x, xerr = line.split()

        if value[2] == 0.0:
            value[2] = float(x)
            statistic[2] = float(xerr)

        photonSFs.append(float(x))
        
    for i in range(1, len(photonSFs), 2):
        up = abs(photonSFs[i] - photonSFs[0])
        down = abs(photonSFs[i+1] - photonSFs[0])
        avg = (up+down)/2
            
        systematic[2] = math.sqrt(systematic[2]*systematic[2] + avg*avg)

    for line in jetFile:
        if 'systematic' in line:
            continue
        if 'scale' in line or 'pdf' in line:
            continue
        
        name, x, xerr = line.split()

        if value[3] == 0.0:
            value[3] = float(x)
            statistic[3] = float(xerr)

        jetSFs.append(float(x))
        
    for i in range(1, len(jetSFs), 2):
        up = abs(jetSFs[i] - jetSFs[0])
        down = abs(jetSFs[i+1] - jetSFs[0])
        avg = (up+down)/2
            
        systematic[3] = math.sqrt(systematic[3]*systematic[3] + avg*avg)

    channel = 'e'
    if 'muon' in purityPath:
        channel = '$\\mu$'

    plusMinus = ' $\\pm$ '
    columnBreak = ' & '

    parts = []
    for i in range(0, 4):
        precision = 3 if (i < 2) else 2
        x = round(value[i], precision)
        y = round(statistic[i], precision)
        z = round(systematic[i], precision)
        parts.append(str(x) + plusMinus + str(y) + plusMinus + str(z))

    print channel, columnBreak, parts[0], columnBreak, parts[1], columnBreak, parts[2], columnBreak, parts[3], ' \\\\'


    purityFile.close()
    photonFile.close()
    jetFile.close()

go('../scaleFactors/ttbarSF_M3_ele_bjj.txt')
go('../scaleFactors/ttbarSF_M3_muon_bjj.txt')
go('../scaleFactors/wjetsSF_ele_bjj.txt')
go('../scaleFactors/wjetsSF_muon_bjj.txt')

purityTable('../fit/photonPurity_chHadIso_ele_bjj.txt', '../fit/photonSigSF_chHadIso_ele_bjj.txt', '../fit/jetBkgSF_chHadIso_ele_bjj.txt')
purityTable('../fit/photonPurity_chHadIso_muon_bjj.txt', '../fit/photonSigSF_chHadIso_muon_bjj.txt', '../fit/jetBkgSF_chHadIso_muon_bjj.txt')

purityTable('../fit/photonPurity_chHadIso_metCut_50_ele_bjj.txt', '../fit/photonSigSF_chHadIso_metCut_50_ele_bjj.txt', '../fit/jetBkgSF_chHadIso_metCut_50_ele_bjj.txt')
purityTable('../fit/photonPurity_chHadIso_metCut_50_muon_bjj.txt', '../fit/photonSigSF_chHadIso_metCut_50_muon_bjj.txt', '../fit/jetBkgSF_chHadIso_metCut_50_muon_bjj.txt')

#purityTable('../fit/photonPurity_sigma_ele_bjj.txt', '../fit/photonSigSF_sigma_ele_bjj.txt', '../fit/jetBkgSF_sigma_ele_bjj.txt')
#purityTable('../fit/photonPurity_sigma_muon_bjj.txt', '../fit/photonSigSF_sigma_muon_bjj.txt', '../fit/jetBkgSF_sigma_muon_bjj.txt')

