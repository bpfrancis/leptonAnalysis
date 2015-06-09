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

def goPhoton(path):

    systMC = 0.0
    systData = 0.0
    
    mc = []
    data = []

    f = open(path, 'r')
    
    for line in f:
        if 'systematic' in line:
            continue
        if 'scale' in line or 'pdf' in line:
            continue
        name, x, xerr, y, yerr = line.split()
        mc.append(float(x))
        data.append(float(y))
        
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
            
    print path, ' -- purityMC(', systMC, ') -- purityData(', systData, ')'

    f.close()

go('../scaleFactors/ttbarSF_M3_ele_bjj.txt')
go('../scaleFactors/ttbarSF_M3_muon_bjj.txt')
go('../scaleFactors/wjetsSF_ele_bjj.txt')
go('../scaleFactors/wjetsSF_muon_bjj.txt')

goPhoton('../fit/photonPurity_chHadIso_ele_bjj.txt')
goPhoton('../fit/photonPurity_chHadIso_muon_bjj.txt')
goPhoton('../fit/photonPurity_sigma_ele_bjj.txt')
goPhoton('../fit/photonPurity_sigma_muon_bjj.txt')
