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

go('../scaleFactors/ttbarSF_M3_ele_bjj.txt')
go('../scaleFactors/ttbarSF_M3_muon_bjj.txt')
go('../scaleFactors/wjetsSF_ele_bjj.txt')
go('../scaleFactors/wjetsSF_muon_bjj.txt')
