from utils import *

def doFit(channel, controlRegion, systematic, output, xlo, xhi):

    input = '../histograms_'+channel+'_'+controlRegion+'.root'

    varName = 'z_mass'

    dataHist = get1DHist(input, varName+'_gg_'+channel)

    zHist = get1DHist(input, varName+'_dy1JetsToLL_'+channel+systematic)
    zHist.Add(get1DHist(input, varName+'_dy2JetsToLL_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_dy3JetsToLL_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_dy4JetsToLL_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_ZGToLLG_'+channel+systematic))
    zHist.Add(get1DHist(input, varName+'_WGToLNuG_'+channel+systematic))

    bkgHist = get1DHist(input, varName+'_TBar_s_'+channel+systematic)
    bkgHist.Add(get1DHist(input, varName+'_TBar_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TBar_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_s_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_t_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_T_tW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_WW_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_WZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ZZ_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTWJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTZJets_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsHadronic_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsFullLep_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_ttJetsSemiLep_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_TTGamma_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_W3JetsToLNu_'+channel+systematic))
    bkgHist.Add(get1DHist(input, varName+'_W4JetsToLNu_'+channel+systematic))

    (fitFrac, fitFracErr) = makeFit(varName, xlo, xhi, zHist, bkgHist, dataHist)

    lowbin = dataHist.FindBin(xlo)
    highbin = dataHist.FindBin(xhi) - 1
    
    dataInt = dataHist.Integral(lowbin, highbin)
    zInt = zHist.Integral(lowbin, highbin)
    bkgInt = bkgHist.Integral(lowbin, highbin)

    zSF = fitFrac * dataInt / zInt
    zSFerror = fitFracErr * dataInt / zInt

    bkgSF = (1.0-fitFrac) * dataInt / bkgInt
    bkgSFerror = fitFracErr * dataInt / bkgInt

    output_z.write(systematic+'\t'+
                   str(zSF)+'\t'+
                   str(zSFerror)+'\n')

    drawPlots(dataHist, zHist, zSF, 'ttbar', bkgHist, bkgSF, 'bkg', xlo, xhi, varName+'_'+channel+systematic)

    return (zSF, zSFerror)
