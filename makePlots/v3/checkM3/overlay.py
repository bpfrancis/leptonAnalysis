import ROOT
from ROOT import gROOT
from ROOT import gStyle

import math

openfiles = {}

def get1DHist(filename, histname):
    if filename not in openfiles:
        openfiles[filename] = ROOT.TFile(filename,'READ')
    file = openfiles[filename]

    hist = ROOT.TH1D()
    hist = file.Get(histname)
    hist.SetDirectory(0)
    return hist

def compareQCDNormalization(channel):

    inputCR0 = '../histograms_'+channel+'_CR0.root'
    inputCR1 = '../histograms_'+channel+'_CR1.root'
    inputSR1 = '../histograms_'+channel+'_SR1.root'
    inputSR2 = '../histograms_'+channel+'_SR2.root'
    inputCR2 = '../histograms_'+channel+'_CR2.root'

    controlRegions = ['CR0', 'CR1', 'SR1', 'SR2', 'CR2', 'CR2a']

    fitSF = 0.329381930897
    fitSFerror = 0.00388688459937
    if 'muon' in channel:
        fitSF = 0.0145930197407
        fitSFerror = 0.000532915435935

    print 'Fit results in CR0: ', fitSF, ' +- ', fitSFerror

    for region in controlRegions:
        filename = '../histograms_'+channel+'_'+region+'.root'

        data = get1DHist(filename, 'pfMET_gg_'+channel)
        qcd = get1DHist(filename, 'pfMET_qcd_'+channel)

        mc = get1DHist(filename, 'pfMET_ttJetsHadronic_'+channel)
        mc.Add(get1DHist(filename, 'pfMET_ttJetsSemiLep_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_ttJetsFullLep_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_W1JetsToLNu_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_W2JetsToLNu_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_W3JetsToLNu_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_W4JetsToLNu_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_TBar_s_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_TBar_t_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_TBar_tW_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_T_s_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_T_t_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_T_tW_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_dy1JetsToLL_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_dy2JetsToLL_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_dy3JetsToLL_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_dy4JetsToLL_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_ttA_2to5_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_WW_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_WZ_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_ZZ_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_TTWJets_'+channel))
        mc.Add(get1DHist(filename, 'pfMET_TTZJets_'+channel))

        lowbin = 1
        highbin = data.GetXaxis().FindBin(20.0) - 1

        ndata = data.Integral(lowbin, highbin)
        nqcd = qcd.Integral(lowbin, highbin)
        nmc = mc.Integral(lowbin, highbin)

        sigmadata = 0
        sigmaqcd = 0
        sigmamc = 0

        for ibin in range(lowbin, highbin):
            sigmadata += data.GetBinError(ibin) * data.GetBinError(ibin)
            sigmaqcd += qcd.GetBinError(ibin) * qcd.GetBinError(ibin)
            sigmamc += mc.GetBinError(ibin) * mc.GetBinError(ibin)

        sigmadata = math.sqrt(sigmadata)
        sigmaqcd = math.sqrt(sigmaqcd)
        sigmamc = math.sqrt(sigmamc)

        if nqcd == 0:
            print 'No QCD to normalize in ', region
            continue

        scale = (ndata - nmc) / nqcd
        scaleError = sigmadata*sigmadata + sigmamc*sigmamc
        scaleError /= (ndata - nmc) * (ndata - nmc)
        scaleError += sigmaqcd*sigmaqcd / nqcd / nqcd
        scaleError = scale * math.sqrt(scaleError)

        print 'Normalizing in ', region, ': ', scale, ' +- ', scaleError

###########################################################
def overlayWJets(channel):

    inputCR0 = '../histograms_'+channel+'_CR0.root'
    inputCR1 = '../histograms_'+channel+'_CR1.root'
    inputSR1 = '../histograms_'+channel+'_SR1.root'
    inputSR2 = '../histograms_'+channel+'_SR2.root'
    inputCR2 = '../histograms_'+channel+'_CR2.root'

    hCR0 = get1DHist(inputCR0, 'm3_W1JetsToLNu_'+channel)
    hCR0.Add(get1DHist(inputCR0, 'm3_W2JetsToLNu_'+channel))
    hCR0.Add(get1DHist(inputCR0, 'm3_W3JetsToLNu_'+channel))
    hCR0.Add(get1DHist(inputCR0, 'm3_W4JetsToLNu_'+channel))
    hCR0.Rebin(4)
    hCR0.Scale(1.0 / hCR0.Integral())
    hCR0.SetLineColor(ROOT.kBlue)
    hCR0.SetLineWidth(3)

    hCR1 = get1DHist(inputCR1, 'm3_W1JetsToLNu_'+channel)
    hCR1.Add(get1DHist(inputCR1, 'm3_W2JetsToLNu_'+channel))
    hCR1.Add(get1DHist(inputCR1, 'm3_W3JetsToLNu_'+channel))
    hCR1.Add(get1DHist(inputCR1, 'm3_W4JetsToLNu_'+channel))
    hCR1.Rebin(4)
    hCR1.Scale(1.0 / hCR1.Integral())
    hCR1.SetLineColor(ROOT.kRed)
    hCR1.SetLineWidth(3)

    hSR1 = get1DHist(inputSR1, 'm3_W1JetsToLNu_'+channel)
    hSR1.Add(get1DHist(inputSR1, 'm3_W2JetsToLNu_'+channel))
    hSR1.Add(get1DHist(inputSR1, 'm3_W3JetsToLNu_'+channel))
    hSR1.Add(get1DHist(inputSR1, 'm3_W4JetsToLNu_'+channel))
    hSR1.Rebin(4)
    hSR1.Scale(1.0 / hSR1.Integral())
    hSR1.SetLineColor(8)
    hSR1.SetLineWidth(3)

    hSR2 = get1DHist(inputSR2, 'm3_W1JetsToLNu_'+channel)
    hSR2.Add(get1DHist(inputSR2, 'm3_W2JetsToLNu_'+channel))
    hSR2.Add(get1DHist(inputSR2, 'm3_W3JetsToLNu_'+channel))
    hSR2.Add(get1DHist(inputSR2, 'm3_W4JetsToLNu_'+channel))
    hSR2.Rebin(4)
    if hSR2.Integral() > 0:
        hSR2.Scale(1.0 / hSR2.Integral())
    hSR2.SetLineWidth(3)

    hCR2 = get1DHist(inputCR2, 'm3_W1JetsToLNu_'+channel)
    hCR2.Add(get1DHist(inputCR2, 'm3_W2JetsToLNu_'+channel))
    hCR2.Add(get1DHist(inputCR2, 'm3_W3JetsToLNu_'+channel))
    hCR2.Add(get1DHist(inputCR2, 'm3_W4JetsToLNu_'+channel))
    hCR2.Rebin(4)
    if hCR2.Integral() > 0:
        hCR2.Scale(1.0 / hCR2.Integral())
    hCR2.SetLineColor(ROOT.kMagenta)
    hCR2.SetLineWidth(3)

    leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.95, '', 'brNDC')
    leg.AddEntry(hCR0, 'CR0', 'L')
    leg.AddEntry(hCR1, 'CR1', 'L')
    leg.AddEntry(hSR1, 'SR1', 'L')
    leg.AddEntry(hSR2, 'SR2', 'L')
    leg.AddEntry(hCR2, 'CR2', 'L')

    can = ROOT.TCanvas()
    #can.SetLogy(True)

    hCR1.GetYaxis().SetRangeUser(0, 0.3)
    hCR1.Draw('e1')
    hCR0.Draw('hist same')
    hSR1.Draw('hist same')
    hSR2.Draw('hist same')
    hCR2.Draw('hist same')
    leg.Draw()
    can.SaveAs('m3_wjets_'+channel+'.png')

def overlayTTJets(channel):

    inputCR0 = '../histograms_'+channel+'_CR0.root'
    inputCR1 = '../histograms_'+channel+'_CR1.root'
    inputSR1 = '../histograms_'+channel+'_SR1.root'
    inputSR2 = '../histograms_'+channel+'_SR2.root'
    inputCR2 = '../histograms_'+channel+'_CR2.root'

    hCR0 = get1DHist(inputCR0, 'm3_ttJetsSemiLep_'+channel)
    hCR0.Add(get1DHist(inputCR0, 'm3_ttJetsFullLep_'+channel))
    hCR0.Add(get1DHist(inputCR0, 'm3_ttJetsHadronic_'+channel))
    hCR0.Add(get1DHist(inputCR0, 'm3_TTWJets_'+channel))
    hCR0.Add(get1DHist(inputCR0, 'm3_TTZJets_'+channel))
    hCR0.Add(get1DHist(inputCR0, 'm3_ttA_2to5_'+channel))
    hCR0.Rebin(4)
    hCR0.Scale(1.0 / hCR0.Integral())
    hCR0.SetLineColor(ROOT.kBlue)
    hCR0.SetLineWidth(3)

    hCR1 = get1DHist(inputCR1, 'm3_ttJetsSemiLep_'+channel)
    hCR1.Add(get1DHist(inputCR1, 'm3_ttJetsFullLep_'+channel))
    hCR1.Add(get1DHist(inputCR1, 'm3_ttJetsHadronic_'+channel))
    hCR1.Add(get1DHist(inputCR1, 'm3_TTWJets_'+channel))
    hCR1.Add(get1DHist(inputCR1, 'm3_TTZJets_'+channel))
    hCR1.Add(get1DHist(inputCR1, 'm3_ttA_2to5_'+channel))
    hCR1.Rebin(4)
    hCR1.Scale(1.0 / hCR1.Integral())
    hCR1.SetLineColor(ROOT.kRed)
    hCR1.SetLineWidth(3)

    hSR1 = get1DHist(inputSR1, 'm3_ttJetsSemiLep_'+channel)
    hSR1.Add(get1DHist(inputSR1, 'm3_ttJetsFullLep_'+channel))
    hSR1.Add(get1DHist(inputSR1, 'm3_ttJetsHadronic_'+channel))
    hSR1.Add(get1DHist(inputSR1, 'm3_TTWJets_'+channel))
    hSR1.Add(get1DHist(inputSR1, 'm3_TTZJets_'+channel))
    hSR1.Add(get1DHist(inputSR1, 'm3_ttA_2to5_'+channel))
    hSR1.Rebin(4)
    hSR1.Scale(1.0 / hSR1.Integral())
    hSR1.SetLineColor(8)
    hSR1.SetLineWidth(3)

    hSR2 = get1DHist(inputSR2, 'm3_ttJetsSemiLep_'+channel)
    hSR2.Add(get1DHist(inputSR2, 'm3_ttJetsFullLep_'+channel))
    hSR2.Add(get1DHist(inputSR2, 'm3_ttJetsHadronic_'+channel))
    hSR2.Add(get1DHist(inputSR2, 'm3_TTWJets_'+channel))
    hSR2.Add(get1DHist(inputSR2, 'm3_TTZJets_'+channel))
    hSR2.Add(get1DHist(inputSR2, 'm3_ttA_2to5_'+channel))
    hSR2.Rebin(4)
    if hSR2.Integral() > 0:
        hSR2.Scale(1.0 / hSR2.Integral())
    hSR2.SetLineWidth(3)

    hCR2 = get1DHist(inputCR2, 'm3_ttJetsSemiLep_'+channel)
    hCR2.Add(get1DHist(inputCR2, 'm3_ttJetsFullLep_'+channel))
    hCR2.Add(get1DHist(inputCR2, 'm3_ttJetsHadronic_'+channel))
    hCR2.Add(get1DHist(inputCR2, 'm3_TTWJets_'+channel))
    hCR2.Add(get1DHist(inputCR2, 'm3_TTZJets_'+channel))
    hCR2.Add(get1DHist(inputCR2, 'm3_ttA_2to5_'+channel))
    hCR2.Rebin(4)
    if hCR2.Integral() > 0:
        hCR2.Scale(1.0 / hCR2.Integral())
    hCR2.SetLineColor(ROOT.kMagenta)
    hCR2.SetLineWidth(3)

    leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.95, '', 'brNDC')
    leg.AddEntry(hCR0, 'CR0', 'L')
    leg.AddEntry(hCR1, 'CR1', 'L')
    leg.AddEntry(hSR1, 'SR1', 'L')
    leg.AddEntry(hSR2, 'SR2', 'L')
    leg.AddEntry(hCR2, 'CR2', 'L')

    can = ROOT.TCanvas()
    #can.SetLogy(True)

    hCR1.GetYaxis().SetRangeUser(0, 0.3)
    hCR1.Draw('e1')
    hCR0.Draw('hist same')
    hSR1.Draw('hist same')
    hSR2.Draw('hist same')
    hCR2.Draw('hist same')
    leg.Draw()
    can.SaveAs('m3_ttjets_'+channel+'.png')
