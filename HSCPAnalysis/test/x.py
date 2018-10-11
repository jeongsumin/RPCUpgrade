#!/usr/bin/env python
from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

import os
def buildLegend(x1, y1, x2, y2):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetLineWidth(0)
    leg.SetFillStyle(0)
    return leg

#fS = TFile("hist_HSCPppstau_m1599_LGW25.root")
#fS = TFile("test_hist_HSCPppstau_m1599_LGW25.root")
#fB = TFile("test_hist_DYJetsToLL_M-50_noPU.root")
#fS = TFile("hist_fitslope_HSCPppstau_m1599_LGW25.root")
#fB = TFile("hist_DYJetsToLL_M-50_PU200.root")
#fS = TFile("HSCP_1599_noPU_test.root")
#fS = TFile("hist_HSCPppstau_M_1599_14TeV_PU200.root")
#fS = TFile("RangeTest_hist_DYJetsToLL_M-50_noPU.root")
#fS_PU200 = TFile("hist_fitslope_HSCPppstau_M_1599_14TeV_PU200.root")
#fS_PU200 = TFile("hist_HSCPppstau_M_1599_14TeV_PU200.root")
#fS_PU200 = TFile("test_hist_HSCPppstau_M_1599_14TeV_PU200.root")
#fS_PU200 = TFile("RangeTest_hist_HSCPppstau_M_1599_14TeV_PU200.root")
#outfile = TFile('Effplot.root', 'RECREATE')
#tB = fB.Get("tree")

fS = TFile("test_hist_DYJetsToLL_M-50_noPU.root")
fS_PU200 = TFile("test_hist_DYJetsToLL_M-50_PU200.root")
#fS = TFile("hist_fitslope_DYJetsToLL_M-50_noPU.root")
#fS_PU200 = TFile("hist_fitslope_DYJetsToLL_M-50_PU200.root")

tS = fS.Get("tree")
tS_PU200 = fS_PU200.Get("tree")

c_beta_Eff = TCanvas("c_beta_Eff", "c_beta_Eff", 1000, 1000)
c_beta_Eff.SetGridx()
c_beta_Eff.SetGridy()

#hS_beta_genTot = TH1D("hS_beta_genTot", "hS_beta_genTot;Total Generated #beta", 100, -1, 2)
#hS_beta_gen = TH1D("hS_beta_gen", "hS_beta_gen;Generated #beta", 100, -1, 2)
hS_beta_genTot = TH1D("hS_beta_genTot", "hS_beta_genTot;Total Generated #beta", 25, 0., 1.)
hS_beta_gen = TH1D("hS_beta_gen", "hS_beta_gen;Generated #beta", 25, 0., 1.)
hS_PU200_beta_genTot = TH1D("hS_PU200_beta_genTot", "hS_PU200_beta_genTot;Total Generated #beta", 25, 0., 1.)
hS_PU200_beta_gen = TH1D("hS_PU200_beta_gen", "hS_PU200_beta_gen;Generated #beta", 25, 0., 1.)
hS_noPUhit = TH1D("hS_noPUhit", "", 10, 0., 10.)
hS_PU200hit = TH1D("hS_PU200hit", "", 10, 0., 10.)

#tS.Draw("genBeta>>hS_beta_gen", "TMath::Abs(genEta) < 2.4 && TMath::Abs(rpcBetaErr/rpcBeta) < 100 && rpcBeta < 0.7 && rpcHits_n >= 3", "goff")
#tS.Draw("genBeta>>hS_beta_gen", "TMath::Abs(genEta) < 2.4 && rpcBetaErr < 10e9 && rpcBeta < 0.7 && rpcHits_n >= 3 && rpcBetaErr/rpcBeta < 0.3", "goff")
#tS.Draw("genBeta>>hS_beta_gen", "TMath::Abs(genEta) < 1.6 && rpcBetaErr < 1e9  && rpcHits_n >= 3", "goff")
#####high eta region########
'''
tS.Draw("gen1_p4.Beta()>>hS_beta_gen", "TMath::Abs(gen1_p4.Eta()) < 1.8  && fit_beta1 > 0 && fit_beta1 < 0.7 && fit_nhit1 >= 2", "goff")
tS.Draw("gen2_p4.Beta()>>+hS_beta_gen", "TMath::Abs(gen2_p4.Eta()) < 1.8 && fit_beta2 > 0 && fit_beta2 < 0.7 && fit_nhit2 >= 2", "goff")
tS_PU200.Draw("gen1_p4.Beta()>>hS_PU200_beta_gen", "TMath::Abs(gen1_p4.Eta()) < 1.8 && fit_beta1 > 0 && fit_beta1 < 0.7 && fit_nhit1 >= 2", "goff")
tS_PU200.Draw("gen2_p4.Beta()>>+hS_PU200_beta_gen", "TMath::Abs(gen2_p4.Eta()) < 1.8 && fit_beta2 > 0 && fit_beta2 < 0.7 && fit_nhit2 >= 2", "goff")
'''

######bxconstrained#########
tS.Draw("gen1_p4.Beta()>>hS_beta_gen", "TMath::Abs(gen1_p4.Eta()) >= 1.8 && TMath::Abs(gen1_p4.Eta()) < 2.4  && fit_beta1 > 0 && fit_beta1 < 0.7 && fit_nhit1 >= 2", "goff")
tS.Draw("gen2_p4.Beta()>>+hS_beta_gen", "TMath::Abs(gen2_p4.Eta()) >= 1.8 && TMath::Abs(gen2_p4.Eta()) < 2.4  && fit_beta2 > 0 && fit_beta2 < 0.7 && fit_nhit2 >= 2", "goff")
tS_PU200.Draw("gen1_p4.Beta()>>hS_PU200_beta_gen", "TMath::Abs(gen1_p4.Eta()) >= 1.8 && TMath::Abs(gen1_p4.Eta()) < 2.4  && fit_beta1 > 0 && fit_beta1 < 0.7 && fit_nhit1 >= 2", "goff")
tS_PU200.Draw("gen2_p4.Beta()>>+hS_PU200_beta_gen", "TMath::Abs(gen2_p4.Eta()) >= 1.8 && TMath::Abs(gen2_p4.Eta()) < 2.4 && fit_beta2 > 0 && fit_beta2 < 0.7 && fit_nhit2 >= 2", "goff")

'''
#######fitslope##########
tS.Draw("gen1_p4.Beta()>>hS_beta_gen", "TMath::Abs(gen1_p4.Eta()) < 1.8 && fit_qual1/fit_beta1 < 0.3 && fit_beta1 > 0 && fit_beta1 < 0.7 && fit_nhit1 >= 3", "goff")
tS.Draw("gen2_p4.Beta()>>+hS_beta_gen", "TMath::Abs(gen2_p4.Eta()) < 1.8 && fit_qual2/fit_beta2 < 0.3 && fit_beta2 > 0 && fit_beta2 < 0.7 && fit_nhit2 >= 3", "goff")
tS_PU200.Draw("gen1_p4.Beta()>>hS_PU200_beta_gen", "TMath::Abs(gen1_p4.Eta()) < 1.8 && fit_qual1/fit_beta1 < 0.3 && fit_beta1 > 0 && fit_beta1 < 0.7 && fit_nhit1 >= 3", "goff")
tS_PU200.Draw("gen2_p4.Beta()>>+hS_PU200_beta_gen", "TMath::Abs(gen2_p4.Eta()) < 1.8 && fit_qual2/fit_beta2 < 0.3 && fit_beta2 > 0 && fit_beta2 < 0.7 && fit_nhit2 >= 3", "goff")
'''

#tS.Draw("gen1_p4.Beta()>>hS_beta_gen", "fit_nhit1 >= 3", "goff")
#tS.Draw("gen2_p4.Beta()>>+hS_beta_gen", "fit_nhit2 >= 3", "goff")
#tS.Draw("genBeta>>hS_beta_gen", " TMath::Abs(genEta) < 2.4")
#tB.Draw("gen1_p4.Beta()>>hB_beta_gen", "gen1_pdgId!=0 && fit_qual1 < 1e9 && fit_nhit1 >= 3 && fit_qual1/fit_beta1 < 0.3 && fit_beta1 > 0 && fit_beta1 < 0.7", "goff"

#tS.Draw("gen1_p4.Beta()>>hS_beta_genTot", "TMath::Abs(gen1_p4.Eta()) < 1.6", "goff")
#tS.Draw("gen2_p4.Beta()>>+hS_beta_genTot", "TMath::Abs(gen2_p4.Eta()) < 1.6", "goff")

tS.Draw("gen1_p4.Beta()>>hS_beta_genTot", "TMath::Abs(gen1_p4.Eta()) >= 1.8 && TMath::Abs(gen1_p4.Eta()) < 2.4", "goff")
tS.Draw("gen2_p4.Beta()>>+hS_beta_genTot", "TMath::Abs(gen2_p4.Eta()) >= 1.8 && TMath::Abs(gen2_p4.Eta()) < 2.4", "goff")
tS_PU200.Draw("gen1_p4.Beta()>>hS_PU200_beta_genTot", "TMath::Abs(gen1_p4.Eta()) >= 1.8 && TMath::Abs(gen1_p4.Eta()) < 2.4", "goff")
tS_PU200.Draw("gen2_p4.Beta()>>+hS_PU200_beta_genTot", "TMath::Abs(gen2_p4.Eta()) >= 1.8 && TMath::Abs(gen2_p4.Eta()) < 2.4", "goff")
'''
tS.Draw("gen1_p4.Beta()>>hS_beta_genTot", "TMath::Abs(gen1_p4.Eta()) < 1.8", "goff")
tS.Draw("gen2_p4.Beta()>>+hS_beta_genTot", "TMath::Abs(gen2_p4.Eta()) < 1.8", "goff")
tS_PU200.Draw("gen1_p4.Beta()>>hS_PU200_beta_genTot", "TMath::Abs(gen1_p4.Eta()) < 1.8", "goff")
tS_PU200.Draw("gen2_p4.Beta()>>+hS_PU200_beta_genTot", "TMath::Abs(gen2_p4.Eta()) < 1.8", "goff")
'''
print "(noPU) total iRPC region efficiency = ", hS_beta_gen.Integral()/hS_beta_genTot.Integral()
print "(PU200) total iRPC region efficiency = ", hS_PU200_beta_gen.Integral()/hS_PU200_beta_genTot.Integral()

#tS.Draw("genBeta>>hS_beta_genTot")

tS.Draw("fit_nhit1>>hS_noPUhit","TMath::Abs(gen1_p4.Eta())<2.4","goff")
tS.Draw("fit_nhit2>>+hS_noPUhit","TMath::Abs(gen2_p4.Eta())<2.4","goff")
tS_PU200.Draw("fit_nhit1>>hS_PU200hit","TMath::Abs(gen1_p4.Eta())<2.4","goff")
tS_PU200.Draw("fit_nhit2>>hS_PU200hit","TMath::Abs(gen2_p4.Eta())<2.4","goff")

Sig_betaEffCurve = TEfficiency(hS_beta_gen, hS_beta_genTot)
Sig_PU200_betaEffCurve = TEfficiency(hS_PU200_beta_gen, hS_PU200_beta_genTot)
Sig_betaEffCurve.SetTitle(";generated #beta;Efficiency")
Sig_PU200_betaEffCurve.SetMarkerColor(kRed+2)
Sig_PU200_betaEffCurve.SetLineColor(kRed+2)
Sig_betaEffCurve.Draw()
Sig_PU200_betaEffCurve.Draw("SAME")
gPad.Update()
Sig_betaEffCurve.GetPaintedGraph().GetXaxis().SetRangeUser(0.,1.)
#Sig_betaEffCurve.GetPaintedGraph().SetMinimum(0.)
Sig_betaEffCurve.GetPaintedGraph().SetMinimum(0.)
Sig_betaEffCurve.GetPaintedGraph().SetMaximum(1.2)
gPad.Update()

leg = TLegend(0.45, 0.78 , 0.85, 0.88)
leg.SetTextSize(0.03)
leg.AddEntry(Sig_betaEffCurve,"RPC-HSCP trigger (PU=0)", "lp")
leg.AddEntry(Sig_PU200_betaEffCurve,"RPC-HSCP trigger (PU=200)", "lp")
leg.SetBorderSize(0)
leg.Draw()

text = TText(0.35, 0.92, "Work in progress");
#TText = TText(0.83, 0.9, "Preliminary");
text.SetTextFont(52)
text.SetTextAngle(0)
text.SetTextColor(kBlack)
text.SetTextSize(0.04)
text.SetTextAlign(22)
text.SetNDC(kTRUE) #IMPORTANT! visible Text;True or not;False
text.Draw()

#c_beta_Eff.Print("noPU_beta_Eff.png")
#outfile.Write()
#outfile.Close()
'''
cc = TCanvas()
hS_beta_genTot.Draw("")
hS_beta_gen.Draw("sametext")
hS_PU200_beta_genTot.SetLineColor(kBlue)
hS_PU200_beta_genTot.Draw("same")
hS_PU200_beta_gen.Draw("sametext")
'''
noPUhit = hS_noPUhit.GetMean()
PU200hit = hS_PU200hit.GetMean()
print "noPU num. of hit = ", noPUhit
print "PU200 num. of hit = ", PU200hit

####### 2d histogram beta VS |eta| #######
c2D = TCanvas("betaVS|eta|","betaVS|eta|", 700, 700)
hS_noPU_genBetaVSgenEta_Tot = TH2D("hS_noPU_genBetaVSgenEta_Tot", "", 25, 0., 1., 75, 0., 2.5)
hS_noPU_genBetaVSgenEta = TH2D("hS_noPU_genBetaVSgenEta", "", 25, 0., 1., 75, 0., 2.5)

tS.Draw("gen1_p4.Eta():gen1_p4.Beta()>>hS_noPU_genBetaVSgenEta_Tot", "TMath::Abs(gen1_p4.Eta() < 2.4)")
tS.Draw("gen2_p4.Eta():gen2_p4.Beta()>>+hS_noPU_genBetaVSgenEta_Tot", "TMath::Abs(gen2_p4.Eta() < 2.4)")
tS.Draw("gen1_p4.Eta():gen1_p4.Beta()>>hS_noPU_genBetaVSgenEta", "TMath::Abs(gen1_p4.Eta()) < 2.4 && fit_qual1 <= 0.01 && fit_beta1 > 0 && fit_beta1 < 0.7 && fit_nhit1 >= 1")
tS.Draw("gen2_p4.Eta():gen2_p4.Beta()>>+hS_noPU_genBetaVSgenEta", "TMath::Abs(gen2_p4.Eta()) < 2.4 && fit_qual2 <= 0.01 && fit_beta2 > 0 && fit_beta2 < 0.7 && fit_nhit2 >= 1")

hS_noPU_genBetaVSgenEta_Tot.Sumw2()
hS_noPU_genBetaVSgenEta.Sumw2()
hS_noPU_genBetaVSgenEta.Divide(hS_noPU_genBetaVSgenEta_Tot)
hS_noPU_genBetaVSgenEta.GetXaxis().SetTitle("#beta")
hS_noPU_genBetaVSgenEta.GetYaxis().SetTitle("|#eta|")
hS_noPU_genBetaVSgenEta.GetYaxis().SetTitleOffset(1)
hS_noPU_genBetaVSgenEta.GetZaxis().SetTitle("Efficiency")
hS_noPU_genBetaVSgenEta.GetZaxis().SetTitleOffset(0.97)
hS_noPU_genBetaVSgenEta.SetMaximum(1.0)
gStyle.SetPalette(kDeepSea)
gPad.SetRightMargin(0.1)
hS_noPU_genBetaVSgenEta.Draw("colz")


####### noPU, PU200 fit_qual1 & 2 #######
c2 = TCanvas()
hS_noPU_qual = TH1D("hS_noPU_qual", "quality(dRerr)", 100, -0.001, 0.015)
hS_PU200_qual = TH1D("hS_PU200_qual", "", 100, -0.001, 0.015)
hS_noPU_qual.SetLineColor(kBlack)
hS_PU200_qual.SetLineColor(kRed)
tS.Draw("fit_qual1>>hS_noPU_qual","TMath::Abs(gen1_p4.Eta()) >= 1.8 && TMath::Abs(gen1_p4.Eta()) < 2.4","")
tS.Draw("fit_qual2>>+hS_noPU_qual","TMath::Abs(gen2_p4.Eta()) >= 1.8 && TMath::Abs(gen2_p4.Eta()) < 2.4","")
tS_PU200.Draw("fit_qual1>>hS_PU200_qual","TMath::Abs(gen1_p4.Eta()) >= 1.8 && TMath::Abs(gen1_p4.Eta()) < 2.4","same")
tS_PU200.Draw("fit_qual2>>+hS_PU200_qual","TMath::Abs(gen2_p4.Eta()) >= 1.8 && TMath::Abs(gen2_p4.Eta()) < 2.4","same")
leg2 = TLegend(0.4, 0.78 , 0.9, 0.88)
leg2.SetTextSize(0.03)
leg2.AddEntry(hS_noPU_qual,"PU=0", "l")
leg2.AddEntry(hS_PU200_qual,"PU=200", "l")
leg2.SetBorderSize(0)
leg2.Draw()
'''
####### dphi, deta (hscp generated particle, simDigi(already only hscp)) ########
c3 = TCanvas()
hS_PU200_genPhi = TH1D("hS_PU200_genPhi", "gen particle phi", 100, -3.14, 3.14)
hS_PU200_simDigiPhi = TH1D("hS_PU200_simDigiPhi", "simDigi particle phi", 100, -3.14, 3.14)
tS_PU200.Draw("gen1_p4.Phi()>>hS_PU200_genPhi")
tS_PU200.Draw("gen2_p4.Phi()>>+hS_PU200_genPhi")

c4 = TCanvas()
hS_PU200_genEta = TH1D("hS_PU200_genEta", "gen particle eta", 100, -3, 3)
hS_PU200_simDigiPhi = TH1D("hS_PU200_simDigiEta", "simDigi particle eta", 100, -3, 3)
tS_PU200.Draw("gen1_p4.Eta()>>hS_PU200_genEta")
tS_PU200.Draw("gen2_p4.Eta()>>+hS_PU200_genEta")
'''

