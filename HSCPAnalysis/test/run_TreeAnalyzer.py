#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

import os
objs = []
#mass = 651
mass = 1599
samples = [
    "/data/users/jhgoh/HSCP/20170804_1/DYJetsToLL_M-50_noPU.root",
    #"/data/users/jhgoh/HSCP/20170627_1/HSCPppstau_m%d_LGW25.root" % mass,
    #"/data/users/sumin/HSCPAnalysis/hist_HSCP_m1599_LGW25_input80th.root",
    #"/data/users/jhgoh/HSCP/20170804_1/DYJetsToLL_M-50_PU200*.root",
    #"/data/users/jhgoh/HSCP/20171204_1/HSCPppstau_M_%d_PU200*.root" % mass,
    #"/data/users/jhgoh/HSCP/20171204_1/HSCPppstau_M_1599_14TeV_PU200.root",
]

gSystem.CompileMacro("TreeAnalyzer.C", "k");

for sample in samples:
    foutName = "RangeTest_hist_" + os.path.basename(sample.replace('*',''))
    #foutName = "hist_" + os.path.basename(sample.replace('*',''))

    tree = TChain("HSCPTree/tree")
    tree.Add(sample)
    ana = TreeAnalyzer(tree)
    fout = TFile(foutName, "RECREATE")
    ana.Loop(fout)
    ana = None

