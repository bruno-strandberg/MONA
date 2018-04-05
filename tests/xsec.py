#!/usr/bin/python -i
from ROOT import *
gSystem.Load("$NMHDIR/common_software/libnmhsoft.so")

xsec = NuXsec()

flavors      = [0,1,2] #3 flavors
interactions = [0,1]   #nc/cc
isnubar      = [0,1]   #particle/antiparticle

graphs = []

for flav in flavors:
    for inter in interactions:
        for nubar in isnubar:

            xsec.SelectInteraction(flav,inter,nubar)
            graphs.append( TGraph() )
            name = "f_{0}_inter_{1}_nubar_{2}".format(flav, inter, nubar)
            title = "flavour {0}, is CC {1}, is nubar {2}".format(flav, inter, nubar)
            graphs[-1].SetNameTitle(name,title)

            E = 1
            while (E < 100):
                graphs[-1].SetPoint( graphs[-1].GetN(), E, xsec.GetXsec(E)/E )
                E += 1

c1 = TCanvas("c1","c1",1)
c1.Divide(4,4)

for i, g in enumerate(graphs):
    c1.cd(i+1)
    c1.GetPad(i+1).SetLogx()
    g.Draw()
    g.GetXaxis().SetRangeUser(1,100)
    g.GetYaxis().SetRangeUser(0,1)
    gPad.Modified()
    
