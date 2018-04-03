#!/usr/bin/python -i
from ROOT import *
gSystem.Load("$NMHDIR/common_software/libnmhsoft.so")

# get martijn's histogram, create an equivalent histogram for my calculation

fmj = TFile("mjongen_plots.root","READ");
hmj = fmj.Get("NH/Atmospheric fluxes/hAtmosphericFlux_num");

hme = hmj.Clone();
hme.Reset();

# init the flux, set flux for my histogram

f = AtmFlux()

for xbin in range(1, hme.GetXaxis().GetNbins()+1):
    for ybin in range(1, hme.GetYaxis().GetNbins()+1):

        E    = hme.GetXaxis().GetBinCenter(xbin)
        cosz = hme.GetYaxis().GetBinCenter(ybin)
        
        flux_me = ( f.Flux_dE_dcosz(1, 0, E, cosz) *
	            hme.GetXaxis().GetBinWidth(xbin) * hme.GetYaxis().GetBinWidth(ybin) )

        hme.SetBinContent(xbin, ybin, flux_me);
        
# create 1D projections for comparison
        
h_proj_me = hme.ProjectionX("proj_me", 5, 5)
h_proj_mj = hmj.ProjectionX("proj_mj", 5, 5)

# subtract martijn's plot from mine
        
hme.Add(hmj, -1);

c1 = TCanvas("c1","c1",1);
c1.Divide(2,2);
c1.cd(1);
hmj.Draw();
c1.cd(2);
hme.Draw();
c1.cd(3);
h_proj_mj.Draw();
c1.cd(4);
h_proj_me.Draw();
h_proj_mj.SetLineColor(kRed);
h_proj_mj.Draw("same");
