#!/usr/bin/python -i
from ROOT import *
import collections
gSystem.Load("../common_software/libcommonsoft.so")

#----------------------------------------------------------------------------
#This creates comparison graphs between the AtmFlux calculator value and actual data points
#----------------------------------------------------------------------------

def comp(f_dict, Energy = 1.0000E+00, flavor = 1, is_nubar = 0):

    # init graphs, load library
    
    flux = AtmFlux(h_frj_min)
    g_inpol = TGraph()
    g_data  = TGraph()

    # add points to data and interpolated graphs
    
    cosz   = -0.95      # for the dictionary data this needs to be between [-0.95, 0.95] +- 0.1
    while cosz <= 0.95:
        cosz = round(cosz, 3)
        g_inpol.SetPoint( g_inpol.GetN(), cosz, flux.GetHondaFlux(flavor, is_nubar, Energy, cosz) )
        g_data.SetPoint (  g_data.GetN(), cosz, f_dict[Energy][cosz])
        cosz += 0.1

    # add extra points away from data to check interpolation is OK

    cosz = -1
    while cosz <= 1:
        g_inpol.SetPoint( g_inpol.GetN(), cosz, flux.GetHondaFlux(flavor, is_nubar, Energy, cosz) )
        cosz += 0.1
        
    g_inpol.SetMarkerStyle(20)
    g_inpol.SetMarkerColor(kBlack)
    g_inpol.SetLineColor(kBlack)
    g_data.SetMarkerStyle(21)
    g_data.SetMarkerColor(kRed)
    g_data.SetLineColor(kRed)
            
    return g_inpol, g_data

#----------------------------------------------------------------------------
#This function reads the the data from a honda flux file for comparison
#----------------------------------------------------------------------------
def read_data():

    f = open("../data/honda_flux/frj-ally-20-01-solmin.d")

    cosz = 2

    # 2d dictionaries
    flux_mu  = collections.defaultdict(dict)
    flux_mub = collections.defaultdict(dict)
    flux_e   = collections.defaultdict(dict)
    flux_eb  = collections.defaultdict(dict)
        
    for line in f:

        #get the cosz
        if ("average" in line):
            cosz = round(float(line[ line.index("=")+1 : line.index("--") ]) + 0.05, 3)
            continue

        if ("Enu" in line): continue

        E, f_mu, f_mub, f_e, f_eb = map(double, line.split())
        flux_mu [E][cosz] = f_mu
        flux_mub[E][cosz] = f_mub
        flux_e  [E][cosz] = f_e
        flux_eb [E][cosz] = f_eb

    # end loop
    f.close()

    return flux_mu, flux_mub, flux_e, flux_eb
    
if __name__=="__main__":

    flux_mu, flux_mub, flux_e, flux_eb = read_data()
    g1, g2 = comp(flux_mu)
    g3, g4 = comp(flux_mu, 1.0000E+01)
    g5, g6 = comp(flux_e, 1.0000E+00, 0, 0)
    g7, g8 = comp(flux_e, 1.0000E+01, 0, 0)

    c1 = TCanvas("c1","c1",1)
    c1.Divide(2,2)

    c1.cd(1)
    g1.Draw()
    g2.Draw("Psame")

    c1.cd(2)
    g3.Draw()
    g4.Draw("Psame")

    c1.cd(3)
    g5.Draw()
    g6.Draw("Psame")

    c1.cd(4)
    g7.Draw()
    g8.Draw("Psame")
