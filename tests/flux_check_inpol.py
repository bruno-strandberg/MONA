#!/usr/bin/python -i
from ROOT import *
import collections
import math
gSystem.Load("$NMHDIR/common_software/libnmhsoft.so")

#----------------------------------------------------------------------------
#This creates comparison graphs between the AtmFlux calculator value and actual data points
#along the cosz axis
#----------------------------------------------------------------------------

def comp_cosz(f_dict, Energy = 1.0000E+00, flavor = 1, is_nubar = 0):

    # init graphs, load library
    
    flux = AtmFlux(h_frj_min)
    g_inpol = TGraph()
    g_data  = TGraph()

    # add points to data and interpolated graphs
    
    cosz   = -0.95      # for the dictionary data this needs to be between [-0.95, 0.95] +- 0.1
    while cosz <= 0.95:
        cosz = round(cosz, 3)
        g_inpol.SetPoint( g_inpol.GetN(),
                          cosz, flux.Flux_dE_dcosz(flavor, is_nubar, Energy, cosz)/(2*math.pi) )
        g_data.SetPoint (  g_data.GetN(), cosz, f_dict[Energy][cosz])
        cosz += 0.1

    # add extra points away from data to check interpolation is OK

    cosz = -1
    while cosz <= 1:
        g_inpol.SetPoint( g_inpol.GetN(),
                          cosz, flux.Flux_dE_dcosz(flavor, is_nubar, Energy, cosz)/(2*math.pi) )
        cosz += 0.1
        
    g_inpol.SetMarkerStyle(20)
    g_inpol.SetMarkerColor(kBlack)
    g_inpol.SetLineColor(kBlack)
    g_data.SetMarkerStyle(21)
    g_data.SetMarkerColor(kRed)
    g_data.SetLineColor(kRed)
            
    return g_inpol, g_data

#----------------------------------------------------------------------------
#This creates comparison graphs between the AtmFlux calculator value and actual data points
#along the energy axis
#----------------------------------------------------------------------------

def comp_E(f_dict, cosz = -0.95, flavor = 1, is_nubar = 0):

    # init graphs, load library
    
    flux = AtmFlux(h_frj_min)
    g_inpol = TGraph()
    g_data  = TGraph()

    # add points to data and interpolated graphs

    for E, angles in f_dict.iteritems():
        cosz = round(cosz, 3)
        g_inpol.SetPoint( g_inpol.GetN(),
                          E, flux.Flux_dE_dcosz(flavor, is_nubar, E, cosz)/(2*math.pi) )
        g_data.SetPoint (  g_data.GetN(), E, f_dict[E][cosz])

    # add extra points away from data to check interpolation is OK

    E = 50
    while E <= 10000:
        g_inpol.SetPoint( g_inpol.GetN(),
                          E, flux.Flux_dE_dcosz(flavor, is_nubar, E, cosz)/(2*math.pi) )
        E += 100
        
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

#----------------------------------------------------------------------------
#main
#----------------------------------------------------------------------------
if __name__=="__main__":

    flux_mu, flux_mub, flux_e, flux_eb = read_data()
    g1, g2 = comp_cosz(flux_mu)
    g3, g4 = comp_cosz(flux_mu, 1.0000E+01)
    g5, g6 = comp_cosz(flux_e, 1.0000E+00, 0, 0)
    g7, g8 = comp_cosz(flux_e, 1.0000E+01, 0, 0)
    g9, g10= comp_E(flux_mu)

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

    c2 = TCanvas("c1","c1",1)
    g9.Draw("AP")
    g10.Draw("Psame")
