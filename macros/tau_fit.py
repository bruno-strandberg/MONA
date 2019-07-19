#!/usr/bin/python
"""

This script can be used to create a profile-likelhood scan in tau. Note that the profile-likelihood entails many fits, this script can take a few hours to run. The exclusion-bands for tau normalisation can be read directly from the TGraph in the output file.

Usage:
   tau_fit.py -s SUMMARYFILE -e EFFMFILE -t TIME -o OUTFILE
   tau_fit.py -h

Option:
   -s SUMMARYFILE  MONA-format ORCA MC summary data
   -e EFFMFILE     Corresponding effective mass file
   -t TIME         Run time in months
   -o OUTFILE      Output file where the TGraph with the profile LLH curve is written
   -h --help       Show this screen

"""

from docopt import docopt
args = docopt(__doc__)

from ROOT import *

#===========================================================================================
# initialise responses & fill
#===========================================================================================

ebins  = 40
ctbins = 40
bybins = 1

# response for clearly track-like events
trkR = DetResponse(DetResponse.track, "trkR", ebins, 1, 100, ctbins, -1, 1, bybins, 0, 1 )
trkR.AddCut("SummaryEvent.Get_track_ql0"      , ">", 0.5 , True)
trkR.AddCut("SummaryEvent.Get_track_ql2"      , ">", 0.5 , True)
trkR.AddCut("SummaryEvent.Get_RDF_track_score", ">", 0.7 , True)
trkR.AddCut("SummaryEvent.Get_RDF_muon_score" , "<", 0.01, True)
trkR.AddCut("SummaryEvent.Get_RDF_noise_score", "<", 0.01, True)

# response for clearly shower-like events
shwR = DetResponse(DetResponse.shower, "shwR", ebins, 1, 100, ctbins, -1, 1, bybins, 0, 1 )
shwR.AddCut("SummaryEvent.Get_shower_ql0"     , ">", 0.5 , True)
shwR.AddCut("SummaryEvent.Get_shower_ql2"     , ">", 0.5 , True)
shwR.AddCut("SummaryEvent.Get_RDF_track_score", "<", 0.3 , True)
shwR.AddCut("SummaryEvent.Get_RDF_muon_score" , "<", 0.01, True)
shwR.AddCut("SummaryEvent.Get_RDF_noise_score", "<", 0.01, True)

# response for the middle sample (using shower reco)
midR = DetResponse(DetResponse.shower, "midR", ebins, 1, 100, ctbins, -1, 1, bybins, 0, 1)
midR.AddCut("SummaryEvent.Get_shower_ql0"     , ">" , 0.5 , True)
midR.AddCut("SummaryEvent.Get_shower_ql2"     , ">" , 0.5 , True)
midR.AddCut("SummaryEvent.Get_RDF_track_score", ">=", 0.3 , True)
midR.AddCut("SummaryEvent.Get_RDF_track_score", "<" , 0.7 , True)
midR.AddCut("SummaryEvent.Get_RDF_muon_score" , "<" , 0.01, True)
midR.AddCut("SummaryEvent.Get_RDF_noise_score", "<" , 0.01, True)

sp = SummaryParser( args['-s'] )
for i in range( 0, sp.GetTree().GetEntries() ):
    if ( i % 250000 == 0): print("Event: {}/{}".format(i, sp.GetTree().GetEntries()) )
    trkR.Fill( sp.GetEvt(i) )
    shwR.Fill( sp.GetEvt(i) )
    midR.Fill( sp.GetEvt(i) )

#===========================================================================================
# init FitUtil, pdf's and simultaneous fitting
#===========================================================================================

time = float(args['-t'])/12

# main worker class where fit parameters are defined; run time is the first argument
fu = FitUtilWsyst(time, trkR.GetHist3DTrue(), trkR.GetHist3DReco(), 3, 70, -1, -1e-5, 0, 1, args['-e'])

#roofit-style probability density functions for track and shower
trkPdf = FitPDF("trkPdf","trkPdf", fu, trkR)
shwPdf = FitPDF("shwPdf","shwPdf", fu, shwR)
midPdf = FitPDF("midpdf","midpdf", fu, midR)

# categories to link data <--> pdf
categs = RooCategory("categs", "data categories")
categs.defineType("TRK")
categs.defineType("SHW")
categs.defineType("MID")

# simultaneous pdf
simPdf = RooSimultaneous("simPdf","simPdf", categs)
simPdf.addPdf( trkPdf, "TRK")
simPdf.addPdf( shwPdf, "SHW")
simPdf.addPdf( midPdf, "MID")

#===========================================================================================
# create simultaneous dataset
#===========================================================================================

# just to illustrate that we find some sort of a minimum around 1.1
fu.GetVar("Tau_norm").setVal(1.1)

# create trk-channel data
trkData = trkPdf.GetExpValHist()                                                             # create expectation value data (pars and time configured in FitUtil)
trkData.SetNameTitle("trkdata","trkdata")                                                    # set name and title (needs to be unique for ROOT to work)
rf_trkData = RooDataHist("rf_trkdata", "rf_trkdata", fu.GetObs(), RooFit.Import(trkData) )   # import data to RooFIt

# same story for showers
shwData = shwPdf.GetExpValHist()
shwData.SetNameTitle("shwdata","shwdata")
rf_shwData = RooDataHist("rf_shwdata", "rf_shwdata", fu.GetObs(), RooFit.Import(shwData) )

# same story for middle sample
midData = midPdf.GetExpValHist()
midData.SetNameTitle("middata","middata")
rf_midData = RooDataHist("rf_middata", "rf_middata", fu.GetObs(), RooFit.Import(midData) )

# combine the datasets in RooFit
combD = RooDataHist("combd","combd", fu.GetObs(), RooFit.Index(categs), RooFit.Import("TRK", rf_trkData), RooFit.Import("SHW", rf_shwData), RooFit.Import("MID", rf_midData) )

#===========================================================================================
# fix some parameters are create priors for others - do `monafitpars -f FitUtilWsyst` in the
# terminal to see all fit parameters
#===========================================================================================

# fix solar parameters and e-scale (has little effect as currently implemented)
fu.GetVar("SinsqTh12").setConstant(True)
fu.GetVar("Dm21").setConstant(True)
fu.GetVar("E_scale").setConstant(True)

# add some priors for the systematics
p_mu_amu = RooGaussian("mu_amu_prior", "mu_amu_prior", fu.GetVar("skew_mu_amu"), RooFit.RooConst(0.), RooFit.RooConst(0.1) )
p_e_ae   = RooGaussian("e_ae_prior"  , "e_ae_prior"  , fu.GetVar("skew_e_ae")  , RooFit.RooConst(0.), RooFit.RooConst(0.1) )
p_mu_e   = RooGaussian("mu_e_prior"  , "mu_e_prior"  , fu.GetVar("skew_mu_e")  , RooFit.RooConst(0.), RooFit.RooConst(0.1) )
pncnorm  = RooGaussian("ncnorm_prior", "ncnorm_prior", fu.GetVar("NC_norm")    , RooFit.RooConst(1.), RooFit.RooConst(0.1) )

# create a set of priors
priors = RooArgSet( p_mu_amu, p_e_ae, p_mu_e, pncnorm )

#===========================================================================================
# do a manual scan in tau-variable to calculate the profile llh
#===========================================================================================

# first fit the dataset with free tau normalisation to get the likelihood minimum
result  = simPdf.fitTo( combD, RooFit.SumW2Error(kFALSE), RooFit.ExternalConstraints( priors ), RooFit.Save() )
llh_min = result.minNll()

# now scan in tau normalisation in 10 steps
graph = TGraph()
tau_norm = 0.5
while tau_norm <= 1.5:

    print("NOTICE: minimizing at tau value {}".format(tau_norm) )
    
    fu.GetVar("Tau_norm").setVal(tau_norm)
    fu.GetVar("Tau_norm").setConstant(True)
    result = simPdf.fitTo( combD, RooFit.SumW2Error(kFALSE), RooFit.ExternalConstraints( priors ), RooFit.Save() )
    llh = result.minNll()
    graph.SetPoint( graph.GetN(), tau_norm, llh - llh_min )
    tau_norm += 0.1
    
# draw the graph and write the graph to file
graph.SetNameTitle("profile_llh", "profile_llh")
graph.Draw()

fout = TFile(args['-o'], "RECREATE")
graph.Write()
fout.Close()

hang = raw_input("Press enter to exit")
