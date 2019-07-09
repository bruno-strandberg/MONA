#!/usr/bin/python
"""
This python script performs the same calculation as `llh_scan.C`

Usage:
   llh_scan.py [-s SUMMARYFILE] [-e EFFMASSFILE]
   llh_scan.py -h

Option:
  -s SUMMARYFILE  MC summary data in MONA format [default: ../data/ORCA_MCsummary_SEv2_ORCA115_20x9m_ECAP190222.root]
  -e EFFMASSFILE  Corresponding MONA effective mass file [default: ../data/eff_mass/EffMass_ORCA115_20x9m_ECAP190222.root]
  -h --help       Show this screen

"""
from docopt import docopt
args = docopt(__doc__)
from ROOT import *

#initialise a detector response, the numbers configure the binning
DR = DetResponse(DetResponse.track, "DR", 40, 1, 100, 40, -1, 1, 1, 0, 1);

#add selection cuts
DR.AddCut( "SummaryEvent.Get_RDF_track_score", ">", 0.7, True ); #ask PID track score to be larger than > 0.7
DR.AddCut( "SummaryEvent.Get_track_ql2"      , ">", 0.5, True ); #ask for ql2 (gandalf_loose_is_selected), see `MONA/apps/data_sorting/ECAP190222_20m_to_MONA.C`

#fill the response with MC events
sp = SummaryParser( args['-s'] );
for i in range( 0, sp.GetTree().GetEntries() ):
    if (i%500000 == 0): print( "Event: {}/{}".format(i, sp.GetTree().GetEntries()) )
    DR.Fill( sp.GetEvt(i) );
    
#initialise fitter utility and a PDF with the above response; the first number is operation time in years,
#the rest define the fit range. FitUtil is the worker class where calculations are done and parameters are defined
fu = FitUtil(3, DR.GetHist3DTrue(), DR.GetHist3DReco(), 3, 80, -1, -1e-5, 0, 1, args['-e'])

#init the probability density function; this is just a wrapper class that uses FitUtil
#and DetResponse to create a RooFit probability density function
pdf = FitPDF("pdf","pdf", fu, DR)

#create a pseudo-experiment, 3 years of data taking (run time determined at FitUtil construction)
D = pdf.SimplePseudoExp("Data");

#Import data to RooFit
RFD = RooDataHist("RFD","RFD", fu.GetObs(), RooFit.Import(D) );

#create a likelihood function between the pdf and the data
nll = pdf.createNLL( RFD );

#change th23 and evaluate the likelihood
fu.GetVar("SinsqTh23").setVal(0.6);
print( "NOTICE llh_scan() log-likelihood: {}".format( nll.getVal() ) )

#let RooFit to scan the likelihood
frame = fu.GetVar("SinsqTh23").frame();
nll.plotOn( frame, RooFit.ShiftToZero() );

c1 = TCanvas("c1","c1",1);
D.Project3D("yx").Draw("colz");
  
c2 = TCanvas("c2","c2",1);
frame.Draw();

hang = raw_input("Press enter to exit")
