
#include "DetResponse.h"
#include "SummaryEvent.h"
#include "SummaryParser.h"
#include "NMHUtils.h"

#include "FitUtil.h"
#include "FitPDF.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooMinimizer.h"

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF3.h"
#include "TRandom3.h"

#include<stdexcept>

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

/** Namespace that stores structures and functions for the `FitPseudoExperiment` app*/
namespace PSEUDOFIT {

  /** class to store fit parameter instances*/
  struct FitPars {

    Double_t sinsqth12;
    Double_t sinsqth13;
    Double_t sinsqth23;
    Double_t dcp;
    Double_t dm21;
    Double_t dm31;

  };

  /** class to store fit trial results*/
  struct FitTrial {

    FitPars true_vals;    //!< true values
    FitPars roofit_vals;  //!< roofit results
    FitPars roofit_errs;  //!< roofit errors
    FitPars root_vals;    //!< root results
    FitPars root_errs;    //!< root errors

  };

  void ConfigureTree(TTree&, FitTrial&);

};


using namespace RooFit;
using namespace std;
using namespace PSEUDOFIT;

/** This example application demonstrates the generation of a pseudo-experiment in the track channel and fitting the pseudo-experiment using ROOT and RooFit.

    The example uses the `track` and `shower` event selections as were optimised for ORCA115_23x9m MC production; for a different production the selections of this example are probably not optimal.

    This example demonstrates:
    1. the creation of a compiled program with `NMO` libraries and `Jpp` headers;
    2. use of the `DetResponse`, `EventSelection`, `FitUtil`, `FitPDF` classes and `RooFit`;
    3. creation of a simple pseudo-experiment.

*/
int main(const int argc, const char **argv) {

  TString rfd = NMHUtils::Getcwd() + "/rootfiles/";

  //----------------------------------------------------------
  // parse command line arguments with Jpp
  //----------------------------------------------------------
  
  TString       simdata_file;
  bool          refill_response;
  TString       outputfile;
  int           nfits;
  TString       meff_file;
  int           ncpu;

  try {

    JParser<> zap("This example application demonstrates the generation of a pseudo-experiment in the track channel and fitting the pseudo-experiment using ROOT and RooFit. See the documentation in the application for further information.");

    zap['s'] = make_field(simdata_file, "File with all summary data") =
      (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root";

    zap['r'] = make_field(refill_response, "Flag to request re-filling of the detector responses");
    zap['o'] = make_field(outputfile, "File where output histograms are written") = rfd + "pseudoexpfit.root";
    zap['n'] = make_field(nfits, "Number of fits to be performed") = 1;
    zap['M'] = make_field(meff_file, "Effective mass file created by using `EffMass` class") = 
      (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
    zap['N'] = make_field(ncpu, "Number of CPUs when fitting with RooFit") = 1;

    if ( zap.read(argc, argv)!= 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }
  
  Int_t sysret = system("mkdir -p " + rfd);
  if (sysret != 0) { throw std::invalid_argument( "ERROR! ExpctFit1H cannot create dir " + (string)rfd ); }

  //----------------------------------------------------------
  // set up the detector response binning
  //----------------------------------------------------------
  
  Int_t fENB   =  24;
  Int_t fEmin  =   1;
  Int_t fEmax  = 100;
  Int_t fCtNB  =  40;
  Int_t fCtmin =  -1;
  Int_t fCtmax =   1;
  Int_t fByNB  =   1;
  Int_t fBymin =   0;
  Int_t fBymax =   1;
  
  DetResponse track_resp (DetResponse::track, "track_resp", 
			  fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  track_resp.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_resp.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );    

  TString track_resp_name  = rfd + "fitpseudoexperiment_response.root";
  
  if ( !NMHUtils::FileExists(track_resp_name) || refill_response ) {
    cout << "NOTICE ExpctFit1H (Re)filling response" << endl;

    SummaryParser sp(simdata_file);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i % 100000 == 0) cout << "Event: " << i << endl;
      track_resp.Fill( sp.GetEvt(i) );
    }
    
    track_resp.WriteToFile(track_resp_name);
  }
  else {
    cout << "NOTICE ExpctFit1H Reading in response" << endl;
    track_resp.ReadFromFile(track_resp_name);
  }

  cout << "NOTICE ExpctFit1H Response ready" << endl;

  //----------------------------------------------------------
  // set up the PDF and oscillation parameters limits for normal ordering,
  // get roofit variable pointers; fix th12 and dm21
  //----------------------------------------------------------

  // the fit limits in the observables are configured through `FitUtil` initialisation
  Int_t    runtime   = 3;
  Double_t fit_emin  = 1;
  Double_t fit_emax  = 100;
  Double_t fit_ctmin = -1;
  Double_t fit_ctmax = 0;
  Double_t fit_bymin = 0;
  Double_t fit_bymax = 1;
  
  FitUtil *fitutil = new FitUtil(runtime, track_resp.GetHist3D(), fit_emin, fit_emax,
				 fit_ctmin, fit_ctmax, fit_bymin, fit_bymax, meff_file);

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_resp);  

  fitutil->SetNOlims();
  fitutil->SetNOcentvals();

  RooRealVar* rf_sinsqth12 = fitutil->GetVar("SinsqTh12");
  RooRealVar* rf_sinsqth13 = fitutil->GetVar("SinsqTh13");
  RooRealVar* rf_sinsqth23 = fitutil->GetVar("SinsqTh23");
  RooRealVar* rf_dcp       = fitutil->GetVar("dcp");
  RooRealVar* rf_dm21      = fitutil->GetVar("Dm21");
  RooRealVar* rf_dm31      = fitutil->GetVar("Dm31");

  rf_sinsqth12->setConstant(kTRUE);
  rf_dm21->setConstant(kTRUE);

  //----------------------------------------------------------
  // set-up the output tree
  //----------------------------------------------------------
  
  TFile fout((TString)outputfile, "RECREATE");
  TTree tout("fittrials","fit trial results");
  FitTrial trial;
  ConfigureTree(tout, trial);
  
  //----------------------------------------------------------
  // pereform the requested number of fits
  //----------------------------------------------------------

  TRandom3 frand(0);
  RooRandom::randomGenerator()->SetSeed( frand.Uniform(1,1e6) );
  
  for (Int_t fitnr = 0; fitnr < nfits; fitnr++) {

    // randomize the fitted oscillation parameters
    rf_sinsqth13->randomize();
    rf_sinsqth23->randomize();
    rf_dcp->randomize();
    rf_dm31->randomize();
    
    // get the expectation value
    TH3D *tracks_exp = pdf_tracks.SimplePseudoExp("tracks_exp",kTRUE);
    
    // store the true values of the expectation
    trial.true_vals.sinsqth12 = rf_sinsqth12->getVal();
    trial.true_vals.sinsqth13 = rf_sinsqth13->getVal();
    trial.true_vals.sinsqth23 = rf_sinsqth23->getVal();
    trial.true_vals.dcp       = rf_dcp->getVal();
    trial.true_vals.dm21      = rf_dm21->getVal();
    trial.true_vals.dm31      = rf_dm31->getVal();

    cout << "*************************************************************************" << endl;
    cout << "Started fitting with RooFit" << endl;
    cout << "*************************************************************************" << endl;

    // reset to NO central values
    fitutil->SetNOcentvals();

    // fit
    RooDataHist rf_hist("rf_hist", "rf_hist", fitutil->GetObs(), Import(*tracks_exp) );
    RooFitResult *fitres = pdf_tracks.fitTo( rf_hist, Save(kTRUE), Offset(kTRUE), NumCPU(ncpu) );
    RooArgSet result( fitres->floatParsFinal() );

    // store roofit results and errors
    trial.roofit_vals.sinsqth12 = rf_sinsqth12->getVal();
    trial.roofit_vals.sinsqth13 = ((RooRealVar*)result.find("SinsqTh13"))->getVal();
    trial.roofit_vals.sinsqth23 = ((RooRealVar*)result.find("SinsqTh23"))->getVal();
    trial.roofit_vals.dcp       = ((RooRealVar*)result.find("dcp"))->getVal();
    trial.roofit_vals.dm21      = rf_dm21->getVal();
    trial.roofit_vals.dm31      = ((RooRealVar*)result.find("Dm31"))->getVal();
    
    trial.roofit_errs.sinsqth12 = 0;
    trial.roofit_errs.sinsqth13 = ((RooRealVar*)result.find("SinsqTh13"))->getError();
    trial.roofit_errs.sinsqth23 = ((RooRealVar*)result.find("SinsqTh23"))->getError();
    trial.roofit_errs.dcp       = ((RooRealVar*)result.find("dcp"))->getError();
    trial.roofit_errs.dm21      = 0;
    trial.roofit_errs.dm31      = ((RooRealVar*)result.find("Dm31"))->getError();

    cout << "*************************************************************************" << endl;
    cout << "Started fitting with ROOT" << endl;
    cout << "*************************************************************************" << endl;

    // reset to NO central values
    fitutil->SetNOcentvals();    

    // configure fit and fit
    TF3 *func = new TF3("fitfunc", pdf_tracks, fit_emin, fit_emax, fit_ctmin, fit_ctmax, fit_bymin, fit_bymax, 6);
    func->SetParameters( rf_sinsqth12->getVal(), rf_sinsqth13->getVal(), rf_sinsqth23->getVal(),
			 rf_dcp->getVal(), rf_dm21->getVal(), rf_dm31->getVal() );
    func->SetParNames("sinsq12","sinsq13","sinsq23","dcp","dm21","dm31");
    
    func->FixParameter( 0, rf_sinsqth12->getVal() );
    func->SetParLimits( 1, rf_sinsqth13->getMin(), rf_sinsqth13->getMax() );
    func->SetParLimits( 2, rf_sinsqth23->getMin(), rf_sinsqth23->getMax() );
    func->SetParLimits( 3, rf_dcp->getMin(), rf_dcp->getMax() );
    func->FixParameter( 4, rf_dm21->getVal() );
    func->SetParLimits( 5, rf_dm31->getMin(), rf_dm31->getMax() );

    TFitResultPtr rootres = tracks_exp->Fit("fitfunc", "RSLV");

    // store ROOT results
    trial.root_vals.sinsqth12 = rf_sinsqth12->getVal();
    trial.root_vals.sinsqth13 = rootres->Parameter(1);
    trial.root_vals.sinsqth23 = rootres->Parameter(2);
    trial.root_vals.dcp       = rootres->Parameter(3);
    trial.root_vals.dm21      = rf_dm21->getVal();
    trial.root_vals.dm31      = rootres->Parameter(5);
    
    trial.root_errs.sinsqth12 = 0;
    trial.root_errs.sinsqth13 = rootres->ParError(1);
    trial.root_errs.sinsqth23 = rootres->ParError(2);
    trial.root_errs.dcp       = rootres->ParError(3);
    trial.root_errs.dm21      = 0;
    trial.root_errs.dm31      = rootres->ParError(5);

    tout.Fill();

    if (tracks_exp) delete tracks_exp;
    
  } // end loop over trials  

  cout << "NOTICE ExpctFit1H finished loop" << endl;

  tout.Write();
  fout.Close();
  
}

//***************************************************************************************

/** Inline function that configures the output tree with fit results. */
void PSEUDOFIT::ConfigureTree(TTree& tout, FitTrial& trial) {

  tout.Branch("true_vals_sinsqth12"  , &trial.true_vals.sinsqth12);
  tout.Branch("true_vals_sinsqth13"  , &trial.true_vals.sinsqth13);
  tout.Branch("true_vals_sinsqth23"  , &trial.true_vals.sinsqth23);
  tout.Branch("true_vals_dcp"        , &trial.true_vals.dcp);
  tout.Branch("true_vals_dm21"       , &trial.true_vals.dm21);
  tout.Branch("true_vals_dm31"       , &trial.true_vals.dm31);

  tout.Branch("roofit_vals_sinsqth12", &trial.roofit_vals.sinsqth12);
  tout.Branch("roofit_vals_sinsqth13", &trial.roofit_vals.sinsqth13);
  tout.Branch("roofit_vals_sinsqth23", &trial.roofit_vals.sinsqth23);
  tout.Branch("roofit_vals_dcp"      , &trial.roofit_vals.dcp);
  tout.Branch("roofit_vals_dm21"     , &trial.roofit_vals.dm21);
  tout.Branch("roofit_vals_dm31"     , &trial.roofit_vals.dm31);
  
  tout.Branch("roofit_errs_sinsqth12", &trial.roofit_errs.sinsqth12);
  tout.Branch("roofit_errs_sinsqth13", &trial.roofit_errs.sinsqth13);
  tout.Branch("roofit_errs_sinsqth23", &trial.roofit_errs.sinsqth23);
  tout.Branch("roofit_errs_dcp"      , &trial.roofit_errs.dcp);
  tout.Branch("roofit_errs_dm21"     , &trial.roofit_errs.dm21);
  tout.Branch("roofit_errs_dm31"     , &trial.roofit_errs.dm31);

  tout.Branch("root_vals_sinsqth12"  , &trial.root_vals.sinsqth12);
  tout.Branch("root_vals_sinsqth13"  , &trial.root_vals.sinsqth13);
  tout.Branch("root_vals_sinsqth23"  , &trial.root_vals.sinsqth23);
  tout.Branch("root_vals_dcp"        , &trial.root_vals.dcp);
  tout.Branch("root_vals_dm21"       , &trial.root_vals.dm21);
  tout.Branch("root_vals_dm31"       , &trial.root_vals.dm31);
  
  tout.Branch("root_errs_sinsqth12"  , &trial.root_errs.sinsqth12);
  tout.Branch("root_errs_sinsqth13"  , &trial.root_errs.sinsqth13);
  tout.Branch("root_errs_sinsqth23"  , &trial.root_errs.sinsqth23);
  tout.Branch("root_errs_dcp"        , &trial.root_errs.dcp);
  tout.Branch("root_errs_dm21"       , &trial.root_errs.dm21);
  tout.Branch("root_errs_dm31"       , &trial.root_errs.dm31);

}
