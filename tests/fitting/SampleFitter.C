// fitter_software headers
#include "FitUtil.h"
#include "FitPDF.h"

// common_software headers
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "NMHUtils.h"
#include "FileHeader.h"
#include "DetResponse.h"
#include "EventSelection.h"

// root headers
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"

// roofit headers
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"

// cpp headers
#include <iostream>
#include <stdexcept>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

using namespace RooFit;
using namespace std;

/** Namespace that collects functions and variables for the `SampleFitter` application */
namespace FITTER {

  // functions
  void InitRespsAndSels();
  void FillRespsAndSels(TString simdata_file, TString expdata_file, Bool_t refill_response);
  void FillSelections();
  void Cleanup();
  
  // binning variables
  static const Int_t    fENB    =  40;
  static const Double_t fEmin   =   1;
  static const Double_t fEmax   = 100;
  static const Int_t    fCtNB   =  40;
  static const Double_t fCtmin  =  -1;
  static const Double_t fCtmax  =   1;
  static const Int_t    fByNB   =   1;
  static const Double_t fBymin  =   0;
  static const Double_t fBymax  =   1;

  // variables for FitUtil configuration
  static const Double_t fRunTime  = 3;   // 3 years running time
  static const Double_t fFitEmin  = 3;   // minimum fitted energy range
  static const Double_t fFitEmax  = 80;  // maximum fitted energy range
  static const Double_t fFitCTmin = -1;  // minimum fitted cos-theta range
  static const Double_t fFitCTmax =  0;  // maximum fitted cos-theta range
  static const Double_t fFitBYmin =  0;  // minimum fitted bjorken-y range
  static const Double_t fFitBYmax =  1;  // maximum fitted bjorken-y range
  
  // member selections and responses
  DetResponse    *fTRres; // detector response for tracks
  DetResponse    *fSHres; // detector response for showers
  EventSelection *fTRsel; // event selection for tracks
  EventSelection *fSHsel; // event selection for showers

};

/** This example application can be used to fit an event-by-event experiment sample, created with applications in `apps/evt_sampler`. 

    The example uses the `track` and `shower` event selections as were optimised for ORCA115_23x9m MC production; for a different production the selections of this example are probably not optimal.

    This example demonstrates:
    1. the creation of a compiled program with `NMO` libraries and `Jpp` headers;
    2. use of the `DetResponse`, `EventSelection`, `FitUtil`, `FitPDF` classes and `RooFit`;
    3. simultaneous fitting of more than one event selections.

 
*/
int main(const int argc, const char **argv) {

  using namespace FITTER;

  //----------------------------------------------------------
  // parse command line arguments with Jpp
  //----------------------------------------------------------
  
  string        simdata_file;
  string        expdata_file;
  bool          refill_response;
  TString       meff_file;
  Int_t         numcpu;

  try {

    JParser<> zap("This application can be used to fit an event-by-event experiment sample, created with applications in `apps/evt_sampler`. See documentation in the source code for more info.");

    zap['s'] = make_field(simdata_file, "File with all summary data") =
      (string)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root";

    zap['e'] = make_field(expdata_file, "File with experimental data sample");

    zap['r'] = make_field(refill_response, "Flag to request re-filling of the detector responses");
    
    zap['M'] = make_field(meff_file, "Effective mass file created by using `EffMass` class") = 
      (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
    zap['n'] = make_field(numcpu, "Number of CPUs for fitting") = 1;

    if ( zap.read(argc, argv)!= 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  //----------------------------------------------------------
  // call methods to initialise and fill responses
  //----------------------------------------------------------
  InitRespsAndSels();
  FillRespsAndSels(simdata_file, expdata_file, refill_response);

  //----------------------------------------------------------
  // set up the PDFs for fitting
  //----------------------------------------------------------
  FitUtil *fitutil = new FitUtil(fRunTime, fTRres->GetHist3D(), fFitEmin, fFitEmax,
				 fFitCTmin, fFitCTmax, fFitBYmin, fFitBYmax, meff_file);

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, fTRres);  
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, fSHres);

  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------
  std::map<string, TH1* > hist_map = { { (string)fTRsel->Get_SelName() , fTRsel->Get_h_E_costh_by() }, 
				       { (string)fSHsel->Get_SelName() , fSHsel->Get_h_E_costh_by() } };
  
  // create categories
  RooCategory categs("categs","data categories");
  categs.defineType( fTRsel->Get_SelName() );
  categs.defineType( fSHsel->Get_SelName() );
  
  // create combined data set
  RooDataHist combData("combData", "combined data", fitutil->GetObs(), categs, hist_map);

  // create simultaneous pdf
  RooSimultaneous simPdf("simPdf", "simultaneous Pdf", categs);
  simPdf.addPdf(pdf_tracks , fTRsel->Get_SelName() );
  simPdf.addPdf(pdf_showers, fSHsel->Get_SelName() );

  cout << "NOTICE Fitter started fitting" << endl;
  
  TStopwatch timer;
  RooFitResult *fitres = simPdf.fitTo( combData, Save(kTRUE), NumCPU(numcpu) );
  cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;
  
  //----------------------------------------------------------
  // print comparison
  //----------------------------------------------------------
  RooArgSet result ( fitres->floatParsFinal() );
  FileHeader h("Fitter");
  h.ReadHeader(expdata_file);

  cout << "*********Fit result comparison****************************" << endl;
  cout << "dm21x1e5 actual  : " << h.GetParameter("dm21_1e5")   << "\t" << " fitted: " << ((RooRealVar*)result.find("Dm21"))->getVal()*1e5 << endl;
  cout << "dm31x1e5 actual  : " << h.GetParameter("dm31_1e5")   << "\t" << " fitted: " << ((RooRealVar*)result.find("Dm31"))->getVal()*1e5 << endl;
  cout << "sinsq_th12 actual: " << h.GetParameter("sinsq_th12") << "\t" << " fitted: " << ((RooRealVar*)result.find("SinsqTh12"))->getVal() << endl;
  cout << "sinsq_th13 actual: " << h.GetParameter("sinsq_th13") << "\t" << " fitted: " << ((RooRealVar*)result.find("SinsqTh13"))->getVal() << endl;
  cout << "sinsq_th23 actual: " << h.GetParameter("sinsq_th23") << "\t" << " fitted: " << ((RooRealVar*)result.find("SinsqTh23"))->getVal() << endl;
  cout << "dcp        actual: " << h.GetParameter("dcp")        << "\t" << " fitted: " << ((RooRealVar*)result.find("dcp"))->getVal() << endl;
  cout << "*********Fit result comparison****************************" << endl;
  
  //----------------------------------------------------------
  // clean-up
  //----------------------------------------------------------
  
  Cleanup();
  if (fitutil) delete fitutil;
 
}

//**********************************************************************************

/** Internal function to initialise detector responses and event selections */
void FITTER::InitRespsAndSels() {

  using namespace FITTER;

  //----------------------------------------------------------
  // initialise event selection and response for tracks; the event selection
  // and the corresponding response need to have the same selection criteria
  //----------------------------------------------------------

  fTRsel = new EventSelection(EventSelection::track, "track_sel", NULL, 
			      fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  fTRres = new DetResponse(DetResponse::track, "track_resp", 
			   fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  vector< EventFilter* > T = {fTRsel, fTRres};

  for (auto n: T) {
    n->AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
    n->AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    n->AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );    
  }

  //----------------------------------------------------------
  // initialise event selection and response for showers; the event selection
  // and the corresponding response need to have the same selection criteria
  //----------------------------------------------------------

  fSHsel = new EventSelection(EventSelection::shower, "shower_sel", NULL,
			      fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);
  fSHres = new DetResponse(DetResponse::shower, "shower_resp", 
			   fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  vector< EventFilter* > S = {fSHsel, fSHres};

  for (auto n: S) {
    n->AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
    n->AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    n->AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );
  }

}

//**********************************************************************************

/** Internal function to fill the event selections and responses. */
void FITTER::FillRespsAndSels(TString simdata_file, TString expdata_file, Bool_t refill_response) {

  using namespace FITTER;

  if ( !fTRres || !fSHres || !fTRsel || !fSHsel ) {
    throw std::logic_error("ERROR! FITTER::FillRespsAndSels() DetResponse's and EventSelection's are not yet initialized!");
  }

  //----------------------------------------------------------
  // fill the responses
  //----------------------------------------------------------

  TString track_resp_name  = NMHUtils::Getcwd() + "/rootfiles/samplefitter_track_response.root";
  TString shower_resp_name = NMHUtils::Getcwd() + "/rootfiles/samplefitter_shower_response.root";

  if ( !NMHUtils::FileExists(track_resp_name) || !NMHUtils::FileExists(shower_resp_name) || refill_response ) {

    cout << "NOTICE FITTER::FillRespsAndSels() (Re)filling responses" << endl;

    SummaryParser sp(simdata_file);

    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

      if (i % 100000 == 0) cout << "Event: " << i << endl;

      SummaryEvent *evt = sp.GetEvt(i);
      fTRres->Fill(evt);
      fSHres->Fill(evt);

    }

    fTRres->WriteToFile(track_resp_name);
    fSHres->WriteToFile(shower_resp_name);

  }
  else {

    cout << "NOTICE FITTER::FillRespsAndSels() Reading in responses" << endl;

    fTRres->ReadFromFile(track_resp_name);
    fSHres->ReadFromFile(shower_resp_name);
  }

  cout << "NOTICE FITTER::FillRespsAndSels() Responses ready" << endl;

  //----------------------------------------------------------
  // fill the selections
  //----------------------------------------------------------
  SummaryParser sp(expdata_file);

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

    SummaryEvent *evt = sp.GetEvt(i);
    fTRsel->Fill(evt);
    fSHsel->Fill(evt);

  }

  cout << "NOTICE FITTER::FillRespsAndSels() Selections filled" << endl;
}

//**********************************************************************************

/** Interal function to clean dynamic memory. */
void FITTER::Cleanup() {

  using namespace FITTER;

  if (fTRres) delete fTRres;
  if (fSHres) delete fSHres;
  if (fTRsel) delete fTRsel;
  if (fSHsel) delete fSHsel;

}

