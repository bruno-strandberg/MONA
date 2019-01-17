//#include "Fitter.h"

// fitter headers
#include "FitUtil.h"
#include "FitPDF.h"

// NMH headers
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
#include "RooDataSet.h"

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

/** Namespace that collects functions and variables for the `Fitter` application */
namespace FITTER {

  // functions
  void InitRespsAndSels();
  void FillRespsAndSels(TString simdata_file, TString expdata_file, Bool_t refill_response, FitUtil *futil);
  void FillSelections();
  void Cleanup();
  
  // binning variables
  static const Int_t    fENB    =  40;
  static const Double_t fEmin   =   1;
  static const Double_t fEmax   = 100;
  static const Int_t    fCtNB   =  80;
  static const Double_t fCtmin  =  -1;
  static const Double_t fCtmax  =   1;
  static const Int_t    fByNB   =   1;
  static const Double_t fBymin  =   0;
  static const Double_t fBymax  =   1;

  // member selections and responses
  DetResponse    *fTRres;
  DetResponse    *fSHres;
  EventSelection *fTRsel;
  EventSelection *fSHsel;
  RooDataSet     *fTRset;
  RooDataSet     *fSHset;

};

/** Main function of the program */
int main(const int argc, const char **argv) {

  using namespace FITTER;

  //----------------------------------------------------------
  // parse command line arguments with Jpp
  //----------------------------------------------------------
  
  string        simdata_file;
  string        expdata_file;
  bool          refill_response;
  bool          unbinned;
  TString       meff_file;

  try {

    JParser<> zap("Program to fit the experimental data with standard track-shower separation.");

    zap['s'] = make_field(simdata_file, "File with all summary data") =
      (string)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root";

    zap['e'] = make_field(expdata_file, "File with experimental data sample");

    zap['r'] = make_field(refill_response, "Flag to request re-filling of the detector responses");

    zap['u'] = make_field(unbinned, "Flag to request unbinned fitting");
    
    zap['M'] = make_field(meff_file, "Effective mass file created by using `EffMass` class") = 
      (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";

    if ( zap.read(argc, argv)!= 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  //----------------------------------------------------------
  // call methods to initialise and fill responses
  //----------------------------------------------------------

  InitRespsAndSels();
  
  FitUtil *fitutil = new FitUtil(3, fTRres->GetHist3D(), 1, 100, -1, 0, 0, 1, meff_file);
  fTRset = new RooDataSet("tracks_unbinned","tracks_unbinned"  , fitutil->GetObs() );
  fSHset = new RooDataSet("showers_unbinned","showers_unbinned", fitutil->GetObs() );
  
  FillRespsAndSels(simdata_file, expdata_file, refill_response, fitutil);
  
  //----------------------------------------------------------
  // set up the PDFs for fitting
  //----------------------------------------------------------

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, fTRres);  
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, fSHres);

  //----------------------------------------------------------
  // set up data for simultaneous fitting for binned data
  //----------------------------------------------------------
  std::map<string, TH1* > hist_map = { { (string)fTRsel->Get_SelName() , fTRsel->Get_h_E_costh_by() }, 
				       { (string)fSHsel->Get_SelName() , fSHsel->Get_h_E_costh_by() } };
  
  // create categories
  RooCategory categs("categs","data categories");
  categs.defineType( fTRsel->Get_SelName() );
  categs.defineType( fSHsel->Get_SelName() );
  
  // create combined data set for the histogram
  RooDataHist combData("combData", "combined data", fitutil->GetObs(), categs, hist_map);
  
  // create simultaneous pdf
  RooSimultaneous simPdf("simPdf", "simultaneous Pdf", categs);
  simPdf.addPdf(pdf_tracks , fTRsel->Get_SelName() );
  simPdf.addPdf(pdf_showers, fSHsel->Get_SelName() );

  //----------------------------------------------------------
  // set up data for simultaneous fitting for un-binned data
  //----------------------------------------------------------

  RooCategory categsUB("categsUB", "data categories, unbinned");
  categs.defineType( fTRset->GetName() );
  categs.defineType( fSHset->GetName() );
  
  RooDataSet combDataUB("combDataUB","combined data, unbinned", fitutil->GetObs(),
			Index(categsUB), Import(fTRset->GetName(), *fTRset), Import(fSHset->GetName(), *fSHset) );

  RooSimultaneous simPdfUB("simPdf","simultaneous pdf, unbinned", categsUB);
  simPdfUB.addPdf( pdf_tracks , fTRset->GetName() );
  simPdfUB.addPdf( pdf_showers, fSHset->GetName() );
  
  //----------------------------------------------------------
  // do the fitting, either binned or un-binned
  //----------------------------------------------------------
    
  cout << "NOTICE Fitter started fitting" << endl;
  
  TStopwatch timer;

  RooAbsReal *nll = NULL;

  if (unbinned) {
    nll = simPdfUB.createNLL(combDataUB);
  }
  else {
    nll = simPdf.createNLL(combData);
  }
 
  RooMinimizer m(*nll);
  const char* minimizer = "Minuit";
  const char* minalg    = "minuit";
  m.setMinimizerType(minimizer); // select minimizer
  m.optimizeConst(kFALSE);       // turn off optimisation
  m.setOffsetting(kTRUE);        // offset the likelihood, helps with large LLH values
  m.setEps(1e-5);                // increase MIGRAD precision
  m.hesse();                     // do initial error evaluation
  m.minimize(minimizer, minalg);
  m.hesse();
  RooFitResult *fitres = m.save();

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

void FITTER::InitRespsAndSels() {

  using namespace FITTER;

  //----------------------------------------------------------
  //initialise event selection and response for tracks
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
  //initialise event selection and response for showers
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

void FITTER::FillRespsAndSels(TString simdata_file, TString expdata_file,
			      Bool_t refill_response, FitUtil *fitutil) {

  using namespace FITTER;

  if ( !fTRres || !fSHres || !fTRsel || !fSHsel ) {
    throw std::logic_error("ERROR! FITTER::FillRespsAndSels() DetResponse's and EventSelection's are not yet initialized!");
  }
  
  //----------------------------------------------------------
  // fill the responses
  //----------------------------------------------------------

  TString track_resp_name  = NMHUtils::Getcwd() + "/rootfiles/track_response.root";
  TString shower_resp_name = NMHUtils::Getcwd() + "/rootfiles/shower_response.root";

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

  RooRealVar *e  = fitutil->GetVar("E_reco");
  RooRealVar *ct = fitutil->GetVar("ct_reco");
  RooRealVar *by = fitutil->GetVar("by_reco");
  
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

    SummaryEvent *evt = sp.GetEvt(i);
    fTRsel->Fill(evt);
    fSHsel->Fill(evt);

    // fill the sets
    if ( fTRsel->PassesCuts(evt) ) {
      e ->setVal(  fTRsel->fEnergy );
      ct->setVal( -fTRsel->fDir.z() );
      by->setVal(  fTRsel->fBy );
      fTRset->add( fitutil->GetObs() );
    }

    if ( fSHsel->PassesCuts(evt) ) {
      e ->setVal(  fSHsel->fEnergy );
      ct->setVal( -fSHsel->fDir.z() );
      by->setVal(  fSHsel->fBy );
      fSHset->add( fitutil->GetObs() );
    }

  }

  cout << "NOTICE FITTER::FillRespsAndSels() Selections filled" << endl;
}

//**********************************************************************************

void FITTER::Cleanup() {

  using namespace FITTER;

  if (fTRres) delete fTRres;
  if (fSHres) delete fSHres;
  if (fTRsel) delete fTRsel;
  if (fSHsel) delete fSHsel;
  // if (fTRset) delete fTRset;
  // if (fSHsel) delete fSHsel;
  
}

