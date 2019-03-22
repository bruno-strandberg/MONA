
// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "FitPDF.h"

// ROOT headers
#include "TRandom.h"

// RooFit headers
#include "RooDataHist.h"
#include "RooFitResult.h"

// cpp headers
#include <iostream>

using namespace std;
using namespace RooFit;

//******************************************************************************************

/** This application tests that the overall normalisation of the data does not affect the fit*/
int main(const int argc, const char **argv) {
  
  //-------------------------------------------------------------
  // fill response with pseudo-data
  //-------------------------------------------------------------
  
  // create responses
  DetResponse trkR(DetResponse::track, "trkresp" , 40, 1, 100, 40, -1, 1, 1, 0, 1);
  trkR.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  trkR.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );

  SummaryEvent evt;
  for (Int_t i = 0; i < 1e6; i++) {
    evt.FillPseudoData();
    trkR.Fill(&evt);
  }

  cout << "NOTICE TestNorm response ready" << endl;
  
  //-------------------------------------------------------------
  // init fitutil and fitpdf's
  //-------------------------------------------------------------

  Double_t op_time =  3.;
  Double_t emin    =  2.;
  Double_t emax    =  75.;
  Double_t ctmin   = -1;
  Double_t ctmax   =  0;
  Double_t bymin   =  0;
  Double_t bymax   =  1;

  FitUtil *futil = new FitUtil(op_time, trkR.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);

  FitPDF trkpdf("trkpdf","trkpdf", futil, &trkR);
    
  //-------------------------------------------------------------
  // fix all oscillation parameters except dm31; set dm31 and normalisation
  // constants to random values and create data for fitting
  //-------------------------------------------------------------

  TRandom3 rand(0);
  
  // set parameters to NO central values and fix
  futil->SetNOcentvals();
  futil->GetVar("SinsqTh12")->setConstant(kTRUE);
  futil->GetVar("SinsqTh13")->setConstant(kTRUE);
  futil->GetVar("SinsqTh23")->setConstant(kTRUE);
  futil->GetVar("dcp")      ->setConstant(kTRUE);
  futil->GetVar("Dm21")     ->setConstant(kTRUE);

  // randomize dm31, so that there is something to fit
  Double_t dm31_start = futil->GetVar("Dm31")->getVal();
  futil->GetVar("Dm31")->setVal( rand.Uniform(2.5, 2.7)*1e-3 );
  Double_t dm31_true  = futil->GetVar("Dm31")->getVal();
  Double_t norm       = rand.Uniform(1.3, 2.); // random normalisation for the data
  
  // create data with and without scaling
  TH3D* T = trkpdf.GetExpValHist();
  T->SetNameTitle("T","T");
  TH3D* N = (TH3D*)T->Clone("N");
  N->SetNameTitle("N", "N");
  N->Scale( norm );

  //-------------------------------------------------------------
  // fit the normalised and unnormalised data
  //-------------------------------------------------------------

  RooDataHist rf_T ("rf_T", "rf_T", futil->GetObs(), Import( *T ) );
  RooDataHist rf_N ("rf_N", "rf_N", futil->GetObs(), Import( *N ) );

  futil->GetVar( "Dm31" )->setVal( dm31_start );
  RooFitResult *res_T   = trkpdf.fitTo( rf_T, Save(), SumW2Error(kFALSE) );
  
  futil->GetVar( "Dm31" )->setVal( dm31_start );
  RooFitResult *res_N = trkpdf.fitTo( rf_N, Save(), SumW2Error(kFALSE) );

  // check that input data are different by norm
  if ( TMath::Abs( rf_T.sumEntries() - rf_N.sumEntries()/norm ) > 1e-5 ) {
    cout << "NOTICE TestNorm failed, normalisation not applied" << endl;
    return 1;
  }
  
  // check that start values are the same
  if ( ((RooRealVar*)res_T->floatParsInit().find("Dm31"))->getVal() !=
       ((RooRealVar*)res_N->floatParsInit().find("Dm31"))->getVal() ) {
    throw std::logic_error("ERROR! TestNorm starting points for the fit should be the same");
  }
  
  // check that end values are the same 
  if ( TMath::Abs( ((RooRealVar*)res_T->floatParsFinal().find("Dm31"))->getVal() -
		   ((RooRealVar*)res_N->floatParsFinal().find("Dm31"))->getVal() ) > 1e-5 ) {
    cout << "NOTICE TestNorm failed, different fit results" << endl;
    return 1;
  }
 
  // check that both close to truth
  if ( TMath::Abs( ((RooRealVar*)res_T->floatParsFinal().find("Dm31"))->getVal() - dm31_true ) > 1e-5 ) {
    cout << "NOTICE TestNorm failed, result away from truth" << endl;
    return 1;
  }
  
  // check that likelihood is scaled
  Double_t nll_T = res_T->minNll();
  Double_t nll_N = res_N->minNll()/norm;
  if ( TMath::Abs( nll_T - nll_N ) > 1e-5 ) {
    cout << "NOTICE TestNorm failed, likelihoods differ" << endl;
    return 1;
  }

  cout << "NOTICE TestNorm passed " << endl;
  return 0;
  
}
