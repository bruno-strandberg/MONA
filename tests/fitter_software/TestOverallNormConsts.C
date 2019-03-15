
// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "FitUtilWsyst.h"
#include "FitPDF.h"

#include "TRandom.h"

// cpp headers
#include <iostream>

using namespace std;
using namespace RooFit;

//******************************************************************************************

/** This application tests the over-all normalisation constants*/
int main(const int argc, const char **argv) {
  
  //-------------------------------------------------------------
  // fill response with pseudo-data
  //-------------------------------------------------------------
  
  // create responses
  DetResponse trkR(DetResponse::track, "trkresp" , 40, 1, 100, 40, -1, 1, 1, 0, 1);
  trkR.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  trkR.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );

  DetResponse shwR(DetResponse::shower, "shwresp" , 40, 1, 100, 40, -1, 1, 1, 0, 1);
  shwR.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   , 0.5, true );
  shwR.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), 0.6, true );

  SummaryEvent evt;
  for (Int_t i = 0; i < 1e6; i++) {
    evt.FillPseudoData();
    trkR.Fill(&evt);
    shwR.Fill(&evt);
  }

  cout << "NOTICE TestOverallNormConsts response ready" << endl;
  
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

  FitPDF trkpdf("trkpdf","trkpdf", futil, &trkR, kFALSE);
  FitPDF shwpdf("shwpdf","shwpdf", futil, &shwR, kFALSE);
    
  // test 1 -  see that an error is thrown if I init another pdf with an already existing name
  //--------------------------------------------------------------------------------
  Bool_t caught = kFALSE;
  try {
    FitPDF("trkpdf","trkpdf", futil, &trkR, kFALSE);
  }
  catch (const std::invalid_argument &ia) {
    caught = kTRUE;
  }

  if ( !caught ) {
    cout << "NOTICE TestOverallNormConsts failed, `FitPDF` does not generate an error when name is already in use" << endl;
    return 1;
  }

  // test 2 -  check that a normalisation constant exists for both pdf's
  //--------------------------------------------------------------------------------
  RooArgList normPars( futil->GetNormSet() );

  if ( normPars.getSize() != 2 ) {
    cout << "NOTICE TestOverallNormConsts failed, wrong number of normalisation parameters" << endl;
    return 1;
  }

  try {
    futil->GetVar( trkpdf.GetNormName() );
    futil->GetVar( shwpdf.GetNormName() );
  }
  catch (const std::invalid_argument &ia) {
    cout << "NOTICE TestOverallNormConsts failed, cannot access normalisation parameters through FitUtil::GetVar() " << endl;
    return 1;
  }
  
  // loop over trials
  //--------------------------------------------------------------------------------
  
  Int_t ntrial = 10;
  TRandom3 rand(0);
  
  for (Int_t i = 0; i < ntrial; i++) {

    // select some reco variables
    Double_t E_reco  = rand.Uniform( futil->GetVar("E_reco")->getMin(), futil->GetVar("E_reco")->getMax() );
    Double_t Ct_reco = rand.Uniform(-1, 0);
    Double_t By_reco = rand.Uniform(0, 1);

    // set normalisation to 1 and get values
    futil->GetVar( trkpdf.GetNormName() )->setVal( 1. );
    futil->GetVar( shwpdf.GetNormName() )->setVal( 1. );
    Double_t RE_trk_1 = futil->RecoEvts(E_reco, Ct_reco, By_reco, &trkR, trkpdf.GetProxyMap(), trkpdf.GetNorm()).first;
    Double_t RE_shw_1 = futil->RecoEvts(E_reco, Ct_reco, By_reco, &shwR, shwpdf.GetProxyMap(), shwpdf.GetNorm()).first;
    auto EXP_trk_1 = trkpdf.GetExpValHist();
    auto EXP_shw_1 = shwpdf.GetExpValHist();
    Double_t EXP_trkI_1 = EXP_trk_1->Integral();
    Double_t EXP_shwI_1 = EXP_shw_1->Integral();
    delete EXP_trk_1;
    delete EXP_shw_1;
    
    // set some random normalisation
    Double_t trkNorm = rand.Uniform(0.5, 1.5);
    Double_t shwNorm = rand.Uniform(0.5, 1.5);
    futil->GetVar( trkpdf.GetNormName() )->setVal( trkNorm );
    futil->GetVar( shwpdf.GetNormName() )->setVal( shwNorm );
    Double_t RE_trk_N = futil->RecoEvts(E_reco, Ct_reco, By_reco, &trkR, trkpdf.GetProxyMap(), trkpdf.GetNorm()).first;
    Double_t RE_shw_N = futil->RecoEvts(E_reco, Ct_reco, By_reco, &shwR, shwpdf.GetProxyMap(), shwpdf.GetNorm()).first;
    auto EXP_trk_N = trkpdf.GetExpValHist();
    auto EXP_shw_N = shwpdf.GetExpValHist();
    Double_t EXP_trkI_N = EXP_trk_N->Integral();
    Double_t EXP_shwI_N = EXP_shw_N->Integral();    
    delete EXP_trk_N;
    delete EXP_shw_N;

    // test 3 -  check that the return of `RecoEvts` behaves as expected
    //--------------------------------------------------------------------------------

    if ( RE_trk_1 * trkNorm != RE_trk_N ) {
      cout << "NOTICE TestOverallNormConsts test fails at track normalisation" << endl;
      return 1;
    }

    if ( RE_shw_1 * shwNorm != RE_shw_N ) {
      cout << "NOTICE TestOverallNormConsts test fails at shower normalisation" << endl;
      return 1;
    }

    // test 4 -  check that the the integral is properly normalised
    //--------------------------------------------------------------------------------

    if ( TMath::Abs( EXP_trkI_1 * trkNorm - EXP_trkI_N ) > 1e-5 ) {
      cout << "NOTICE TestOverallNormConsts test fails at track integral normalisation" << endl;
      return 1;
    }

    if ( TMath::Abs( EXP_shwI_1 * shwNorm - EXP_shwI_N ) > 1e-5 ) {
      cout << "NOTICE TestOverallNormConsts test fails at shower integral normalisation" << endl;
      return 1;
    }
    
  }

  cout << "NOTICE TestOverallNormConsts test passed" << endl;
  return 0;
    
}
