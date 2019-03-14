
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

/** This application tests that the energy scale systematic affects `FitUtilWsyst::RecoEvts` as expected.*/
int main(const int argc, const char **argv) {
  
  //-------------------------------------------------------------
  // fill response with pseudo-data
  //-------------------------------------------------------------
  
  // create a response
  DetResponse resp(DetResponse::track, "trkresp" , 40, 1, 100, 40, -1, 1, 1, 0, 1);
  resp.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  resp.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );

  SummaryEvent evt;
  for (Int_t i = 0; i < 1e6; i++) {
    evt.FillPseudoData();
    resp.Fill(&evt);
  }

  cout << "NOTICE TestEscale response ready" << endl;
  
  //-------------------------------------------------------------
  // init fitutil and fitpdf
  //-------------------------------------------------------------

  Double_t op_time =  3.;
  Double_t emin    =  2.;
  Double_t emax    =  75.;
  Double_t ctmin   = -1;
  Double_t ctmax   =  0;
  Double_t bymin   =  0;
  Double_t bymax   =  1;

  //-------------------------------------------------------------
  // init the FitUtil and FitUtilWsyst and the corresponding pdf classes for access to proxymap_t
  //-----------------------------------------------------------------------------
  FitUtil *futil = new FitUtil(op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);
  FitUtilWsyst *futilws = new FitUtilWsyst(op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);

  FitPDF pdf  ("pdf"   ,"pdf"  , futil            , &resp);
  FitPDF pdfws("pdfws" ,"pdfws", (FitUtil*)futilws, &resp);

  //-------------------------------------------------------------
  // run test trials
  //-----------------------------------------------------------------------------
  Int_t trials = 1000;
  TRandom3 rand(0);

  for (Int_t N = 0; N < trials; N++) {

    // select random energy, costheta and bjorken-y
    Double_t ereco  = rand.Uniform( futil->GetVar("E_reco")->getMin(), futil->GetVar("E_reco")->getMax() );
    Double_t ctreco = rand.Uniform(-1,0);
    Double_t byreco = rand.Uniform(0,1);

    // test values when the e-scale is 0; have to be identical
    //------------------------------------------------------------------------
    
    ((FitUtil*)futilws)->GetVar("E_scale")->setVal( 0. );

    Double_t recoevts   = futil->RecoEvts(ereco, ctreco, byreco, &resp, pdf.GetProxyMap()).first;
    Double_t recoevts_s = ((FitUtil*)futilws)->RecoEvts(ereco, ctreco, byreco, &resp, pdfws.GetProxyMap()).first;

    if ( recoevts != recoevts_s ) {
      cout << "NOTICE TestEscale failed, reco events differ when no energy scale applied" << endl;
      return 1;
    }

    // test values when the e-scale is non-zero; this here re-produces the calculation
    // that is performed in `FitUtilWsyst::RecoEvts`
    //------------------------------------------------------------------------

    Double_t scale = rand.Uniform( -0.2, 0.2 );
    ((FitUtil*)futilws)->GetVar("E_scale")->setVal( scale );

    TAxis *EX = resp.GetHist3D()->GetXaxis();
    
    Int_t ebin    = EX->FindBin( ereco );
    Double_t elo  = EX->GetBinLowEdge( ebin ) * (1. + scale);
    Double_t ehi  = EX->GetBinUpEdge( ebin ) * (1. + scale);

    Int_t ebin_lo = EX->FindBin( elo );
    Int_t ebin_hi = EX->FindBin( ehi );

    if ( (ebin_hi - ebin_lo) > 1 ) {
      cout << "NOTICE TestEscale failed, bin range extraction not correct" << endl;
      return 1;
    }

    Double_t binW_lo = ( EX->GetBinUpEdge( ebin_lo ) - elo )/EX->GetBinWidth( ebin_lo );
    Double_t binW_hi = ( ehi - EX->GetBinLowEdge( ebin_hi ) )/EX->GetBinWidth( ebin_hi );

    recoevts = ( futil->RecoEvts( elo, ctreco, byreco, &resp, pdf.GetProxyMap() ).first * binW_lo +
		 futil->RecoEvts( ehi, ctreco, byreco, &resp, pdf.GetProxyMap() ).first * binW_hi );

    recoevts_s = ((FitUtil*)futilws)->RecoEvts( ereco, ctreco, byreco, &resp, pdfws.GetProxyMap() ).first;

    if ( recoevts != recoevts_s ) {
      cout << "NOTICE TestEscale failed, cannot reproduce the energy scale calculation." << endl;
      return 1;
    }

  }

  cout << "NOTICE TestEscale test passed" << endl;
  return 0;
  
}
