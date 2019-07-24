
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
    
    futilws->GetVar("E_scale")->setVal( 1.0 );
    futilws->GetVar("E_scale")->setConstant( kTRUE );

    Double_t recoevts   = futil->RecoEvts(ereco, ctreco, byreco, &resp, pdf.GetProxyMap()).first;
    Double_t recoevts_s = futilws->RecoEvts(ereco, ctreco, byreco, &resp, pdfws.GetProxyMap()).first;

    if ( recoevts != recoevts_s ) {
      cout << "NOTICE TestEscale failed, reco events differ when no energy scale applied, values: " << recoevts << "\t" << recoevts_s << endl;
      return 1;
    }

    // test values when the e-scale is non-zero; this here re-produces the calculation
    // that is performed in `FitUtilWsyst::TrueEvts`
    //------------------------------------------------------------------------

    // set a random value to e-scale
    Double_t scale = rand.Uniform( 0.8, 1.2 );
    futilws->GetVar("E_scale")->setVal( scale );
    futilws->GetVar("E_scale")->setConstant( kFALSE );

    // get a random true bin from the response
    vector<TrueB> true_bins;
    Bool_t bin_found = kFALSE;
    
    while( !bin_found ) {

      true_bins = resp.GetBinWeights( rand.Uniform(5, 50), rand.Uniform(-1,0), rand.Uniform(0,1) );

      // continue if nothing contributes to this reco bin (should happen rarely anyway, but make certain)
      if (true_bins.size() == 0) continue;

      // ignore the first and last couple of bins, as FitUtilWsyst will also pretend they are not there
      if (true_bins[0].fE_true_bin > 5 && true_bins[0].fE_true_bin < 35 ) bin_found = kTRUE;
      
    }

    auto tb = true_bins[0];
    
    // find the bin fractions
    TAxis *EX      = resp.GetHist3D()->GetXaxis();
    Double_t elo   = EX->GetBinLowEdge( tb.fE_true_bin ) * scale;
    Double_t ehi   = EX->GetBinUpEdge ( tb.fE_true_bin ) * scale;
    
    Int_t ebin_lo = EX->FindBin( elo );
    Int_t ebin_hi = EX->FindBin( ehi );
    
    if ( (ebin_hi - ebin_lo) > 1 ) {
      cout << "NOTICE TestEscale failed, bin range extraction not correct" << endl;
      return 1;
    }

    Double_t binW_lo = ( EX->GetBinUpEdge( ebin_lo ) - elo )/EX->GetBinWidth( ebin_lo );
    Double_t binW_hi = ( ehi - EX->GetBinLowEdge( ebin_hi ) )/EX->GetBinWidth( ebin_hi );
    
    TrueB tb_lo( tb );
    TrueB tb_hi( tb );
    tb_lo.fE_true_bin = ebin_lo;
    tb_hi.fE_true_bin = ebin_hi;

    // compare the calculation
    auto trueevts   = futil->TrueEvts( tb_lo, pdf.GetProxyMap() ).first * binW_lo + futil->TrueEvts( tb_hi, pdf.GetProxyMap() ).first * binW_hi;
    auto trueevts_s = futilws->TrueEvts( tb, pdfws.GetProxyMap() ).first;
    
    if ( trueevts != trueevts_s ) {
      cout << "NOTICE TestEscale failed, cannot reproduce the energy scale calculation, values: " << trueevts << "\t" << trueevts_s << ", true bin: " << tb << endl;
      return 1;
    }

  }

  cout << "NOTICE TestEscale test passed" << endl;
  return 0;
  
}
