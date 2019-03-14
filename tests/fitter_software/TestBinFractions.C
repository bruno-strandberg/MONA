
// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "FitUtilWsyst.h"
#include "FitPDF.h"

#include "TH1.h"

// cpp headers
#include <iostream>

using namespace std;
using namespace RooFit;

//******************************************************************************************
// create a class that inherits from `FitUtilWsyst` to access the protected functions for testing
class FU : protected FitUtilWsyst {

public:

  Int_t test_result;

  // perform the tests in the constructor
  FU(Double_t op_time, TH3 *h_template,
     Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax,
     Double_t bymin, Double_t bymax, TString meff_file) :
    FitUtilWsyst(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {

    test_result = 0;    
    TH1D htest("htest","htest", 10, 0, 10);

    // check that wrong input throws an error
    //------------------------------------------------------
    Int_t errors = 0;
    
    try { GetBinFractions( 5, 4, htest.GetXaxis() ); }
    catch (const std::invalid_argument& ia) { errors++; }

    try { GetBinFractions( -1, 4, htest.GetXaxis() ); }
    catch (const std::invalid_argument& ia) { errors++; }

    try { GetBinFractions( 1, 14, htest.GetXaxis() ); }
    catch (const std::invalid_argument& ia) { errors++; }

    if (errors != 3) {
      cout << "NOTICE TestBinFractions test failed at input checks" << endl;
      test_result = 1;
      return;
    }

    // check that correct bin range is found
    //------------------------------------------------------
    auto bins = GetBinFractions( 0.5, 1.5, htest.GetXaxis() );

    if ( bins.size() != 2 || bins[0].first != 1 || bins[0].second != 0.5 ||
	 bins[1].first != 2 || bins[0].second != 0.5 ) {
      cout << "NOTICE TestBinFractions test failed at range with 2 bins" << endl;
      test_result = 1;
      return;      
    }

    bins = GetBinFractions(0.7, 2.5, htest.GetXaxis());

    if (bins.size() != 3   ||
	bins[0].first != 1 || TMath::Abs(bins[0].second - 0.3) > 1e-10 ||
	bins[1].first != 2 || TMath::Abs(bins[1].second - 1. ) > 1e-10 ||
	bins[2].first != 3 || TMath::Abs(bins[2].second - 0.5) > 1e-10 ) {

      cout << "NOTICE TestBinFractions test failed at range with 3 bins" << endl;
      test_result = 1;
      return;
      
    }
    
  };

  virtual ~FU() {};
  
};


//******************************************************************************************

/** Application to test the flux normalisation */
int main(const int argc, const char **argv) {
    
  // create a response
  DetResponse resp(DetResponse::track, "trkresp" , 40, 1, 100, 40, -1, 1, 1, 0, 1);
  resp.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  resp.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );

  FU fu(3, resp.GetHist3D(), 2, 75, -1, 0, 0, 1, EffMass::DUMMYFILE);

  if (fu.test_result == 0) cout << "NOTICE TestBinFractions test passed" << endl;
  return fu.test_result;
  
}
