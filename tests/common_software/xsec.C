#include "NuXsec.h"
#include <iostream>
#include "TRandom.h"
#include "TMath.h"

using namespace std;

/** 
    This applications tests some base functionality of the NuXsec class.

    \return 0 if test worked, 1 otherwise.
*/
int main(const int argc, const char** argv) {

  //---------------------------------------------------------------------
  // test that a random binning in bjorken-y cannot be requested
  //---------------------------------------------------------------------

  Bool_t error_caught = kFALSE;

  try {
    NuXsec(15);
  }
  catch (const std::invalid_argument& ia) {
    error_caught = true;
  }

  if (!error_caught) {
    cout << "NOTICE xsec: test failed at bjorken-y binning" << endl;
    return 1;
  }

  //---------------------------------------------------------------------
  // test consistency between different gSeaGen versions
  //---------------------------------------------------------------------

  TString xsec_v4r1 = (TString)getenv("MONADIR") + "/data/xsec_data/xsec_gSeaGen_v4r1.root";
  TString by_v4r1   = (TString)getenv("MONADIR") + "/data/xsec_data/by_dists_gSeaGen_v4r1.root";
  NuXsec v4r1(2, xsec_v4r1, by_v4r1);

  TString xsec_v5r1 = (TString)getenv("MONADIR") + "/data/xsec_data/xsec_gSeaGen_v5r1.root";
  TString by_v5r1   = (TString)getenv("MONADIR") + "/data/xsec_data/by_dists_gSeaGen_v5r1.root";
  NuXsec v5r1(2, xsec_v5r1, by_v5r1);

  TRandom rand(0);
 
  Double_t xsec_tolerance = 1;
  Double_t by_tolerance   = 5;

  for (Int_t N = 0; N < 10000; N++) {

    Int_t flav  = rand.Integer(3);
    Int_t iscc  = rand.Integer(2);
    Int_t isnb  = rand.Integer(2);
    Double_t E  = rand.Uniform(7, 50); // choose regions with "reasonable" stats for bjornen-y
    Double_t by = rand.Uniform(0,1);
    
    Double_t Xv4 = v4r1.GetXsec(flav, iscc, isnb, E);
    Double_t Xv5 = v5r1.GetXsec(flav, iscc, isnb, E);
    Double_t Xd  = TMath::Abs( (Xv4 - Xv5)/Xv4*100 );

    Double_t Bv4 = v4r1.GetBYfrac(flav, iscc, isnb, E, by);
    Double_t Bv5 = v5r1.GetBYfrac(flav, iscc, isnb, E, by);
    Double_t Bd  = TMath::Abs( (Bv4 - Bv5)/Bv4*100 );

    // xsec tolerance is always tested
    if ( Xd > xsec_tolerance ) {
      cout << "NOTICE xsec: test failed at xsec tolerance, neutrino (" << flav << " " << iscc << " " << isnb << "), xsec's and diff " << Xv4 << "\t" << Xv5 << "\t" << Xd << endl;
      return 1;
    }

    // bjorken-y tolerance is tested only in bins where bjorken-y fraction is larger than 0.01; this is
    // to avoid comparisons based on low statistics
    if ( Bd > by_tolerance && Bv4 > 0.01 ) {
      cout << "NOTICE xsec: test failed at bjorken-y tolerance, neutrino (" << flav << " " << iscc << " " << isnb << "), bjorken-y's and diff " << Bv4 << "\t" << Bv5 << "\t" << Bd << endl;
      return 1;
    }
    
  }

  //---------------------------------------------------------------------
  // test that bjorken-y fraction fails outside the MC range
  //---------------------------------------------------------------------

  error_caught = kFALSE;

  try {
    v5r1.GetBYfrac(1,1,0, 120, 0.49);
  }
  catch (const std::invalid_argument& ia) {
    error_caught = kTRUE;
  }

  if (!error_caught) {
    cout << "NOTICE xsec: test failed at catching an error for bjorken-y calculation outside MC E range" << endl;
    return 1;
  }

  //---------------------------------------------------------------------
  // test that bjorken-y fraction is 1 in any range if 1 bin is used, xsec available up to 1000
  //---------------------------------------------------------------------

  NuXsec xsec(1, xsec_v5r1, by_v5r1);

  for (Int_t i = 0; i < 1000; i++) {

    Int_t flav  = rand.Integer(3);
    Int_t iscc  = rand.Integer(2);
    Int_t isnb  = rand.Integer(2);
    Double_t E  = rand.Uniform(1, 1000);
    Double_t by = rand.Uniform(0,1);

    try {
      xsec.GetXsec(flav, iscc, isnb, E);
    }
    catch (const std::invalid_argument& ia) {
      cout << "NOTICE xsec: test failed, caught an error at calculating the cross-section at E = " << E << endl;
      return 1;
    }

    try {
      Double_t byf = xsec.GetBYfrac(flav, iscc, isnb, E, by);

      if (  byf != 1.0 ) {
	cout << "NOTICE xsec: test failed, bjorken-y equal to " << byf << ", expected 1.0" << endl;
	return 1;
      }

    }
    catch (const std::invalid_argument& ia) {
      cout << "NOTICE xsec: test failed, caught an error at calculating the bjorken-y dist at E = " << E << endl;
      return 1;
    }

  }
  
  return 0;

}
