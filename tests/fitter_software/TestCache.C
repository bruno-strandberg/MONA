
// ROOT headers
#include "TRandom.h"

// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "EffMass.h"
#include "FitPDF.h"
#include "AtmFlux.h"
#include "NuXsec.h"
#include "PMNS_Fast.h"
#include "PremModel.h"

// cpp headers
#include <iostream>
using namespace std;
using namespace OscProb;

/** This application tests the cache functionality in the `FitUtil` class. */
int main(const int argc, const char **argv) {
  
  Int_t testpoints = 100;
  
  //-------------------------------------------------------------
  // init the necessary classes
  //-------------------------------------------------------------
  
  // create an empty response
  DetResponse resp(DetResponse::track, "dummy_trk", 40, 1, 100, 40, -1, 1, 1, 0, 2);
  TH3D *hb = resp.GetHist3D();
  
  // init fitutil with dummy effective mass
  FitUtil futil(3, hb, 3, 80, -1, 0, 0, 1, EffMass::DUMMYFILE);
  futil.SetNOlims(); 
  
  // init a FitPDF to get access to the parameter map with the empty response - as long as we're not
  // calculating true events the response can be empty
  FitPDF pdf("pdf","pdf", &futil, &resp);

  // get a reference to the proxies for calling cache functions
  proxymap_t proxies = pdf.GetProxyMap();

  // get pointers to the flux, osc and xsec calculators for testing
  AtmFlux   *flux = futil.GetFluxCalculator();
  PMNS_Fast *osc  = futil.GetOscCalculator();
  PremModel *prem = futil.GetEarthModel();
  NuXsec    *xsec = futil.GetXsecCalculator();
  EffMass   *meff = futil.GetEffMassCalculator();

  //seconds in a tropical year
  double sec_per_y   = 365.2421897 * 24 * 60 * 60; 
  
  TRandom3 rand(0);

  //-------------------------------------------------------------
  // perform test trials
  //-------------------------------------------------------------
  
  for (Int_t trial = 0; trial < testpoints; trial++) {

    Int_t flav_in = rand.Integer(2); // flav elec(0), muon(1)

    //-------------------------------------------------------------
    // create true bin object to test the caches
    //-------------------------------------------------------------
    Int_t flav = rand.Integer(3); // flav elec(0), muon(1) or tau(2)
    Int_t iscc = rand.Integer(2); // false(0) or true(1)
    Int_t isnb = rand.Integer(2); // false(0) or true(1)
    Int_t ebin  = rand.Integer( hb->GetXaxis()->GetNbins() ) + 1;
    Int_t ctbin = rand.Integer( hb->GetYaxis()->GetNbins() ) + 1;
    Int_t bybin = rand.Integer( hb->GetZaxis()->GetNbins() ) + 1;

    TrueB tb( flav, iscc, isnb, ebin, ctbin, bybin );

    //-------------------------------------------------------------
    // get bin center values and bin widths
    //-------------------------------------------------------------
    Double_t E    = hb->GetXaxis()->GetBinCenter( ebin );
    Double_t ct   = hb->GetYaxis()->GetBinCenter( ctbin );
    Double_t by   = hb->GetZaxis()->GetBinCenter( bybin );

    Double_t e_w  = hb->GetXaxis()->GetBinWidth( ebin );
    Double_t ct_w = hb->GetYaxis()->GetBinWidth( ctbin );

    //-------------------------------------------------------------
    // calculate flux difference
    //-------------------------------------------------------------
    Double_t flux_conv = e_w * ct_w * sec_per_y * 3;
    Double_t flux_diff = futil.GetCachedFlux(flav, isnb, ebin, ctbin) - flux->Flux_dE_dcosz(flav, isnb, E, ct) * flux_conv;

    //-------------------------------------------------------------
    // calcuate osc difference
    //-------------------------------------------------------------

    futil.GetVar("SinsqTh12")->randomize();
    futil.GetVar("SinsqTh13")->randomize();
    futil.GetVar("SinsqTh23")->randomize();
	    
    Double_t osc_cache = futil.GetCachedOsc(flav_in, tb, proxies); 
    
    // get the parameter values from the proxies
    Double_t SinsqTh12 = *( proxies.at( "SinsqTh12" ) );
    Double_t SinsqTh13 = *( proxies.at( "SinsqTh13" ) );
    Double_t SinsqTh23 = *( proxies.at( "SinsqTh23" ) );
    Double_t Dcp       = *( proxies.at( "dcp" ) );
    Double_t Dm21      = *( proxies.at( "Dm21" ) );
    Double_t Dm31      = *( proxies.at( "Dm31" ) );

    osc->SetAngle(1, 2, TMath::ASin( TMath::Sqrt( SinsqTh12 ) ) );
    osc->SetAngle(1, 3, TMath::ASin( TMath::Sqrt( SinsqTh13 ) ) );
    osc->SetAngle(2, 3, TMath::ASin( TMath::Sqrt( SinsqTh23 ) ) );
    osc->SetDelta(1, 3, Dcp * TMath::Pi()  );
    osc->SetDm(2, Dm21);
    osc->SetDm(3, Dm31);

    // set oscillation path and is-nu-bar flag
    prem->FillPath( ct );
    osc->SetPath ( prem->GetNuPath() );
    osc->SetIsNuBar( isnb );

    Double_t osc_calc = osc->Prob(flav_in, tb.fFlav, E);
    
    Double_t osc_diff = osc_cache - osc_calc;

    //-------------------------------------------------------------
    // calcuate xsec difference
    //-------------------------------------------------------------

    Double_t xsec_diff = futil.GetCachedXsec(tb)*1e42 - xsec->GetXsec(tb.fFlav, tb.fIsCC, tb.fIsNB, E)*1e42;

    //-------------------------------------------------------------
    // calcuate by fraction difference difference
    //-------------------------------------------------------------

    Double_t byfrac_diff = futil.GetCachedBYfrac(tb) - xsec->GetBYfrac(tb.fFlav, tb.fIsCC, tb.fIsNB, E, by);

    //-------------------------------------------------------------
    // calcuate effective mass difference
    //-------------------------------------------------------------

    Double_t meff_diff = futil.GetCachedMeff(tb) - meff->GetMeff(tb.fFlav, tb.fIsCC, tb.fIsNB, E, ct, by);

    //-------------------------------------------------------------
    // test
    //-------------------------------------------------------------
    
    if ( flux_diff != 0. || osc_diff != 0. || xsec_diff != 0. || byfrac_diff != 0. || meff_diff != 0.) {
      cout << tb << endl;
      cout << "Flux   : " << flux_diff << endl;
      cout << "Osc    : " << osc_diff << endl;
      cout << "Xsec   : " << xsec_diff << endl;
      cout << "BY frac: " << byfrac_diff << endl;
      cout << "Meff   : " << meff_diff << endl;
      cout << "NOTICE TestCache failed" << endl;
      return 1;
    }
        
  } // end loop over trials

  cout << "NOTICE TestCache passed" << endl;
  return 0;
  
}
