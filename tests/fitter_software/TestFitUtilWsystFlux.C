
// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "FitUtilWsyst.h"
#include "FitPDF.h"

// ROOT headers
#include "TRandom.h"

// cpp headers
#include <iostream>

using namespace std;
using namespace RooFit;


//******************************************************************************************

/** Application to test the flux normalisation */
int main(const int argc, const char **argv) {
  
  //-------------------------------------------------------------
  // fill response with pseudo-data
  //-------------------------------------------------------------
  
  // create a response; not dealing with `RecoEvts`, response can be empty
  DetResponse resp(DetResponse::track, "truthresp", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  resp.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  resp.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );
  
  //-------------------------------------------------------------
  // init fitutil and fitpdf
  //-------------------------------------------------------------

  Double_t op_time =  3.;
  Double_t emin    =  3.;
  Double_t emax    =  70.;
  Double_t ctmin   = -1;
  Double_t ctmax   =  0;
  Double_t bymin   =  0;
  Double_t bymax   =  1;
  
  // init the FitUtil and FitUtilWsyst and the corresponding pdf classes for access to proxymap_t
  //-----------------------------------------------------------------------------
  FitUtil *futil = new FitUtil(op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);
  FitUtilWsyst *futilws = new FitUtilWsyst(op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);

  FitPDF pdf  ("pdf"   ,"pdf"  , futil  , &resp);
  FitPDF pdfws("pdfws" ,"pdfws", futilws, &resp);

  // start testing
  //-----------------------------------------------------------------------------
  Int_t ntrials = 1000;
  TRandom3 rand(0);
  enum flavs { ELEC = 0, MUON, TAU };
  enum pols { NU = 0, NUB };
  
  /* numbers that have to be conserved
     1. tilt systematics alone need to preserve the flux for a specific neutrino type in full E, ct plane
     2. mu_amu skew needs to preserve muon+antimuon flux
     3. e_ae skew needs to preserve elec+antielec flux
     4. e_mu skew needs to preserve mu+amu+e+ae
     5. the overall flux over flavors and nu/nub needs to be conserved
   */

  for (Int_t i = 0; i < ntrials; i++) {

    // test 1: set random tilt parameters and check that the overall flux for a specific flavor is preserved
    // the skew parameters are set to 1, because they do not preserv the integral for e.g. muon-nu only
    //-------------------------------------------------------------------------------------------------------

    futilws->GetVar("E_tilt")      ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("ct_tilt")     ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("skew_mu_amu") ->setVal( 0. );
    futilws->GetVar("skew_e_ae")   ->setVal( 0. );
    futilws->GetVar("skew_mu_e")   ->setVal( 0. );

    UInt_t flav = rand.Integer(3); // select random flavor
    UInt_t isnb = rand.Integer(2); // select random  nu/nub
    
    Double_t fluxint   = 0.; // flux integral over (E,ct) without tilt systematics
    Double_t fluxintws = 0.; // flux integral over (E,ct) with tilt systematics

    // loop over bins
    for (Int_t ebin = 1; ebin <= resp.GetHist3D()->GetXaxis()->GetNbins(); ebin++) {
      for (Int_t ctbin = 1; ctbin <= resp.GetHist3D()->GetYaxis()->GetNbins(); ctbin++) {
	
	Double_t flux   = futil->GetCachedFlux( flav, isnb, ebin, ctbin );
	Double_t fluxws = futilws->GetFluxWsyst( flav, isnb, ebin, ctbin, pdfws.GetProxyMap() );

	fluxint   += flux;
	fluxintws += fluxws;

      }
    }

    if ( TMath::Abs(fluxint - fluxintws) > 1e-3 ) {
      cout << "NOTICE TestFitUtilWsystFlux test failed at flux integral preservation with tilt parameters" << endl;
      return 1;
    }

    // test 2: set tilt parameters to 0, set random skew parameters; check, in random bin, conservation of
    // i) mu+amu, ii) e+ae. If these are conserved for one bin, they are conserve the overall integral.
    // The mu_e skew is 0, because it does not preserve mu-e ratio (obviously)
    //-------------------------------------------------------------------------------------------------------

    futilws->GetVar("E_tilt")      ->setVal( 0. );
    futilws->GetVar("ct_tilt")     ->setVal( 0. );
    futilws->GetVar("skew_mu_amu") ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("skew_e_ae")   ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("skew_mu_e")   ->setVal( 0. );

    Int_t ebin  = 1 + rand.Integer( resp.GetHist3D()->GetXaxis()->GetNbins() );
    Int_t ctbin = 1 + rand.Integer( resp.GetHist3D()->GetYaxis()->GetNbins() );

    // no systematics
    Double_t flux_mu = futil->GetCachedFlux( MUON, NU, ebin, ctbin )
      + futil->GetCachedFlux( MUON, NUB, ebin, ctbin );
    
    Double_t flux_e = futil->GetCachedFlux( ELEC, NU, ebin, ctbin )
      + futil->GetCachedFlux( ELEC, NUB, ebin, ctbin );

    // with systematics
    Double_t fluxws_mu = futilws->GetFluxWsyst( MUON, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
      futilws->GetFluxWsyst( MUON, NUB, ebin, ctbin, pdfws.GetProxyMap() );

    Double_t fluxws_e = futilws->GetFluxWsyst( ELEC, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
      futilws->GetFluxWsyst( ELEC, NUB, ebin, ctbin, pdfws.GetProxyMap() );

    if ( TMath::Abs(flux_mu-fluxws_mu) > 1e-5 ) {
      cout << "NOTICE TestFitUtilWsystFlux test failed at muon flux conservation" << endl;
      return 1;
    }

    if ( TMath::Abs(flux_e-fluxws_e) > 1e-5 ) {
      cout << "NOTICE TestFitUtilWsystFlux test failed at elec flux conservation" << endl;
      return 1;
    }

    // check ratios
    Double_t R_mu_amu = futil->GetCachedFlux( MUON, NU, ebin, ctbin )/
      futil->GetCachedFlux( MUON, NUB, ebin, ctbin ) * (1. + futilws->GetVar("skew_mu_amu")->getVal());

    Double_t Rws_mu_amu = futilws->GetFluxWsyst( MUON, NU, ebin, ctbin, pdfws.GetProxyMap() )/
      futilws->GetFluxWsyst( MUON, NUB, ebin, ctbin, pdfws.GetProxyMap() );

    Double_t R_e_ae = futil->GetCachedFlux( ELEC, NU, ebin, ctbin )/
      futil->GetCachedFlux( ELEC, NUB, ebin, ctbin ) * (1. + futilws->GetVar("skew_e_ae")->getVal());

    Double_t Rws_e_ae = futilws->GetFluxWsyst( ELEC, NU, ebin, ctbin, pdfws.GetProxyMap() )/
      futilws->GetFluxWsyst( ELEC, NUB, ebin, ctbin, pdfws.GetProxyMap() );

    if ( TMath::Abs(R_mu_amu - Rws_mu_amu) > 1e-5 ) {
      cout << "NOTICE TestFitUtilWsystFlux test failed, mu-amu ratio is not modified as expected" << endl;
      return 1;
    }

    if ( TMath::Abs(R_e_ae - Rws_e_ae) > 1e-5 ) {
      cout << "NOTICE TestFitUtilWsystFlux test failed, e-ae ratio is not modified as expected" << endl;
      return 1;
    }

    // test 3: check that mu_e skew preserves muon + elec count in a random bin
    //-------------------------------------------------------------------------------------------------------

    futilws->GetVar("E_tilt")      ->setVal( 0. );
    futilws->GetVar("ct_tilt")     ->setVal( 0. );
    futilws->GetVar("skew_mu_amu") ->setVal( 0. );
    futilws->GetVar("skew_e_ae")   ->setVal( 0. );
    futilws->GetVar("skew_mu_e")   ->setVal( rand.Uniform(-0.5, 0.5) );

    Double_t flux_mu_e = ( futil->GetCachedFlux( MUON, NU, ebin, ctbin ) +
			   futil->GetCachedFlux( ELEC, NU, ebin, ctbin ) +
			   futil->GetCachedFlux( MUON, NUB, ebin, ctbin ) +
			   futil->GetCachedFlux( ELEC, NUB, ebin, ctbin ) );

    Double_t fluxws_mu_e = ( futilws->GetFluxWsyst( MUON, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
			     futilws->GetFluxWsyst( ELEC, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
			     futilws->GetFluxWsyst( MUON, NUB, ebin, ctbin, pdfws.GetProxyMap() ) +
			     futilws->GetFluxWsyst( ELEC, NUB, ebin, ctbin, pdfws.GetProxyMap() ) );

    if ( TMath::Abs(flux_mu_e-fluxws_mu_e) > 1e-5 ) {
      cout << "NOTICE TestFitUtilWsystFlux test failed at muon+elec flux conservation" << endl;
      return 1;
    }
    
    // test 4: check that the flux over all flav nu/nub is always conserved
    //-------------------------------------------------------------------------------------------------------
    futilws->GetVar("E_tilt")      ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("ct_tilt")     ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("skew_mu_amu") ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("skew_e_ae")   ->setVal( rand.Uniform(-0.5, 0.5) );
    futilws->GetVar("skew_mu_e")   ->setVal( rand.Uniform(-0.5, 0.5) );

    fluxint   = 0.;
    fluxintws = 0.;

    // loop over nu types
    for (Int_t flav = ELEC; flav <= TAU; flav++) {
      for (Int_t isnb = NU; isnb <= NUB; isnb++) {
	
	// loop over bins
	for (Int_t ebin = 1; ebin <= resp.GetHist3D()->GetXaxis()->GetNbins(); ebin++) {
	  for (Int_t ctbin = 1; ctbin <= resp.GetHist3D()->GetYaxis()->GetNbins(); ctbin++) {
	
	    Double_t flux   = futil->GetCachedFlux( flav, isnb, ebin, ctbin );
	    Double_t fluxws = futilws->GetFluxWsyst( flav, isnb, ebin, ctbin, pdfws.GetProxyMap() );

	    fluxint   += flux;
	    fluxintws += fluxws;

	  }
	}
	
      }
    }
    
    if ( TMath::Abs(fluxint - fluxintws) > 1e-2 ) {
      cout << "NOTICE TestFitUtilWsystFlux test failed at overall flux conservation" << endl;
      return 1;
    }
    
  } // end loop over trials


  cout << "NOTICE TestFitUtilWsystFlux test passed" << endl;
  return 0;
  
}
