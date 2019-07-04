
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

class FU : public FitUtilWsyst {

public:
  FU(Double_t op_time, TH3 *h_template,
     Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
     TString meff_file) :
    FitUtilWsyst(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {};

  virtual ~FU() {};

  Int_t TestFluxNorm(DetResponse *resp) {

    FitPDF pdfws("pdfws" ,"pdfws", this, resp);
    
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

      GetVar("E_tilt")      ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("ct_tilt")     ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("skew_mu_amu") ->setVal( 0. );
      GetVar("skew_e_ae")   ->setVal( 0. );
      GetVar("skew_mu_e")   ->setVal( 0. );

      UInt_t flav = rand.Integer(3); // select random flavor
      UInt_t isnb = rand.Integer(2); // select random  nu/nub
    
      Double_t fluxint   = 0.; // flux integral over (E,ct) without tilt systematics
      Double_t fluxintws = 0.; // flux integral over (E,ct) with tilt systematics

      // loop over bins
      for (Int_t ebin = 1; ebin <= resp->GetHist3D()->GetXaxis()->GetNbins(); ebin++) {
	for (Int_t ctbin = 1; ctbin <= resp->GetHist3D()->GetYaxis()->GetNbins(); ctbin++) {
	
	  Double_t flux   = GetCachedFlux( flav, isnb, ebin, ctbin );
	  Double_t fluxws = GetFluxWsyst( flav, isnb, ebin, ctbin, pdfws.GetProxyMap() );

	  fluxint   += flux;
	  fluxintws += fluxws;

	}
      }

      if ( TMath::Abs(fluxint - fluxintws) > 1e-2 ) {
	cout << "NOTICE TestFitUtilWsystFlux test failed at flux integral preservation with tilt parameters, fluxes and diff: " << fluxint << "\t" << fluxintws << "\t" << fluxint - fluxintws << endl;
	this->PrintParameters();
	return 1;
      }

      // test 2: set tilt parameters to 0, set random skew parameters; check, in random bin, conservation of
      // i) mu+amu, ii) e+ae. If these are conserved for one bin, they are conserve the overall integral.
      // The mu_e skew is 0, because it does not preserve mu-e ratio (obviously)
      //-------------------------------------------------------------------------------------------------------

      GetVar("E_tilt")      ->setVal( 0. );
      GetVar("ct_tilt")     ->setVal( 0. );
      GetVar("skew_mu_amu") ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("skew_e_ae")   ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("skew_mu_e")   ->setVal( 0. );

      Int_t ebin  = 1 + rand.Integer( resp->GetHist3D()->GetXaxis()->GetNbins() );
      Int_t ctbin = 1 + rand.Integer( resp->GetHist3D()->GetYaxis()->GetNbins() );

      // no systematics
      Double_t flux_mu = GetCachedFlux( MUON, NU, ebin, ctbin )
	+ GetCachedFlux( MUON, NUB, ebin, ctbin );
    
      Double_t flux_e = GetCachedFlux( ELEC, NU, ebin, ctbin )
	+ GetCachedFlux( ELEC, NUB, ebin, ctbin );

      // with systematics
      Double_t fluxws_mu = GetFluxWsyst( MUON, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
	GetFluxWsyst( MUON, NUB, ebin, ctbin, pdfws.GetProxyMap() );

      Double_t fluxws_e = GetFluxWsyst( ELEC, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
	GetFluxWsyst( ELEC, NUB, ebin, ctbin, pdfws.GetProxyMap() );

      if ( TMath::Abs(flux_mu-fluxws_mu) > 1e-5 ) {
	cout << "NOTICE TestFitUtilWsystFlux test failed at muon flux conservation, fluxes and diff: " << flux_mu << "\t" << fluxws_mu << "\t" << flux_mu - fluxws_mu << endl;
	this->PrintParameters();
	return 1;
      }

      if ( TMath::Abs(flux_e-fluxws_e) > 1e-5 ) {
	cout << "NOTICE TestFitUtilWsystFlux test failed at elec flux conservation, fluxes and diff: " << flux_e << "\t" << fluxws_e << "\t" << flux_e - fluxws_e << endl;
	this->PrintParameters();
	return 1;
      }

      // check ratios
      Double_t R_mu_amu = GetCachedFlux( MUON, NU, ebin, ctbin )/
	GetCachedFlux( MUON, NUB, ebin, ctbin ) * (1. + GetVar("skew_mu_amu")->getVal());

      Double_t Rws_mu_amu = GetFluxWsyst( MUON, NU, ebin, ctbin, pdfws.GetProxyMap() )/
	GetFluxWsyst( MUON, NUB, ebin, ctbin, pdfws.GetProxyMap() );

      Double_t R_e_ae = GetCachedFlux( ELEC, NU, ebin, ctbin )/
	GetCachedFlux( ELEC, NUB, ebin, ctbin ) * (1. + GetVar("skew_e_ae")->getVal());

      Double_t Rws_e_ae = GetFluxWsyst( ELEC, NU, ebin, ctbin, pdfws.GetProxyMap() )/
	GetFluxWsyst( ELEC, NUB, ebin, ctbin, pdfws.GetProxyMap() );

      if ( TMath::Abs(R_mu_amu - Rws_mu_amu) > 1e-5 ) {
	cout << "NOTICE TestFitUtilWsystFlux test failed, mu-amu ratio is not modified as expected, ratios and diff: " << R_mu_amu << "\t" << Rws_mu_amu << "\t" << R_mu_amu - Rws_mu_amu << endl;
	this->PrintParameters();
	return 1;
      }

      if ( TMath::Abs(R_e_ae - Rws_e_ae) > 1e-5 ) {
	cout << "NOTICE TestFitUtilWsystFlux test failed, e-ae ratio is not modified as expected, ratios and diff: " << R_e_ae << "\t" << Rws_e_ae << "\t" << R_e_ae - Rws_e_ae << endl;
	this->PrintParameters();
	return 1;
      }

      // test 3: check that mu_e skew preserves muon + elec count in a random bin
      //-------------------------------------------------------------------------------------------------------

      GetVar("E_tilt")      ->setVal( 0. );
      GetVar("ct_tilt")     ->setVal( 0. );
      GetVar("skew_mu_amu") ->setVal( 0. );
      GetVar("skew_e_ae")   ->setVal( 0. );
      GetVar("skew_mu_e")   ->setVal( rand.Uniform(-0.5, 0.5) );

      Double_t flux_mu_e = ( GetCachedFlux( MUON, NU, ebin, ctbin ) +
			     GetCachedFlux( ELEC, NU, ebin, ctbin ) +
			     GetCachedFlux( MUON, NUB, ebin, ctbin ) +
			     GetCachedFlux( ELEC, NUB, ebin, ctbin ) );

      Double_t fluxws_mu_e = ( GetFluxWsyst( MUON, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
			       GetFluxWsyst( ELEC, NU, ebin, ctbin, pdfws.GetProxyMap() ) +
			       GetFluxWsyst( MUON, NUB, ebin, ctbin, pdfws.GetProxyMap() ) +
			       GetFluxWsyst( ELEC, NUB, ebin, ctbin, pdfws.GetProxyMap() ) );

      if ( TMath::Abs(flux_mu_e-fluxws_mu_e) > 1e-5 ) {
	cout << "NOTICE TestFitUtilWsystFlux test failed at muon+elec flux conservation, fluxes and diff: " << flux_mu_e << "\t" << fluxws_mu_e << "\t" << flux_mu_e - fluxws_mu_e << endl;
	this->PrintParameters();
	return 1;
      }
    
      // test 4: check that the flux over all flav nu/nub is always conserved
      //-------------------------------------------------------------------------------------------------------
      GetVar("E_tilt")      ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("ct_tilt")     ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("skew_mu_amu") ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("skew_e_ae")   ->setVal( rand.Uniform(-0.5, 0.5) );
      GetVar("skew_mu_e")   ->setVal( rand.Uniform(-0.5, 0.5) );

      fluxint   = 0.;
      fluxintws = 0.;

      // loop over nu types
      for (Int_t flav = ELEC; flav <= TAU; flav++) {
	for (Int_t isnb = NU; isnb <= NUB; isnb++) {
	
	  // loop over bins
	  for (Int_t ebin = 1; ebin <= resp->GetHist3D()->GetXaxis()->GetNbins(); ebin++) {
	    for (Int_t ctbin = 1; ctbin <= resp->GetHist3D()->GetYaxis()->GetNbins(); ctbin++) {
	
	      Double_t flux   = GetCachedFlux( flav, isnb, ebin, ctbin );
	      Double_t fluxws = GetFluxWsyst( flav, isnb, ebin, ctbin, pdfws.GetProxyMap() );

	      fluxint   += flux;
	      fluxintws += fluxws;

	    }
	  }
	
	}
      }
    
      if ( TMath::Abs(fluxint - fluxintws) > 1e-2 ) {
	cout << "NOTICE TestFitUtilWsystFlux test failed at overall flux conservation, fluxes and diff: " << fluxint << "\t" << fluxintws << "\t" << fluxint - fluxintws << endl;
	this->PrintParameters();
	return 1;
      }
    
    } // end loop over trials


    cout << "NOTICE TestFitUtilWsystFlux test passed" << endl;
    return 0;
    
  }
  
};

//******************************************************************************************

/** Application to test the flux normalisation */
int main(const int argc, const char **argv) {
    
  DetResponse resp(DetResponse::track, "truthresp", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  FU fu(3, resp.GetHist3D(), 3, 70, -1, 0, 0, 1, EffMass::DUMMYFILE);
  return fu.TestFluxNorm(&resp);
  
}
