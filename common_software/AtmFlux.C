#include "AtmFlux.h"
#include <fstream>
#include "TMath.h"
#include <sstream>
#include <stdexcept>

//*************************************************************************

/** Constructor.
 * \param opt   Determines which flux to use, options defined through enum AtmFluxOpt.
 * \param debug Data reading from file writes out everything for visual comparison.
 */
AtmFlux::AtmFlux(UInt_t opt, Bool_t debug) {

  fEmin  = 0;
  fEmax  = 0;
  fCtmin = 0;
  fCtmax = 0;
  
  //======================================================
  //init the 2D graphs that hold the flux data
  //======================================================
  
  for (Int_t f = 0; f < fNflavs; f++) {
    for (Int_t p = 0; p < fPtypes; p++) {
      fGraphs[f][p] = new TGraph2D();
    }
  }
  
  //======================================================
  //init the map options vs filenames
  //======================================================
  
  TString nmhdir = getenv("NMHDIR");

  if ( nmhdir == "" ) {
    throw std::invalid_argument("ERROR! AtmFlux::AtmFlux() $NMHDIR not set (source setenv.sh), init failed.");
  }
  
  TString datadir = nmhdir + "/data/honda_flux/";
  
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_frj_min    , datadir + "frj-ally-20-01-solmin.d"} );
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_frj_max    , datadir + "frj-ally-20-01-solmax.d"} );
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_frj_mtn_min, datadir + "frj-ally-20-01-mtn-solmin.d"} );
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_frj_mtn_max, datadir + "frj-ally-20-01-mtn-solmax.d"} );
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_grn_min    , datadir + "grn-ally-20-01-solmin.d"} );
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_grn_max    , datadir + "grn-ally-20-01-solmax.d"} );
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_grn_mtn_min, datadir + "grn-ally-20-01-mtn-solmin.d"} );
  fFileMap.insert( pair<UInt_t, TString> {AtmFluxOpt::h_grn_mtn_max, datadir + "grn-ally-20-01-mtn-solmax.d"} );
  
  //======================================================
  //check that option exists
  //======================================================
  
  if ( fFileMap.find(opt) == fFileMap.end() ) {
    throw std::invalid_argument("ERROR! AtmFlux::AtmFlux() option " + to_string(opt) + 
				" not supported, init failed.");
  }

  //======================================================
  //read in the data
  //======================================================
  ReadHondaFlux( fFileMap[opt], debug );

}

//*************************************************************************

/**
 * Destructor.
 */
AtmFlux::~AtmFlux() {

  //delete the 2D graphs that hold the flux data

  for (Int_t f = 0; f < fNflavs; f++) {
    for (Int_t p = 0; p < fPtypes; p++) {
      if (fGraphs[f][p]) delete fGraphs[f][p];
    }
  }
  
}

//*************************************************************************

/** Function to read azimuth-averaged honda flux tables.
 *
 *  The function assumes a very specific file format, exemplified in 
 *  NMH/data/frj-ally-20-01-mtn-solmax.d, downloaded in March 2018. If new azimuth-averaged
 *  files are added, checks are required to confirm that the input honda tables are in the
 *  same format.
 *
 *  Other Honda formats (azimuth dependent or all-direction averaged) are not supported.
 * *
 *  \param fname Name of the azimuth-averaged Honda flux file.
 *  \param debug Debug flag (default false). In debug mode it prints out 
 *               E, cosz, flux_mu, flux_mub, flux_e, flux_eb.
 */
void AtmFlux::ReadHondaFlux(TString fname, Bool_t debug) {
  
  ifstream in;             //create in-stream
  in.open(fname);          //load the file to in-stream
  Int_t nlines = 0;        //lines counter
  string line;             //line from file
  Double_t cosz = 2.;      //cosz value, updated from specific lines in file

  vector<Double_t> cos_range; // vector to help determine the cos-theta range
  vector<Double_t> e_range;   // vector to help determine the energy range
  
  //===================================================
  //check if file open failed
  //===================================================  
  if ( !in.is_open() ) {
    throw std::invalid_argument("ERROR! AtmFlux::ReadHondaFlux() could not open file " + (string)fname);
  }
  //===================================================
  //file open, read the data
  //===================================================  
  else {
    
    while ( getline (in, line) ) {

      nlines++;

      //===================================================
      //Skip the header line, extract the cosz from the line starting with 'average'
      //===================================================  
      
      if ( line.find("Enu") != string::npos ) { continue; }
      
      if ( line.find("average") != string::npos ) {

	// because string parsing in cpp is great, right. In the line that states the cosz range
	// fetch the string between the first "=" and "--", convert this to double, add bin width

	cosz = stod( line.substr( line.find("=") + 1, line.find("--") - line.find("=") - 1 ) ) + fCt_bw/2;
	cos_range.push_back(cosz);
	
	continue;
      }

      //===================================================
      //Get the energy and the fluxes, fill 2D graphs
      //===================================================  
      
      Double_t E_nu, f_mu, f_mub, f_e, f_eb;
      stringstream(line) >> E_nu >> f_mu >> f_mub >> f_e >> f_eb;

      if ( cos_range.size() == 1 ) e_range.push_back( E_nu ); //count the energies for the first cos-theta value
      
      fGraphs[0][0]->SetPoint( fGraphs[0][0]->GetN(), E_nu, cosz, f_e  * TMath::Power(E_nu, fPower) );
      fGraphs[0][1]->SetPoint( fGraphs[0][1]->GetN(), E_nu, cosz, f_eb * TMath::Power(E_nu, fPower) );
      fGraphs[1][0]->SetPoint( fGraphs[1][0]->GetN(), E_nu, cosz, f_mu * TMath::Power(E_nu, fPower) );
      fGraphs[1][1]->SetPoint( fGraphs[1][1]->GetN(), E_nu, cosz, f_mub* TMath::Power(E_nu, fPower) );
      
      if (debug) {
	cout << "E_nu, cosz, mu, mubar, e, ebar: " <<  E_nu << "\t" << cosz << "\t" << f_mu
	     << "\t" << f_mub << "\t" << f_e << "\t" << f_eb << endl;
      }
      
    } //end while

    //close file

    in.close();

    if (debug) {
      printf("AtmFlux::ReadHondaFlux(): found %d points in file %s\n",nlines, (const Char_t*)fname);
    }
    
  } //end fin is open

  fEmin  = *std::min_element( e_range.begin(), e_range.end() );
  fEmax  = *std::max_element( e_range.begin(), e_range.end() );
  fCtmin = *std::min_element( cos_range.begin(), cos_range.end() );
  fCtmax = *std::max_element( cos_range.begin(), cos_range.end() );

  // created from strings and have values like 0.1 + 1e-16; round them to 5 digits
  fEmin  = round(fEmin * 1e5)/1e5;
  fEmax  = round(fEmax * 1e5)/1e5;
  fCtmin = round(fCtmin * 1e5)/1e5 - fCt_bw/2;
  fCtmax = round(fCtmax * 1e5)/1e5 + fCt_bw/2; // extend range from 0.95 --> 1
  
  if (debug) {
    printf("AtmFlux::ReadHondaFlux(): energy range=[%f, %f], cos-theta range=[%f, %f]\n",
	   fEmin, fEmax, fCtmin, fCtmax );
  }
}

//*************************************************************************

/** Function that returns the differential atmospheric neutrino in units (m^2 sec rad GeV)^-1.
*
*   If you wish to calculate the flux in a certain bin of a (Energy, cosz) histogram, you need to 
*   mutliply the return value of this function by 
*   \f${\rm bin\_width}_{E} \times {\rm bin\_width}_{cosz}\f$.
*   The factor \f$2\pi\f$ from averaging over \f$\phi\f$ (azimuth) is taken into account.
*
*   A honda flux table points have been read to a TGraph2D, the TGraph2D method Interpolate()
*   is employed to estimate the flux at the input energy and angle.
*
*   \param nu_flavor neutrino flavor (0 - e, 1 - mu, 2 - tau)
*   \param is_nubar  0 - particle, 1 - antiparticle
*   \param E         neutrino energy (linear scale, not logE)
*   \param cosz      neutrino direction cosz
*
*   \return          flux/(dE * d(cosz) * m^2 * dt).
*
*/
Double_t AtmFlux::Flux_dE_dcosz(UInt_t nu_flavor, Bool_t is_nubar, Double_t E, Double_t cosz) {

  //----------------------------------------------------------
  // safety checks
  //----------------------------------------------------------

  if (nu_flavor == 2) return 0.; //tau flux is 0
  
  if (nu_flavor > 2 || (UInt_t)is_nubar > 1) {
    throw std::invalid_argument("ERROR! AtmFlux:Flux_dE_dcosz() wrong flavour (supported 0-2): " + to_string(nu_flavor) + " or is_nubar (supported 0, 1): " + to_string(is_nubar));
  }

  if ( (E <= fEmin) || (E >= fEmax) ) {
    throw std::invalid_argument("ERROR! AtmFlux::Flux_dE_dcosz() Energy " + to_string(E) + " outside the supported range of " + to_string(fEmin) + " - " + to_string(fEmax));
  }

  if ( cosz < fCtmin || cosz > fCtmax) {
    throw std::invalid_argument("ERROR! AtmFlux::Flux_dE_dcosz() cos-theta " + to_string(cosz) + " outside the supported range of " + to_string(fCtmin) + " - " + to_string(fCtmax));
  }
  
  //-----------------------------------------------------------------------
  // the below hokus-pokus is to deal with values in the cosz range 0.95-1. The last bin center in a honda table
  // is 0.95 and hence the TGraph cannot interpolate if cosz >= 0.95. However it is desired to be able to get
  // a value also for e.g. 0.975, for which the below logic returns the value corresponding to 0.95.
  //-----------------------------------------------------------------------
  
  if ( cosz <= (fCtmin + fCt_bw/2) ) {
    cosz = fCtmin + fCt_bw/2 + 1e-5; //interpolation works for cosz < -0.95, hence +1e-5
  }

  if ( cosz >= (fCtmax - fCt_bw/2) ) {
    cosz = fCtmax - fCt_bw/2 - 1e-5; //interpolation works for cosz > +0.95, hence -1e-5
  }
   
  return fGraphs[nu_flavor][(UInt_t)is_nubar]->Interpolate(E, cosz) * 2 * TMath::Pi() * TMath::Power(E, -fPower);
  
}
