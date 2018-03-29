#include "AtmFlux.h"
#include <fstream>
#include "TMath.h"
#include<sstream>

//*************************************************************************

/** Constructor.
 * \param opt   Determines which flux to use, options defined through enum AtmFluxOpt.
 * \param debug Data reading from file writes out everything for visual comparison.
 */
AtmFlux::AtmFlux(UInt_t opt, Bool_t debug) {

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

  Bool_t InitOK = true;
  
  TString nmhdir = getenv("NMHDIR");

  if ( nmhdir == "" ) {
    cout << "ERROR! AtmFlux::AtmFlux() $NMHDIR not set (source setenv.sh), init failed." << endl;
    InitOK = false;
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
    cout << "ERROR! AtmFlux::AtmFlux() option " << opt
	 << " not supported, init failed." << endl;
    InitOK = false;
  }

  //======================================================
  //read in the data
  //======================================================

  if (InitOK) {
    if ( !ReadHondaFlux( fFileMap[opt], debug ) ) {
      cout << "ERROR! AtmFlux::AtmFlux() reading fluxes from file failed." << endl;
    }    
  }
  
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
 *
 *  Honda flux bin centers are from -0.95 -- 0.95. We also wish to have flux at -1 and 1.
 *  To accommodate this, the flux corresponding to E, cosz = +-0.95 is also set for points 
 *  E, cosz = +-1 (otherwise TGraph2D->Interpolate() does not work for abs(cosz) > 0.95).
 *
 *  \param fname Name of the azimuth-averaged Honda flux file.
 *  \param debug Debug flag (default false). In debug mode it prints out 
 *               E, cosz, flux_mu, flux_mub, flux_e, flux_eb.
 */
Bool_t AtmFlux::ReadHondaFlux(TString fname, Bool_t debug) {

  Bool_t ReadOK = false;
  
  ifstream in;             //create in-stream
  in.open(fname);          //load the file to in-stream
  Int_t nlines = 0;        //lines counter
  string line;             //line from file
  Double_t cosz     = 2.;  //cosz value, updated from specific lines in file
  Double_t cosz_bin = 0.1; //cosz bin width

  
  //===================================================
  //check if file open failed
  //===================================================  
  if ( !in.is_open() ) {
    cout << "ERROR! AtmFlux::ReadHondaFlux() could not open file " << fname
	 << ", data read failed." << endl;
    ReadOK = false;
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

	cosz = stod( line.substr( line.find("=") + 1, line.find("--") - line.find("=") - 1 ) )
	  + cosz_bin/2;
	
	continue;
      }

      //===================================================
      //Get the energy and the fluxes, fill 2D graphs
      //===================================================  
      
      Double_t E_nu, f_mu, f_mub, f_e, f_eb;
      stringstream(line) >> E_nu >> f_mu >> f_mub >> f_e >> f_eb;

      fGraphs[0][0]->SetPoint( fGraphs[0][0]->GetN(), E_nu, cosz, f_e  );
      fGraphs[0][1]->SetPoint( fGraphs[0][1]->GetN(), E_nu, cosz, f_eb );
      fGraphs[1][0]->SetPoint( fGraphs[1][0]->GetN(), E_nu, cosz, f_mu );
      fGraphs[1][1]->SetPoint( fGraphs[1][1]->GetN(), E_nu, cosz, f_mub);

      //Honda flux last bin is cosz +- 0.95. Add flux(E, cosz=+-0.95) value also to points
      //cosz=+-1. Otherwise TGraph2D->Interpolate() will give zeroes for abs(cosz) > 0.95.
            
      if ( TMath::Abs(cosz) > 0.9 ) {

	Double_t cosz_end = 2;
	cosz > 0 ? cosz_end = 1. : cosz_end = -1. ;
	
	fGraphs[0][0]->SetPoint( fGraphs[0][0]->GetN(), E_nu, cosz_end, f_e  );
	fGraphs[0][1]->SetPoint( fGraphs[0][1]->GetN(), E_nu, cosz_end, f_eb );
	fGraphs[1][0]->SetPoint( fGraphs[1][0]->GetN(), E_nu, cosz_end, f_mu );
	fGraphs[1][1]->SetPoint( fGraphs[1][1]->GetN(), E_nu, cosz_end, f_mub);
	
      }
      
      if (debug) {
	cout << "E_nu, cosz, mu, mubar, e, ebar: " <<  E_nu << "\t" << cosz << "\t" << f_mu
	     << "\t" << f_mub << "\t" << f_e << "\t" << f_eb << endl;
      }
      
    } //end while

    //close file

    ReadOK = true;
    in.close();

    if (debug) {
      printf("AtmFlux::ReadHondaFlux(): found %d points in file %s\n",nlines, (const Char_t*)fname);
    }
    
  } //end fin is open

  return ReadOK;
  
}

//*************************************************************************

/** Function that returns the atmospheric neutrino in (m^2 sec sr GeV)^-1.
*
*   A honda flux table points have been read to a TGraph2D, the TGraph2D method Interpolate()
*   is employed to estimate the flux at the input energy and angle.
*
*   \param nu_flavor neutrino flavor (0 - e, 1 - mu, 2 - tau)
*   \param is_nubar  0 - particle, 1 - antiparticle
*   \param E         neutrino energy
*   \param cosz      neutrino direction cosz
*
*/
Double_t AtmFlux::GetHondaFlux(UInt_t nu_flavor, Bool_t is_nubar, Double_t E, Double_t cosz) {

  // safety checks
  
  if (nu_flavor == 2) return 0.; //tau flux is 0

  if (nu_flavor > 2 || (UInt_t)is_nubar > 1) {
    cout << "ERROR! AtmFlux::GetHondaFlux() wrong flavour (supported 0-2): " << nu_flavor
	 << " or is_nubar (supported true, false): " << is_nubar << ", returning 0." << endl;
    return 0.;
  }
  
  if ( (E < 0.1) || (E > 1.0000E+04) || TMath::Abs(cosz) > 1. ) {
    cout << "ERROR! AtmFlux::GetHondaFlux() wrong Energy (supported 0.1 - 10000 GeV): " << E
	 << " or cosz (supported -1 - 1): " << cosz << ", returning 0." << endl;
    return 0.;
  }

  return fGraphs[nu_flavor][(UInt_t)is_nubar]->Interpolate(E, cosz);
  
}
