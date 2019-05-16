#include "SummaryEvent.h"

SummaryEvent::SummaryEvent() {

  fMC_runID        = 0.;       
  fMC_evtID        = 0.;       
  fMC_w2           = 0.;          
  fMC_w1y          = 0.;
  fMC_w2denom      = 0.;         
  fMC_erange_start = 0.;
  fMC_erange_stop  = 0.;
  fMC_is_CC        = 0.;       
  fMC_is_neutrino  = 0.; 
  fMC_type         = 0.;
  fMC_ichan        = 0.;        
  fMC_dir_x        = 0.;
  fMC_dir_y        = 0.;
  fMC_dir_z        = 0.;
  fMC_pos_x        = 0.;
  fMC_pos_y        = 0.;
  fMC_pos_z        = 0.;
  fMC_energy       = 0.;
  fMC_bjorkeny     = 0.;

  fTrack_dir_x     = 0.;
  fTrack_dir_y     = 0.;
  fTrack_dir_z     = 0.;
  fTrack_pos_x     = 0.;
  fTrack_pos_y     = 0.;
  fTrack_pos_z     = 0.;
  fTrack_energy    = 0.;   
  fTrack_bjorkeny  = 0.; 
  fTrack_ql0       = 0.;      
  fTrack_ql1       = 0.;      
  fTrack_ql2       = 0.;
  fTrack_ql3       = 0.;      

  fShower_dir_x    = 0.;
  fShower_dir_y    = 0.;
  fShower_dir_z    = 0.;
  fShower_pos_x    = 0.;
  fShower_pos_y    = 0.;
  fShower_pos_z    = 0.;
  fShower_energy   = 0.;   
  fShower_bjorkeny = 0.; 
  fShower_ql0      = 0.;      
  fShower_ql1      = 0.;
  fShower_ql2      = 0.;
  fShower_ql3      = 0.;      

  fRDF_muon_score  = 0.;  
  fRDF_track_score = 0.; 
  fRDF_noise_score = 0.; 

}

//***************************************************************************************

SummaryEvent::~SummaryEvent() {}

//***************************************************************************************

/** This function writes simplistic pseudodata to the member variables and is mainly useful for testing.

    For example, the pseudo-events can be processed through event filters/selections/detector responses.

    \param logE If true, energy sampled from \f$ 10^{{\rm Uniform}(0,2)}\f$,otherwise from \f$ {\rm Uniform}(1,100)\f$
*/
void SummaryEvent::FillPseudoData(Bool_t logE) {

  //-----------------------------------------------------------------------------
  // some hardcoded limits that may be configured through the parameter list in the future
  //-----------------------------------------------------------------------------

  Double_t emin   =    1;  //min energy
  Double_t emax   =  100;  //max energy
  Double_t posmin = -150;  //dimensions of the position cube that contains the pseudoevents in meters
  Double_t posmax =  150;

  Double_t trkEres   = 50;    //FWHM of a gaussian energy resolution for tracks in %
  Double_t trkCTres  = 10;    //FWHM of a gaussian angle resolution for tracks in %
  Double_t trkPOSres = 30;    //FWHM of a gaussian position resolution for tracks in %
  Double_t trkBYres  = 50;    //FWHM of a gaussian bjorken-y resultion for tracks in %

  Double_t shwEres   = 30;    //FWHM of a gaussian energy resolution for tracks in %
  Double_t shwCTres  = 30;    //FWHM of a gaussian angle resolution for tracks in %
  Double_t shwPOSres = 20;    //FWHM of a gaussian position resolution for tracks in %
  Double_t shwBYres  = 40;    //FWHM of a gaussian bjorken-y resultion for tracks in %

  Double_t conv = 2.35;       //to convert a FWHM resolution to sigmas; e.g. 50% fwhm --> sigma = 21.2

  //-----------------------------------------------------------------------------
  //not using variables that are used to identify files and weight variables
  //-----------------------------------------------------------------------------
  fMC_runID        = 0.;
  fMC_evtID        = 0.;       
  fMC_w2           = 0.;          
  fMC_w1y          = 0.;
  fMC_w2denom      = 0.;         
  fMC_erange_start = 0.;
  fMC_erange_stop  = 0.;
  fMC_ichan        = 0.;
  
  //-----------------------------------------------------------------------------
  //set random MC-truth variables
  //-----------------------------------------------------------------------------
  fMC_is_CC        = (Double_t)fRand.Integer(2);                    // either nc or cc
  fMC_is_neutrino  = 1 ;                                            // pseudodata only neutrinos
  fMC_type         = fFlavs[ fRand.Integer( fFlavs.size() ) ];      // flavor number from flavor vector
  fRand.Sphere(fMC_dir_x, fMC_dir_y, fMC_dir_z, 1);                 // create a random unit vector
  fMC_pos_x        = fRand.Uniform(posmin, posmax);                 // select within a random cube for position
  fMC_pos_y        = fRand.Uniform(posmin, posmax);
  fMC_pos_z        = fRand.Uniform(posmin, posmax);
  if (logE) { fMC_energy = TMath::Power( 10, fRand.Uniform(0,2) ); } // log falloff in range 1-100 GeV
  else      { fMC_energy = fRand.Uniform(emin,emax); }               // uniform 1-100 GeV
  fMC_bjorkeny     = fRand.Uniform(0,1);                             // bjorken-y range 0 to 1

  //-----------------------------------------------------------------------------
  //set track reco variables, depending on the resolution parameters
  //-----------------------------------------------------------------------------
  TVector3 trkdir( fMC_dir_x + fRand.Gaus(0, trkCTres/conv)/100  * fMC_dir_x,
		   fMC_dir_y + fRand.Gaus(0, trkCTres/conv)/100  * fMC_dir_y,
		   fMC_dir_z + fRand.Gaus(0, trkCTres/conv)/100  * fMC_dir_z);
  fTrack_dir_x     = trkdir.Unit().X();
  fTrack_dir_y     = trkdir.Unit().Y();
  fTrack_dir_z     = trkdir.Unit().Z();
  fTrack_pos_x     = fMC_pos_x    + fRand.Gaus(0, trkPOSres/conv)/100 * fMC_pos_x;
  fTrack_pos_y     = fMC_pos_y    + fRand.Gaus(0, trkPOSres/conv)/100 * fMC_pos_y;
  fTrack_pos_z     = fMC_pos_z    + fRand.Gaus(0, trkPOSres/conv)/100 * fMC_pos_z;
  fTrack_energy    = fMC_energy   + fRand.Gaus(0, trkEres/conv)/100   * fMC_energy;
  fTrack_bjorkeny  = fMC_bjorkeny + fRand.Gaus(0, trkBYres/conv)/100  * fMC_bjorkeny;

  if (fTrack_bjorkeny > 1.) fTrack_bjorkeny = 1 - 1e-5; // to avoid overflow in by
  if (fTrack_bjorkeny < 0.) fTrack_bjorkeny = 0;

  fTrack_ql0       = 1.;      //quality level 0 always passed
  fTrack_ql1       = 0.;      
  fTrack_ql2       = 0.;      

  //-----------------------------------------------------------------------------
  //set shower reco variables, depending on the resolution parameters
  //-----------------------------------------------------------------------------
  TVector3 shwdir(fMC_dir_x + fRand.Gaus(0, shwCTres/conv)/100  * fMC_dir_x,
		  fMC_dir_y + fRand.Gaus(0, shwCTres/conv)/100  * fMC_dir_y,
		  fMC_dir_z + fRand.Gaus(0, shwCTres/conv)/100  * fMC_dir_z);
  fShower_dir_x    = shwdir.Unit().X();
  fShower_dir_y    = shwdir.Unit().Y();
  fShower_dir_z    = shwdir.Unit().Z();
  fShower_pos_x    = fMC_pos_x    + fRand.Gaus(0, shwPOSres/conv)/100 * fMC_pos_x;
  fShower_pos_y    = fMC_pos_y    + fRand.Gaus(0, shwPOSres/conv)/100 * fMC_pos_y;
  fShower_pos_z    = fMC_pos_z    + fRand.Gaus(0, shwPOSres/conv)/100 * fMC_pos_z;
  fShower_energy   = fMC_energy   + fRand.Gaus(0, shwEres/conv)/100   * fMC_energy;
  fShower_bjorkeny = fMC_bjorkeny + fRand.Gaus(0, shwBYres/conv)/100  * fMC_bjorkeny;

  if (fShower_bjorkeny > 1.) fShower_bjorkeny = 1 - 1e-5; // to avoid overflow in by
  if (fShower_bjorkeny < 0.) fShower_bjorkeny = 0;

  fShower_ql0       = 1.;      //quality level 0 always passed
  fShower_ql1       = 0.;      
  fShower_ql2       = 0.;      

  //-----------------------------------------------------------------------------
  //finally, set some PID score; nothing for noise and muon score
  //-----------------------------------------------------------------------------

  if (fMC_is_CC < 0.5) {
    fRDF_track_score = fRand.Uniform(0., 0.5); // NC events score always below 0.5
  }
  else {
    if ( TMath::Abs(fMC_type) == 14 ) {
      fRDF_track_score = fRand.Uniform(0.5, 1); // muon-CC events score always above 0.5
    }
    else {
      fRDF_track_score = fRand.Uniform(0., 1.); // others get a random number in the interval
    }
  }

  fRDF_muon_score  = 0.;  
  fRDF_noise_score = 0.; 
  
}
