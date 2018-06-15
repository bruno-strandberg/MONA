#include "SummaryEvent.h"

SummaryEvent::SummaryEvent() {

  fMC_runID        = 0.;       
  fMC_evtID        = 0.;       
  fMC_w2           = 0.;          
  fMC_w1y          = 0.;         
  fMC_erange_start = 0.;
  fMC_is_CC        = 0.;       
  fMC_is_neutrino  = 0.; 
  fMC_type         = 0.;        
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

  fRDF_muon_score  = 0.;  
  fRDF_track_score = 0.; 
  fRDF_noise_score = 0.; 

}

SummaryEvent::~SummaryEvent() {}
