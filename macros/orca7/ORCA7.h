#ifndef ORCA7_h
#define ORCA7_h

#include<vector>
using namespace std;

#include "NMHUtils.h"
#include "SummaryEvent.h"
#include "DetResponse.h"
#include "SummaryParser.h"
#include "FitUtilWsyst.h"
#include "FitPDF.h"

#include "TVector3.h"
#include "RooArgSet.h"

//===================================================================================================
// A namespace that stores some functions and structures necessary in the ORCA7 class
//===================================================================================================
namespace O7 {

  // functions to  be able to use the shower energy and track direction for high-purity tracks
  //---------------------------------------------------------------------------------------
  Double_t CustomEnergy(SummaryEvent* evt) {

    if ( evt->Get_shower_ql0() > 0.5 ) { return evt->Get_shower_energy(); }
    else                               { return evt->Get_track_energy();  }

  }

  TVector3 CustomDir(SummaryEvent *evt) { return evt->Get_track_dir();      }
  TVector3 CustomPos(SummaryEvent *evt) { return evt->Get_track_pos();      }
  Double_t CustomBY (SummaryEvent *evt) { return evt->Get_track_bjorkeny(); }

  // structure to store PID bin configurations
  //---------------------------------------------------------------------------------------
  struct PidBinConf {

    Double_t pid_min;
    Double_t pid_max;
    Double_t muon_cut;
    Double_t noise_cut;
    TString  name;

  PidBinConf() : pid_min(0), pid_max(0), muon_cut(0), noise_cut(0), name("") {}
    PidBinConf(Double_t _pid_min, Double_t _pid_max, Double_t _muon_cut, Double_t _noise_cut, TString _name) {
      pid_min   = _pid_min;
      pid_max   = _pid_max;
      muon_cut  = _muon_cut;
      noise_cut = _noise_cut;
      name = _name;
    }

  };

};

using namespace O7;

//===================================================================================================
/** Class for ORCA 7-line analysis with some variables/members that are used throughout several macros*/
//===================================================================================================
struct ORCA7 {

  /** Constructor */
  ORCA7(Bool_t ReadResponses);
  ~ORCA7();

  //*********************************************************************************************
  //*********************************************************************************************

  TString fDataF  = (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root";
  TString fEffmF  = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root";

  // detector response binning configuration
  Int_t f_R_ebins    = 20;     
  Int_t f_R_ctbins   = 40;     
  Int_t f_R_bybins   = 1;

  // detector response limits
  Double_t f_R_emin  = 1.0;
  Double_t f_R_emax  = 100.0;
  Double_t f_R_ctmin = -1.0;
  Double_t f_R_ctmax =  1.0;
  Double_t f_R_bymin = 0.0;
  Double_t f_R_bymax = 1.0;

  vector< PidBinConf > fPidBins;
  std::map< TString, DetResponse*> fResps;
  std::map< TString, FitPDF* > fPdfs;

  FitUtilWsyst *fFitUtil;
  RooArgSet fPriors;

  // fitutil input parameters
  Double_t f_F_runtime = 1.;
  Double_t f_F_emin    = 3;
  Double_t f_F_emax    = 60;
  Double_t f_F_ctmin   = -1;
  Double_t f_F_ctmax   = -1e-3;
  Double_t f_F_bymin   = 0;
  Double_t f_F_bymax   = 1;

};

#endif
