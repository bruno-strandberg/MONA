#ifndef ORCA7_h
#define ORCA7_h

#include<vector>
using namespace std;

#include "NMHUtils.h"
#include "SummaryEvent.h"
#include "DetResponse.h"
#include "SummaryParser.h"

//===================================================================================================
// functions to  be able to use the shower energy and track direction for high-purity tracks
//===================================================================================================
Double_t CustomEnergy(SummaryEvent* evt) {

  if ( evt->Get_shower_ql0() > 0.5 ) { return evt->Get_shower_energy(); }
  else                               { return evt->Get_track_energy();  }

}

TVector3 CustomDir(SummaryEvent *evt) { return evt->Get_track_dir();      }
TVector3 CustomPos(SummaryEvent *evt) { return evt->Get_track_pos();      }
Double_t CustomBY (SummaryEvent *evt) { return evt->Get_track_bjorkeny(); }


//===================================================================================================
/** Class for ORCA 7-line analysis with some variables/members that are used throughout several macros*/
//===================================================================================================
struct ORCA7 {

  /** Constructor */
  ORCA7(Bool_t ReadResponses) {

    // set the PID bins
    //-------------------------------------------------------------------------------
    fPidBins.push_back( PidBinConf(0.0, 0.3, 0.05, 0.0, "shw") );
    fPidBins.push_back( PidBinConf(0.3, 0.7, 0.01, 0.0, "mid") );
    fPidBins.push_back( PidBinConf(0.7, 1.0+1e-5, 0.05, 0.1, "trk") );

    // create a response for each PID bin
    //-------------------------------------------------------------------------------
    for (auto PB: fPidBins) {
      
      if ( PB.pid_min < 0.7 ) {
	fResps.push_back( new DetResponse(DetResponse::shower, PB.name+"R", f_R_ebins, f_R_emin, f_R_emax, f_R_ctbins, f_R_ctmin, f_R_ctmax, f_R_bybins, f_R_bymin, f_R_bymax) );
	fResps.back()->AddCut( &SummaryEvent::Get_shower_ql0, std::greater<double>(), 0.5, true);
      }
      else {
	fResps.push_back( new DetResponse(DetResponse::customreco, PB.name+"R", f_R_ebins, f_R_emin, f_R_emax, f_R_ctbins, f_R_ctmin, f_R_ctmax, f_R_bybins, f_R_bymin, f_R_bymax) );
	fResps.back()->SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );	
	fResps.back()->AddCut( &SummaryEvent::Get_track_ql0, std::greater<double>(), 0.5, true);
      }

      fResps.back()->AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>()   , PB.noise_cut, true);
      fResps.back()->AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>()   , PB.muon_cut , true);
      fResps.back()->AddCut( &SummaryEvent::Get_RDF_track_score , std::greater_equal<double>(), PB.pid_min  , true);
      fResps.back()->AddCut( &SummaryEvent::Get_RDF_track_score , std::less<double>()         , PB.pid_max  , true);

    }

    // create an output name for each response
    //-------------------------------------------------------------------------------

    TString out_dir = NMHUtils::Getcwd() + "/rootfiles/";   // directory where responses are stored
    system("mkdir -p " + out_dir);

    vector< std::pair<TString, DetResponse*> > resp_names;  // vector that stores the outname and the response pairs
    Bool_t ReadFromFile = ReadResponses;                    // flag to indicate whether responses should be read

    for (auto R: fResps) {
      TString outname = out_dir + "resp_" + R->Get_RespName() + ".root";
      resp_names.push_back( std::make_pair(outname, R) );
      ReadFromFile = ReadFromFile && NMHUtils::FileExists(outname);
    }

    // read the responses from file or fill from summary data
    //-------------------------------------------------------------------------------
    if ( ReadFromFile ) {
      for (auto p: resp_names) p.second->ReadFromFile( p.first );
    }
    else {
      
      SummaryParser sp(fDataF);
      for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
	if (i % 100000 == 0) cout << "NOTICE ORCA7::ORCA7() event " << i << " of " << sp.GetTree()->GetEntries() << endl;
	for (auto R: fResps) R->Fill( sp.GetEvt(i) );
      }
      
      for (auto p: resp_names) p.second->WriteToFile( p.first );

    }

  }

  TString fDataF  = (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root";
  TString fEffmF  = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root";

  // structure to store PID bin configurations
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

  // vector with PID bin edges, each has an individual noise and muon suppression score
  vector< PidBinConf > fPidBins;

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

  // responses
  vector< DetResponse* > fResps;

  // fitutil input parameters
  Double_t f_F_runtime = 1.;
  Double_t f_F_emin    = 3;
  Double_t f_F_emax    = 75;
  Double_t f_F_ctmin   = -1;
  Double_t f_F_ctmax   = 0;
  Double_t f_F_bymin   = 0;
  Double_t f_F_bymax   = 1;

};

#endif
