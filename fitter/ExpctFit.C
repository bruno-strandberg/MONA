#include "DetResponse.h"
#include "SummaryEvent.h"
#include "SummaryParser.h"
#include "NMHUtils.h"

#include "FitUtil.h"
#include "FitPDF.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"

#include "TFile.h"

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

using namespace RooFit;
using namespace std;

int main(int argc, char **argv) {

  //----------------------------------------------------------
  // parse command line arguments with Jpp
  //----------------------------------------------------------
  
  string        simdata_file;
  bool          refill_response;
  string        outputfile;
  string        effmh_elecCC;
  string        effmh_muonCC;
  string        effmh_tauCC;
  string        effmh_elecNC;
  
  try {

    JParser<> zap("Program to test RooFit fitter by fitting the expectation value");

    zap['s'] = make_field(simdata_file, "File with all summary data") =
      (string)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";

    zap['r'] = make_field(refill_response, "Flag to request re-filling of the detector responses");
    zap['o'] = make_field(outputfile, "File where output histograms are written") = "rootfiles/asymmetry.root";

    zap['w'] = make_field(effmh_elecCC, "Eff mass histograms for elec-CC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_elec_CC.root";

    zap['x'] = make_field(effmh_muonCC, "Eff mass histograms for muon-CC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_muon_CC.root";

    zap['y'] = make_field(effmh_tauCC , "Eff mass histograms for tau-CC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_tau_CC.root";

    zap['z'] = make_field(effmh_elecNC, "Eff mass histograms for elec-NC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_elec_NC.root";    
    zap(argc, argv);
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }
  
  //----------------------------------------------------------
  // set up the detector response
  //----------------------------------------------------------
  
  Int_t fENB   =  40;
  Int_t fEmin  =   1;
  Int_t fEmax  = 100;
  Int_t fCtNB  =  40;
  Int_t fCtmin =  -1;
  Int_t fCtmax =   1;
  Int_t fByNB  =   1;
  Int_t fBymin =   0;
  Int_t fBymax =   1;
  
  DetResponse track_resp (DetResponse::track, "track_resp", 
			  fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  track_resp.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_resp.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );    

  TString track_resp_name  = "rootfiles/track_response.root";
  
  if ( !NMHUtils::FileExists(track_resp_name) || refill_response ) {
    cout << "NOTICE ExpctFit() (Re)filling response" << endl;

    SummaryParser sp(simdata_file);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i % 100000 == 0) cout << "Event: " << i << endl;
      track_resp.Fill( sp.GetEvt(i) );
    }
    
    track_resp.WriteToFile(track_resp_name);
  }
  else {
    cout << "NOTICE ExpctFit() Reading in response" << endl;
    track_resp.ReadFromFile(track_resp_name);
  }

  cout << "NOTICE ExpctFit() Response ready" << endl;

  //----------------------------------------------------------
  // set up the PDF and get expectation value
  //----------------------------------------------------------
  
  FitUtil *fitutil = new FitUtil(3, track_resp.GetHist3D(),
  				 1, 100, -1, 0, 0, 1, effmh_elecCC, effmh_muonCC, effmh_tauCC, effmh_elecNC);

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_resp);  

  Double_t sinsqth12 = TMath::Power( TMath::Sin( 33.4 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth13 = TMath::Power( TMath::Sin( 8.42 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth23 = TMath::Power( TMath::Sin( 45   * TMath::Pi()/180. ), 2 );
  Double_t dcp       = 0.;
  Double_t dm32      = 2.44e-3;
  Double_t dm21      = 7.53e-5;
  Double_t DM        = dm32 + 0.5*dm21;
  Double_t dm31      = DM + 0.5*dm21;

  RooRealVar* rf_sinsqth12 = (RooRealVar*)fitutil->GetSet().find("SinsqTh12");
  RooRealVar* rf_sinsqth13 = (RooRealVar*)fitutil->GetSet().find("SinsqTh13");
  RooRealVar* rf_sinsqth23 = (RooRealVar*)fitutil->GetSet().find("SinsqTh23");
  RooRealVar* rf_dcp  = (RooRealVar*)fitutil->GetSet().find("dcp");
  RooRealVar* rf_dm21 = (RooRealVar*)fitutil->GetSet().find("Dm21");
  RooRealVar* rf_dm31 = (RooRealVar*)fitutil->GetSet().find("Dm31");

  rf_sinsqth12->setVal(sinsqth12);
  rf_sinsqth13->setVal(sinsqth13);
  rf_sinsqth23->setVal(sinsqth23);
  rf_dcp->setVal(dcp);
  rf_dm21->setVal(dm21);
  rf_dm31->setVal(dm31);
  
  TH3D *tracks_expct  = pdf_tracks.GetExpValHist();
  tracks_expct->SetNameTitle("tracks_expected", "tracks_expected");

  //----------------------------------------------------------
  // randomize starting point
  //----------------------------------------------------------
  rf_sinsqth12->randomize();
  rf_sinsqth13->randomize();
  rf_sinsqth23->randomize();
  rf_dcp->randomize();
  rf_dm21->randomize();
  rf_dm31->randomize();
  
  //----------------------------------------------------------
  // import data to RooFit and fit
  //----------------------------------------------------------

  RooDataHist rf_hist("rf_hist", "rf_hist", fitutil->GetObs(), Import(*tracks_expct) );
  RooFitResult *fitres = pdf_tracks.fitTo( rf_hist, Save(kTRUE) );
  RooArgSet result( fitres->floatParsFinal() );

  cout << "*********Fit result comparison****************************" << endl;
  cout << "dm21      actual: " << dm21      << "\t" << " fitted: " << ((RooRealVar*)result.find("Dm21"))->getVal() << endl;
  cout << "dm31      actual: " << dm31      << "\t" << " fitted: " << ((RooRealVar*)result.find("Dm31"))->getVal() << endl;
  cout << "sinsqth12 actual: " << sinsqth12 << "\t" << " fitted: " << ((RooRealVar*)result.find("SinsqTh12"))->getVal() << endl;
  cout << "sinsqth13 actual: " << sinsqth13 << "\t" << " fitted: " << ((RooRealVar*)result.find("SinsqTh13"))->getVal() << endl;
  cout << "sinsqth23 actual: " << sinsqth23 << "\t" << " fitted: " << ((RooRealVar*)result.find("SinsqTh23"))->getVal() << endl;
  cout << "dcp       actual: " << dcp       << "\t" << " fitted: " << ((RooRealVar*)result.find("dcp"))->getVal() << endl;
  cout << "*********Fit result comparison****************************" << endl;
  
}

