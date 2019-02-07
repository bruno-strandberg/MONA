#include "HelperFunctions.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TString.h"

#include "DetResponse.h"
#include "EventSelection.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include "FitUtil.h"
#include "FitPDF.h"

#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooRealVar.h"


#include <iostream>
using namespace std;
using namespace RooFit;

/* Script that fits seperately track and shower expectation values of both orderings
 * to the wrong ordering. In particular the track and shower pdfs are not linked, ie.
 * this is not a simultaneous fit, meaning the final output parameters are not identical.
 * This means that technically the fits are these two points in phase space can note be
 * compared, however this script does it. THIS SCRIPT is meant as an EXAMPLE, not an 
 * actual analysis result.
 */


void AsimovFitSeperate() {

  TString filefolder = "./default_detres/";

  // DetRes and EvSel input values
  Int_t EBins = 40;
  Int_t EMin = 1;
  Int_t EMax = 100;
  Int_t ctBins = 40;
  Int_t ctMin = -1;
  Int_t ctMax = 1;
  Int_t byBins = 1;
  Int_t byMin = 0;
  Int_t byMax = 1;

  
  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  DetResponse track_response(DetResponse::track, "track_response", EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
  track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  DetResponse shower_response(DetResponse::shower, "shower_response", EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
  shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------
  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);

  TString track_file = "track_response.root";
  TString shower_file = "shower_response.root";

  if ( !NMHUtils::FileExists(filefolder + track_file) or !NMHUtils::FileExists(filefolder + shower_file) ) {

    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
      SummaryEvent *evt = sp.GetEvt(i);
      track_response.Fill(evt);
      shower_response.Fill(evt);
    }

    track_response.WriteToFile("track_response.root");
    shower_response.WriteToFile("shower_response.root");

    cout << "NOTICE: Finished filling response" << endl;
  } 
  else {
    cout << "NOTICE: Reading responses from disk" << endl;
    track_response.ReadFromFile(filefolder + track_file);
    shower_response.ReadFromFile(filefolder + shower_file);
  }

  //----------------------------------------------------------
  // set up the PDFs and static oscillation parameters
  //----------------------------------------------------------

  auto meff_file = (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(), 2, 80, -1, 0, 0, 1, meff_file);

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_response);
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response);

 
  // Set values to NO
  fitutil->SetNOlims();
  fitutil->SetNOcentvals();

  Double_t dm31 = 0; // placeholder for later
  Double_t sinSqTh23 = 0; // placeholder for later

  TH3D* tracks_no   = (TH3D*)pdf_tracks.GetExpValHist();
  TH3D* showers_no  = (TH3D*)pdf_showers.GetExpValHist();
  tracks_no->SetName("tracks_expval_no");
  showers_no->SetName("showers_expval_no");


  // Set values to IO
  fitutil->SetIOlims();
  fitutil->SetIOcentvals();

  TH3D* tracks_io   = (TH3D*)pdf_tracks.GetExpValHist();
  TH3D* showers_io  = (TH3D*)pdf_showers.GetExpValHist();
  tracks_io->SetName("tracks_expval_io");
  showers_io->SetName("showers_expval_io");

  fitutil->GetVar("SinsqTh12")->setConstant(kTRUE);
  fitutil->GetVar("SinsqTh13")->setConstant(kTRUE);
  fitutil->GetVar("dcp")->setConstant(kTRUE);
  fitutil->GetVar("Dm21")->setConstant(kTRUE);
  
  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------
  cout << "NOTICE Fitter started fitting" << endl;
  
  TStopwatch timer;

  // Fit under IO model, NO data
  RooDataHist track_data_no("track_data", "track data", fitutil->GetObs(), tracks_no);
  RooDataHist shower_data_no("shower_data", "shower data", fitutil->GetObs(), showers_no);
  RooFitResult *fr_track_no  = pdf_tracks.fitTo( track_data_no, Save(kTRUE) );
  RooFitResult *fr_shower_no = pdf_showers.fitTo( shower_data_no, Save(kTRUE) );
  cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;
  
  RooArgSet result_track_no ( fr_track_no->floatParsFinal() );
  RooArgSet result_shower_no ( fr_shower_no->floatParsFinal() );

  // Fit under NO model, IO data
  fitutil->SetNOlims();
  fitutil->SetNOcentvals();

  RooDataHist track_data_io("track_data", "track data", fitutil->GetObs(), tracks_io);
  RooDataHist shower_data_io("shower_data", "shower data", fitutil->GetObs(), showers_io);
  RooFitResult *fr_track_io  = pdf_tracks.fitTo( track_data_io, Save(kTRUE) );
  RooFitResult *fr_shower_io = pdf_showers.fitTo( shower_data_io, Save(kTRUE) );
  cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;

  RooArgSet result_track_io ( fr_track_io->floatParsFinal() );
  RooArgSet result_shower_io ( fr_shower_io->floatParsFinal() );

  cout << "*********Fit result comparison: data NO, model IO ********" << endl;
  cout << "Track results:" << endl;
  cout << "dm31       fitted: " << ((RooRealVar*)result_track_no.find("Dm31"))->getVal() << endl;
  cout << "sinsq_th23 fitted: " << ((RooRealVar*)result_track_no.find("SinsqTh23"))->getVal() << endl;
  cout << "Shower results:" << endl;
  cout << "dm31       fitted: " << ((RooRealVar*)result_shower_no.find("Dm31"))->getVal() << endl;
  cout << "sinsq_th23 fitted: " << ((RooRealVar*)result_shower_no.find("SinsqTh23"))->getVal() << endl;
  cout << "*********Fit result comparison****************************" << endl;

  cout << "*********Fit result comparison: data IO, model NO ********" << endl;
  cout << "Track results:" << endl;
  cout << "dm31       fitted: " << ((RooRealVar*)result_track_io.find("Dm31"))->getVal() << endl;
  cout << "sinsq_th23 fitted: " << ((RooRealVar*)result_track_io.find("SinsqTh23"))->getVal() << endl;
  cout << "Shower results:" << endl;
  cout << "dm31       fitted: " << ((RooRealVar*)result_shower_io.find("Dm31"))->getVal() << endl;
  cout << "sinsq_th23 fitted: " << ((RooRealVar*)result_shower_io.find("SinsqTh23"))->getVal() << endl;
  cout << "*********Fit result comparison****************************" << endl;
  //----------------------------------------------------------
  // set hierarchy to fitted values
  //----------------------------------------------------------

  // WARNING: THE NAMES NO AND IO IN THIS SECTION ARE TO SEPERATE THE VARIABLES, NOT TO STATE UNDER
  // WHICH ORDERING THE OBJECTS ARE USED/EVALUATED/GENERATED.
  dm31      = ((RooRealVar*)result_track_no.find("Dm31"))->getVal();
  sinSqTh23 = ((RooRealVar*)result_track_no.find("SinsqTh23"))->getVal();
  fitutil->GetVar("Dm31")->setVal( dm31 );
  fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );
  TH3D *tracks_fitted_no  = (TH3D*)pdf_tracks.GetExpValHist();
  tracks_fitted_no->SetName("tracks_fitted_no");

  dm31      = ((RooRealVar*)result_shower_no.find("Dm31"))->getVal();
  sinSqTh23 = ((RooRealVar*)result_shower_no.find("SinsqTh23"))->getVal();
  fitutil->GetVar("Dm31")->setVal( dm31 );
  fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );
  TH3D *showers_fitted_no = (TH3D*)pdf_showers.GetExpValHist();
  showers_fitted_no->SetName("showers_fitted_no");

  Double_t n_chi2tr_no = HistoChi2Test(tracks_no, tracks_fitted_no, 2, 80, -1, 0);
  Double_t n_chi2sh_no = HistoChi2Test(showers_no, showers_fitted_no, 2, 80, -1, 0);

  cout << "NMHUtils: Chi2 between tracks  NO and tracks  fitted on IO is: " << n_chi2tr_no << endl;
  cout << "NMHUtils: Chi2 between showers NO and showers fitted on IO is: " << n_chi2sh_no << endl;
  cout << "Squared sum is : " << std::sqrt(std::pow(n_chi2tr_no, 2) + std::pow(n_chi2sh_no, 2)) << endl;

  dm31      = ((RooRealVar*)result_track_io.find("Dm31"))->getVal();
  sinSqTh23 = ((RooRealVar*)result_track_io.find("SinsqTh23"))->getVal();
  fitutil->GetVar("Dm31")->setVal( dm31 );
  fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );
  TH3D *tracks_fitted_io = (TH3D*)pdf_tracks.GetExpValHist();
  tracks_fitted_io->SetName("tracks_fitted_io");

  dm31      = ((RooRealVar*)result_shower_io.find("Dm31"))->getVal();
  sinSqTh23 = ((RooRealVar*)result_shower_io.find("SinsqTh23"))->getVal();
  fitutil->GetVar("Dm31")->setVal( dm31 );
  fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );
  TH3D *showers_fitted_io = (TH3D*)pdf_showers.GetExpValHist();
  showers_fitted_io->SetName("showers_fitted_io");

  Double_t n_chi2tr_io = HistoChi2Test(tracks_io, tracks_fitted_io, 2, 80, -1, 0);
  Double_t n_chi2sh_io = HistoChi2Test(showers_io, showers_fitted_io, 2, 80, -1, 0);

  cout << "NMHUtils: Chi2 between tracks  IO and tracks  fitted on NO is: " << n_chi2tr_io << endl;
  cout << "NMHUtils: Chi2 between showers IO and showers fitted on NO is: " << n_chi2sh_io << endl;
  cout << "Squared sum is : " << std::sqrt(std::pow(n_chi2tr_io, 2) + std::pow(n_chi2sh_io, 2)) << endl;

  //----------------------------------------------------------
  // print asymmetry 
  //----------------------------------------------------------

  auto asym_track  = NMHUtils::Asymmetry( (TH2D*)tracks_no ->Project3D("yx"), (TH2D*)tracks_io ->Project3D("yx"), "asymmetry_track" , 2, 80, -1, 0);
  auto asym_shower = NMHUtils::Asymmetry( (TH2D*)showers_no->Project3D("yx"), (TH2D*)showers_io->Project3D("yx"), "asymmetry_shower", 2, 80, -1, 0);
  auto asym_val_track  = std::get<1>(asym_track);
  auto asym_val_shower = std::get<1>(asym_shower);
  cout << "Asym track : " << asym_val_track << endl;
  cout << "Asym shower: " << asym_val_shower << endl;
  cout << "Squared sum is : " << std::sqrt(std::pow(asym_val_track, 2) + std::pow(asym_val_shower, 2)) << endl;
}
