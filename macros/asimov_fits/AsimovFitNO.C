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

void AsimovFitNO() {

  TString filefolder = "./default_detres/";

  // DetRes input values
  Int_t EBins = 40;
  Int_t EMin = 1;
  Int_t EMax = 100;
  Int_t ctBins = 40;
  Int_t ctMin = -1;
  Int_t ctMax = 1;
  Int_t byBins = 1;
  Int_t byMin = 0;
  Int_t byMax = 1;

  // Fitter ranges
  Int_t fitEMin  = 2;
  Int_t fitEMax  = 80;
  Int_t fitctMin = -1;
  Int_t fitctMax = 0;

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

    track_response.WriteToFile(filefolder + track_file);
    shower_response.WriteToFile(filefolder + shower_file);

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
  FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_response);
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response);

  // Set values to IO
  fitutil->SetIOcentvals();

  TH3D* tracks_io   = (TH3D*)pdf_tracks.GetExpValHist();
  TH3D* showers_io  = (TH3D*)pdf_showers.GetExpValHist();
  tracks_io->SetName("tracks_expval_IO");
  showers_io->SetName("showers_expval_IO");

  fitutil->GetVar("SinsqTh12")->setConstant(kTRUE);
  fitutil->GetVar("SinsqTh13")->setConstant(kTRUE);
  fitutil->GetVar("dcp")->setConstant(kTRUE);
  fitutil->GetVar("Dm21")->setConstant(kTRUE);

  // create quantiles for theta23
  fitutil->GetVar("SinsqTh23")->setRange("firstq" , 0., 0.5);
  fitutil->GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.);
  
  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------
  cout << "NOTICE Fitter started fitting" << endl;
  
  TStopwatch timer;

  // Fit under NO model, IO data
  std::map<string, TH1*> hist_map_io = { {(string)tracks_io->GetName(),  tracks_io },
                                         {(string)showers_io->GetName(), showers_io }};

  SetNOlimsChi2Fit(fitutil); // True for free Th23
  fitutil->SetNOcentvals();

  RooCategory cats_io("categories_io","data categories");
  cats_io.defineType( tracks_io->GetName() );
  cats_io.defineType( showers_io->GetName() );

  RooSimultaneous simPdf_io("simPdf_io", "simultaneous Pdf for IO", cats_io);
  simPdf_io.addPdf(pdf_tracks,  tracks_io->GetName() );
  simPdf_io.addPdf(pdf_showers, showers_io->GetName() );

  RooDataHist data_hists_io("data_hists_io", "track and shower data", fitutil->GetObs(), cats_io, hist_map_io);

  // Fit in both quadrants to find the real minimum of Th23.
  ResetToCentral(*fitutil);
  fitutil->GetVar("SinsqTh23")->setVal(0.4);
  RooFitResult *fitres_1q_io = simPdf_io.chi2FitTo( data_hists_io, Save(), Range("firstq"), DataError(RooAbsData::Expected) );
  RooArgSet result_1q_io ( fitres_1q_io->floatParsFinal() );

  ResetToCentral(*fitutil);
  fitutil->GetVar("SinsqTh23")->setVal(0.6);
  RooFitResult *fitres_2q_io = simPdf_io.chi2FitTo( data_hists_io, Save(), Range("secondq"), DataError(RooAbsData::Expected) );
  RooArgSet result_2q_io ( fitres_2q_io->floatParsFinal() );

  RooArgSet *result_io;
  Double_t fitChi2_1q = fitres_1q_io->minNll();
  Double_t fitChi2_2q = fitres_2q_io->minNll();
  cout << "first q" << fitChi2_1q << endl;
  cout << "second q" << fitChi2_2q << endl;
  if (fitChi2_1q == fitChi2_2q) cout << "NOTICE: Minimizer found same minimum for both quadrants." << endl;
  if (fitChi2_1q < fitChi2_2q) result_io = &result_1q_io;
  else                         result_io = &result_2q_io;

  cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;

  cout << "*********Fit result comparison****************************" << endl;
  cout << "dm31       fitted: " << ((RooRealVar*)result_io->find("Dm31"))->getVal() << endl;
  cout << "sinsq_th23 fitted: " << ((RooRealVar*)result_io->find("SinsqTh23"))->getVal() << endl;
  cout << "*********Fit result comparison****************************" << endl;

  //----------------------------------------------------------
  // set hierarchy to fitted values
  //----------------------------------------------------------

  Double_t dm31      = ((RooRealVar*)result_io->find("Dm31"))->getVal();
  Double_t sinSqTh23 = ((RooRealVar*)result_io->find("SinsqTh23"))->getVal();
  fitutil->GetVar("Dm31")->setVal( dm31 );
  fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );
  TH3D *tracks_fitted_io  = (TH3D*)pdf_tracks.GetExpValHist();
  tracks_fitted_io->SetName("tracks_fitted_io");
  TH3D *showers_fitted_io = (TH3D*)pdf_showers.GetExpValHist();
  showers_fitted_io->SetName("showers_fitted_io");

  std::tuple<TH1*, Double_t, Double_t> n_chi2tr_io = NMHUtils::Asymmetry(tracks_io, tracks_fitted_io, "sensitivity_track",
                                                        fitEMin, fitEMax, fitctMin, fitctMax);
  std::tuple<TH1*, Double_t, Double_t> n_chi2sh_io = NMHUtils::Asymmetry(showers_io, showers_fitted_io, "sensitivity_shower",
                                                        fitEMin, fitEMax, fitctMin, fitctMax);

  Double_t chi2tr_io = std::get<1>(n_chi2tr_io);
  Double_t chi2sh_io = std::get<1>(n_chi2sh_io);

  cout << "NMHUtils: Chi2 between tracks  IO and tracks  fitted on NO is: " << chi2tr_io << endl;
  cout << "NMHUtils: Chi2 between showers IO and showers fitted on NO is: " << chi2sh_io << endl;
  cout << "Squared sum is : " << std::sqrt(std::pow(chi2tr_io, 2) + std::pow(chi2sh_io, 2)) << endl;

}
