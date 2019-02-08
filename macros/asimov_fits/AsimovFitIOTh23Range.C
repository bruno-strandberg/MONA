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

void AsimovFitIOTh23Range() {

  TString filefolder = "./default_detres/";
  TString s_outputfile = "./AsimovFitIOTh23Range.txt";

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

  // Open output stream to save sensitivity values
  ofstream outputfile(s_outputfile);
  outputfile << "th23,sinSqTh23,n_chi2tr_io,n_chi2sh_io" << endl;

  for (Int_t i = 0; i < 11; i++) {
    FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(), 1, 100, -1, 0, 0, 1, meff_file);

    FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_response);
    FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response);
    Double_t th23 = 40 + i;
    Double_t sinSqTh23_true = TMath::Power(TMath::Sin(th23 * TMath::Pi()/180.), 2);
 
    // Set values to NO
    fitutil->SetNOlims();
    fitutil->SetNOcentvals();
    fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23_true );

    TH3D* tracks_no  = (TH3D*)pdf_tracks.GetExpValHist();
    TH3D* showers_no = (TH3D*)pdf_showers.GetExpValHist();
    tracks_no->SetName("tracks_expval_NO");
    showers_no->SetName("showers_expval_NO");

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
    std::map<string, TH1*> hist_map_no = { {(string)tracks_no->GetName(),  tracks_no },
                                           {(string)showers_no->GetName(), showers_no }};

    fitutil->SetIOlims();
    fitutil->SetIOcentvals();

    RooCategory cats_no("categories","data categories"); // I love cats :3
    cats_no.defineType( tracks_no->GetName() );
    cats_no.defineType( showers_no->GetName() );

    RooSimultaneous simPdf_no("simPdf_no", "simultaneous Pdf for NO", cats_no);
    simPdf_no.addPdf(pdf_tracks,  tracks_no->GetName() );
    simPdf_no.addPdf(pdf_showers, showers_no->GetName() );

    RooDataHist data_hists_no("data_hists", "track and shower data", fitutil->GetObs(), cats_no, hist_map_no);
    RooFitResult *fitres_no = simPdf_no.fitTo( data_hists_no, Save(kTRUE) );
    cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;

    RooArgSet result_no ( fitres_no->floatParsFinal() );

    cout << "*********Fit result comparison****************************" << endl;
    cout << "dm31       fitted: " << ((RooRealVar*)result_no.find("Dm31"))->getVal() << endl;
    cout << "sinsq_th23 fitted: " << ((RooRealVar*)result_no.find("SinsqTh23"))->getVal() << endl;
    cout << "*********Fit result comparison****************************" << endl;

    //----------------------------------------------------------
    // set hierarchy to fitted values
    //----------------------------------------------------------

    Double_t dm31      = ((RooRealVar*)result_no.find("Dm31"))->getVal();
    Double_t sinSqTh23 = ((RooRealVar*)result_no.find("SinsqTh23"))->getVal();
    fitutil->GetVar("Dm31")->setVal( dm31 );
    TH3D *tracks_fitted_no  = (TH3D*)pdf_tracks.GetExpValHist();
    TH3D *showers_fitted_no = (TH3D*)pdf_showers.GetExpValHist();
    tracks_fitted_no->SetName("tracks_fitted_no");
    showers_fitted_no->SetName("showers_fitted_no");

    Double_t n_chi2tr_no = HistoChi2Test(tracks_no, tracks_fitted_no, 2, 80, -1, 0);
    Double_t n_chi2sh_no = HistoChi2Test(showers_no, showers_fitted_no, 2, 80, -1, 0);

    cout << "NMHUtils: Chi2 between tracks  NO and tracks  fitted on IO is: " << n_chi2tr_no << endl;
    cout << "NMHUtils: Chi2 between showers NO and showers fitted on IO is: " << n_chi2sh_no << endl;
    cout << "Squared sum is : " << std::sqrt(std::pow(n_chi2tr_no, 2) + std::pow(n_chi2sh_no, 2)) << endl;

    //----------------------------------------------------------
    // save fit results to file
    //----------------------------------------------------------
    outputfile << th23 << "," << sinSqTh23_true << "," << n_chi2tr_no << "," << n_chi2sh_no << endl;
  }
  outputfile.close();
}
