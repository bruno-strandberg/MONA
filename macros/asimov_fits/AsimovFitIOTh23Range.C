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

/* Script to calculate the asimov sensitivity at the PDG central values under the assumption
 * that Nature is IO. The script moves `\Theta_{23}` over the range [40, 50] in steps of 1
 * and saves results in csv and root files.
 */

void AsimovFitIOTh23Range() {

  TString filefolder   = "./default_detres/";
  TString s_outputfile = "output/csv/AsimovFitIOTh23Range.csv";
  TString s_rootfile   = "output/root/AsimovFitIOTh23Range.root";

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

  // Open root file to save histograms
  TFile fout(s_rootfile, "RECREATE");

  // Open output stream to save sensitivity values
  ofstream outputfile(s_outputfile);
  outputfile << "th23,sinSqTh23,n_chi2tr_io,n_chi2sh_io," << endl;

  for (Int_t i = 0; i < 11; i++) {
    FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);

    FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_response);
    FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response);
    Double_t th23 = 40 + i;
    Double_t sinSqTh23_true = TMath::Power(TMath::Sin(th23 * TMath::Pi()/180.), 2);
 
    // Set values to NO
    fitutil->SetNOcentvals();
    fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23_true );

    TH3D* tracks_true  = (TH3D*)pdf_tracks.GetExpValHist();
    TH3D* showers_true = (TH3D*)pdf_showers.GetExpValHist();
    tracks_true->SetName("tracks_expval_true");
    showers_true->SetName("showers_expval_true");

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

    // Fit under IO model, NO data
    std::map<string, TH1*> hist_map = { { (string)tracks_true->GetName(),  tracks_true },
                                        { (string)showers_true->GetName(), showers_true }};

    SetIOlimsChi2Fit(fitutil);
    fitutil->SetIOcentvals();

    RooCategory cats("categories","data categories");
    cats.defineType( tracks_true->GetName() );
    cats.defineType( showers_true->GetName() );

    RooSimultaneous simPdf("simPdf", "simultaneous Pdf for NO", cats);
    simPdf.addPdf(pdf_tracks,  tracks_true->GetName() );
    simPdf.addPdf(pdf_showers, showers_true->GetName() );

    RooDataHist data_hists("data_hists", "track and shower data", fitutil->GetObs(), cats, hist_map);

    // Fit in both quadrants to find the real minimum of Th23.
    fitutil->SetIOcentvals();
    fitutil->GetVar("SinsqTh23")->setVal(0.4);
    RooFitResult *fitres_1q = simPdf.chi2FitTo( data_hists, Save(), Range("firstq"), DataError(RooAbsData::Poisson) );
    RooArgSet result_1q ( fitres_1q->floatParsFinal() );

    fitutil->SetIOcentvals();
    fitutil->GetVar("SinsqTh23")->setVal(0.6);
    RooFitResult *fitres_2q = simPdf.chi2FitTo( data_hists, Save(), Range("secondq"), DataError(RooAbsData::Poisson) );
    RooArgSet result_2q ( fitres_2q->floatParsFinal() );

    RooArgSet *result;
    Double_t fitChi2_1q = fitres_1q->minNll();
    Double_t fitChi2_2q = fitres_2q->minNll();
    cout << "first q " << TMath::Sqrt( fitChi2_1q ) << endl;
    cout << "second q " << TMath::Sqrt( fitChi2_2q ) << endl;
    if (fitChi2_1q == fitChi2_2q) cout << "NOTICE: Minimizer found same minimum for both quadrants." << endl;
    if (fitChi2_1q < fitChi2_2q) result = &result_1q;
    else                         result = &result_2q;

    cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;

    cout << "*********Fit result comparison****************************" << endl;
    cout << "dm31       fitted: " << ((RooRealVar*)result->find("Dm31"))->getVal() << endl;
    cout << "sinsq_th23 fitted: " << ((RooRealVar*)result->find("SinsqTh23"))->getVal() << endl;
    cout << "*********Fit result comparison****************************" << endl;

    //----------------------------------------------------------
    // set hierarchy to fitted values
    //----------------------------------------------------------

    Double_t dm31      = ((RooRealVar*)result->find("Dm31"))->getVal();
    Double_t sinSqTh23 = ((RooRealVar*)result->find("SinsqTh23"))->getVal();
    fitutil->GetVar("Dm31")->setVal( dm31 );
    fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );
    TH3D *tracks_fitted  = (TH3D*)pdf_tracks.GetExpValHist();
    tracks_fitted->SetName("tracks_fitted");
    TH3D *showers_fitted = (TH3D*)pdf_showers.GetExpValHist();
    showers_fitted->SetName("showers_fitted");

    std::tuple<TH1*, Double_t, Double_t> n_chi2tr = NMHUtils::Asymmetry(tracks_true, tracks_fitted, "sensitivity_track",
                                                        fitEMin, fitEMax, fitctMin, fitctMax);
    std::tuple<TH1*, Double_t, Double_t> n_chi2sh = NMHUtils::Asymmetry(showers_true, showers_fitted, "sensitivity_shower",
                                                        fitEMin, fitEMax, fitctMin, fitctMax);

    // Write the histograms containing the expectation values and sensitivity for tracks and showers
    TH2D* h_track  = (TH2D*)std::get<0>(n_chi2tr);
    TH2D* h_shower = (TH2D*)std::get<0>(n_chi2sh);
    
    fout.cd();
    tracks_true->Write(    Form("track_expval_true_%.0f", th23));
    showers_true->Write(   Form("shower_expval_true_%.0f", th23));
    tracks_fitted->Write(  Form("track_expval_fitted_%.0f", th23));
    showers_fitted->Write( Form("shower_expval_fitted_%.0f", th23));
    h_track->Write(        Form("sensitivity_track_%.0f", th23));
    h_shower->Write(       Form("sensitivity_shower_%.0f", th23));

    // Get the total values of sensitivities
    Double_t chi2tr = std::get<1>(n_chi2tr);
    Double_t chi2sh = std::get<1>(n_chi2sh);

    cout << "NMHUtils: Chi2 between tracks  NO and tracks  fitted on IO is: " << chi2tr << endl;
    cout << "NMHUtils: Chi2 between showers NO and showers fitted on IO is: " << chi2sh << endl;
    cout << "Squared sum is : " << std::sqrt(std::pow(chi2tr, 2) + std::pow(chi2sh, 2)) << endl;

    //----------------------------------------------------------
    // save fit results to file
    //----------------------------------------------------------
    outputfile << th23 << "," << sinSqTh23_true << "," << chi2tr << "," << chi2sh << "," << endl;
  }
  outputfile.close();
  fout.Close();
}
