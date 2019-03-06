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
 * that Nature is NO. The script uses 5*n PID bins, where the bins are evenly spaced in both channels:
 * q = (0, 0.2) to q = (0.2, 0.4), etc...
 * The script moves `\Theta_{23}` over the range [40, 50] in steps of 1 and saves results in 
 * csv and root files.
 */

void AsimovFitNBinsNOTh23Range() {

  const int N_PID_CLASSES = 5;
  const Double_t PID_CUT = 0.6;
  const Double_t PID_STEP = 1 / float(N_PID_CLASSES);
  const Double_t PID_EDGE = PID_CUT * N_PID_CLASSES;


  std::map<Int_t, Double_t> pid_map;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    pid_map.insert(std::make_pair(i, i*PID_STEP));
  }
  TString filefolder = TString::Format("./pid_detres/pid_binning_%i/", N_PID_CLASSES);
  TString s_outputfile = "output/csv/AsimovFitNBinsNOTh23Range.txt";
  TString s_rootfile   = "output/root/AsimovFitNBinsNOTh23Range.root";

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
  std::vector<DetResponse> track_response_vector;
  std::vector<DetResponse> shower_response_vector;

  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }

    DetResponse track_response(DetResponse::track, Form("track_response_%.2f", pid_map[i]), 
                               EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , pid_map[i]  , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), pid_map[i+1], true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
    track_response_vector.push_back(track_response);

    DetResponse shower_response(DetResponse::shower, Form("shower_response_%.2f", pid_map[i]), 
                                EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , pid_map[i]  , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), pid_map[i+1], true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
    shower_response_vector.push_back(shower_response);
  }

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------

  // Start at true and if any one of the files is missing, become false
  Bool_t files_exist = kTRUE; 
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    TString track_file  = Form("track_response_%.2f.root" , pid_map[i]);
    TString shower_file = Form("shower_response_%.2f.root", pid_map[i]);
    Bool_t track_exists  = NMHUtils::FileExists(filefolder + track_file);
    Bool_t shower_exists = NMHUtils::FileExists(filefolder + shower_file);
    files_exist = ((files_exist and track_exists) and shower_exists);
  }

  if (not files_exist) {

    auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
    SummaryParser sp(summary_file);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
      SummaryEvent *evt = sp.GetEvt(i);
      for (Int_t i = 0; i < N_PID_CLASSES; i++) {
        track_response_vector[i].Fill(evt);
        shower_response_vector[i].Fill(evt);
      }
    }

    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i].WriteToFile(  filefolder +  Form("track_response_%.2f.root" , pid_map[i]) );
      shower_response_vector[i].WriteToFile( filefolder + Form("shower_response_%.2f.root", pid_map[i]) );
    }
    cout << "NOTICE: Finished filling response" << endl;
  }
  else {
    cout << "NOTICE: Reading responses from disk" << endl;
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i].ReadFromFile(  filefolder + Form("track_response_%.2f.root" , pid_map[i]) );
      shower_response_vector[i].ReadFromFile( filefolder + Form("shower_response_%.2f.root", pid_map[i]) );
    }
  }
  
  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // set up the PDFs and static oscillation parameters
  //----------------------------------------------------------

  auto meff_file = (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";

  // Open root file to save histograms
  TFile fout(s_rootfile, "RECREATE");

  // Open output stream to save sensitivity values
  ofstream outputfile(s_outputfile);
  outputfile << "th23,sinSqTh23,n_chi2_no," << endl;

  for (Int_t i = 0; i < 11; i++) {

    FitUtil *fitutil = new FitUtil(3, track_response_vector[0].GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);

    std::vector<TH3D*> track_vector_true;
    std::vector<TH3D*> shower_vector_true;

    std::vector<FitPDF> pdf_tracks_vector;
    std::vector<FitPDF> pdf_showers_vector;

    Double_t th23 = 40 + i;
    Double_t sinSqTh23_true = TMath::Power(TMath::Sin(th23 * TMath::Pi()/180.), 2);

    for (int i = 0; i < N_PID_CLASSES; i++) {
      FitPDF pdf_tracks(  Form("pdf_tracks_%.2f", pid_map[i]),  "pdf_tracks",  fitutil, &track_response_vector[i]);
      FitPDF pdf_showers( Form("pdf_showers_%.2f", pid_map[i]), "pdf_showers", fitutil, &shower_response_vector[i]);
      pdf_tracks_vector.push_back( pdf_tracks );
      pdf_showers_vector.push_back( pdf_showers );

      // Set IO values and make expectation histograms
      fitutil->SetIOcentvals();
      fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23_true );

      track_vector_true.push_back(  (TH3D*)pdf_tracks.GetExpValHist() );
      shower_vector_true.push_back( (TH3D*)pdf_showers.GetExpValHist() );
      TString track_name_true  = Form("tracks_expval_true_%.2f", pid_map[i]);
      TString shower_name_true = Form("showers_expval_true_%.2f", pid_map[i]);
      track_vector_true[i] ->SetName(track_name_true);
      shower_vector_true[i]->SetName(shower_name_true);

      // Set parameters to constant
      fitutil->GetVar("SinsqTh12")->setConstant(kTRUE);
      fitutil->GetVar("SinsqTh13")->setConstant(kTRUE);
      fitutil->GetVar("dcp")->setConstant(kTRUE);
      fitutil->GetVar("Dm21")->setConstant(kTRUE);
    }

    // create quantiles for theta23
    fitutil->GetVar("SinsqTh23")->setRange("firstq" , 0., 0.5);
    fitutil->GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.);

    //----------------------------------------------------------
    // set up data for simultaneous fitting and fit
    //----------------------------------------------------------
    cout << "NOTICE Fitter started fitting" << endl;

    TStopwatch timer;

    // Fit under NO model, IO data
    SetNOlimsChi2Fit(fitutil); // This turned out to be necessary for more than 2 bins
    fitutil->SetNOcentvals();

    std::map<string, TH1*> hist_map;
    RooCategory cats("categories", "data categories");
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      if (pid_map[i] < PID_CUT) {
        hist_map.insert( { (string)shower_vector_true[i]->GetName(), shower_vector_true[i] } );
        cats.defineType( shower_vector_true[i]->GetName() );
        cout << "NOTICE: Added hist and cat to shower" << endl;
      }
      else {
        hist_map.insert( { (string)track_vector_true[i]->GetName(), track_vector_true[i] } );
        cats.defineType( track_vector_true[i]->GetName() );
        cout << "NOTICE: Added hist and cat to track" << endl;
      }
    }

    RooSimultaneous simPdf("simPdf", "simultaneous Pdf for IO", cats);
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      if (pid_map[i] < PID_CUT) {
        simPdf.addPdf(pdf_showers_vector[i], shower_vector_true[i]->GetName() );
        cout << "NOTICE: Added simpdf to shower" << endl;
      }
      else {
        simPdf.addPdf(pdf_tracks_vector[i],  track_vector_true[i]->GetName() );
        cout << "NOTICE: Added simpdf to track" << endl;
      }
    }

    RooDataHist data_hists("data_hists", "track and shower data", fitutil->GetObs(), cats, hist_map);

    // Fit in both quadrants to find the real minimum of Th23.
    fitutil->SetNOcentvals();
    fitutil->GetVar("SinsqTh23")->setVal(0.4);
    RooFitResult *fitres_1q = simPdf.chi2FitTo( data_hists, Save(), Range("firstq"), DataError(RooAbsData::Poisson) );
    RooArgSet result_1q ( fitres_1q->floatParsFinal() );

    fitutil->SetNOcentvals();
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

    std::vector<TH3D*> fitted;
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      if (pid_map[i] < PID_CUT) {
        TH3D *showers_fitted = (TH3D*)pdf_showers_vector[i].GetExpValHist();
        showers_fitted->SetName( Form("showers_fitted_io_%.2f", pid_map[i]) );
        fitted.push_back( showers_fitted );
      }
      else {
        TH3D *tracks_fitted = (TH3D*)pdf_tracks_vector[i].GetExpValHist();
        tracks_fitted->SetName( Form("tracks_fitted_io_%.2f", pid_map[i]) );
        fitted.push_back( tracks_fitted );
      }
    }

    std::vector< std::tuple<TH1*, Double_t, Double_t> > chi2;
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      if (pid_map[i] < PID_CUT) chi2.push_back( NMHUtils::Asymmetry( shower_vector_true[i], fitted[i], Form("sensitivity_shower_%i", i),
                                                   fitEMin, fitEMax, fitctMin, fitctMax) );
      else                      chi2.push_back( NMHUtils::Asymmetry( track_vector_true[i],  fitted[i], Form("sensitivity_track_%i", i),
                                                   fitEMin, fitEMax, fitctMin, fitctMax) );
    }

    Double_t chi2_tot = 0;
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      Double_t chi2_i = std::get<1>(chi2[i]);

      // Write the histograms containing the expectation values and sensitivity for tracks and showers
      TH2D* h_sensitivity = (TH2D*)std::get<0>(chi2[i]);
      
      fout.cd();
      if (pid_map[i] < PID_CUT) {
        shower_vector_true[i]->Write( Form("shower_expval_true_%i_%.0f", i, th23) );
      }
      else {
        track_vector_true[i]->Write(  Form("track_expval_true_%i_%.0f", i, th23) );
      }
      fitted[i]->Write(     Form("fitted_expval_%i_%.0f", i, h23) );
      h_sensitivity->Write( Form("sensitivity_%i_%.0f", i, th23) );

      cout << "NMHUtils: Chi2 between events IO and events fitted on NO is: " << chi2_i << endl;
      chi2_tot += chi2_i * chi2_i;
    }
    cout << "Squared sum is : " << std::sqrt( chi2_tot ) << endl;

    outputfile << th23 << "," << sinSqTh23_true << "," << std::sqrt( chi2_tot ) << "," << endl;
  }
  outputfile.close();
  fout.Close();
}
