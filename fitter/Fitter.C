#include "Fitter.h"

// fitter headers
#include "FitUtil.h"
#include "FitPDF.h"

// NMH headers
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "NMHUtils.h"

// root headers
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"

// roofit headers
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

// cpp headers
#include <iostream>
#include <stdexcept>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace RooFit;
using namespace std;

int main(int argc, char **argv) {

  using namespace Fitter;

  auto pars = CommandLineParser(argc, argv);
  InitRespsAndSels();
  FillRespsAndSels(pars.simdata_file, pars.expdata_file, pars.refill_response);

  //----------------------------------------------------------
  // set up the PDFs for fitting
  //----------------------------------------------------------
  FitUtil *fitutil = new FitUtil(3, fTRres->GetHist3D(), 
				 (TString)getenv("NMHDIR") + "/data/eff_mass/EffMhists_elec_CC.root", 
				 (TString)getenv("NMHDIR") + "/data/eff_mass/EffMhists_muon_CC.root", 
				 (TString)getenv("NMHDIR") + "/data/eff_mass/EffMhists_tau_CC.root", 
				 (TString)getenv("NMHDIR") + "/data/eff_mass/EffMhists_elec_NC.root");

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, fTRres);  
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, fSHres);

  //----------------------------------------------------------
  // set up data for simultaneous fitting
  //----------------------------------------------------------
  std::map<string, TH1* > hist_map = { { (string)fTRsel->Get_SelName() , fTRsel->Get_h_E_costh_by() }, 
				       { (string)fSHsel->Get_SelName() , fSHsel->Get_h_E_costh_by() } };
  
  // create categories
  RooCategory categs("categs","data categories");
  categs.defineType( fTRsel->Get_SelName() );
  categs.defineType( fSHsel->Get_SelName() );
  
  // create combined data set
  RooDataHist combData("combData", "combined data", fitutil->GetObs(), categs, hist_map);

  // create simultaneous pdf
  RooSimultaneous simPdf("simPdf", "simultaneous Pdf", categs);
  simPdf.addPdf(pdf_tracks , fTRsel->Get_SelName() );
  simPdf.addPdf(pdf_showers, fSHsel->Get_SelName() );

  cout << "NOTICE main() started fitting" << endl;
  TStopwatch timer;

  // limit the range in cos-theta and... fit
  ( (RooRealVar*)fitutil->GetSet().find("ct_reco") )->setRange("R1", -1, 0);
  simPdf.fitTo( combData, Range("R1") );

  cout << "NOTICE main() finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;

  //FillExpectationValues(&pdf_tracks, &pdf_showers);
  Cleanup();
  if (fitutil) delete fitutil;
 
}

//**********************************************************************************

Fitter::cmdpars Fitter::CommandLineParser(int argc, char **argv) {

  using namespace Fitter;

  Fitter::cmdpars pars;

  pars.refill_response = kFALSE;
  pars.simdata_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  pars.expdata_file = (TString)getenv("NMHDIR") + "/evt_sampler/output/Experiments/Experiment_oscpars1_sample_0_NH.root";

  int c;
  while ( (c = getopt (argc, argv, "s:e:r")) != -1 ) {

    switch (c)
      {

      case 'r':
	pars.refill_response = kTRUE;
	break;

      case 'e':
	pars.expdata_file = optarg;
	break;

      case 's':
	pars.simdata_file = optarg;
	break;

      case '?':
      
	if (optopt == 'e' || optopt == 's') {
	  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	}
	else if (isprint (optopt)) {
	  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	}
	else {
	  fprintf (stderr,
		   "Unknown option character `\\x%x'.\n",
		   optopt);
	}
	throw std::invalid_argument("Command-line argument parsing failed.");
      default:
	abort ();
      }

  } // end while loop

  // note that order matters here!
  return pars;
  
}

//**********************************************************************************

void Fitter::InitRespsAndSels() {

  using namespace Fitter;

  //----------------------------------------------------------
  //initialise event selection and response for tracks
  //----------------------------------------------------------

  fTRsel = new EventSelection(EventSelection::track, "track_sel", NULL, 
			      fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  fTRres = new DetResponse(DetResponse::track, "track_resp", 
			   fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  vector< EventFilter* > T = {fTRsel, fTRres};

  for (auto n: T) {
    n->AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
    n->AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    n->AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );    
  }

  //----------------------------------------------------------
  //initialise event selection and response for showers
  //----------------------------------------------------------

  fSHsel = new EventSelection(EventSelection::shower, "shower_sel", NULL,
			      fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);
  fSHres = new DetResponse(DetResponse::shower, "shower_resp", 
			   fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  vector< EventFilter* > S = {fSHsel, fSHres};

  for (auto n: S) {
    n->AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
    n->AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
    n->AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    n->AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );
  }

}

//**********************************************************************************

void Fitter::FillRespsAndSels(TString simdata_file, TString expdata_file, Bool_t refill_response) {

  using namespace Fitter;

  if ( !fTRres || !fSHres || !fTRsel || !fSHsel ) {
    throw std::logic_error("ERROR! Fitter::FillRespsAndSels() DetResponse's and EventSelection's are not yet initialized!");
  }

  //----------------------------------------------------------
  // fill the responses
  //----------------------------------------------------------

  TString track_resp_name  = "track_response.root";
  TString shower_resp_name = "shower_response.root";

  if ( !NMHUtils::FileExists(track_resp_name) || !NMHUtils::FileExists(shower_resp_name) || refill_response ) {

    cout << "NOTICE Fitter::FillRespsAndSels() (Re)filling responses" << endl;

    SummaryParser sp(simdata_file);

    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

      if (i % 100000 == 0) cout << "Event: " << i << endl;

      SummaryEvent *evt = sp.GetEvt(i);
      fTRres->Fill(evt);
      fSHres->Fill(evt);

    }

    fTRres->WriteToFile(track_resp_name);
    fSHres->WriteToFile(shower_resp_name);

  }
  else {

    cout << "NOTICE Fitter::FillRespsAndSels() Reading in responses" << endl;

    fTRres->ReadFromFile(track_resp_name);
    fSHres->ReadFromFile(shower_resp_name);
  }

  cout << "NOTICE Fitter::FillRespsAndSels() Responses ready" << endl;

  //----------------------------------------------------------
  // fill the selections
  //----------------------------------------------------------
  SummaryParser sp(expdata_file);

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

    SummaryEvent *evt = sp.GetEvt(i);
    fTRsel->Fill(evt);
    fSHsel->Fill(evt);

  }

  cout << "NOTICE Fitter::FillRespsAndSels() Selections filled" << endl;
}

//**********************************************************************************

void Fitter::FillExpectationValues(FitPDF *track_pdf, FitPDF *shower_pdf) {

  TH3D *hdet_tracks  = (TH3D*)track_pdf->GetResponse()->GetHist3D()->Clone("detected_tracks");
  TH3D *hdet_showers = (TH3D*)shower_pdf->GetResponse()->GetHist3D()->Clone("detected_showers");
  hdet_tracks->Reset();
  hdet_showers->Reset();

  //----------------------------------------------------------
  // fill the track and shower histograms - should be EXACTLY as before..
  //----------------------------------------------------------

  double p[] = {0.297, 0.0215, 0.425, 1.38, 7.37e-5, 2.56e-3};

  TStopwatch timer;

  for (Int_t ebin = 1; ebin <= hdet_tracks->GetXaxis()->GetNbins(); ebin++) {
    for (Int_t ctbin = 1; ctbin <= hdet_tracks->GetYaxis()->GetNbins(); ctbin++) {
      for (Int_t bybin = 1; bybin <= hdet_tracks->GetZaxis()->GetNbins(); bybin++) {

        Double_t E  = hdet_tracks->GetXaxis()->GetBinCenter( ebin );
        Double_t ct = hdet_tracks->GetYaxis()->GetBinCenter( ctbin );
        Double_t by = hdet_tracks->GetZaxis()->GetBinCenter( bybin );

        double x[] = {E, ct, by};

        hdet_tracks->SetBinContent( ebin, ctbin, bybin, track_pdf->operator()(x, p) );
        hdet_showers->SetBinContent( ebin, ctbin, bybin, shower_pdf->operator()(x, p) );
      }
    }
  }

  cout << "NOTICE main() time for filling histograms " << (Double_t)timer.RealTime() << endl;

  TFile fout("expectation_values.root","RECREATE");
  hdet_tracks->Write();
  hdet_showers->Write();
  fout.Close();

}

//**********************************************************************************

void Fitter::Cleanup() {

  using namespace Fitter;

  if (fTRres) delete fTRres;
  if (fSHres) delete fSHres;
  if (fTRsel) delete fTRsel;
  if (fSHsel) delete fSHsel;

}

