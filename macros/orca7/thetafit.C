#include "NMHUtils.h"

//******************************************************************************************
// functions to  be able to use the shower energy and track direction for high-purity tracks
//******************************************************************************************
Double_t CustomEnergy(SummaryEvent* evt) {

  if ( evt->Get_shower_ql0() > 0.5 ) { return evt->Get_shower_energy(); }
  else                               { return evt->Get_track_energy();  }

}

TVector3 CustomDir(SummaryEvent *evt) { return evt->Get_track_dir();      }
TVector3 CustomPos(SummaryEvent *evt) { return evt->Get_track_pos();      }
Double_t CustomBY (SummaryEvent *evt) { return evt->Get_track_bjorkeny(); }

//******************************************************************************************

using namespace RooFit;

/* 
   This macro fits a pseudo-experiment for theta-23. In creates:
   1. the expectation value plot
   2. the plot of the pseudo-experiment
   3. a plot depicting a likelihood scan in theta-23
   4. plots that compare the model and the pseudo-data at different theta-23 value. It creates many 1D plots in energy for different cos-theta values.

*/
void thetafit() {

  TString dataf     = (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root";
  TString effmf     = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root";
  Bool_t   redoresp = kFALSE;

  system("mkdir -p rootfiles");
  TString respn  = "rootfiles/thetafit_resp.root";

  if (redoresp) {
    system("rm " + respn);
  }

  Double_t muoncut  = 0.05;
  Double_t noisecut = 0.01;
  Double_t trackcut = 0.8;   // try to select a relatively pure sample

  Int_t    runtime  = 1;  // runtime 1 year
  Double_t nebins   = 20; // number of energy bins in the range 1-100
  Double_t nctbins  = 40; // number of cos-theta bins in the range -1 to 1
  Double_t nbybins  =  1; // number of bjorken-y bins

  //------------------------------------------------------------
  // initialise the response for tracks with custom reco to use shower energy
  //------------------------------------------------------------

  DetResponse track(DetResponse::customreco, "track", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1);
  track.SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );

  track.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,      0.5, true);
  track.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), noisecut, true);
  track.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(),  muoncut, true);
  track.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , trackcut, true);

  //------------------------------------------------------------
  // fill the response
  //------------------------------------------------------------

  if ( !NMHUtils::FileExists( respn ) ) {

    SummaryParser sp(dataf);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      track.Fill( sp.GetEvt(i) );
    }

    track.WriteToFile( respn );

  }
  else {
    track.ReadFromFile( respn );
  }

  cout << "NOTICE thetafit() response ready" << endl;

  //------------------------------------------------------------
  // create fitutil and fitpdf; fix all parameters except th23 to NO central values
  //------------------------------------------------------------
  
  FitUtil futil(runtime, track.GetHist3D(), 3, 80, -1, 0, 0, 1, effmf);
  FitPDF  trkpdf("fitpdf","fitpdf", &futil, &track);

  futil.GetVar("SinsqTh12")->setConstant(kTRUE);
  futil.GetVar("SinsqTh13")->setConstant(kTRUE);
  futil.GetVar("dcp")->setConstant(kTRUE);
  futil.GetVar("Dm21")->setConstant(kTRUE);
  futil.GetVar("Dm31")->setConstant(kTRUE);

  //------------------------------------------------------------
  // randomize th23, create pseudoexperiment
  //------------------------------------------------------------

  TRandom3 frand(0);
  Double_t th23 = frand.Uniform(0.25, 0.75);
  futil.GetVar("SinsqTh23")->setVal( th23 );

  // experiment
  TH3D* experiment  = trkpdf.SimplePseudoExp("experiment",kTRUE);
  TH3D* expectation = (TH3D*)trkpdf.GetExpValHist()->Clone("expectation");
  expectation->SetDirectory(0);

  //------------------------------------------------------------
  // manipulate theta-23 to be able to fit in a wider range and specifically
  // in first and second quarter.
  //------------------------------------------------------------
  futil.GetVar("SinsqTh23")->setMin(0);
  futil.GetVar("SinsqTh23")->setMax(1);
  futil.GetVar("SinsqTh23")->setRange("firstq" , 0. , 0.5  );
  futil.GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.   );

  //------------------------------------------------------------
  // import data to RooFit and fit in first quarter and second quarter
  //------------------------------------------------------------
  
  RooDataHist rfdata("rfdata","rfdata", futil.GetObs(), Import(*experiment) );

  futil.GetVar("SinsqTh23")->setVal(0.4);
  RooFitResult *fitres_q1 = trkpdf.fitTo( rfdata, Save(kTRUE), SumW2Error(kFALSE), Range("firstq" ) );
  futil.GetVar("SinsqTh23")->setVal(0.6);
  RooFitResult *fitres_q2 = trkpdf.fitTo( rfdata, Save(kTRUE), SumW2Error(kFALSE), Range("secondq") );

  RooFitResult *fitres = fitres_q1;
  if ( fitres_q2->minNll() < fitres_q1->minNll() ) fitres = fitres_q2;
  
  cout << "*****************************************************************************************" << endl;
  cout << "Theta-23, experiment and fitted: " << th23 << "\t" 
       << ((RooRealVar*)fitres->floatParsFinal().find("SinsqTh23"))->getVal() << endl;

  cout << "Min LLH and best fit in first quarter : " << fitres_q1->minNll() << "'\t"
       << ((RooRealVar*)fitres_q1->floatParsFinal().find("SinsqTh23"))->getVal() << endl;

  cout << "Min LLH and best fit in second quarter: " << fitres_q2->minNll() << "'\t"
       << ((RooRealVar*)fitres_q2->floatParsFinal().find("SinsqTh23"))->getVal() << endl;
  cout << "*****************************************************************************************" << endl;

  //------------------------------------------------------------
  // create likelihood profile
  //------------------------------------------------------------
  RooNLLVar nll("nll","nll",trkpdf,rfdata, NumCPU(5));

  RooPlot* frame = futil.GetVar("SinsqTh23")->frame( Range(0.3, 0.7), Title("-log(L) scan vs sinsqth23") );
  nll.plotOn(frame, LineColor(kRed), ShiftToZero() ) ;

  //------------------------------------------------------------
  // draw the profile; the experiment (at sampled theta value) data; the expectation at central theta value
  // (fit start point)
  //------------------------------------------------------------

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->Divide(2,2);
  c1->cd(1);
  expectation->GetYaxis()->SetRangeUser(-1,0);
  expectation->Project3D("yx")->Draw("colz");
  c1->cd(2);
  experiment->GetYaxis()->SetRangeUser(-1,0);
  experiment->Project3D("yx")->Draw("colz");
  c1->cd(3);
  frame->Draw();

  //------------------------------------------------------------
  // for visualisation purposes, compare the pseudo-experiment with a model with different th23
  //------------------------------------------------------------
  if ( th23 <= 0.5 ) { futil.GetVar("SinsqTh23")->setVal( frand.Uniform(0.55, 0.75) ); }
  else               { futil.GetVar("SinsqTh23")->setVal( frand.Uniform(0.25, 0.45) ); }

  TH2D* expected_2D   = (TH2D*)trkpdf.GetExpValHist()->Project3D("yx");
  TH2D* experiment_2D = (TH2D*)experiment->Project3D("yx");

  // for some reason ROOT does not calculate the error when I call sumw2 on experiment_2D, so I do this manually
  for (Int_t ebin = 1; ebin <=experiment_2D->GetXaxis()->GetNbins(); ebin++ ) {
    for (Int_t ctbin = 1; ctbin <=experiment_2D->GetYaxis()->GetNbins(); ctbin++ ) {
      experiment_2D->SetBinError(ebin, ctbin, TMath::Sqrt( experiment_2D->GetBinContent(ebin, ctbin) ) );
    }
  }  

  vector<TH1D*> proj_expected, proj_experiment;

  for (Int_t ctbin = 1; ctbin <=experiment_2D->GetYaxis()->GetNbins(); ctbin++ ) {

    Double_t ct = experiment_2D->GetYaxis()->GetBinCenter(ctbin);
    if ( ct >= 0 ) break;

    TString sfx = "ct=" + (TString)to_string(ct);
    TString title = "#theta_{23}^{model}=" + (TString)to_string( futil.GetVar("SinsqTh23")->getVal() ) + 
      "__#theta_{23}^{exp}=" + (TString)to_string( th23 ) + "__" + sfx; 

    proj_expected.push_back( expected_2D->ProjectionX("model_" + sfx, ctbin, ctbin) );
    proj_experiment.push_back( experiment_2D->ProjectionX("experiment_" + sfx, ctbin, ctbin) );
    proj_expected.back()->SetTitle( title );
    proj_experiment.back()->SetTitle( title );

  }

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->DivideSquare( proj_expected.size() );
  
  for (Int_t i = 0; i < proj_expected.size(); i++) {
    c2->cd(i+1);
    

    proj_experiment[i]->SetLineColor(kBlue);
    proj_experiment[i]->SetMarkerColor(kBlue);
    proj_experiment[i]->SetMarkerStyle(7);

    proj_expected[i]->SetLineColor(kRed);
    proj_expected[i]->SetMarkerColor(kRed);
    proj_expected[i]->SetMarkerStyle(4);

    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->AddEntry(proj_experiment[i], "experiment", "lep");
    leg->AddEntry(proj_expected[i], "model", "lep");

    proj_experiment[i]->Draw();
    proj_expected[i]->Draw("sameHISTE");
    leg->Draw();
  }

}
