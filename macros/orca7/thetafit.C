#include "NMHUtils.h"
#include "ORCA7.h"

using namespace RooFit;

/* 
   This macro fits a pseudo-experiment for theta-23. In creates:
   1. the expectation value plot
   2. the plot of the pseudo-experiment
   3. a plot depicting a likelihood scan in theta-23
   4. plots that compare the model and the pseudo-data at different theta-23 value. It creates many 1D plots in energy for different cos-theta values.

*/
void thetafit() {

  //------------------------------------------------------------
  // initialise ORCA7 to get the track response
  //------------------------------------------------------------
  ORCA7 o7(kTRUE);
  DetResponse *track = o7.fResps.back();

  //------------------------------------------------------------
  // create fitutil and fitpdf; fix all parameters except th23 to NO central values
  //------------------------------------------------------------
  
  FitUtil futil(o7.f_F_runtime, track->GetHist3D(), o7.f_F_emin, o7.f_F_emax, o7.f_F_ctmin, o7.f_F_ctmax, o7.f_F_bymin, o7.f_F_bymax, o7.fEffmF);
  FitPDF  trkpdf("fitpdf","fitpdf", &futil, track);

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
  TH3D* experiment  = trkpdf.SimplePseudoExp("experiment",kFALSE);
  TH3D* expectation = (TH3D*)trkpdf.GetExpValHist()->Clone("expectation");
  expectation->SetDirectory(0);

  //------------------------------------------------------------
  // manipulate theta-23 to be able to fit in a wider range and specifically
  // in first and second quarter.
  //------------------------------------------------------------
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

  RooPlot* frame = futil.GetVar("SinsqTh23")->frame( Range(0.25, 0.75), Title("-log(L) scan vs sinsqth23") );
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
