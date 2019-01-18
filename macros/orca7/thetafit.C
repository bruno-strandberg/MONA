#include "NMHUtils.h"


Double_t CustomEnergy(SummaryEvent* evt) {

  if ( evt->Get_shower_ql0() > 0.5 ) { return evt->Get_shower_energy(); }
  else                               { return evt->Get_track_energy();  }

}

TVector3 CustomDir(SummaryEvent *evt) { return evt->Get_track_dir();      }
TVector3 CustomPos(SummaryEvent *evt) { return evt->Get_track_pos();      }
Double_t CustomBY (SummaryEvent *evt) { return evt->Get_track_bjorkeny(); }

//******************************************************************************************

using namespace RooFit;

void thetafit(TString dataf  = "../../data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root", 
	      TString effmf  = "../../data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root",
	      Bool_t redoresp = kFALSE) {

  TString respn  = "thetafit_resp.root";

  if (redoresp) {
    system("rm " + respn);
  }

  Double_t muoncut  = 0.05;
  Double_t noisecut = 0.01;
  Double_t trackcut = 0.8;   // try to select a relatively pure sample

  //------------------------------------------------------------
  // initialise the response for tracks with custom reco to use shower energy
  //------------------------------------------------------------

  DetResponse track(DetResponse::customreco, "track", 20, 1, 100, 40, -1, 1, 1, 0, 1);
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
  // create fitutil and fitpdf; fix all parameters except th23 and th13
  //------------------------------------------------------------
  
  FitUtil futil(1, track.GetHist3D(), 3, 80, -1, 0, 0, 1, effmf);
  FitPDF  trkpdf("fitpdf","fitpdf", &futil, &track);
  futil.SetNOlims();

  futil.GetVar("SinsqTh12")->setConstant(kTRUE);
  futil.GetVar("SinsqTh13")->setConstant(kTRUE);
  futil.GetVar("dcp")->setConstant(kTRUE);
  futil.GetVar("Dm21")->setConstant(kTRUE);
  futil.GetVar("Dm31")->setConstant(kTRUE);

  //------------------------------------------------------------
  // randomize th13 and th23, create pseudoexperiment, revert to central values
  //------------------------------------------------------------

  TRandom3 frand(0);
  RooRandom::randomGenerator()->SetSeed( frand.Uniform(1,1e6) );
  Double_t th23;

  futil.GetVar("SinsqTh23")->randomize();
  th23 = futil.GetVar("SinsqTh23")->getVal();

  TH3D* experiment = (TH3D*)trkpdf.GetExpValHist()->Clone("experiment");
  experiment->SetDirectory(0);

  futil.SetNOcentvals(); // set back to central values and expand th23 limits

  TH3D* expectation = (TH3D*)trkpdf.GetExpValHist()->Clone("expectation");
  expectation->SetDirectory(0);

  //------------------------------------------------------------
  // import data to RooFit and fit
  //------------------------------------------------------------
  
  RooDataHist rfdata("rfdata","rfdata", futil.GetObs(), Import(*experiment) );
  RooFitResult *fitres = trkpdf.fitTo( rfdata, Save(kTRUE), SumW2Error(kFALSE) );
  cout << "Theta-23, experiment and fitted: " << th23 << "\t" 
       << ((RooRealVar*)fitres->floatParsFinal().find("SinsqTh23"))->getVal() << endl;

  cout << "*****************************************************************************************" << endl;
  cout << "                                BEZERK MODE ACTIVATED " << endl;
  cout << "*****************************************************************************************" << endl;

  // calculate the asymmetry between the experiment and the model
  auto asym = NMHUtils::Asymmetry( (TH2D*)experiment->Project3D("yx")->Clone(), 
				   (TH2D*)expectation->Project3D("yx")->Clone(), "asym" );

  //------------------------------------------------------------
  // create likelihood profile
  //------------------------------------------------------------
  RooNLLVar nll("nll","nll",trkpdf,rfdata, NumCPU(5));

  // TGraph *g1 = new TGraph();
  // for (Double_t th = 0.2; th < 0.8; th = th + 0.1) {

  //   futil.GetVar("SinsqTh23")->setVal(th23);
  //   g1->SetPoint( g1->GetN(), th, nll.getVal() );
  //   cout << th << "\t" << nll.getVal() << endl;

  // }

  // g1->Draw();
  // cout << "Experiment at theta-23 value: " << th23 << endl;

  // return;

  Double_t min = futil.GetVar("SinsqTh23")->getMin();
  Double_t max = futil.GetVar("SinsqTh23")->getMax();

  RooPlot* frame = futil.GetVar("SinsqTh23")->frame( Range(min, max),Title("-log(L) scan vs sinsqth23") );
  nll.plotOn(frame, LineColor(kRed), ShiftToZero() ) ;


  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->Divide(2,2);
  c1->cd(1);
  frame->Draw();
  c1->cd(2);
  experiment->Project3D("yx")->Draw("colz");
  c1->cd(3);
  expectation->Project3D("yx")->Draw("colz");
  c1->cd(4);
  std::get<0>(asym)->Draw("colz");



}
