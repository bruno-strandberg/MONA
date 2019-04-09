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
   This macro creates and fits pseudoexperiments to extract theta-23, dm31. The fit results are written to a ROOT file. The fits are performed simultaneously in three event selections (tracks, showers, middle). For each pseudo-experiment, the output will contain a `RooFitResult` and three histograms of the pseudo-experiment. 
*/
void th23dm31(TString jobname = "job0", Int_t nexps = 5,
	      TString dataf  = (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root", 
	      TString effmf  = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root",
	      Bool_t redoresp = kFALSE) {

  system("mkdir -p rootfiles");
  TString trkrespn = NMHUtils::Getcwd() + "/rootfiles/th23dm31_trkresp.root";
  TString shwrespn = NMHUtils::Getcwd() + "/rootfiles/th23dm31_shwresp.root";
  TString midrespn = NMHUtils::Getcwd() + "/rootfiles/th23dm31_midresp.root";

  if (redoresp) {
    system("rm " + trkrespn + " " + shwrespn + " " + midrespn);
  }

  Double_t muoncut  = 0.05;
  Double_t noisecut = 0.01;
  Double_t trackcut  = 0.7;   // try to select a relatively pure sample
  Double_t showercut = 0.3;   // very shower-like events

  Int_t    runtime  = 1;     // runtime 1 year
  Double_t nebins   = 20;    // number of energy bins in the range 1-100
  Double_t nctbins  = 40;    // number of cos-theta bins in the range -1 to 1
  Double_t nbybins  =  1;    // number of bjorken-y bins

  //------------------------------------------------------------
  // initialise the responses; tracks use custom reco to use shower energy
  //------------------------------------------------------------

  DetResponse track(DetResponse::customreco, "track", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1);
  track.SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );
  track.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,      0.5, true);
  track.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , trackcut, true);

  DetResponse shower(DetResponse::shower, "shower", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1);
  shower.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,      0.5, true);
  shower.AddCut( &SummaryEvent::Get_RDF_track_score , std::less<double>()      ,showercut, true);

  // middle response uses shower reco
  DetResponse middle(DetResponse::shower, "middle", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1);
  middle.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()      ,      0.5, true);
  middle.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>()   , trackcut, true);
  middle.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater_equal<double>(),showercut, true);

  // add the same noise suppression cuts to all selections
  vector<DetResponse*> responses = { &track, &shower, &middle };
  for (auto r: responses) {
    cout << "NOTICE th23dm31() adding noise suppression cuts to " << r->GetRespName() << endl;
    r->AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), noisecut, true);
    r->AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(),  muoncut, true);
  }
  
  //------------------------------------------------------------
  // fill the responses
  //------------------------------------------------------------

  if ( !NMHUtils::FileExists( trkrespn ) || !NMHUtils::FileExists( shwrespn ) || 
       !NMHUtils::FileExists( midrespn ) ) {

    SummaryParser sp(dataf);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if ( i % 100000 == 0 ) cout << "Event: " << i << endl;
      for (auto r: responses) { r->Fill( sp.GetEvt(i) ); }
    }

    track.WriteToFile ( trkrespn );
    shower.WriteToFile( shwrespn );
    middle.WriteToFile( midrespn );

  }
  else {
    track.ReadFromFile ( trkrespn );
    shower.ReadFromFile( shwrespn );
    middle.ReadFromFile( midrespn );
  }

  cout << "NOTICE th23dm31() response ready" << endl;

  //------------------------------------------------------------
  // create fitutil and fitpdf and create pseudo-experiments;
  // select random parameter values for the pseudoexperiments
  //------------------------------------------------------------
  
  FitUtil futil(runtime, track.GetHist3D(), 3, 80, -1, 0, 0, 1, effmf);
  FitPDF  trkpdf("trkpdf","trkpdf", &futil, &track);
  FitPDF  shwpdf("shwpdf","shwpdf", &futil, &shower);
  FitPDF  midpdf("midpdf","midpdf", &futil, &middle);

  RooRandom::randomGenerator()->SetSeed(416); // this seed controls the randomisation of fit parameters
  trkpdf.SetSeed(0);                          // this seed controls pseudo-experiments; 0 means always different

  futil.SetNOlims();
  futil.GetVar("SinsqTh12")->randomize();
  futil.GetVar("SinsqTh13")->randomize();
  futil.GetVar("SinsqTh23")->randomize();
  futil.GetVar("dcp")->randomize();
  futil.GetVar("Dm21")->randomize();
  futil.GetVar("Dm31")->randomize();
  futil.FreeParLims();

  enum { EXPNAME, EXPTRK, EXPSHW, EXPMID };
  vector< std::tuple<TString, TH3D*, TH3D*, TH3D*> > experiments;

  for (Int_t i = 0; i < nexps; i++) {
    TString expname = jobname + "_pseudoexp_" + (TString)to_string(i);
    TH3D* trkexp = trkpdf.SimplePseudoExp("trk_" + expname, kTRUE);
    TH3D* shwexp = shwpdf.SimplePseudoExp("shw_" + expname, kTRUE);
    TH3D* midexp = midpdf.SimplePseudoExp("mid_" + expname, kTRUE);
    experiments.push_back( std::make_tuple(expname, trkexp, shwexp, midexp) );
  }

  // save the pseudo-experiment parameters in the header
  FileHeader head("th23dm31");
  head.AddParameter( "SinsqTh12", (TString)to_string( futil.GetVar("SinsqTh12")->getVal() ) );
  head.AddParameter( "SinsqTh13", (TString)to_string( futil.GetVar("SinsqTh13")->getVal() ) );
  head.AddParameter( "SinsqTh23", (TString)to_string( futil.GetVar("SinsqTh23")->getVal() ) );
  head.AddParameter( "dcp"      , (TString)to_string( futil.GetVar("dcp")->getVal() ) );
  head.AddParameter( "Dm21"     , (TString)to_string( futil.GetVar("Dm21")->getVal() ) );
  head.AddParameter( "Dm31"     , (TString)to_string( futil.GetVar("Dm31")->getVal() ) );

  //------------------------------------------------------------
  // manipulate theta-23 to be able to fit in a wider range and specifically
  // in first and second quarter; release dm31 range; fix other parameters to NO central values
  //------------------------------------------------------------
  futil.SetNOcentvals();

  futil.GetVar("SinsqTh23")->setRange("firstq" , 0. , 0.5  );
  futil.GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.   );

  futil.GetVar("SinsqTh12")->setConstant(kTRUE);
  futil.GetVar("SinsqTh13")->setConstant(kTRUE);
  futil.GetVar("dcp")->setConstant(kTRUE);
  futil.GetVar("Dm21")->setConstant(kTRUE);

  //------------------------------------------------------------
  // fit each pseudo-experiment; write the experiments and results to file
  //------------------------------------------------------------
  TString outname = NMHUtils::Getcwd() + "/rootfiles/" + jobname + "_th23dm31_out.root";
  TFile fout(outname,"RECREATE");
  head.WriteHeader(&fout);
  
  for (auto exp: experiments) {

    // get the experiment histograms from the tuple and set things up for simultaneous fitting
    TString expname = std::get<EXPNAME>(exp);
    TH1 *trk_exp    = std::get<EXPTRK> (exp);
    TH1 *shw_exp    = std::get<EXPSHW> (exp);
    TH1 *mid_exp    = std::get<EXPMID> (exp);

    std::map< string, TH1* > hist_map = { {trk_exp->GetName(), trk_exp},
					  {shw_exp->GetName(), shw_exp},
					  {mid_exp->GetName(), mid_exp} };

    RooCategory categ("categ","data categories, trk, shw and mid");
    categ.defineType( trk_exp->GetName() );
    categ.defineType( shw_exp->GetName() );
    categ.defineType( mid_exp->GetName() );

    RooDataHist comb("comb","combined trk, shw and mid", futil.GetObs(), categ, hist_map);
  
    RooSimultaneous simPdf("simPdf","simultaneous pdf", categ);
    simPdf.addPdf( trkpdf, trk_exp->GetName() );
    simPdf.addPdf( shwpdf, shw_exp->GetName() );
    simPdf.addPdf( midpdf, mid_exp->GetName() );

    // fitting
    futil.GetVar("SinsqTh23")->setVal(0.4);
    RooFitResult *fitres_q1 = simPdf.fitTo( comb, Save(kTRUE), SumW2Error(kFALSE), Range("firstq" ) );

    futil.GetVar("SinsqTh23")->setVal(0.6);
    RooFitResult *fitres_q2 = simPdf.fitTo( comb, Save(kTRUE), SumW2Error(kFALSE), Range("secondq") );

    // choose between two theta-23 quarters the result with smaller likelihood value
    RooFitResult *fitres = fitres_q1;
    if ( fitres_q2->minNll() < fitres_q1->minNll() ) fitres = fitres_q2;

    trk_exp->Write();
    shw_exp->Write();
    mid_exp->Write();
    fitres->Write("result_" + expname);

  }

  fout.Close();

}
