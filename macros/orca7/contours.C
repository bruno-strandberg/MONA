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
  This macro compares different event selections in PID variable for theta-23, dm31 sensitivity. MC events are split to very track-like, very shower-like and a middle class. For each selection a model (`FitPDF`) is created, which is in turn used to create an expectation value plot for each selection at a given set of oscillation parameters. Then the likelihood of the model with respect to the expectation value is scanned in theta-23, dm31. 

  A plot is drawn that depicts the likelihood profile in theta-23 if the fit was performed only using tracks, only using showers, only using the middle range, and if the fit was performed simultaneously to tracks and showers and tracks and showers and middle.

  Finally, a contour plot is created in (theta-23, dm31), using only tracks and using all three event selections.

*/

void contours(TString jobname = "contours",
	      TString dataf  = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root", 
	      TString effmf  = (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root",
	      Bool_t redoresp = kFALSE) {
  
  system("mkdir -p rootfiles");

  TString trackresp  = NMHUtils::Getcwd() + "/rootfiles/contours_trackresp.root";
  TString showerresp = NMHUtils::Getcwd() + "/rootfiles/contours_showerresp.root";
  TString midresp    = NMHUtils::Getcwd() + "/rootfiles/contours_midresp.root";

  if (redoresp) {
    system("rm " + trackresp + " " + showerresp + " " + midresp);
  }

  Double_t muoncut   = 0.05;
  Double_t noisecut  = 0.01;
  Double_t trackcut  = 0.7;   // try to select a relatively pure sample of tracks
  Double_t showercut = 0.3;   // same for showers, try to get a clean sample

  Int_t    runtime  = 1;     // runtime 1 year
  Double_t nebins   = 20;    // number of energy bins in the range 1-100
  Double_t nctbins  = 40;    // number of cos-theta bins in the range -1 to 1
  Double_t nbybins  =  1;    // number of bjorken-y bins

  //------------------------------------------------------------
  // initialise the responses
  //------------------------------------------------------------

  DetResponse track(DetResponse::customreco, "track", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1);
  track.SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );

  track.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,      0.5, true);
  track.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), noisecut, true);
  track.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(),  muoncut, true);
  track.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , trackcut, true);

  DetResponse shower(DetResponse::shower, "shower", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1);
  shower.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,      0.5, true);
  shower.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), noisecut, true);
  shower.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(),  muoncut, true);
  shower.AddCut( &SummaryEvent::Get_RDF_track_score , std::less<double>()      ,showercut, true);

  // middle response uses shower reco
  DetResponse middle(DetResponse::shower, "middle", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1);
  middle.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()      ,      0.5, true);
  middle.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>()   , noisecut, true);
  middle.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>()   ,  muoncut, true);
  middle.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>()   , trackcut, true);
  middle.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater_equal<double>(),showercut, true);

  //------------------------------------------------------------
  // fill the responses
  //------------------------------------------------------------

  if ( !NMHUtils::FileExists( trackresp ) || !NMHUtils::FileExists( showerresp ) || 
       !NMHUtils::FileExists( midresp ) ) {

    SummaryParser sp(dataf);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      track.Fill( sp.GetEvt(i) );
      shower.Fill( sp.GetEvt(i) );
      middle.Fill( sp.GetEvt(i) );
    }

    track.WriteToFile( trackresp );
    shower.WriteToFile( showerresp );
    middle.WriteToFile( midresp );

  }
  else {
    track.ReadFromFile( trackresp );
    shower.ReadFromFile( showerresp );
    middle.ReadFromFile( midresp );
  }

  cout << "NOTICE contours() responses ready" << endl;

  //------------------------------------------------------------
  // create fitutil and fitpdf and create expectation value plots
  //------------------------------------------------------------
  FitUtil futil(runtime, track.GetHist3D(), 3, 80, -1, 0, 0, 1, effmf);
  FitPDF  trkpdf("trkpdf","trkpdf", &futil, &track);
  FitPDF  shwpdf("shwpdf","shwpdf", &futil, &shower);
  FitPDF  midpdf("midpdf","midpdf", &futil, &middle);

  RooRandom::randomGenerator()->SetSeed(416); // this seed controls the randomisation of osc parameters

  futil.GetVar("SinsqTh12")->randomize();
  futil.GetVar("SinsqTh13")->randomize();
  futil.GetVar("SinsqTh23")->randomize();
  futil.GetVar("dcp")->randomize();
  futil.GetVar("Dm21")->randomize();
  futil.GetVar("Dm31")->randomize();

  cout << "**********************************************************************************" << endl;
  cout << "True th23, dm31: " << futil.GetVar("SinsqTh23")->getVal() << "\t" 
       << futil.GetVar("Dm31")->getVal() << endl;
  cout << "**********************************************************************************" << endl;

  TH3D* trk_exp = trkpdf.GetExpValHist();
  TH3D* shw_exp = shwpdf.GetExpValHist();
  TH3D* mid_exp = midpdf.GetExpValHist();
  trk_exp->SetNameTitle("trkexp","trkexp");
  shw_exp->SetNameTitle("shwexp","shwexp");
  mid_exp->SetNameTitle("midexp","midexp");

  // save the pseudo-experiment parameters in the header
  FileHeader head("contours");
  head.AddParameter( "SinsqTh12", (TString)to_string( futil.GetVar("SinsqTh12")->getVal() ) );
  head.AddParameter( "SinsqTh13", (TString)to_string( futil.GetVar("SinsqTh13")->getVal() ) );
  head.AddParameter( "SinsqTh23", (TString)to_string( futil.GetVar("SinsqTh23")->getVal() ) );
  head.AddParameter( "dcp"      , (TString)to_string( futil.GetVar("dcp")->getVal() ) );
  head.AddParameter( "Dm21"     , (TString)to_string( futil.GetVar("Dm21")->getVal() ) );
  head.AddParameter( "Dm31"     , (TString)to_string( futil.GetVar("Dm31")->getVal() ) );

  //------------------------------------------------------------
  // manipulate theta-23 and dm31 to be able to create contours in a wide range
  // other parameters are fixed to NO central values
  //------------------------------------------------------------
  futil.SetNOcentvals();

  futil.GetVar("SinsqTh23")->setMin(0);
  futil.GetVar("SinsqTh23")->setMax(1);
  futil.GetVar("Dm31")->setMin(1e-3);
  futil.GetVar("Dm31")->setMax(5e-3);

  futil.GetVar("SinsqTh12")->setConstant(kTRUE);
  futil.GetVar("SinsqTh13")->setConstant(kTRUE);
  futil.GetVar("dcp")->setConstant(kTRUE);
  futil.GetVar("Dm21")->setConstant(kTRUE);

  //------------------------------------------------------------
  // setup simultaneous fitting using 1) only trk and shw 2) trk, shw and mid
  //------------------------------------------------------------

  // this just imports individual selections for LLH profile comparisons
  RooDataHist rftrk("rftrk","trk data in roofit", futil.GetObs(), Import(*trk_exp) );
  RooDataHist rfshw("rfshw","shw data in roofit", futil.GetObs(), Import(*shw_exp) );
  RooDataHist rfmid("rfmid","mid data in roofit", futil.GetObs(), Import(*mid_exp) );
  
  // for simfit: in this case I don't use the middle one
  std::map< string, TH1* > hist_map_1 = { {trk_exp->GetName(), trk_exp},
					  {shw_exp->GetName(), shw_exp} };

  // for simfit: here I also use the middle one
  std::map< string, TH1* > hist_map_2 = { {trk_exp->GetName(), trk_exp},
					  {shw_exp->GetName(), shw_exp},
					  {mid_exp->GetName(), mid_exp} };

  RooCategory categ1("categ1","data categories, trk and shw");
  categ1.defineType( trk_exp->GetName() );
  categ1.defineType( shw_exp->GetName() );

  RooCategory categ2("categ2","data categories, trk, shw and mid");
  categ2.defineType( trk_exp->GetName() );
  categ2.defineType( shw_exp->GetName() );
  categ2.defineType( mid_exp->GetName() );

  RooDataHist comb1("comb1","combined trk and shw"     , futil.GetObs(), categ1, hist_map_1);
  RooDataHist comb2("comb2","combined trk, shw and mid", futil.GetObs(), categ2, hist_map_2);

  RooSimultaneous simPdf1("simPdf1","simultaneous pdf 1", categ1);
  simPdf1.addPdf( trkpdf, trk_exp->GetName() );
  simPdf1.addPdf( shwpdf, shw_exp->GetName() );

  RooSimultaneous simPdf2("simPdf2","simultaneous pdf 2", categ2);
  simPdf2.addPdf( trkpdf, trk_exp->GetName() );
  simPdf2.addPdf( shwpdf, shw_exp->GetName() );
  simPdf2.addPdf( midpdf, mid_exp->GetName() );

  //------------------------------------------------------------
  // create likelihood scans
  //------------------------------------------------------------

  RooPlot* frame = futil.GetVar("SinsqTh23")->frame( Range(0.25, 0.75), Title("-log(L) scan vs sinsqth23") );

  RooNLLVar nll_trk("nll_trk", "trk",  trkpdf, rftrk, NumCPU(5) );
  RooNLLVar nll_shw("nll_shw", "shw",  shwpdf, rfshw, NumCPU(5) );
  RooNLLVar nll_mid("nll_mid", "mid",  midpdf, rfmid, NumCPU(5) );
  RooNLLVar nll_1("nll_1", "trk and shw"  , simPdf1, comb1, NumCPU(5) );
  RooNLLVar nll_2("nll_2", "trk, shw, mid", simPdf2, comb2, NumCPU(5) );

  nll_trk.plotOn(frame, LineColor(kRed)   , ShiftToZero(), Name("trks") ) ;
  nll_shw.plotOn(frame, LineColor(kBlue)  , ShiftToZero(), Name("shws") ) ;
  nll_mid.plotOn(frame, LineColor(kGreen) , ShiftToZero(), Name("mid") ) ;
  nll_1.plotOn(frame  , LineColor(kBlack) , ShiftToZero(), LineStyle(9), Name("trkshw") );
  nll_2.plotOn(frame  , LineColor(kYellow), ShiftToZero(), LineStyle(9), Name("trkshwmid") );

  TLegend *leg1 = new TLegend(0.3, 0.6, 0.7, 0.9);
  leg1->AddEntry( frame->findObject("trks")     , "trk"   , "l" );
  leg1->AddEntry( frame->findObject("shws")     , "shw"   , "l" );
  leg1->AddEntry( frame->findObject("mid")      , "mid"   , "l" );
  leg1->AddEntry( frame->findObject("trkshw")   , "trkshw", "l" );
  leg1->AddEntry( frame->findObject("trkshwmid"), "all"   , "l" );
  leg1->SetLineWidth(0);
  leg1->SetFillStyle(0);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  frame->Draw();
  leg1->Draw();
  
  //------------------------------------------------------------
  // Create contour plots. I noticed that several RooMinuit instances do not work together, I am guessing there
  // are some static variables for the FORTRAN interface. This is solved by dynamic memory allocation.
  //------------------------------------------------------------

  RooMinuit *min0 = new RooMinuit(nll_trk);
  RooPlot* contplot_0 = min0->contour( *(futil.GetVar("SinsqTh23")), *(futil.GetVar("Dm31")), 2);
  delete min0;

  RooMinuit *min2 = new RooMinuit(nll_2);
  RooPlot* contplot_2 = min2->contour( *(futil.GetVar("SinsqTh23")), *(futil.GetVar("Dm31")), 2);
  delete min2;

  // clone the graphs from RooPlot for easier draw option manipulation
  TGraph *g_cont0 = (TGraph*)contplot_0->getObject(1)->Clone("g_cont0");
  TGraph *g_cont2 = (TGraph*)contplot_2->getObject(1)->Clone("g_cont2");

  g_cont0->SetLineColor(kRed);
  g_cont0->SetMarkerColor(kRed);
  g_cont0->SetLineWidth(2);

  g_cont2->SetLineColor(kBlue);
  g_cont2->SetMarkerColor(kBlue);
  g_cont2->SetLineWidth(2);
  g_cont2->SetLineStyle(9);

  g_cont0->GetXaxis()->SetTitle("sin^2#theta_{23}");
  g_cont0->GetYaxis()->SetTitle("#Delta m_{31}^2");
  g_cont0->GetYaxis()->SetRangeUser(0.5*1e-3, 4.5*1e-3);

  TLegend *leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg2->AddEntry(g_cont0, "tracks", "l");
  leg2->AddEntry(g_cont2, "all", "l");
  leg2->SetLineWidth(0);
  leg2->SetFillStyle(0);

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->SetTicks();
  g_cont0->Draw("AL");
  g_cont2->Draw("sameL");
  leg2->Draw();

  //------------------------------------------------------------
  // write the plots to output for further manipulation
  //------------------------------------------------------------
  TFile fout("rootfiles/contours.root","RECREATE");
  head.WriteHeader(&fout);
  g_cont0->Write();
  g_cont2->Write();
  frame->Write("llhplot");
  fout.Close();

}
