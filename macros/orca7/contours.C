#include "NMHUtils.h"
#include "ORCA7.h"

using namespace RooFit;

/*
  This macro compares different event selections in PID variable for theta-23, dm31 sensitivity. MC events are split to very track-like, very shower-like and a middle class. For each selection a model (`FitPDF`) is created, which is in turn used to create an expectation value plot for each selection at a given set of oscillation parameters. Then the likelihood of the model with respect to the expectation value is scanned in theta-23, dm31. 

  A plot is drawn that depicts the likelihood profile in theta-23 if the fit was performed only using tracks, only using showers, only using the middle range, and if the fit was performed simultaneously to tracks and showers and tracks and showers and middle.

  Finally, a contour plot is created in (theta-23, dm31), using only tracks and using all three event selections.

*/

void contours(TString jobname = "contours",
	      Bool_t redoresp = kFALSE) {

  ORCA7 o7( !redoresp );

  //------------------------------------------------------------
  // create fitutil, fitpdfs and expectation value plots
  //------------------------------------------------------------
  FitUtil futil(o7.f_F_runtime, o7.fResps.back()->GetHist3D(), o7.f_F_emin, o7.f_F_emax, o7.f_F_ctmin, o7.f_F_ctmax, o7.f_F_bymin, o7.f_F_bymax, o7.fEffmF);

  // randomise oscillation parameters
  RooRandom::randomGenerator()->SetSeed(418); // this seed controls the randomisation of osc parameters
  futil.SetNOlims();
  futil.GetVar("SinsqTh12")->randomize();
  futil.GetVar("SinsqTh13")->randomize();
  futil.GetVar("SinsqTh23")->randomize();
  futil.GetVar("dcp")->randomize();
  futil.GetVar("Dm21")->randomize();
  futil.GetVar("Dm31")->randomize();
  futil.FreeParLims();

  // initialise pdfs and create expectation value plots
  vector< std::pair<TH3D*, FitPDF*> > exp_pdf_pairs;
  for (auto R: o7.fResps) {
    TString pdfname = "pdf_" + R->Get_RespName();
    TString hname   = "exp_" + R->Get_RespName();
    FitPDF *pdf = new FitPDF(pdfname, pdfname, &futil, R);
    TH3D   *exp = pdf->GetExpValHist();
    exp->SetNameTitle( hname, hname );
    exp_pdf_pairs.push_back( make_pair(exp, pdf) );
  }
  
  // save the pseudo-experiment parameters in the header
  FileHeader head("contours");
  head.AddParameter( "SinsqTh12", (TString)to_string( futil.GetVar("SinsqTh12")->getVal() ) );
  head.AddParameter( "SinsqTh13", (TString)to_string( futil.GetVar("SinsqTh13")->getVal() ) );
  head.AddParameter( "SinsqTh23", (TString)to_string( futil.GetVar("SinsqTh23")->getVal() ) );
  head.AddParameter( "dcp"      , (TString)to_string( futil.GetVar("dcp")->getVal() ) );
  head.AddParameter( "Dm21"     , (TString)to_string( futil.GetVar("Dm21")->getVal() ) );
  head.AddParameter( "Dm31"     , (TString)to_string( futil.GetVar("Dm31")->getVal() ) );

  cout << "**********************************************************************************" << endl;
  cout << "True th23, dm31: " << futil.GetVar("SinsqTh23")->getVal() << "\t" 
       << futil.GetVar("Dm31")->getVal() << endl;
  cout << "**********************************************************************************" << endl;

  //------------------------------------------------------------
  // manipulate theta-23 and dm31 to be able to create contours in a wide range
  // other parameters are fixed to NO central values
  //------------------------------------------------------------
  futil.SetNOcentvals();
  futil.GetVar("SinsqTh12")->setConstant(kTRUE);
  futil.GetVar("SinsqTh13")->setConstant(kTRUE);
  futil.GetVar("dcp")->setConstant(kTRUE);
  futil.GetVar("Dm21")->setConstant(kTRUE);

  //------------------------------------------------------------
  // setup simultaneous fitting using 1) only trk and shw 2) trk, shw and mid
  //------------------------------------------------------------

  // for simfit 1: in this case I don't use the middle one
  //----------------------------------------------------------
  RooCategory categ1("categ1","data categories, trk and shw");
  categ1.defineType( exp_pdf_pairs.front().first->GetName() );
  categ1.defineType( exp_pdf_pairs.back().first->GetName() );

  std::map< string, TH1* > hist_map_1 = { {exp_pdf_pairs.front().first->GetName(), exp_pdf_pairs.front().first },
					  {exp_pdf_pairs.back().first ->GetName(), exp_pdf_pairs.back().first  } };

  RooDataHist comb1("comb1","combined trk and shw", futil.GetObs(), categ1, hist_map_1);  

  RooSimultaneous simPdf1("simPdf1","simultaneous pdf 1", categ1);
  simPdf1.addPdf( *(exp_pdf_pairs.front().second), exp_pdf_pairs.front().first->GetName() );
  simPdf1.addPdf( *(exp_pdf_pairs.back().second) , exp_pdf_pairs.back().first->GetName()  );

  // for simfit 2: here I also use the middle one
  //----------------------------------------------------------
  RooCategory categ2("categ2","data categories, trk, shw and mid");
  std::map< string, TH1* > hist_map_2;
  for (auto p: exp_pdf_pairs) {
    hist_map_2.insert( std::make_pair( p.first->GetName(), p.first ) );
    categ2.defineType( p.first->GetName() );
  }

  RooDataHist comb2("comb2","combined trk, shw and mid", futil.GetObs(), categ2, hist_map_2);
  
  RooSimultaneous simPdf2("simPdf2","simultaneous pdf 2", categ2);
  for (auto p: exp_pdf_pairs) {
    simPdf2.addPdf( *(p.second), p.first->GetName() );
  }

  //------------------------------------------------------------
  // create likelihood scans
  //------------------------------------------------------------

  RooPlot* frame = futil.GetVar("SinsqTh23")->frame( Range(0.25, 0.75), Title("-log(L) scan vs sinsqth23") );

  vector<RooDataHist*> rfdhs;
  vector<RooNLLVar*> nlls;
  for (auto p: exp_pdf_pairs) {
    TString hname   = "rf" + (TString)p.first->GetName();
    TString nllname = "nll_" + (TString)p.second->GetName();
    RooDataHist *rfh = new RooDataHist( hname, hname, futil.GetObs(), Import(*p.first) );
    rfdhs.push_back( rfh );
    nlls.push_back( new RooNLLVar(nllname, nllname, *p.second, *rfh, NumCPU(4)) );
  }

  RooNLLVar nll_1("nll_1", "trk and shw"  , simPdf1, comb1, NumCPU(5) );
  RooNLLVar nll_2("nll_2", "trk, shw, mid", simPdf2, comb2, NumCPU(5) );

  Int_t i = 0;
  for (auto nll: nlls) {
    nll->plotOn( frame, LineColor(2+i), ShiftToZero(), LineStyle(1+i), Name(nll->GetName()) );
    i++;
  }
  nll_1.plotOn(frame  , LineColor(kBlack) , ShiftToZero(), LineStyle(9), Name("trkshw") );
  nll_2.plotOn(frame  , LineColor(kYellow), ShiftToZero(), LineStyle(9), Name("trkshwmid") );

  TLegend *leg1 = new TLegend(0.3, 0.6, 0.7, 0.9);
  for (auto nll: nlls) {
    leg1->AddEntry( frame->findObject(nll->GetName()), nll->GetName(), "l" );
  }
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

  RooMinuit *min0 = new RooMinuit( *( nlls.back() ) );
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
