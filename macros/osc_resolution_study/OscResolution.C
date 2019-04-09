#include "DetResponse.h"
#include "SummaryParser.h"
#include "TFile.h"
#include "TH3.h"
#include "FitUtil.h"
#include "FitPDF.h"
#include "TGraph.h"
#include "NMHUtils.h"

#include <iostream>
using namespace std;

/** This macro creates expectation values plots at identical oscillation parameter values, using a different value for oscillation sampling. n=1 means that bin center is used to estimate the oscillation probability, n=10 means that the oscillation probablitity is an average of 10^2=100 samples (10 in E x 10 in ct) inside the bin.*/

void OscResolution(TString datafile    = "../../data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root", 
		   TString effmassfile = "../../data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root",
		   TString output = "OscResolution.root") {

  //================================================================================
  // initialise responses
  //================================================================================

  DetResponse trk  (DetResponse::track   , "trk"  , 24, 1, 100, 40, -1, 1, 1, 0, 1);
  DetResponse shw  (DetResponse::shower  , "shw"  , 24, 1, 100, 40, -1, 1, 1, 0, 1);
  DetResponse trkMC(DetResponse::mc_truth, "trkMC", 24, 1, 100, 40, -1, 1, 1, 0, 1);
  DetResponse shwMC(DetResponse::mc_truth, "shwMC", 24, 1, 100, 40, -1, 1, 1, 0, 1);

  // semi-standard track response and shower response
  //-----------------------------------------------------
  trk.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  trk.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );

  shw.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5, true );
  shw.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), 0.6, true );


  // ideal track and shower selections
  //-----------------------------------------------------
  trkMC.AddCut( &SummaryEvent::Get_MC_type, std::equal_to<double>(),  14., false );
  trkMC.AddCut( &SummaryEvent::Get_MC_type, std::equal_to<double>(), -14., false );

  shwMC.AddCut( &SummaryEvent::Get_MC_type, std::not_equal_to<double>(),  14., true );
  shwMC.AddCut( &SummaryEvent::Get_MC_type, std::not_equal_to<double>(), -14., true );

  // for all selections only use neutrino events, noise and atm muons not sensitive to oscillation probabilities
  vector<DetResponse*> resps = { &trkMC, &trk, &shwMC, &shw };
  for (auto R: resps) R->AddCut( &SummaryEvent::Get_MC_is_neutrino, std::greater<double>(), 0.5, true );

  SummaryParser sp(datafile);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % 200000 == 0) cout << "Event: " << i << "/" << sp.GetTree()->GetEntries() << endl;
    for (auto R: resps) R->Fill( sp.GetEvt(i) );
  }

  //================================================================================
  // init fitutil and pdf's, create expectation value histograms at different oscillation sampling
  //================================================================================
  FitUtil futil(3, trk.GetHist3D(), 3, 75, -1, 0, 0, 1, effmassfile);

  FitPDF pdf_trkMC("pdftrkMC","pdftrkMC", &futil, &trkMC);
  FitPDF pdf_shwMC("pdfshwMC","pdfshwMC", &futil, &shwMC);
  FitPDF pdf_trk  ("pdftrk"  ,"pdftrk"  , &futil, &trk);
  FitPDF pdf_shw  ("pdfshw"  ,"pdfshw"  , &futil, &shw);

  vector<TH3D*> hTrkMC, hShwMC, hTrk, hShw, hTrkIO, hShwIO;

  Int_t max_samples = 10;

  for (Int_t n = 1; n <= max_samples; n++) {

    cout << "NOTICE N osc samples: " << n*n << endl;

    TString sfx = "_n" + (TString)to_string(n);
    futil.SetOscSamplesN( n );
    
    // expectation values at given n, ideal MC, normal ordering
    hTrkMC.push_back( pdf_trkMC.GetExpValHist() );
    hTrkMC.back()->SetNameTitle("trkMC"+sfx, "trkMC"+sfx);

    hShwMC.push_back( pdf_shwMC.GetExpValHist() );
    hShwMC.back()->SetNameTitle("shwMC"+sfx, "shwMC"+sfx);

    // expectation values at given n, realistic response, normal ordering
    hTrk.push_back( pdf_trk.GetExpValHist() );
    hTrk.back()->SetNameTitle("trk"+sfx, "trk"+sfx);

    hShw.push_back( pdf_shw.GetExpValHist() );
    hShw.back()->SetNameTitle("shw"+sfx, "shw"+sfx);

    // expectation values at given n, realistic response, inverted ordering
    futil.SetIOcentvals();
    
    hTrkIO.push_back( pdf_trk.GetExpValHist() );
    hTrkIO.back()->SetNameTitle("trkIO"+sfx, "trkIO"+sfx);

    hShwIO.push_back( pdf_shw.GetExpValHist() );
    hShwIO.back()->SetNameTitle("shwIO"+sfx, "shwIO"+sfx);

    futil.SetNOcentvals();
  }

  //================================================================================
  // perform analysis
  //================================================================================
  TGraph g_trkMC, g_shwMC, g_trk, g_shw, g_trkMC_m1, g_shwMC_m1, g_trk_m1, g_shw_m1, g_asym_trk, g_asym_shw;

  vector<TH1*> asym_hists;

  for (UInt_t i = 0; i < hTrkMC.size(); i++) {

    // calculate the asym (=sqrt(X2)) between histograms with samples = 1 and samples = n
    auto A_trkMC = NMHUtils::Asymmetry( hTrkMC.front(), hTrkMC[i], "asym_trkMC_n" + (TString)to_string(i+1) );
    auto A_shwMC = NMHUtils::Asymmetry( hShwMC.front(), hShwMC[i], "asym_shwMC_n" + (TString)to_string(i+1) );
    auto A_trk   = NMHUtils::Asymmetry( hTrk.front(), hTrk[i], "asym_trk_n" + (TString)to_string(i+1) );
    auto A_shw   = NMHUtils::Asymmetry( hShw.front(), hShw[i], "asym_shw_n" + (TString)to_string(i+1) );

    g_trkMC.SetPoint( g_trkMC.GetN(), i+1, std::get<1>(A_trkMC) );
    g_shwMC.SetPoint( g_shwMC.GetN(), i+1, std::get<1>(A_shwMC) );
    g_trk.SetPoint( g_trk.GetN(), i+1, std::get<1>(A_trk) );
    g_shw.SetPoint( g_shw.GetN(), i+1, std::get<1>(A_shw) );

    asym_hists.push_back( std::get<0>(A_trkMC) );
    asym_hists.push_back( std::get<0>(A_shwMC) );
    asym_hists.push_back( std::get<0>(A_trk) );
    asym_hists.push_back( std::get<0>(A_shw) );

    // calculate the asym (=sqrt(X2)) between histograms with samples = n-1 and samples = n
    if (i >= 1) {
      auto A_trkMC_m1 = NMHUtils::Asymmetry(hTrkMC[i-1], hTrkMC[i], "asym_trkMC_m1_n" + (TString)to_string(i+1) );
      auto A_shwMC_m1 = NMHUtils::Asymmetry(hShwMC[i-1], hShwMC[i], "asym_shwMC_m1_n" + (TString)to_string(i+1) );
      auto A_trk_m1   = NMHUtils::Asymmetry(hTrk[i-1], hTrk[i], "asym_trk_m1_n" + (TString)to_string(i+1) );
      auto A_shw_m1   = NMHUtils::Asymmetry(hShw[i-1], hShw[i], "asym_shw_m1_n" + (TString)to_string(i+1) );

      g_trkMC_m1.SetPoint( g_trkMC_m1.GetN(), i+1, std::get<1>(A_trkMC_m1) );
      g_shwMC_m1.SetPoint( g_shwMC_m1.GetN(), i+1, std::get<1>(A_shwMC_m1) );
      g_trk_m1.SetPoint( g_trk_m1.GetN(), i+1, std::get<1>(A_trk_m1) );
      g_shw_m1.SetPoint( g_shw_m1.GetN(), i+1, std::get<1>(A_shw_m1) );

      asym_hists.push_back( std::get<0>(A_trkMC_m1) );
      asym_hists.push_back( std::get<0>(A_shwMC_m1) );
      asym_hists.push_back( std::get<0>(A_trk_m1) );
      asym_hists.push_back( std::get<0>(A_shw_m1) );

    }

    // calculate the asymmetry between NO and IO at given n
    auto asym_trk = NMHUtils::Asymmetry(hTrk[i], hTrkIO[i], "trk_NOvsIO_n" + (TString)to_string(i+1) );
    auto asym_shw = NMHUtils::Asymmetry(hShw[i], hShwIO[i], "shw_NOvsIO_n" + (TString)to_string(i+1) );

    asym_hists.push_back( std::get<0>(asym_trk) );
    asym_hists.push_back( std::get<0>(asym_shw) );

    g_asym_trk.SetPoint( g_asym_trk.GetN(), i+1, std::get<1>(asym_trk) );
    g_asym_shw.SetPoint( g_asym_shw.GetN(), i+1, std::get<1>(asym_shw) );

  }

  //================================================================================
  // write data to output
  //================================================================================

  TFile fout(output, "RECREATE");
  for (auto h: hTrkMC) h->Write();
  for (auto h: hShwMC) h->Write();
  for (auto h: hTrk) h->Write();
  for (auto h: hShw) h->Write();
  g_trkMC    .Write("g_trkMC");
  g_shwMC    .Write("g_shwMC");
  g_trk	     .Write("g_trk");
  g_shw	     .Write("g_shw");
  g_trkMC_m1 .Write("g_trkMC_m1");
  g_shwMC_m1 .Write("g_shwMC_m1");
  g_trk_m1   .Write("g_trk_m1");
  g_shw_m1   .Write("g_shw_m1");
  g_asym_trk .Write("g_asym_trk");
  g_asym_shw .Write("g_asym_shw");
  for (auto h: asym_hists) h->Write();
  fout.Close();

  // remove the junk from compilation
  Int_t sysret = system("rm OscResolution_C*");
  if (sysret != 0) cout << "NOTICE 'rm OscResoltion_C*' returned " << sysret << endl;
  
}
