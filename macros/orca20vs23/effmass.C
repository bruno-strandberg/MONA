#include "EffMass.h"
#include "DetResponse.h"
#include "SummaryParser.h"
#include "NMHUtils.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include<tuple>
#include<vector>

/* This macro created plots to compare the effective masses between ORCA 23 and ORCA 20 geometries. It relies on ORCA MC chain data in NNMO format and effective mass files created with `apps/effective_mass`. It utilises `DetResponse` functionality to plot effective masses for certain quality criteria that aim to match event selections for ORCA 23 and ORCA 20 geometries. The macro can be run interactively (root effmass.C) or compiled (root effmass.C+), the latter will be substantially quicker in looping over events.*/

void effmass(Int_t nu_pdg = 14, Bool_t iscc = 1, TString outname = "") {

  //------------------------------------------------------------------
  // input parameters
  //------------------------------------------------------------------

  vector<Int_t> supported = {12,14,16,-12,-14,-16};
  
  if ( std::find( supported.begin(), supported.end(), nu_pdg ) == supported.end() ) {
    throw std::invalid_argument("ERROR! effmass() Unknown neutrino flavor, supported are +- 12, 14, 16");
  }


  TString effmf_orca23 = "../../data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  TString effmf_orca20 = "../../data/eff_mass/EffMass_ORCA115_20x9m_ECAP1218.root";

  TString dataf_orca23 = "../../data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root";
  TString dataf_orca20 = "../../data/ORCA_MC_summary_ORCA115_20x9m_ECAP1218.root";

  Int_t nebins  = 40;
  Int_t nctbins = 40;

  //------------------------------------------------------------------
  // init eff masses; configure the responses to create approximately comparable event selections
  //------------------------------------------------------------------

  EffMass *orca23em = new EffMass(effmf_orca23, nebins, nctbins, 1);
  EffMass *orca20em = new EffMass(effmf_orca20, nebins, nctbins, 1);

  DetResponse r_orca23(DetResponse::mc_truth, "resp_orca23", nebins, 1, 100, nctbins, -1, 1, 1, 0, 1);
  r_orca23.AddCut( &SummaryEvent::Get_track_ql1 , std::greater<double>(), 0.5, false ); // this is "gandalf_loose_is_selected", see 'apps/data_sorting/PIDAlphaToSummary'
  r_orca23.AddCut( &SummaryEvent::Get_shower_ql1, std::greater<double>(), 0.5, false ); // this is "dusj_is_selected", see 'apps/data_sorting/PIDAlphaToSummary'

  DetResponse r_orca20(DetResponse::mc_truth, "resp_orca20", nebins, 1, 100, nctbins, -1, 1, 1, 0, 1);
  r_orca20.AddCut( &SummaryEvent::Get_track_ql1 , std::greater<double>(), 0.5, false ); // this is "gandalf_loose_is_selected", see 'apps/data_sorting/PIDGammaToSummary'
  r_orca20.AddCut( &SummaryEvent::Get_shower_ql2, std::greater<double>(), 0.5, false ); // this is dusj selection as in 23x9, see 'apps/data_sorting/PIDGammaToSummary'

  // here I add cuts that are common to both selections to select neutrino events of the type
  // specified at input
  vector<DetResponse*> resps = { &r_orca23, &r_orca20 };
  for (auto r: resps) {
    r->AddCut( &SummaryEvent::Get_MC_is_neutrino, std::greater<double>() ,              0.5, true);
    r->AddCut( &SummaryEvent::Get_MC_type       , std::equal_to<double>(), (Double_t)nu_pdg, true);
    r->AddCut( &SummaryEvent::Get_MC_is_CC      , std::equal_to<double>(),   (Double_t)iscc, true);
  }

  SummaryParser sp_orca23(dataf_orca23);
  for (Int_t i = 0; i < sp_orca23.GetTree()->GetEntries(); i++) {
    if ( i % 1000000 == 0) cout << "Event: " << i << " of " << sp_orca23.GetTree()->GetEntries() << endl;
    r_orca23.Fill( sp_orca23.GetEvt(i) );
  }
  cout << "NOTICE ORCA 23x9 response ready" << endl;

  SummaryParser sp_orca20(dataf_orca20);
  for (Int_t i = 0; i < sp_orca20.GetTree()->GetEntries(); i++) {
    if ( i % 1000000 == 0) cout << "Event: " << i << " of " << sp_orca20.GetTree()->GetEntries() << endl;
    r_orca20.Fill( sp_orca20.GetEvt(i) );
  }
  cout << "NOTICE ORCA 20x9 response ready" << endl;

  //------------------------------------------------------------------
  // create original and corrected 2D histograms with effective masses
  //------------------------------------------------------------------
  
  std::map<Int_t, Int_t> pdgtoflav = { {12, 0}, {14, 1}, {16, 2} };

  TH2D* hmeff_orca23 = (TH2D*)orca23em->GetMeff3DH( pdgtoflav[TMath::Abs(nu_pdg)], iscc, (nu_pdg < 0) )->Project3D("yx")->Clone("meff_orca23");
  TH2D* hmeff_orca20 = (TH2D*)orca20em->GetMeff3DH( pdgtoflav[TMath::Abs(nu_pdg)], iscc, (nu_pdg < 0) )->Project3D("yx")->Clone("meff_orca20");
  hmeff_orca23->SetDirectory(0);
  hmeff_orca20->SetDirectory(0);
  
  TH2D* hmeff_orca23_c = (TH2D*)hmeff_orca23->Clone("meff_orca23_c");
  TH2D* hmeff_orca20_c = (TH2D*)hmeff_orca20->Clone("meff_orca20_c");
  hmeff_orca23_c->SetDirectory(0);
  hmeff_orca20_c->SetDirectory(0);
  hmeff_orca23_c->Reset();
  hmeff_orca20_c->Reset();

  if (! NMHUtils::BinsMatch(hmeff_orca23, hmeff_orca20) ) {
    throw std::logic_error("ERROR! effmass() 2D effective mass histograms have different binning.");
  }

  for (Int_t xbin = 1; xbin <= hmeff_orca23->GetXaxis()->GetNbins(); xbin++) {
    for (Int_t ybin = 1; ybin <= hmeff_orca23->GetYaxis()->GetNbins(); ybin++) {

      Double_t E  = hmeff_orca23->GetXaxis()->GetBinCenter(xbin);
      Double_t ct = hmeff_orca23->GetYaxis()->GetBinCenter(ybin);
      Double_t by = 0.5;  // in this script I always use 1 by bin with center 0.5
      
      Double_t w23 = 0;
      Double_t w20 = 0;

      auto wghs_23 = r_orca23.GetBinWeights(E, ct, by);
      auto wghs_20 = r_orca20.GetBinWeights(E, ct, by);

      if ( wghs_23.size() > 1 || wghs_20.size() > 1 ) {
	throw std::logic_error("ERROR! Something is wrong, each bin should have 1 or 0 associated weights");
      }

      if ( wghs_23.size() == 1) w23 = wghs_23[0].fW;
      if ( wghs_20.size() == 1) w20 = wghs_20[0].fW;

      Double_t bc_23 = hmeff_orca23->GetBinContent(xbin, ybin);
      Double_t be_23 = hmeff_orca23->GetBinError(xbin, ybin);
      Double_t bc_20 = hmeff_orca20->GetBinContent(xbin, ybin);
      Double_t be_20 = hmeff_orca20->GetBinError(xbin, ybin);

      hmeff_orca23_c->SetBinContent(xbin, ybin, bc_23 * w23);
      hmeff_orca23_c->SetBinError(xbin, ybin, be_23 * w23);
      hmeff_orca20_c->SetBinContent(xbin, ybin, bc_20 * w20);
      hmeff_orca20_c->SetBinError(xbin, ybin, be_20 * w20);

    }
  }

  gStyle->SetOptStat(0);

  TCanvas *c0 = new TCanvas("c0","c0",1);
  c0->Divide(2,2);
  c0->cd(1);
  hmeff_orca23->Draw("colz");
  c0->cd(2);
  hmeff_orca20->Draw("colz");
  c0->cd(3);
  hmeff_orca23_c->Draw("colz");
  c0->cd(4);
  hmeff_orca20_c->Draw("colz");

  //------------------------------------------------------------------
  // create slices
  //------------------------------------------------------------------

  vector< std::tuple<TH1*, TH1*, TH1*, TH1*> > slices;

  for (Int_t ybin = 1; ybin <= hmeff_orca23_c->GetYaxis()->GetNbins(); ybin++) {

    Double_t ct = hmeff_orca23_c->GetYaxis()->GetBinCenter(ybin);
    if (ct > 0) break;
    
    slices.push_back( std::make_tuple( hmeff_orca23  ->ProjectionX("o23_ct="  + (TString)to_string(ct), ybin, ybin), 
				       hmeff_orca20  ->ProjectionX("o20_ct="  + (TString)to_string(ct), ybin, ybin),
				       hmeff_orca23_c->ProjectionX("o23C_ct=" + (TString)to_string(ct), ybin, ybin),
				       hmeff_orca20_c->ProjectionX("o20C_ct=" + (TString)to_string(ct), ybin, ybin)
				       ) );
    
    std::get<0>( slices.back() )->SetLineColor(kBlue);
    std::get<1>( slices.back() )->SetLineColor(kRed);
    std::get<2>( slices.back() )->SetLineColor(kBlue);
    std::get<3>( slices.back() )->SetLineColor(kRed);

    std::get<0>( slices.back() )->SetLineWidth(2);
    std::get<1>( slices.back() )->SetLineWidth(2);
    std::get<2>( slices.back() )->SetLineWidth(2);
    std::get<3>( slices.back() )->SetLineWidth(2);

    std::get<0>( slices.back() )->SetTitle( std::get<0>( slices.back() )->GetName() );
    std::get<1>( slices.back() )->SetTitle( std::get<1>( slices.back() )->GetName() );
    std::get<2>( slices.back() )->SetTitle( std::get<2>( slices.back() )->GetName() );
    std::get<3>( slices.back() )->SetTitle( std::get<3>( slices.back() )->GetName() );

    std::get<2>( slices.back() )->SetLineStyle(9);
    std::get<3>( slices.back() )->SetLineStyle(9);
  }

  TLegend *leg = new TLegend(0.1, 0.5, 0.5, 0.9);
  leg->AddEntry( std::get<0>( slices.back() ), "orca23"     , "l" );
  leg->AddEntry( std::get<1>( slices.back() ), "orca20"     , "l" );
  leg->AddEntry( std::get<2>( slices.back() ), "orca23-corr", "l" );
  leg->AddEntry( std::get<3>( slices.back() ), "orca20-corr", "l" );
  leg->SetLineWidth(0);
  leg->SetFillStyle(0);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare( slices.size() );
  Int_t pad = 1;
  for (auto slice: slices) {
    c1->cd(pad);
    std::get<1>(slice)->Draw("HIST");
    std::get<0>(slice)->Draw("HISTsame");
    std::get<2>(slice)->Draw("HISTsame");
    std::get<3>(slice)->Draw("HISTsame");
    leg->Draw();
    pad++;
  }

  if (outname != "") {
    TFile fout(outname, "RECREATE");
    hmeff_orca23->Write();
    hmeff_orca20->Write();
    hmeff_orca23_c->Write();
    hmeff_orca20_c->Write();
    for (auto slice: slices) {
      std::get<1>(slice)->Write();
      std::get<0>(slice)->Write();
      std::get<2>(slice)->Write();
      std::get<3>(slice)->Write();
    }
  }
  
}
