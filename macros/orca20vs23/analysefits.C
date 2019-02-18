#include "AsimovFit.h"
#include "AsimovFit.C"

// return type of GetChi2 and GetAsym
typedef std::tuple<Double_t, Double_t, Double_t> gof_t;

// goodness-of-fit graphs
struct GofGraphs {
  TGraph *g_trk = new TGraph(); // graph for tracks only
  TGraph *g_shw = new TGraph(); // graph for showers only
  TGraph *g_ts  = new TGraph(); // graph for tracks+showers
  TGraph *g_tse = new TGraph(); // graph for tracks+shower+chi2 from external constraint on th13
  TGraph *g_e   = new TGraph(); // graph showing the X2 ONLY from external constraint on th13
  vector<TGraph*> memgs = { g_tse, g_ts, g_trk, g_shw, g_e }; // for looping
};

// plots depicting fitted variables vs th23 true value
struct ParGraphs {
  TGraphErrors *g_th13 = new TGraphErrors();
  TGraphErrors *g_th23 = new TGraphErrors();
  TGraphErrors *g_dm31 = new TGraphErrors();
  TGraphErrors *g_dcp  = new TGraphErrors();

  // these show the true value at which the data was generated
  TGraph *g_th13_true = new TGraph();
  TGraph *g_th23_true = new TGraph();
  TGraph *g_dm31_true = new TGraph();
  TGraph *g_dcp_true  = new TGraph();

  vector< TGraph* > memgs = { g_th13, g_th23, g_dm31, g_dcp, g_th13_true, g_th23_true, g_dm31_true, g_dcp_true };
  vector< TGraph* > memgs_true = { g_th13_true, g_th23_true, g_dm31_true, g_dcp_true };
};

//==========================================================================

Double_t SinsqToDeg(Double_t sinsq) {
  return TMath::ASin( TMath::Sqrt(sinsq) ) * 180./TMath::Pi();
}

//==========================================================================

void AddGofPoint(GofGraphs* g, gof_t data, Double_t th23_true) {

  // calculate combined chi2 without external constraint
  Double_t trkX2 = std::get<1>(data);
  Double_t shwX2 = std::get<2>(data);
  Double_t combchi2 = trkX2 + shwX2;

  g->g_trk->SetPoint( g->g_trk->GetN(), th23_true, trkX2);
  g->g_shw->SetPoint( g->g_shw->GetN(), th23_true, shwX2);
  g->g_ts ->SetPoint( g->g_ts ->GetN(), th23_true, combchi2);
  g->g_tse->SetPoint( g->g_tse->GetN(), th23_true, std::get<0>(data));
  g->g_e  ->SetPoint( g->g_e  ->GetN(), th23_true, std::get<0>(data)-combchi2 );

}

//==========================================================================

void AddParPoint(ParGraphs* g, RooFitResult *fr, RooArgSet* pars_start, Double_t th23_true) {

    //----------------------------------------------------------------
    // get the fitted parameters
    //----------------------------------------------------------------
    RooRealVar* sinsqth13_fit = (RooRealVar*)fr->floatParsFinal().find("SinsqTh13");
    RooRealVar* sinsqth23_fit = (RooRealVar*)fr->floatParsFinal().find("SinsqTh23");
    RooRealVar* dm31_fit      = (RooRealVar*)fr->floatParsFinal().find("Dm31");
    RooRealVar* dcp_fit       = (RooRealVar*)fr->floatParsFinal().find("dcp");
    
    g->g_th13->SetPoint( g->g_th13->GetN(), th23_true, SinsqToDeg( sinsqth13_fit->getVal() ) );
    g->g_th13->SetPointError( g->g_th13->GetN()-1, 0, SinsqToDeg( sinsqth13_fit->getError() ) );

    g->g_th23->SetPoint( g->g_th23->GetN(), th23_true, SinsqToDeg( sinsqth23_fit->getVal() ) );
    g->g_th23->SetPointError( g->g_th23->GetN()-1, 0, SinsqToDeg( sinsqth23_fit->getError() ) );

    g->g_dm31->SetPoint( g->g_dm31->GetN(), th23_true, dm31_fit->getVal() );
    g->g_dm31->SetPointError( g->g_dm31->GetN()-1, 0, dm31_fit->getError() );

    g->g_dcp->SetPoint( g->g_dcp->GetN(), th23_true, dcp_fit->getVal() );
    g->g_dcp->SetPointError( g->g_dcp->GetN()-1, 0, dcp_fit->getError() );

    //----------------------------------------------------------------
    // fill the graphs that show the true value
    //----------------------------------------------------------------
    RooRealVar *sinsqth13_true = (RooRealVar*)pars_start->find("SinsqTh13");
    RooRealVar *sinsqth23_true = (RooRealVar*)pars_start->find("SinsqTh23");
    RooRealVar *dm31_true      = (RooRealVar*)pars_start->find("Dm31");
    RooRealVar *dcp_true       = (RooRealVar*)pars_start->find("dcp");

    g->g_th13_true->SetPoint( g->g_th13_true->GetN(), th23_true, SinsqToDeg( sinsqth13_true->getVal() ) );
    g->g_th23_true->SetPoint( g->g_th23_true->GetN(), th23_true, SinsqToDeg( sinsqth23_true->getVal() ) );
    g->g_dm31_true->SetPoint( g->g_dm31_true->GetN(), th23_true, -dm31_true->getVal() ); // flipped!
    g->g_dcp_true ->SetPoint( g->g_dcp_true ->GetN(), th23_true, dcp_true->getVal() );
    
};

//==========================================================================

void analysefits(AsimovFit::Detector det = AsimovFit::ORCA20, 
		 TString infile = "rootfiles/asimov_th23scan_0.35-0.65_v6_chi2fit.root") {

  AsimovFit AF(det);
  AF.ReadFromFile(infile);
  auto fps = AF.fFPs;

  GofGraphs *g_chi2_best  = new GofGraphs;
  GofGraphs *g_chi2_other = new GofGraphs; 
  GofGraphs *g_asym_best  = new GofGraphs;
  GofGraphs *g_asym_other = new GofGraphs;
  
  ParGraphs *g_pars = new ParGraphs;

  Int_t counter = 0;

  for (auto fp: AF.fFPs) {

    // get the th23 value
    Double_t sinsqth23 = ( (RooRealVar*)fp->fParData->find("SinsqTh23") )->getVal();
    Double_t th23 = SinsqToDeg( sinsqth23 );

    cout << "Analyzing fit at: " << sinsqth23 << " for detector " << fp->fDetString << endl;

    // calculate chi2 data and asym data
    gof_t d_chi2_1q = AF.GetChi2(*fp, kTRUE);
    gof_t d_chi2_2q = AF.GetChi2(*fp, kFALSE);
    gof_t d_asym_1q = AF.GetAsym(*fp, kTRUE);
    gof_t d_asym_2q = AF.GetAsym(*fp, kFALSE);

    gof_t d_chi2_best, d_chi2_other, d_asym_best, d_asym_other;
    RooFitResult *result_best, *result_other;

    cout << "*****************MINIMIZER STATUSES IN FIRST AND SECOND QUADRANT***********" << endl;
    cout << fp->fRes_1q->status() << "\t" << fp->fRes_2q->status() << endl;
    cout << "*****************MINIMIZER STATUSES IN FIRST AND SECOND QUADRANT***********" << endl;

    if ( std::get<0>(d_chi2_1q) < std::get<0>(d_chi2_2q) ) {
      d_chi2_best  = d_chi2_1q;
      d_chi2_other = d_chi2_2q;
      d_asym_best  = d_asym_1q;
      d_asym_other = d_asym_2q;
      result_best  = fp->fRes_1q;
      result_other = fp->fRes_2q;
    }
    else {
      d_chi2_best  = d_chi2_2q;
      d_chi2_other = d_chi2_1q;
      d_asym_best  = d_asym_2q;
      d_asym_other = d_asym_1q;
      result_best  = fp->fRes_2q;
      result_other = fp->fRes_1q;
    }

    // fill goodness-of-fit graphs
    AddGofPoint( g_chi2_best , d_chi2_best , th23 );
    AddGofPoint( g_chi2_other, d_chi2_other, th23 );
    AddGofPoint( g_asym_best , d_asym_best , th23 );
    AddGofPoint( g_asym_other, d_asym_other, th23 );    

    // fill the parameter plots
    AddParPoint(g_pars, result_best, fp->fParData, th23);

    counter++;

  }

  //-----------------------------------------------------------------
  // operations on the graphs
  //-----------------------------------------------------------------
  vector<GofGraphs*> gcolls = { g_chi2_best, g_chi2_other, g_asym_best, g_asym_other };

  for (auto gc: gcolls) { 

    for (auto g: gc->memgs) {
      g->Sort();          // sort points along X values
      g->SetLineWidth(2); // increase line width
    }

    gc->g_trk->SetLineColor(kBlue);  // tracks always blue
    gc->g_shw->SetLineColor(kRed); // showers always red
    gc->g_ts->SetLineColor(kGreen); // tracks & showers combined are green
    gc->g_tse->SetLineColor(kBlack); // tracks & showers & external combined are black
    gc->g_e->SetLineColor(6);        // external X2 lilla

  }

  // operations on asymmetry graphs only
  vector<GofGraphs*> gasym = { g_asym_best, g_asym_other };

  for (auto gc: gasym) { 
    for (auto g: gc->memgs) {
      g->SetLineStyle(9); // asym plots are dashed
    }
  }

  // operations on the parameter plots
  for (auto gp: g_pars->memgs) {
    gp->Sort();
    gp->SetLineWidth(2);
    gp->SetMarkerStyle(20);
    gp->SetMarkerColor(kRed);
  }
  
  for (auto gp: g_pars->memgs_true) {
    gp->SetLineColor(kBlue);
  }

  //-----------------------------------------------------------------
  // drawing
  //-----------------------------------------------------------------

  // 4 canvases; each canvas has combined, trk and shw result; 1 canvas for each of
  // {bestchi2, otherchi2, bestasym, otherasym}
  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare( gcolls.size() );
  
  Int_t pad = 1;
  for (auto gc: gcolls) {

    c1->cd(pad);
    Int_t gcount = 0;
    for (auto g: gc->memgs) {
      if (gcount == 0) { g->Draw(); }
      else g->Draw("same");
      gcount++;
    }

    gcount = 0;
    pad++;

  }

  // top left - combined chi2 best and other, top right same for asym
  // bottom left trk chi2 best and other, bottom right shw chi2 best and other
  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->DivideSquare(4);
  c2->cd(1);
  g_chi2_best->g_ts->Draw();
  g_chi2_other->g_ts->Draw("same");
  c2->cd(2);
  g_asym_best->g_ts->Draw();
  g_asym_other->g_ts->Draw("same");
  c2->cd(3);
  g_chi2_best->g_trk->Draw();
  g_chi2_other->g_trk->Draw("same");  
  c2->cd(4);
  g_chi2_best->g_shw->Draw();
  g_chi2_other->g_shw->Draw("same");  

  // left - comparison of best chi2 and asym, right - chi2 with and without X2 term for external constraint
  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->DivideSquare(4);
  c3->cd(1);
  g_chi2_best->g_ts->Draw();
  g_asym_best->g_ts->Draw("same");
  c3->cd(2);
  g_chi2_best->g_tse->Draw();
  g_chi2_best->g_ts->Draw("same");


  // parameter fitted vs true at various th23
  TCanvas *c4 = new TCanvas("c4","c4",1);
  c4->DivideSquare(4);
  c4->cd(1);
  g_pars->g_th13->Draw();
  g_pars->g_th13_true->Draw("sameL");
  c4->cd(2);
  g_pars->g_th23->Draw();
  g_pars->g_th23_true->Draw("sameL");
  c4->cd(3);
  Double_t x,y;
  g_pars->g_dm31->GetPoint(0, x, y);
  if ( y > 0) { g_pars->g_dm31->GetYaxis()->SetRangeUser(2e-3, 3e-3); }
  else        { g_pars->g_dm31->GetYaxis()->SetRangeUser(-3e-3, -2e-3); }
  g_pars->g_dm31->Draw();
  g_pars->g_dm31_true->Draw("sameL");
  c4->cd(4);
  g_pars->g_dcp->Draw();
  g_pars->g_dcp_true->Draw("sameL");

}
