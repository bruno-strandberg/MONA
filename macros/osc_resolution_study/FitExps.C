#include "DetResponse.h"
#include "SummaryParser.h"
#include "TFile.h"
#include "TH3.h"
#include "FitUtil.h"
#include "FitPDF.h"
#include "TGraph.h"
#include "NMHUtils.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooConstVar.h"

#include <iostream>
#include <stdexcept>
using namespace std;
using namespace RooFit;

/** This macro reads the trk and shower histogram from the input file (output of OscResolution) and fits them simultaneously, using the oscillation sampling N=1 (bin centers). The outputs of this macro indicate the systematic uncertainty introduced by sampling the oscillation probability sparsely.*/

void FitExps(TString filename="OscResolution.root", TString h_trk, TString h_shw, TString outfile,
	     TString datafile    = "../../data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root", 
	     TString effmassfile = "../../data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root") {

  //================================================================================
  // initialise responses
  //================================================================================

  DetResponse trk  (DetResponse::track   , "trk"  , 24, 1, 100, 40, -1, 1, 1, 0, 1);
  DetResponse shw  (DetResponse::shower  , "shw"  , 24, 1, 100, 40, -1, 1, 1, 0, 1);

  // semi-standard track response and shower response
  //-----------------------------------------------------
  trk.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  trk.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );

  shw.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5, true );
  shw.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), 0.6, true );

  // for all selections only use neutrino events, noise and atm muons not sensitive to oscillation probabilities
  vector<DetResponse*> resps = { &trk, &shw };
  for (auto R: resps) R->AddCut( &SummaryEvent::Get_MC_is_neutrino, std::greater<double>(), 0.5, true );

  SummaryParser sp(datafile);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % 200000 == 0) cout << "Event: " << i << "/" << sp.GetTree()->GetEntries() << endl;
    for (auto R: resps) R->Fill( sp.GetEvt(i) );
  }

  //================================================================================
  // init fitutil and pdf's, get histograms, fit and write out
  //================================================================================
  FitUtil futil(3, trk.GetHist3D(), 3, 75, -1, -1e-3, 0, 1, effmassfile);
  FitPDF pdf_trk  ("pdftrk"  ,"pdftrk"  , &futil, &trk);
  FitPDF pdf_shw  ("pdfshw"  ,"pdfshw"  , &futil, &shw);

  TFile fin(filename, "READ");
  TH3D* htrk = (TH3D*)fin.Get(h_trk);
  TH3D* hshw = (TH3D*)fin.Get(h_shw);

  if (htrk == NULL || hshw == NULL) {
    throw std::invalid_argument("ERROR! FitExps() Cannot find histogram(s)");
  }

  std::map<string, TH1*> hist_map = { {"trk", (TH1*)htrk}, {"shw", (TH1*)hshw}};

  RooCategory categ("categs","categs");
  categ.defineType("trk");
  categ.defineType("shw");

  RooDataHist combData("combData","combData", futil.GetObs(), categ, hist_map);
    
  RooSimultaneous simPdf("simPdf","simPdf",categ);
  simPdf.addPdf( pdf_trk, "trk" );
  simPdf.addPdf( pdf_shw, "shw" );

  futil.GetVar("SinsqTh12")->setConstant(kTRUE); //dm21 and th12 are typically fixed
  futil.GetVar("Dm21")     ->setConstant(kTRUE);

  RooFitResult *res = simPdf.fitTo( combData, Save(kTRUE), SumW2Error(kFALSE) );

  TFile fout(outfile, "RECREATE");
  res->Write("fitresult");
  htrk->Write();
  hshw->Write();
  fout.Close();
  
}
