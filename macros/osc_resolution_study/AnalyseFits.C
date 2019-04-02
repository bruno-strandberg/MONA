#include "NMHUtils.h"

/** Macro to analyse outputs of FitExps.C . The input file is a list of files output by FitExps.C */
void AnalyseFits(TString file_list="outlist.dat") {

  auto files = NMHUtils::ReadLines(file_list);

  map<TString, TGraph*> graphs;

  for (auto f: files) {

    TFile fin(f, "READ");
    RooFitResult *res = (RooFitResult*)fin.Get("fitresult");
    
    TString n_str = f( f.Index("_n")+2, f.Length() );
    n_str = n_str( 0, n_str.Index(".root") );
    Int_t N = stoi( (string)n_str );

    RooArgSet par_T = (RooArgSet)res->floatParsInit();
    RooArgSet par_F = (RooArgSet)res->floatParsFinal();
    
    TIterator *it = par_T.createIterator();
    RooRealVar *v = NULL;

    while( ( v = (RooRealVar*)it->Next() ) ) {
      TString vname  = v->GetName();
      Double_t val_T = v->getVal();                                      // true value (same as start)
      Double_t val_F = ( (RooRealVar*)par_F.find( vname ) )->getVal();   // fitted value

      TGraph *g = NULL;
      if ( graphs.find( vname ) == graphs.end() ) {
	graphs.insert( make_pair( vname, new TGraph() ) );
	graphs[vname]->SetNameTitle(vname, vname);
      }
      g = graphs[vname];

      g->SetPoint( g->GetN(), N, val_F );

    }

  }

  map<TString, TLine*> lines;

  for (auto kv: graphs) {
    kv.second->Sort();
    Double_t N, V;
    kv.second->GetPoint(0, N, V);
    if (N != 1.0) throw std::logic_error("Expected N==1, got N=" + to_string(N));
    lines[kv.first] = new TLine(1, V, 10, V);

    kv.second->SetLineColor(kBlue);
    kv.second->SetLineWidth(2);
    kv.second->SetMarkerStyle(20);

    kv.second->GetXaxis()->SetTitle("N");

    lines[kv.first]->SetLineColor(kRed);
    lines[kv.first]->SetLineWidth(2);
    lines[kv.first]->SetLineStyle(9);
  }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare( graphs.size() );

  Int_t pad = 1;
  for (auto kv: graphs) {
    c1->cd(pad);
    kv.second->Draw();
    lines[kv.first]->Draw("same");
    TLegend leg(0.6, 0.2, 0.9, 0.5);
    leg.AddEntry( kv.second      , "Fitted", "l" );
    leg.AddEntry( lines[kv.first], "True"  , "l" );
    leg.SetLineWidth(0);
    leg.SetFillStyle(0);
    leg.DrawClone();
    pad++;
  }

}
