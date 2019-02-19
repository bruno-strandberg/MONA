#include "AsimovFit.C"

Bool_t FitparAtLimit(RooFitResult* fitres) {

  // a fraction of the variable range that defines the parameter limit
  // for example: dcp range is 0-2; if the fit is within 2/1000 of one of the limits,
  // the parameter is assumed to be at limit
  Double_t frac = 1000;

  // prepare looping over fitted variables
  RooArgList fitpars = fitres->floatParsFinal();
  TIterator *iter    = fitpars.createIterator();
  RooRealVar* var;
  
  while ( ( var = (RooRealVar*)iter->Next() ) ) {

    // get the width of the range
    Double_t width = var->getMax() - var->getMin();

    // return true if any fitted variable is at upper or lower limit
    if ( ( var->getVal() - var->getMin() ) < width/frac ) {
      cout << "NOTICE FitparAtLimit() variable " << var->GetName() << " value " << var->getVal()
	   << " is at lower limit " << var->getMin() << ", returning true" << endl;
      return kTRUE;
    }

    if ( ( var->getMax() - var->getVal() ) < width/frac ) {
      cout << "NOTICE FitparAtLimit() variable " << var->GetName() << " value " << var->getVal()
	   << " is at upper limit " << var->getMax() << ", returning true" << endl;
      return kTRUE;
    }
    
  }

  return kFALSE;

}

//*************************************************************************************

/**
   \param ExcludeBadFits  - excludes fits where one of the parameters is at limit
   \param TSCHI2          - use chi2 calculated from track and shower, do not use the external constraint
 */
void FillChi2Graph(AsimovFit &AF, TGraph *g, Bool_t ExcludeBadFits = kFALSE, Bool_t TSCHI2 = kFALSE) {

  for (auto fp: AF.fFPs) {

    cout << "NOTICE FillChi2Graph() analysing point " << g->GetN() << " for detector " << fp->fDetString << endl;

    auto d_chi2_1q = AF.GetChi2(*fp, kTRUE);
    auto d_chi2_2q = AF.GetChi2(*fp, kFALSE);
    Double_t best_chi2;
    RooFitResult *bestfit;

    // use a chi2 based only on tracks and showers, external constraint is not included
    if (TSCHI2) {
     
      if ( ( std::get<1>(d_chi2_1q) + std::get<2>(d_chi2_1q) ) < 
	   ( std::get<1>(d_chi2_2q) + std::get<2>(d_chi2_2q) ) ) {
	best_chi2 = std::get<1>(d_chi2_1q) + std::get<2>(d_chi2_1q);
	bestfit = fp->fRes_1q;
      }
      else {
	best_chi2 = std::get<1>(d_chi2_2q) + std::get<2>(d_chi2_2q);
	bestfit = fp->fRes_2q;
      }
 
    }

    // select best fit based on the combined chi2 which includes the X2 from external constraint on dm31
    else {
      if ( std::get<0>(d_chi2_1q) < std::get<0>(d_chi2_2q) ) { 
	best_chi2 = std::get<0>(d_chi2_1q); 
	bestfit = fp->fRes_1q;
      }
      else { 
	best_chi2 = std::get<0>(d_chi2_2q);
	bestfit = fp->fRes_2q;
      }
    }

    // logic to ingore "bad" points
    if ( ExcludeBadFits ) {
      if ( FitparAtLimit( bestfit ) ) continue;
    }

    // get the th23 value
    Double_t sinsqth23 = ( (RooRealVar*)fp->fParData->find("SinsqTh23") )->getVal();
    Double_t th23 = TMath::ASin( TMath::Sqrt( sinsqth23 ) ) * 180./TMath::Pi();

    g->SetPoint( g->GetN(), th23, TMath::Sqrt(best_chi2) );
    
  }
  
}

//*************************************************************************************

void plot20vs23(TString infile = "rootfiles/asimov_th23scan_0.35-0.65_v8_chi2fit_NOdata.root", 
		TString plottitle = "ORCA Asimov NMO sensitivity for NO data",
		Bool_t excludebadfits = kFALSE, Bool_t TSCHI2 = kFALSE) {

  AsimovFit AF20(AsimovFit::ORCA20);
  AsimovFit AF23(AsimovFit::ORCA23);

  // read the results to both instances
  AF20.ReadFromFile(infile);
  AF23.ReadFromFile(infile);

  TGraph *g_O20 = new TGraph();
  TGraph *g_O23 = new TGraph();

  FillChi2Graph( AF20, g_O20, excludebadfits, TSCHI2 );
  FillChi2Graph( AF23, g_O23, excludebadfits, TSCHI2 );

  vector<TGraph*> gs = { g_O20, g_O23 };

  for (auto g: gs) {

    g->GetXaxis()->SetTitle("#theta_{23} [Deg]");
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetLabelSize(0.04);
    g->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetLabelSize(0.04);
    
    g->GetXaxis()->SetTitleOffset(1.1);
    g->GetYaxis()->SetTitleOffset(1.1);
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();

    g->SetLineWidth(2);

    g->Sort();
    
  }

  g_O20->SetLineColor(kBlack);
  g_O23->SetLineColor(kBlue);
  g_O23->SetLineStyle(9);

  g_O20->GetYaxis()->SetRangeUser(0,7);
  g_O20->SetMinimum(0);
  g_O20->SetMaximum(7);

  TLegend *leg = new TLegend(0.15, 0.6, 0.45, 0.9);
  leg->AddEntry(g_O20, "ORCA 20x9m", "l");
  leg->AddEntry(g_O23, "ORCA 23x9m", "l");
  leg->SetLineWidth(0);
  leg->SetFillStyle(0);

  g_O20->SetTitle(plottitle);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SetTicks();
  c1->SetGrid();
  g_O20->Draw();
  g_O23->Draw("sameL");
  leg->Draw();

}
