#include "ORCA7.C"
#include "NMHUtils.h"
#include "FitUtilWsyst.h"
#include "FitPDF.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

void EoverZ() {
  
  ORCA7 o7( kTRUE );
  FitUtilWsyst *fu = o7.fFitUtil;
  FitPDF *trkpdf = o7.fPdfs["trk"];

  //---------------------------------------------------------------------
  // set to NuFit3.2, free the limits
  //---------------------------------------------------------------------

  o7.Set_NuFit_3p2_NO();
  fu->FreeParLims();

  //---------------------------------------------------------------------
  // define binning for the E/cosz plot
  //---------------------------------------------------------------------

  TH3D* ht = o7.fResps["trk"]->GetHist3DReco();

  Double_t emin  = ht->GetXaxis()->GetBinCenter( ht->GetXaxis()->FindBin( o7.f_F_emin ) );
  Double_t emax  = ht->GetXaxis()->GetBinCenter( ht->GetXaxis()->FindBin( o7.f_F_emax ) );
  Double_t ctmin = ht->GetYaxis()->GetBinCenter( ht->GetYaxis()->FindBin( o7.f_F_ctmin ) );
  Double_t ctmax = ht->GetYaxis()->GetBinCenter( ht->GetYaxis()->FindBin( o7.f_F_ctmax ) );

  Double_t xmin = -emin/ctmin; // something like 3 / -0.95
  Double_t xmax = -emax/ctmax; // something like 60/-0.1

  // template histogram
  auto bins = NMHUtils::GetLogBins(20, xmin, xmax);
  TH1D* hTemplate = new TH1D("template","template", 20, &bins[0]);
  vector< std::pair<TString,TH1D*> > hists;

  //---------------------------------------------------------------------
  // scan in theta-23
  //---------------------------------------------------------------------
  
  for (Double_t th23 = 0; th23 <= 1.0; th23 += 0.25) {

    fu->GetVar("SinsqTh23")->setVal( th23 );
    TH3D *htrk = trkpdf->GetExpValHist();

    TString htitle = "Model, sin^{2}#theta_{23}=0." + to_string( (Int_t)(th23*100) );
    if (th23 == 1.0) htitle = "Model, sin^{2}#theta_{23}=1.0";
    TString hname  = "hist_th23_" + (TString)to_string(th23);
    TH1D* h = (TH1D*)hTemplate->Clone( hname );
    h->SetNameTitle(hname, htitle);

    for (Int_t xbin = 1; xbin <= htrk->GetXaxis()->GetNbins(); xbin++) {
      for (Int_t ybin = 1; ybin <= htrk->GetYaxis()->GetNbins(); ybin++) {
	
	Double_t E  = htrk->GetXaxis()->GetBinCenter( xbin );
	Double_t ct = htrk->GetYaxis()->GetBinCenter( ybin );

	if ( ct >= 0 ) continue;

	Double_t bc = htrk->GetBinContent(xbin, ybin, 1);

	if ( bc == 0. ) continue;

	Int_t xbin = h->GetXaxis()->FindBin( -E/ct );
	h->AddBinContent( xbin, bc );

      }
    }

    delete htrk;
    hists.push_back( std::make_pair(htitle, h) );

  }

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->SetTicks();
  c1->GetPad(0)->SetLeftMargin(0.12);
  c1->GetPad(0)->SetBottomMargin(0.12);

  hists[0].second->SetTitle("-E/cosz dependence on sin^{2}#theta_{23}");

  hists[0].second->GetXaxis()->SetRangeUser(5,500);
  
  hists[0].second->GetXaxis()->SetTitle("-E_{reco}/cosz_{reco}");
  hists[0].second->GetYaxis()->SetTitle("N_{expected}");

  hists[0].second->GetXaxis()->SetTitleSize(0.05);
  hists[0].second->GetYaxis()->SetTitleSize(0.05);
  hists[0].second->GetXaxis()->SetLabelSize(0.04);
  hists[0].second->GetYaxis()->SetLabelSize(0.04);
  hists[0].second->GetXaxis()->SetTitleOffset(0.95);
  hists[0].second->GetYaxis()->SetTitleOffset(0.95);
  hists[0].second->GetXaxis()->CenterTitle();
  hists[0].second->GetYaxis()->CenterTitle();

  TLegend *leg = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  
  Int_t counter = 2;

  for (auto g: hists) {

    g.second->SetLineWidth(2);
    g.second->SetLineColor(counter++);

    if ( g.first == "Model, sin^{2}#theta_{23}=0." + to_string( (Int_t)(0.5*100) ) ) {
      g.second->SetMarkerStyle(21);
      g.second->SetLineColor(kBlack);
      g.second->Draw("Esame");
      leg->AddEntry(g.second, "Pseudo-data at sin^{2}#theta_{23}=0.5", "lep");
    }
    else {
      leg->AddEntry(g.second, g.first, "l");
    }
    
    g.second->Draw("HISTsame");
  }

  leg->Draw();

  TLatex latex;
  latex.DrawLatex(350, 50, "KM3NeT preliminary");

  
}
