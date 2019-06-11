#include "ORCA7.C"
#include "NMHUtils.h"

void EoverZ() {
  
  ORCA7 o7( kTRUE );
  FitUtilWsyst *fu = o7.fFitUtil;
  FitPDF *trkpdf = o7.fPdfs["trk"];

  // find the bin limits
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
  vector< TH1D* > hists, hists_norm;

  for (Double_t th23 = 0; th23 <= 1; th23 += 0.05) {

    fu->GetVar("SinsqTh23")->setVal( th23 );
    TH3D *htrk = trkpdf->GetExpValHist();

    TString hname      = "hist_th23_" + (TString)to_string(th23);
    TString hname_norm = "hist_th23_norm_" + (TString)to_string(th23);

    TH1D* h = (TH1D*)hTemplate->Clone( hname );
    h->SetNameTitle(hname, hname);

    TH1D* h_norm = (TH1D*)hTemplate->Clone( hname_norm );
    h_norm->SetNameTitle(hname_norm, hname_norm);

    for (Int_t xbin = 1; xbin <= htrk->GetXaxis()->GetNbins(); xbin++) {
      for (Int_t ybin = 1; ybin <= htrk->GetYaxis()->GetNbins(); ybin++) {
	
	Double_t E  = htrk->GetXaxis()->GetBinCenter( xbin );
	Double_t ct = htrk->GetYaxis()->GetBinCenter( ybin );

	if ( ct >= 0 ) continue;

	Double_t bc = htrk->GetBinContent(xbin, ybin, 1);

	if ( bc == 0. ) continue;

	Int_t xbin = h->GetXaxis()->FindBin( -E/ct );
	h->AddBinContent( xbin, bc );
	h_norm->AddBinContent( xbin, bc/htrk->Integral() );

      }
    }

    delete htrk;
    hists.push_back( h );
    hists_norm.push_back( h_norm );

  }

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1);

  Int_t counter = 1;

  hists[0]->GetXaxis()->SetTitle("-E_{reco}/cosz_{reco}");
  hists[0]->GetYaxis()->SetTitle("N_{expected}");
  hists[0]->SetTitle("Scan in sin^{2}#theta_{23}");
  hists[0]->GetXaxis()->SetRangeUser(5,500);

  for (auto g: hists) {
    g->SetLineWidth(2);
    g->SetLineColor(counter++);
    g->Draw("same");
  }

  TCanvas *c2 = new TCanvas("c2","c2",1);

  hists_norm[0]->GetXaxis()->SetTitle("-E_{reco}/cosz_{reco}");
  hists_norm[0]->GetYaxis()->SetTitle("Event density");
  hists_norm[0]->SetTitle("Scan in sin^{2}#theta_{23}");
  hists_norm[0]->GetXaxis()->SetRangeUser(5,500);

  counter = 1;
  for (auto g: hists_norm) {
    g->SetLineWidth(2);
    g->SetLineColor(counter++);
    g->Draw("same");
  }

}
