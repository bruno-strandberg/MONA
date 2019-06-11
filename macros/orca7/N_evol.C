#include "ORCA7.C"

void N_evol() {

  ORCA7 o7( kTRUE );
  FitUtilWsyst *fu = o7.fFitUtil;
  FitPDF *trkpdf = o7.fPdfs["trk"];
  FitPDF *midpdf = o7.fPdfs["mid"];
  FitPDF *shwpdf = o7.fPdfs["shw"];

  TH3D *htrk = trkpdf->GetExpValHist();
  TH3D *hshw = shwpdf->GetExpValHist();
  TH3D *hmid = midpdf->GetExpValHist();

  Double_t N_start_trk = htrk->Integral();
  Double_t N_start_shw = hshw->Integral();
  Double_t N_start_mid = hmid->Integral();
  Double_t th23_start  = fu->GetVar("SinsqTh23")->getVal();
  delete htrk;
  delete hshw;
  delete hmid;

  TGraph *gtrk = new TGraph();
  TGraph *gshw = new TGraph();
  TGraph *gmid = new TGraph();
  TGraph *gtot = new TGraph();

  for (Double_t th23 = 0; th23 <= 1; th23 += 0.05) {

    fu->GetVar("SinsqTh23")->setVal( th23 );
    TH3D *htrk = trkpdf->GetExpValHist();
    TH3D *hshw = shwpdf->GetExpValHist();
    TH3D *hmid = midpdf->GetExpValHist();

    gtrk->SetPoint( gtrk->GetN(), th23, htrk->Integral() );
    gshw->SetPoint( gshw->GetN(), th23, hshw->Integral() );
    gmid->SetPoint( gmid->GetN(), th23, hmid->Integral() );
    gtot->SetPoint( gtot->GetN(), th23, htrk->Integral()+hshw->Integral()+hmid->Integral() );

    delete htrk;
    delete hshw;
    delete hmid;

  }

  vector< std::pair<TString,TGraph*> > graphs = { {"PID [0.7,1.0]" , gtrk}, {"PID [0.3, 0.7)", gmid}, 
						  {"PID [0.0, 0.3)", gshw}, {"PID [0.0, 1.0]", gtot} };
  for (auto g: graphs) {
    g.second->SetTitle(g.first);
    g.second->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    g.second->GetYaxis()->SetTitle("Total expected events");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare(4);
  c1->cd(1);
  gtrk->Draw();
  c1->cd(2);
  gmid->Draw();
  c1->cd(3);
  gshw->Draw();
  c1->cd(4);
  gtot->Draw();

}
