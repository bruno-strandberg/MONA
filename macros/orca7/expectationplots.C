#include "NMHUtils.h"
#include "ORCA7.h"

/* 
   In this script ORCA 7-line events are distributed to three PID classes. This script creates plots that depict the number of detected neutrinos, muons and noise in 1 year with ORCA 7-line detector in energy dimension for each PID class. The relevant numbers are also printed to terminal.
*/

void expectationplots() {

  //------------------------------------------------------------
  // init the ORCA class where responses and some other common
  // variables are set
  //------------------------------------------------------------
  ORCA7 o7(kTRUE);

  //------------------------------------------------------------
  // initialise the pdf's
  //------------------------------------------------------------

  FitUtil futil(o7.f_F_runtime, o7.fResps[0]->GetHist3D(), o7.f_F_emin, o7.f_F_emax, o7.f_F_ctmin, o7.f_F_ctmax, o7.f_F_bymin, o7.f_F_bymax, o7.fEffmF);
  vector<FitPDF*> pdfs;
  for (auto R: o7.fResps) {
    TString name = "pdf_" + R->Get_RespName();
    pdfs.push_back( new FitPDF(name, name, &futil, R) );
  }

  //------------------------------------------------------------
  // plot the expected number of neutrinos, muons and noise 
  //------------------------------------------------------------
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare( pdfs.size() );
  Int_t pad = 1;

  for (auto pdf: pdfs) {

    TH1D* nus   = (TH1D*)pdf->GetExpValHist()->Project3D("x");
    TH1D* muons = (TH1D*)pdf->GetResponse()->GetHistAtmMu1y()->Project3D("x");
    TH1D* noise = (TH1D*)pdf->GetResponse()->GetHistNoise1y()->Project3D("x");

    nus->SetLineColor(kRed);
    muons->SetLineColor(kBlack);
    noise->SetLineColor(kGreen);

    nus->GetXaxis()->SetTitle("Reco Energy [GeV]");
    nus->GetYaxis()->SetTitle("Events [1 year]");

    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->AddEntry(nus, "neutrinos", "l");
    leg->AddEntry(muons, "muons", "l");
    leg->AddEntry(noise, "noise", "l");

    c1->cd(pad);
    nus->Draw("HIST");
    noise->Draw("HISTsame");
    muons->Draw("HISTsame");
    leg->Draw();
    pad++;

    cout << "Neutrinos, muons and noise in 1 y for " << pdf->GetResponse()->Get_RespName() << ": " 
	 << nus->Integral() << "\t" << muons->Integral() << "\t" << noise->Integral() << endl;

  }

  //------------------------------------------------------------
  // plot  in 2D
  //------------------------------------------------------------
  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->DivideSquare( pdfs.size() );
  pad = 1;
  for (auto pdf: pdfs) {
    c3->cd(pad);
    pdf->GetExpValHist()->Project3D("yx")->Draw("colz");
    pad++;
  }

  //------------------------------------------------------------
  // plot the MC errors in 2D
  //------------------------------------------------------------

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->DivideSquare( pdfs.size() );
  pad = 1;

  for (auto pdf: pdfs) {
    c2->cd(pad);
    pdf->GetExpValErrHist()->Project3D("yx")->Draw("colz");
    pad++;
  }
 
}
