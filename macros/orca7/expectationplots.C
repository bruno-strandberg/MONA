#include "NMHUtils.h"
#include "ORCA7.h"
#include "ORCA7.C"

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
  // plot the expected number of neutrinos, muons and noise 
  //------------------------------------------------------------
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare( o7.fPdfs.size() );
  Int_t pad = 1;

  for (auto P: o7.fPdfs) {

    auto pdf = P.second;
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

    cout << "Neutrinos, muons and noise in 1 y for " << pdf->GetResponse()->GetRespName() << ": " 
	 << nus->Integral() << "\t" << muons->Integral() << "\t" << noise->Integral() << endl;

  }

  //------------------------------------------------------------
  // plot  in 2D
  //------------------------------------------------------------
  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->DivideSquare( o7.fPdfs.size() );
  pad = 1;
  for (auto P: o7.fPdfs) {
    auto pdf = P.second;
    c3->cd(pad);
    pdf->GetExpValHist()->Project3D("yx")->Draw("colz");
    pad++;
  }

  //------------------------------------------------------------
  // plot the MC errors in 2D
  //------------------------------------------------------------

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->DivideSquare( o7.fPdfs.size() );
  pad = 1;

  for (auto P: o7.fPdfs) {
    auto pdf = P.second;
    c2->cd(pad);
    pdf->GetExpValErrHist()->Project3D("yx")->Draw("colz");
    pad++;
  }
 
}
