#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "DetResponse.h"
#include "FitFunction.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "HelperFunctions.C"

#include <iostream>
using namespace std;

// This script calculates the asymmetry for unbinned data in default_detres/

void Asymmetry_Default() {

  bool plot = true;

  TString filefolder = "./default_detres/";

  TString file_NO = filefolder + "default_expectated_evts_NO.root";
  TString file_IO = filefolder + "default_expectated_evts_IO.root";
  TString output  = filefolder + "asymmertry_default.root";
  
  TFile *f_IO = TFile::Open(file_IO, "READ");

  std::tuple<TH2D*, TH2D*, TH2D*> h_tuple_NO = ReadDetectorResponseFile(file_NO);
  std::tuple<TH2D*, TH2D*, TH2D*> h_tuple_IO = ReadDetectorResponseFile(file_IO);

  TH2D *h_t_NO = std::get<0>(h_tuple_NO);
  TH2D *h_s_NO = std::get<1>(h_tuple_NO);
  TH2D *h_m_NO = std::get<2>(h_tuple_NO);
  TH2D *h_t_IO = std::get<0>(h_tuple_IO);
  TH2D *h_s_IO = std::get<1>(h_tuple_IO);
  TH2D *h_m_IO = std::get<2>(h_tuple_IO);

  h_t_NO->SetName("detected_tracks_NO");
  h_s_NO->SetName("detected_showers_NO");
  h_m_NO->SetName("detected_mc_NO");
  h_t_IO->SetName("detected_tracks_IO");
  h_s_IO->SetName("detected_showers_IO");
  h_m_IO->SetName("detected_mc_IO");

  std::tuple<TH2D*, Double_t, Double_t, Double_t> asym_t = 
  NMHUtils::Asymmetry(h_t_NO, h_t_IO, "asymmetry_track", 2, 80, -1, 0);
  std::tuple<TH2D*, Double_t, Double_t, Double_t> asym_s = 
  NMHUtils::Asymmetry(h_s_NO, h_s_IO, "asymmetry_shower", 2, 80, -1, 0);
  std::tuple<TH2D*, Double_t, Double_t, Double_t> asym_m = 
  NMHUtils::Asymmetry(h_m_NO, h_m_IO, "asymmetry_mc", 2, 80, -1, 0);

  TH2D* h_asym_t = std::get<0>(asym_t);
  TH2D* h_asym_s = std::get<0>(asym_s);
  TH2D* h_asym_m = std::get<0>(asym_m);

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1200); // has to be outside of the if, so that the writing below works.
  gStyle->SetPalette(kBird);
  if (plot) { 
    c1->Divide(3,3);

    vector<TH2D*> plots = {h_t_NO, h_s_NO, h_m_NO, h_t_IO, h_s_IO, h_m_IO, h_asym_t, h_asym_s, h_asym_m};
    int i = 1;
    for (auto plot: plots) {
      c1->cd(i);
      plot->Draw("colz");
      plot->GetXaxis()->SetTitle("Energy [GeV]");
      plot->GetYaxis()->SetTitle("cos(#theta)");
      plot->GetXaxis()->SetRangeUser(3,100);
      plot->GetYaxis()->SetRangeUser(-1,0);
      if (i>6) { plot->GetZaxis()->SetRangeUser(-1,1); }
      c1->cd(i)->SetLogx();

      i++;
    }
  }

  TFile fout(output, "RECREATE");
  c1->Write();
  h_t_NO->Write();
  h_t_IO->Write();
  h_s_NO->Write();
  h_s_IO->Write();
  h_m_NO->Write();
  h_m_IO->Write();
  h_asym_t->Write();
  h_asym_s->Write();
  h_asym_m->Write();
  fout.Write();
  fout.Close();

  c1->cd(1);
  c1->SaveAs(filefolder + "asym_plots.pdf");

  if (plot) { // To keep the memory low (esp. over tunnels) close the canvas.
    c1->Close(); 
    gSystem->ProcessEvents(); 
  }

  double asym_t_value = std::get<1>(asym_t);
  double asym_s_value = std::get<1>(asym_s); 
  double asym_m_value = std::get<1>(asym_m); 
  double asym_t_err   = std::get<2>(asym_t);
  double asym_s_err   = std::get<2>(asym_s);
  double asym_m_err   = std::get<2>(asym_m);

  std::tuple<Double_t, Double_t> sq_sum = NMHUtils::SquaredSumErrorProp({asym_t_value, asym_s_value}, {asym_t_err, asym_s_err});
  double asym_tot_value = std::get<0>(sq_sum);
  double asym_tot_err   = std::get<1>(sq_sum);

  cout << "Asymmetry total for tracks : " << asym_t_value   << " +- " << asym_t_err   << " (" << 100*asym_t_err/asym_t_value     << "%)" << endl;
  cout << "Asymmetry total for showers: " << asym_s_value   << " +- " << asym_s_err   << " (" << 100*asym_s_err/asym_s_value     << "%)" << endl;
  cout << "Asymmetry for mc           : " << asym_m_value   << " +- " << asym_m_err   << " (" << 100*asym_m_err/asym_m_value     << "%)" << endl;
  cout << "Asymmetry total combined   : " << asym_tot_value << " +- " << asym_tot_err << " (" << 100*asym_tot_err/asym_tot_value << "%)" << endl;

}
