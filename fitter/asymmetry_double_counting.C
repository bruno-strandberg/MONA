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

// This script calculates the asymmetry for 10 binned data in pid_detres.

void asymmetry_double_counting() {

  bool plot = true;
  const int N_PID_CLASSES = 2;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  vector<Double_t> asym_ts(N_PID_CLASSES);
  vector<Double_t> asym_ss(N_PID_CLASSES);
  vector<Double_t> asym_ts_err(N_PID_CLASSES);
  vector<Double_t> asym_ss_err(N_PID_CLASSES);

  std::map<int, float> cut_map;
  cut_map.insert(std::make_pair(0, 0.0));
  cut_map.insert(std::make_pair(1, 0.6));
  cut_map.insert(std::make_pair(2, 1.0));

  cout << "Asymmetries per single bin: " << endl;
  for (int i = 0; i < N_PID_CLASSES; i++) {
    TString file_NO = TString::Format("./quality_detres/split_pid_NO_%.1f.root", cut_map[i]);
    TString file_IO = TString::Format("./quality_detres/split_pid_IO_%.1f.root", cut_map[i]);
    TString output  = TString::Format("./quality_detres/asym_pid_%.1f.root", cut_map[i]);
    
    TFile *f_IO = TFile::Open(file_IO, "READ");

    std::tuple<TH2D*, TH2D*> h_tuple_NO = ReadDetectorResponseFile2(file_NO);
    std::tuple<TH2D*, TH2D*> h_tuple_IO = ReadDetectorResponseFile2(file_IO);

    TH2D *h_t_NO = std::get<0>(h_tuple_NO);
    TH2D *h_s_NO = std::get<1>(h_tuple_NO);
    TH2D *h_t_IO = std::get<0>(h_tuple_IO);
    TH2D *h_s_IO = std::get<1>(h_tuple_IO);

    h_t_NO->SetName("detected_tracks_NO");
    h_s_NO->SetName("detected_showers_NO");
    h_t_IO->SetName("detected_tracks_IO");
    h_s_IO->SetName("detected_showers_IO");

    std::tuple<TH2D*, Double_t, Double_t, Double_t> asym_t = 
    NMHUtils::Asymmetry(h_t_NO, h_t_IO, TString::Format("asymmetry_track_%.1f",  cut_map[i]), 2, 80, -1, 0);
    std::tuple<TH2D*, Double_t, Double_t, Double_t> asym_s = 
    NMHUtils::Asymmetry(h_s_NO, h_s_IO, TString::Format("asymmetry_shower_%.1f", cut_map[i]), 2, 80, -1, 0);

    TH2D* h_asym_t = std::get<0>(asym_t);
    TH2D* h_asym_s = std::get<0>(asym_s);

    TCanvas *c1 = new TCanvas(TString::Format("c1_%.1f", cut_map[i]), "c1", 800, 1200); // has to be outside of the if, so that the writing below works.
    gStyle->SetPalette(kBird);
    if (plot) { 
      c1->Divide(2,3);

      vector<TH2D*> plots = {h_t_NO, h_s_NO, h_t_IO, h_s_IO, h_asym_t, h_asym_s};
      int i = 1;
      for (auto plot: plots) {
        c1->cd(i);
        plot->Draw("colz");
        plot->GetXaxis()->SetTitle("Energy [GeV]");
        plot->GetYaxis()->SetTitle("cos(#theta)");
        plot->GetXaxis()->SetRangeUser(3,100);
        plot->GetYaxis()->SetRangeUser(-1,0);
        c1->cd(i)->SetLogx();

        i++;
      }
    }

    asym_ts[i] = std::get<1>(asym_t);
    asym_ss[i] = std::get<1>(asym_s);
    asym_ts_err[i] = std::get<2>(asym_t);
    asym_ss_err[i] = std::get<2>(asym_s);

    cout << "Asymmetry for tracks : " << asym_ts[i] << " +- " << asym_ts_err[i] <<  endl;
    cout << "Asymmetry for showers: " << asym_ss[i] << " +- " << asym_ss_err[i] << endl;
  
    TFile fout(output, "RECREATE");
    c1->Write();
    h_t_NO->Write();
    h_t_IO->Write();
    h_s_NO->Write();
    h_s_IO->Write();
    h_asym_t->Write();
    h_asym_s->Write();
    fout.Write();
    fout.Close();

    if (plot) { // To keep the memory low (esp. over tunnels) close the canvas.
      c1->Close(); 
      gSystem->ProcessEvents(); 
    }
  }

  // DOUBLE COUNTING
  cout << "Double counting: " << endl;

  std::tuple<Double_t, Double_t> track_value_squared  = NMHUtils::SquaredSumErrorProp(asym_ts, asym_ts_err);
  std::tuple<Double_t, Double_t> shower_value_squared = NMHUtils::SquaredSumErrorProp(asym_ss, asym_ss_err);

  std::tuple<Double_t, Double_t> total_value_squared = 
    NMHUtils::SquaredSumErrorProp({std::get<0>(track_value_squared), std::get<0>(shower_value_squared)},
                                  {std::get<1>(track_value_squared), std::get<1>(shower_value_squared)});

  cout << "Asymmetry total for tracks : " << std::get<0>(track_value_squared)  << " +- " << std::get<1>(track_value_squared)  << endl;
  cout << "Asymmetry total for showers: " << std::get<0>(shower_value_squared) << " +- " << std::get<1>(shower_value_squared) << endl;
  cout << "Asymmetry total combined   : " << std::get<0>(total_value_squared)  << " +- " << std::get<1>(total_value_squared) << endl;

}
