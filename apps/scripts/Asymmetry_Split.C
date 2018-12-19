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

void Asymmetry_Split() {

  bool b_plot = false;
  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = TString::Format("./pid_detres/pid_binning_%i/", N_PID_CLASSES);

  vector<Double_t> asym_ts(N_PID_CLASSES);
  vector<Double_t> asym_ss(N_PID_CLASSES);
  vector<Double_t> asym_ms(N_PID_CLASSES);
  vector<Double_t> asym_ts_err(N_PID_CLASSES);
  vector<Double_t> asym_ss_err(N_PID_CLASSES);
  vector<Double_t> asym_ms_err(N_PID_CLASSES);

  cout << "Asymmetries per single bin: " << endl;
  for (int i = 0; i < N_PID_CLASSES; i++) {
    TString file_NO = filefolder + TString::Format("split_expected_evts_NO_%.2f.root", PID_step * i);
    TString file_IO = filefolder + TString::Format("split_expected_evts_IO_%.2f.root", PID_step * i);
    TString output  = filefolder + TString::Format("asymmetry_split_%.2f.root", PID_step * i);
    
    TFile *f_NO  = TFile::Open(file_NO, "READ");
    TFile *f_IO  = TFile::Open(file_IO, "READ");

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
    NMHUtils::Asymmetry(h_t_NO, h_t_IO, TString::Format("asymmetry_track_%.2f",  i * PID_step), 2, 80, -1, 0);
    std::tuple<TH2D*, Double_t, Double_t, Double_t> asym_s = 
    NMHUtils::Asymmetry(h_s_NO, h_s_IO, TString::Format("asymmetry_shower_%.2f", i * PID_step), 2, 80, -1, 0);
    std::tuple<TH2D*, Double_t, Double_t, Double_t> asym_m = 
    NMHUtils::Asymmetry(h_m_NO, h_m_IO, TString::Format("asymmetry_mc_%.2f", i * PID_step), 2, 80, -1, 0);

    TH2D* h_asym_t = std::get<0>(asym_t);
    TH2D* h_asym_s = std::get<0>(asym_s);
    TH2D* h_asym_m = std::get<0>(asym_m);

    TCanvas *c1 = new TCanvas(TString::Format("c1_%.2f", PID_step * i), "c1", 1200, 1200); // has to be outside of the if, so that the writing below works.
    gStyle->SetPalette(kBird);
    if (b_plot) { 
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

    asym_ts[i] = std::get<1>(asym_t);
    asym_ss[i] = std::get<1>(asym_s);
    asym_ms[i] = std::get<1>(asym_m);
    asym_ts_err[i] = std::get<2>(asym_t);
    asym_ss_err[i] = std::get<2>(asym_s);
    asym_ms_err[i] = std::get<2>(asym_m);

    cout << "Asymmetry for tracks : " << asym_ts[i] << " +- " << asym_ts_err[i] << " (" << 100*asym_ts_err[i]/asym_ts[i] << "%)" << endl;
    cout << "Asymmetry for showers: " << asym_ss[i] << " +- " << asym_ss_err[i] << " (" << 100*asym_ss_err[i]/asym_ss[i] << "%)" << endl;
    cout << "Asymmetry for mc     : " << asym_ms[i] << " +- " << asym_ms_err[i] << " (" << 100*asym_ms_err[i]/asym_ms[i] << "%)" << endl;
  
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

    if (b_plot) { // To keep the memory low (esp. over tunnels) close the canvas.
      c1->Close(); 
      gSystem->ProcessEvents(); 
    }
  }

  // SINGLE COUNTING WITH q < 0.6
  cout << "Single counting (q<0.6): " << endl;
  int offset = N_PID_CLASSES * 0.6;

  // Wow, slicing a vector in cpp is hard...
  std::vector<Double_t> asym_ts_q6(asym_ts.begin()+offset, asym_ts.end());
  std::vector<Double_t> asym_ts_err_q6(asym_ts_err.begin()+offset, asym_ts_err.end());
  std::vector<Double_t> asym_ss_q6(asym_ss.begin(), asym_ss.begin()+offset);
  std::vector<Double_t> asym_ss_err_q6(asym_ss_err.begin(), asym_ss_err.begin()+offset);

  std::tuple<Double_t, Double_t> track_value_squared_q6  = 
    NMHUtils::SquaredSumErrorProp(asym_ts_q6, asym_ts_err_q6);
  std::tuple<Double_t, Double_t> shower_value_squared_q6 = 
    NMHUtils::SquaredSumErrorProp(asym_ss_q6, asym_ss_err_q6);
  std::tuple<Double_t, Double_t> total_value_squared_q6  = 
    NMHUtils::SquaredSumErrorProp({std::get<0>(track_value_squared_q6), std::get<0>(shower_value_squared_q6)},
                                  {std::get<1>(track_value_squared_q6), std::get<1>(shower_value_squared_q6)});

  Double_t track_value_q6  = std::get<0>(track_value_squared_q6);
  Double_t shower_value_q6 = std::get<0>(shower_value_squared_q6);
  Double_t total_value_q6  = std::get<0>(total_value_squared_q6);
  Double_t track_error_q6  = std::get<1>(track_value_squared_q6);
  Double_t shower_error_q6 = std::get<1>(shower_value_squared_q6);
  Double_t total_error_q6  = std::get<1>(total_value_squared_q6);

  cout << "Asymmetry total for tracks : " << track_value_q6  << " +- " << track_error_q6  << " (" << 100*track_error_q6/track_value_q6    << "%)" << endl;
  cout << "Asymmetry total for showers: " << shower_value_q6 << " +- " << shower_error_q6 << " (" << 100*shower_error_q6/shower_value_q6  << "%)" << endl;
  cout << "Asymmetry total combined   : " << total_value_q6  << " +- " << total_error_q6  << " (" << 100*total_error_q6/total_value_q6    << "%)" << endl;

  // SINGLE COUNTING WITH q < 0.9
  cout << "Single counting (q<0.9): " << endl;
  offset = 0.9 * N_PID_CLASSES;

  std::vector<Double_t> asym_ts_q9(asym_ts.begin()+offset, asym_ts.end());
  std::vector<Double_t> asym_ts_err_q9(asym_ts_err.begin()+offset, asym_ts_err.end());
  std::vector<Double_t> asym_ss_q9(asym_ss.begin(), asym_ss.begin()+offset);
  std::vector<Double_t> asym_ss_err_q9(asym_ss_err.begin(), asym_ss_err.begin()+offset);

  std::tuple<Double_t, Double_t> track_value_squared_q9  = NMHUtils::SquaredSumErrorProp(asym_ts_q9, asym_ts_err_q9);
  std::tuple<Double_t, Double_t> shower_value_squared_q9 = NMHUtils::SquaredSumErrorProp(asym_ss_q9, asym_ss_err_q9);
  std::tuple<Double_t, Double_t> total_value_squared_q9  = 
    NMHUtils::SquaredSumErrorProp({std::get<0>(track_value_squared_q9), std::get<0>(shower_value_squared_q9)},
                                  {std::get<1>(track_value_squared_q9), std::get<1>(shower_value_squared_q9)});

  Double_t track_value_q9  = std::get<0>(track_value_squared_q9);
  Double_t shower_value_q9 = std::get<0>(shower_value_squared_q9);
  Double_t total_value_q9  = std::get<0>(total_value_squared_q9);
  Double_t track_error_q9  = std::get<1>(track_value_squared_q9);
  Double_t shower_error_q9 = std::get<1>(shower_value_squared_q9);
  Double_t total_error_q9  = std::get<1>(total_value_squared_q9);

  cout << "Asymmetry total for tracks : " << track_value_q9  << " +- " << track_error_q9  << " (" << 100*track_error_q9/track_value_q9    << "%)" << endl;
  cout << "Asymmetry total for showers: " << shower_value_q9 << " +- " << shower_error_q9 << " (" << 100*shower_error_q9/shower_value_q9  << "%)" << endl;
  cout << "Asymmetry total combined   : " << total_value_q9  << " +- " << total_error_q9  << " (" << 100*total_error_q9/total_value_q9    << "%)" << endl;

}
