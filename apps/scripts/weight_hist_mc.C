#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"

#include "DetResponse.h"
#include "FitFunction.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

void weight_hist_mc(Double_t ebin=10, Double_t ctbin=-0.8) {
  TString filefolder = "./default_detres/";

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  DetResponse mc_response(DetResponse::shower, "mc_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  mc_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
  mc_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );

  mc_response.ReadFromFile(filefolder + "mc_response.root");

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------

  Double_t bybin = 0.5; // Theres only one by

  auto effmass_folder = (TString)getenv("NMHDIR") + "/data/eff_mass/";
  FitFunction mfitf(&mc_response, 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  std::vector<TrueB> mc_weights;
  mc_weights  = mc_response.GetBinWeights(ebin, ctbin, bybin);
 

  TH1D h_weights_mc = TH1D("h_w_m_q", "Weights mc",  100, 0, 2); 
  for (int j = 0; j < mc_weights.size(); j++) {
    h_weights_mc.Fill(mc_weights[j].fW);
  }
  
  TString output = TString::Format("weights_mc_E%.1f_ct%.1f.root", ebin, ctbin);
  TFile fout(filefolder + output, "RECREATE");
  h_weights_mc.Write();
  fout.Close();
  cout << "NOTICE: Written " << output << endl;
}
