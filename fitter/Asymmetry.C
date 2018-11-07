#include "DetResponse.h"
#include "SummaryEvent.h"
#include "SummaryParser.h"
#include "NMHUtils.h"

#include "FitUtil.h"
#include "FitPDF.h"

#include "RooRealVar.h"

#include "TFile.h"

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

using namespace RooFit;
using namespace std;

int main(int argc, char **argv) {

  //----------------------------------------------------------
  // parse command line arguments with Jpp
  //----------------------------------------------------------
  
  string        simdata_file;
  bool          refill_response;
  string        outputfile;
  string        effmh_elecCC;
  string        effmh_muonCC;
  string        effmh_tauCC;
  string        effmh_elecNC;
  
  try {

    JParser<> zap("Program to calculate asymmetries with standard track-shower separation.");

    zap['s'] = make_field(simdata_file, "File with all summary data") =
      (string)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";

    zap['r'] = make_field(refill_response, "Flag to request re-filling of the detector responses");
    zap['o'] = make_field(outputfile, "File where output histograms are written") = "rootfiles/asymmetry.root";

    zap['w'] = make_field(effmh_elecCC, "Eff mass histograms for elec-CC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_elec_CC.root";

    zap['x'] = make_field(effmh_muonCC, "Eff mass histograms for muon-CC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_muon_CC.root";

    zap['y'] = make_field(effmh_tauCC , "Eff mass histograms for tau-CC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_tau_CC.root";

    zap['z'] = make_field(effmh_elecNC, "Eff mass histograms for elec-NC") =
      (string)getenv("NMHDIR") + "/data/eff_mass/EffMhists_elec_NC.root";    
    zap(argc, argv);
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }
  
  //----------------------------------------------------------
  // set up the detector responses
  //----------------------------------------------------------
  
  Int_t fENB   =  40;
  Int_t fEmin  =   1;
  Int_t fEmax  = 100;
  Int_t fCtNB  =  40;
  Int_t fCtmin =  -1;
  Int_t fCtmax =   1;
  Int_t fByNB  =   1;
  Int_t fBymin =   0;
  Int_t fBymax =   1;
  
  DetResponse track_resp (DetResponse::track, "track_resp", 
			  fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  DetResponse shower_resp(DetResponse::shower, "shower_resp", 
			  fENB, fEmin, fEmax, fCtNB, fCtmin, fCtmax, fByNB, fBymin, fBymax);

  track_resp.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_resp.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_resp.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );    

  shower_resp.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  shower_resp.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  shower_resp.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  shower_resp.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  shower_resp.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  TString track_resp_name  = "rootfiles/track_response.root";
  TString shower_resp_name = "rootfiles/shower_response.root";
  
  if ( !NMHUtils::FileExists(track_resp_name) || !NMHUtils::FileExists(shower_resp_name) || refill_response ) {

    cout << "NOTICE Asymmetry() (Re)filling responses" << endl;

    SummaryParser sp(simdata_file);

    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

      if (i % 100000 == 0) cout << "Event: " << i << endl;

      SummaryEvent *evt = sp.GetEvt(i);
      track_resp.Fill(evt);
      shower_resp.Fill(evt);

    }

    track_resp.WriteToFile(track_resp_name);
    shower_resp.WriteToFile(shower_resp_name);

  }
  else {

    cout << "NOTICE Asymmetry() Reading in responses" << endl;

    track_resp.ReadFromFile(track_resp_name);
    shower_resp.ReadFromFile(shower_resp_name);
  }

  cout << "NOTICE Asymmetry() Responses ready" << endl;

  //----------------------------------------------------------
  // set up the PDFs
  //----------------------------------------------------------
  
  FitUtil *fitutil = new FitUtil(3, track_resp.GetHist3D(),
  				 1, 100, -1, 0, 0, 1, effmh_elecCC, effmh_muonCC, effmh_tauCC, effmh_elecNC);

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_resp);  
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_resp);

  Double_t sinsqth12 = TMath::Power( TMath::Sin( 33.4 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth13 = TMath::Power( TMath::Sin( 8.42 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth23 = TMath::Power( TMath::Sin( 45   * TMath::Pi()/180. ), 2 );
  Double_t dcp       = 0.;
  Double_t dm32      = 2.44e-3;
  Double_t dm21      = 7.53e-5;
  Double_t DM        = dm32 + 0.5*dm21;

  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh12") )->setVal( sinsqth12 );
  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh13") )->setVal( sinsqth13 );
  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh23") )->setVal( sinsqth23 );
  ( (RooRealVar*)fitutil->GetSet().find("dcp") )->setVal( dcp );
  ( (RooRealVar*)fitutil->GetSet().find("Dm21") )->setVal( dm21 );

  //----------------------------------------------------------
  // set normal hierarchy
  //----------------------------------------------------------
  fitutil->SetNOlims();
  Double_t dm31 = DM + 0.5*dm21;
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setVal( dm31 );
  
  TH2D *tracks_NH  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
  TH2D *showers_NH = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
  tracks_NH->SetNameTitle("tracks_NH","tracks_NH");
  showers_NH->SetNameTitle("showers_NH","showers_NH");

  //----------------------------------------------------------
  // set inverted hierarchy
  //----------------------------------------------------------
  dm31 = -DM + 0.5*dm21; //IH
  fitutil->SetIOlims();
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setVal( dm31 );
  
  TH2D *tracks_IH  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
  TH2D *showers_IH = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
  tracks_IH->SetNameTitle("tracks_IH","tracks_IH");
  showers_IH->SetNameTitle("showers_IH","showers_IH");

  //----------------------------------------------------------
  // asymmetry calculation
  //----------------------------------------------------------
  
  auto trackasym  = NMHUtils::Asymmetry( tracks_NH, tracks_IH, "asym_tracks", 1,100,-1,0);
  auto showerasym = NMHUtils::Asymmetry( showers_NH, showers_IH, "asym_showers", 1,100,-1,0);

  cout << "NOTICE Asymmetry() track  asym: " << std::get<1>(trackasym) << endl;
  cout << "NOTICE Asymmetry() shower asym: " << std::get<1>(showerasym) << endl;
  
  // write stuff out
  TFile fout( (TString)outputfile, "RECREATE");
  tracks_NH->Write("tracks_NH");
  showers_NH->Write("showers_NH");
  tracks_IH->Write("tracks_IH");
  showers_IH->Write("showers_IH");
  std::get<0>(trackasym)->Write();
  std::get<0>(showerasym)->Write();
  fout.Close();

}

