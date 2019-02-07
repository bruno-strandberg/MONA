#include "DetResponse.h"
#include "SummaryEvent.h"
#include "SummaryParser.h"
#include "NMHUtils.h"

#include "FitUtil.h"
#include "FitPDF.h"

#include "RooRealVar.h"

#include "TFile.h"
#include "TGraph.h"

#include <stdexcept>

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"
#include "JTools/JRange.hh"

using namespace RooFit;
using namespace std;
using namespace JTOOLS;

/** This example application calculates the asymmetry between normal and inverted mass ordering at specified theta-23 values.

    The example uses the `track` and `shower` event selections as were optimised for ORCA115_23x9m MC production; for a different production the selections of this example are probably not optimal.

*/
int main(const int argc, const char **argv) {

  //----------------------------------------------------------
  // parse command line arguments with Jpp
  //----------------------------------------------------------
  
  TString        simdata_file;
  bool           refill_response;
  TString        outputfile;
  JRange<double> MC_E_range;
  Int_t          nebins;
  Int_t          nctbins;
  JRange<double> range_th23;
  Int_t          nsteps;
  TString        meff_file;
  
  try {

    JParser<> zap("This example application calculates the asymmetry between normal and inverted mass ordering at specified theta-23 values.");

    zap['s'] = make_field(simdata_file, "File with all summary data") =
      (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root";

    zap['r'] = make_field(refill_response, "Flag to request re-filling of the detector responses");
    zap['o'] = make_field(outputfile, "File where output histograms are written") = NMHUtils::Getcwd() + "/rootfiles/asymmetry.root";
    zap['e'] = make_field(MC_E_range, "Energy range of the MC data") = JRange<double>(1, 100);
    zap['E'] = make_field(nebins, "Number of energy bins in the energy range") = 40;
    zap['C'] = make_field(nctbins, "Number of cos-theta bins in the range -1 to 1") = 80;
    zap['a'] = make_field(range_th23, "sin^2 theta23 range") = JRange<double>(0.5, 0.5);
    zap['n'] = make_field(nsteps, "Number of steps in sin^2 theta23 range") = 1;
    zap['M'] = make_field(meff_file, "Effective mass file created by using `EffMass` class") = 
      (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";

    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  if ( range_th23.getLowerLimit() == range_th23.getUpperLimit() && nsteps != 1 ) {
    throw std::invalid_argument("ERROR! requesting theta-23 steps in range width 0.");
  }

  if ( range_th23.getLowerLimit() != range_th23.getUpperLimit() && nsteps == 1 ) {
    throw std::invalid_argument("ERROR! at least two steps are required in a range with width > 0.");
  }
  
  //----------------------------------------------------------
  // set up the detector responses
  //----------------------------------------------------------
  
  Int_t emin  = MC_E_range.getLowerLimit();
  Int_t emax  = MC_E_range.getUpperLimit();
  Int_t ctmin   = -1;
  Int_t ctmax   =  1;

  // for asym analysis always use 1 bjorken-y bin
  Int_t nbybins =  1;
  Int_t bymin   =  0;
  Int_t bymax   =  1;
  
  DetResponse track_resp (DetResponse::track, "track_resp", 
			  nebins, emin, emax, nctbins, ctmin, ctmax, nbybins, bymin, bymax);

  DetResponse shower_resp(DetResponse::shower, "shower_resp", 
			  nebins, emin, emax, nctbins, ctmin, ctmax, nbybins, bymin, bymax);

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

  if ( system("mkdir -p " + NMHUtils::Getcwd() + "/rootfiles") != 0 ) {
    cout << "WARNING! Asymmetry() creation of rootfiles/ directory returned non-zero" << endl;
  }
  TString track_resp_name  = NMHUtils::Getcwd() + "/rootfiles/asymmetry_track_response.root";
  TString shower_resp_name = NMHUtils::Getcwd() + "/rootfiles/asymmetry_shower_response.root";
  
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
  // set up the PDFs and static oscillation parameters
  //----------------------------------------------------------

  Int_t    runtime   = 3;
  Double_t fit_emin  = 1;
  Double_t fit_emax  = 100;
  Double_t fit_ctmin = -1;
  Double_t fit_ctmax = 0;
  Double_t fit_bymin = 0;
  Double_t fit_bymax = 1;

  FitUtil *fitutil = new FitUtil(runtime, track_resp.GetHist3D(), fit_emin, fit_emax,
				 fit_ctmin, fit_ctmax, fit_bymin, fit_bymax, meff_file);
  
  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_resp);  
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_resp);

  Double_t sinsqth12 = TMath::Power( TMath::Sin( 33.4 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth13 = TMath::Power( TMath::Sin( 8.42 * TMath::Pi()/180. ), 2 );
  Double_t dcp       = 0.;
  Double_t dm32      = 2.44e-3;
  Double_t dm21      = 7.53e-5;
  Double_t DM        = dm32 + 0.5*dm21;

  fitutil->GetVar("SinsqTh12")->setVal( sinsqth12 );
  fitutil->GetVar("SinsqTh13")->setVal( sinsqth13 );
  fitutil->GetVar("dcp")->setVal( dcp );
  fitutil->GetVar("Dm21")->setVal( dm21 );

  //----------------------------------------------------------
  // loop over theta23 range and calculate asymmetries
  //----------------------------------------------------------

  TGraph asym_vs_th23;
  vector<TH2D*> hlist;
  Double_t step_size = 0;

  if (nsteps > 1) {
    step_size = ( range_th23.getUpperLimit() - range_th23.getLowerLimit() ) / (nsteps-1);
  }

  for (Int_t step = 0; step < nsteps; step++) {

    Double_t sinsqth23 = range_th23.getLowerLimit() + step * step_size;
    fitutil->GetVar("SinsqTh23")->setVal( sinsqth23 );
    TString th23str = (TString)to_string(sinsqth23*1e3);
    TString prefix = "th23x1e3_" + (TString)th23str(0,3) + "_";

    //----------------------------------------------------------
    // set normal hierarchy
    //----------------------------------------------------------
    Double_t dm31 = DM + 0.5*dm21;
    fitutil->GetVar("Dm31")->setVal( dm31 );
    
    TH2D *tracks_NH  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
    TH2D *showers_NH = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
    tracks_NH->SetNameTitle(prefix+"tracks_NH",prefix+"tracks_NH");
    showers_NH->SetNameTitle(prefix+"showers_NH",prefix+"showers_NH");
    hlist.push_back(tracks_NH);
    hlist.push_back(showers_NH);

    //----------------------------------------------------------
    // set inverted hierarchy
    //----------------------------------------------------------
    dm31 = -DM + 0.5*dm21; //IH
    fitutil->GetVar("Dm31")->setVal( dm31 );
  
    TH2D *tracks_IH  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
    TH2D *showers_IH = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
    tracks_IH->SetNameTitle(prefix+"tracks_IH",prefix+"tracks_IH");
    showers_IH->SetNameTitle(prefix+"showers_IH",prefix+"showers_IH");
    hlist.push_back(tracks_IH);
    hlist.push_back(showers_IH);

    //----------------------------------------------------------
    // asymmetry calculation
    //----------------------------------------------------------
  
    auto trackasym  = NMHUtils::Asymmetry( tracks_NH, tracks_IH, prefix+"asym_tracks",    -1,1e3,-1,0);
    auto showerasym = NMHUtils::Asymmetry( showers_NH, showers_IH, prefix+"asym_showers", -1,1e3,-1,0);

    Double_t th23 = TMath::ASin( TMath::Sqrt(sinsqth23) ) * 180./TMath::Pi();
    cout << "NOTICE Asymmetry() at th23 = " << th23 << ", tracks: " << std::get<1>(trackasym) << endl;
    cout << "NOTICE Asymmetry() at th23 = " << th23 << ", showers: " << std::get<1>(showerasym) << endl;
    hlist.push_back( (TH2D*)std::get<0>(trackasym)  );
    hlist.push_back( (TH2D*)std::get<0>(showerasym) );    

    Double_t combA = TMath::Sqrt( TMath::Power( std::get<1>(trackasym), 2 ) + 
				  TMath::Power( std::get<1>(showerasym), 2 ) );

    asym_vs_th23.SetPoint( asym_vs_th23.GetN(), th23, combA );    

  }

  asym_vs_th23.SetMarkerStyle(20);
  asym_vs_th23.SetMarkerColor(kBlue);
  asym_vs_th23.SetLineWidth(2);
  asym_vs_th23.SetLineColor(kBlue);
  asym_vs_th23.SetTitle("Asymmetry vs theta-23");

  // write stuff out
  TFile fout( (TString)outputfile, "RECREATE");
  for (auto h: hlist) h->Write();
  asym_vs_th23.Write("asym_vs_th23");
  fout.Close();

  delete fitutil;

}
