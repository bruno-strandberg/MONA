// mona headers
#include "FitUtil.h"
#include "FitPDF.h"
#include "DetResponse.h"
#include "SummaryParser.h"
#include "NMHUtils.h"

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

// ROOT headers
#include "TStopwatch.h"
#include "TFile.h"

#include <iostream>
using namespace std;

int main(const int argc, const char **argv) {

  TString datafile;
  TString effmfile;
  Int_t   nsteps;
  
  try {

    JParser<> zap("Application to test the performance of FitUtil software");

    zap['s'] = make_field(datafile, "File with all summary data") =
      (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_ORCA115_20x9m_ECAP1218.root";

    zap['e'] = make_field(effmfile, "Effective mass file") =
      (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA115_20x9m_ECAP1218.root";

    zap['n'] = make_field(nsteps, "Number of osc configurations explored") = 5;

    if ( zap.read(argc, argv)!= 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  DetResponse trk(DetResponse::track , "trk", 24, 1, 100, 40, -1, 1, 1, 0, 1);
  trk.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, kTRUE);
  trk.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, kTRUE);

  DetResponse shw(DetResponse::shower, "shw", 24, 1, 100, 40, -1, 1, 1, 0, 1);
  shw.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5, kTRUE);
  shw.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), 0.6, kTRUE);

  TString trkname = "performancestudy_trkresp.root";
  TString shwname = "performancestudy_shwresp.root";

  if ( NMHUtils::FileExists(trkname) && NMHUtils::FileExists(shwname) ) {
    trk.ReadFromFile( trkname );
    shw.ReadFromFile( shwname );
  }
  else {
  
    SummaryParser sp( datafile );
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i%200000 == 0) cout << "Event: " << i << "/" << sp.GetTree()->GetEntries() << endl;
      trk.Fill( sp.GetEvt(i) );
      shw.Fill( sp.GetEvt(i) );
    }
    trk.WriteToFile(trkname);
    shw.WriteToFile(shwname);
    
  }
  cout << "NOTICE PerformanceStudy: responses ready" << endl;

  
  FitUtil FU(3, trk.GetHist3D(), 3, 75, -1, -1e-5, 0, 1, effmfile);
  FU.SetNOlims();
  
  FitPDF pdftrk("pdftrk", "pdftrk", &FU, &trk);
  FitPDF pdfshw("pdfshw", "pdfshw", &FU, &shw);

  TH3D* hb = trk.GetHist3D();

  TStopwatch timer;
  
  for (Int_t N = 0; N < nsteps; N++) {

    cout << "NOTICE PerformanceStudy: sample " << N << endl; 
    
    // randomise oscillation parameters
    FU.GetVar("SinsqTh23")->randomize();
    FU.GetVar("Dm31")->randomize();

    for (Int_t xbin = 1; xbin <= hb->GetXaxis()->GetNbins(); xbin++) {
      for (Int_t ybin = 1; ybin <= hb->GetYaxis()->GetNbins(); ybin++) {

	Double_t E  = hb->GetXaxis()->GetBinCenter(xbin);
	Double_t ct = hb->GetYaxis()->GetBinCenter(ybin);
	Double_t by = 0.5;

	FU.RecoEvts( E, ct, by, &trk, pdftrk.GetProxyMap() );
	FU.RecoEvts( E, ct, by, &shw, pdfshw.GetProxyMap() );
	
      }
    }
    
  }

  cout << "NOTICE PerformanceStudy: Total time [s]: " << timer.RealTime() << endl;

  TFile fout("performancestudy.root","RECREATE");
  pdftrk.GetExpValHist()->Project3D("yx")->Write("trk");
  pdfshw.GetExpValHist()->Project3D("yx")->Write("shw");
  fout.Close();
  
}
