#include "DetResponse.h"
#include "SummaryParser.h"
#include "FitUtil.h"
#include "FitPDF.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"

/** This macro demonstrates the usage of simultaneous fits with MONA/RooFit. See README.md to see what is required to run this macro.*/
void sim_fit(TString dataf = (TString)getenv("MONADIR") + (TString)"/data/ORCA_MCsummary_SEv2_ORCA115_20x9m_ECAP190222.root",
	     TString effmf = (TString)getenv("MONADIR") + (TString)"/data/eff_mass/EffMass_ORCA115_20x9m_ECAP190222.root") {

  //==================================================================
  // initialise a detector responses for tracks and showers & fill them with summary data
  //==================================================================

  DetResponse trkR(DetResponse::track, "trkR", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  trkR.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.7, kTRUE ); // ask PID track score to be larger than > 0.7
  trkR.AddCut( &SummaryEvent::Get_track_ql2      , std::greater<double>(), 0.5, kTRUE ); // ask for ql2 (gandalf_loose_is_selected), see `MONA/apps/data_sorting/ECAP190222_20m_to_MONA.C`

  DetResponse shwR(DetResponse::shower, "shwR", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  shwR.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), 0.3, kTRUE ); // ask PID track score to be less than > 0.3
  shwR.AddCut( &SummaryEvent::Get_shower_ql2     , std::greater<double>()   , 0.5, kTRUE ); // ask for ql2 (dusj_is_selected), see `MONA/apps/data_sorting/ECAP190222_20m_to_MONA.C`

  vector<DetResponse*> drs = { &trkR, &shwR };

  SummaryParser sp( dataf );

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i%500000 == 0) cout << "Event: " << i << "/" << sp.GetTree()->GetEntries() << endl;	
    for (auto R: drs) R->Fill( sp.GetEvt(i) );
  }

  //==================================================================
  // initialise pdf's for tracks and showers and set up simultaneous fitting
  //==================================================================

  FitUtil fu(3, trkR.GetHist3DTrue(), trkR.GetHist3DReco(), 3, 80, -1, -1e-5, 0, 1, effmf);
  FitPDF trkpdf("trkpdf", "trkpdf", &fu, &trkR);
  FitPDF shwpdf("shwpdf", "shwpdf", &fu, &shwR);

  // create pseudo-experiments corresponding to 3years of data taking (run time determined at FitUtil construction)
  TH3D* trkdata = trkpdf.SimplePseudoExp("trkdata");
  TH3D* shwdata = shwpdf.SimplePseudoExp("shwdata");

  // this map and the categories helps RooFit with simultaneous fitting
  std::map<string, TH1*> hmap = { {"TRK", (TH1*)trkdata}, {"SHW", (TH1*)shwdata} };

  RooCategory categs("categs","categs");
  categs.defineType("TRK");
  categs.defineType("SHW");

  // combined dataset, this links RooCategories with the histogram in hmap
  RooDataHist combData("combData","combData", fu.GetObs(), categs, hmap);

  // combined pdf, link pdf's with the categories
  RooSimultaneous simPdf("simPdf", "simPdf", categs);
  simPdf.addPdf(trkpdf, "TRK");
  simPdf.addPdf(shwpdf, "SHW");

  //==================================================================
  // perform the fit
  //==================================================================

  // fix some parameters
  fu.GetVar("Dm21")->setConstant(kTRUE);
  fu.GetVar("SinsqTh12")->setConstant(kTRUE);
  
  simPdf.fitTo( combData );
  
}
