#include <iostream>
#include "NMHUtils.h"

void detresponse() {

  SummaryParser sp( (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root");

  //--------------------------------------------------------------------------
  // initialize det response, add some filters and fill
  //--------------------------------------------------------------------------

  DetResponse dr(DetResponse::shower, "shower_response", 40, 1, 100, 40, -1, 1, 4, 0, 1);
  dr.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
  dr.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
  dr.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(),  0.6, true );
  dr.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  dr.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

    if (i % 100000 == 0) cout << "Entry: " << i << endl;
    
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    
    dr.Fill(evt);

  }

  cout << "NOTICE: Finished filling response" << endl;

  //--------------------------------------------------------------------------
  // write to file; clone the existing response; read-in the response written out
  //--------------------------------------------------------------------------

  dr.WriteToFile("./detresponse-out.root");
  cout << "NOTICE: Wrote to file" << endl;

  DetResponse dr_clone(dr);
  DetResponse dr_readin(DetResponse::shower, "shower_readin");
  dr_readin.ReadFromFile("detresponse-out.root");

  //--------------------------------------------------------------------------
  // get the true bins from each of the DetResponse instance, should give identical results
  //--------------------------------------------------------------------------

  auto true_bins1 = dr.GetBinWeights(10, -0.8, 0.5);
  auto true_bins2 = dr_clone.GetBinWeights(10, -0.8, 0.5);
  auto true_bins3 = dr_readin.GetBinWeights(10, -0.8, 0.5);

  if ( ( true_bins1.size() != true_bins2.size() ) || ( true_bins2.size() != true_bins3.size() ) ) {
    cout << "ERROR! Detector responses that should be the same size differ, the sizes are: "
	 << true_bins1.size() << "\t" << true_bins2.size() << "\t" << true_bins3.size() << ", exiting." << endl;
    return;
  }
  else {

    cout << "###############################################################" << endl;
    cout << "The true bin data in different responses" << endl;
    cout << "###############################################################" << endl;
    for (Int_t i = 0; i < true_bins1.size(); i++) {
      cout << true_bins1[i] << endl;
      cout << true_bins2[i] << endl;
      cout << true_bins3[i] << endl;
      if (i > 20) break;
      cout << "***************************************************************" << endl;

    }
  }

  Int_t sysret = system("rm detresponse-out.root");

}
