#include <iostream>
#include "EvtResponse.h"
#include "SummaryParser.h"
#include "TRandom.h"

/** This example/test application tests the IO functionality of the `EvtResponse` class.*/

int main(const int argc, const char** argv) {

  //--------------------------------------------------------------------------
  // create some pseudodata
  //--------------------------------------------------------------------------

  TString fname    = "detresponseIO_data.root";
  TString respname = "detresponseIO_resp.root";

  SummaryParser wr(fname, kFALSE);
  for (Int_t i = 0; i < 100000; i++) {
    wr.GetEvt()->FillPseudoData();
    wr.GetTree()->Fill();
  }
  wr.WriteAndClose();

  //--------------------------------------------------------------------------
  // initialize det response, add some filters and fill with pseudodata
  //--------------------------------------------------------------------------

  EvtResponse dr(EvtResponse::shower, "shower_response", 40, 1, 100, 40, -1, 1, 4, 0, 1);
  dr.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
  dr.AddCut( &SummaryEvent::Get_shower_dir_z    , std::less_equal<double>(),   0., true );
  dr.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(),  0.6, true );

  SummaryParser sp(fname);

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {        
    dr.Fill( sp.GetEvt(i) );
  }

  cout << "NOTICE evtresponseIO() Finished filling response" << endl;

  //--------------------------------------------------------------------------
  // write to file; clone the existing response; read-in the response written out
  //--------------------------------------------------------------------------

  dr.WriteToFile(respname);
  cout << "NOTICE evtresponseIO() Wrote to file" << endl;

  EvtResponse dr_readin(EvtResponse::shower, "shower_readin");
  dr_readin.ReadFromFile(respname);

  //--------------------------------------------------------------------------
  // get the true bins from each of the EvtResponse instance, should give identical results
  //--------------------------------------------------------------------------

  TRandom3 rand(0);

  Double_t E  = rand.Uniform(5,50);
  Double_t ct = rand.Uniform(-1,0);
  
  auto true_evts1 = dr.GetBinEvts( E, ct, 0.5);
  auto true_evts2 = dr_readin.GetBinEvts( E, ct, 0.5);

  Bool_t testOK = kTRUE;

  if ( true_evts1.size() != true_evts2.size() ) {
    
    cout << "NOTICE evtresponseIO() test failing, detector responses that should be the same size differ, the sizes are: "  << true_evts1.size() << "\t" << true_evts2.size() << ", exiting." << endl;
    testOK = kFALSE;
    
  }
  else {

    for (Int_t i = 0; i < (Int_t)true_evts1.size(); i++) {

      testOK = testOK && ( true_evts1[i].GetTrueE() == true_evts2[i].GetTrueE() );

      if (!testOK) {
	cout << "NOTICE evtresponseIO() test failing, different energies for the following events:" << endl;
	cout << true_evts1[i] << endl;
	cout << true_evts2[i] << endl;
      }

    }
  }

  if ( system("rm " + fname + " " + respname) != 0) {
    cout << "WARNING! evtresponseIO() file removal returned an error code" << endl;
  }

  if (testOK) return 0;
  else return 1;

}
