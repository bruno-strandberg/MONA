#include <iostream>
#include "DetResponse.h"
#include "SummaryParser.h"


/** This example/test application tests the IO functionality of the `DetResponse` class.*/

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

  DetResponse dr(DetResponse::shower, "shower_response", 40, 1, 100, 40, -1, 1, 4, 0, 1);
  dr.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
  dr.AddCut( &SummaryEvent::Get_shower_dir_z    , std::less_equal<double>(),   0., true );
  dr.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(),  0.6, true );

  SummaryParser sp(fname);

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {        
    dr.Fill( sp.GetEvt(i) );
  }

  cout << "NOTICE detresponseIO() Finished filling response" << endl;

  //--------------------------------------------------------------------------
  // write to file; clone the existing response; read-in the response written out
  //--------------------------------------------------------------------------

  dr.WriteToFile(respname);
  cout << "NOTICE detresponseIO() Wrote to file" << endl;

  DetResponse dr_clone("clone",dr);
  DetResponse dr_readin(DetResponse::shower, "shower_readin");
  dr_readin.ReadFromFile(respname);

  //--------------------------------------------------------------------------
  // get the true bins from each of the DetResponse instance, should give identical results
  //--------------------------------------------------------------------------

  auto true_bins1 = dr.GetBinWeights(10, -0.8, 0.5);
  auto true_bins2 = dr_clone.GetBinWeights(10, -0.8, 0.5);
  auto true_bins3 = dr_readin.GetBinWeights(10, -0.8, 0.5);

  Bool_t testOK = kTRUE;

  if ( ( true_bins1.size() != true_bins2.size() ) || ( true_bins2.size() != true_bins3.size() ) ) {
    
    cout << "NOTICE detresponseIO() test failing, detector responses that should be the same size differ, the sizes are: "  << true_bins1.size() << "\t" << true_bins2.size() << "\t" << true_bins3.size() << ", exiting." << endl;
    testOK = kFALSE;
  }
  else {

    for (Int_t i = 0; i < (Int_t)true_bins1.size(); i++) {

      testOK = testOK && (true_bins1[i].fW == true_bins2[i].fW && true_bins1[i].fW == true_bins3[i].fW);

      if (!testOK) {
	cout << "NOTICE detresponseIO() test failing, different bin weights for the following bins:" << endl;
	cout << true_bins1[i] << endl;
	cout << true_bins2[i] << endl;
	cout << true_bins3[i] << endl;
      }

    }
  }

  if ( system("rm " + fname + " " + respname) != 0) {
    cout << "WARNING! detresponseIO() file removal returned an error code" << endl;
  }

  if (testOK) return 0;
  else return 1;

}
