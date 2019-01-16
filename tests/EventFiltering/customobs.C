#include "EventSelection.h"
#include "SummaryParser.h"
#include "TVector3.h"
#include "TMath.h"

/** 
    Namespace for functions of `customobs.C` test/example application. It holds the functions that need to be defined in `EventFilter` or inheriting classes are used with `customreco` reconstruction type. These functions need to adhere to a certain format - all of them take a SummaryEvent as an argument; energy and b-y function return a `Double_t`, direction and position function return a `TVector3`.
*/
namespace CUSTOMOBS {

  /** example function to return reconstructed energy.
      \param evt  Pointer to a summary event instance
      \return     reconstructed energy
   */
  Double_t CustomEnergy(SummaryEvent *evt) {

    // for example, use the energy that is closer to the MC truth
    if ( TMath::Abs( evt->Get_MC_energy() - evt->Get_track_energy() ) < 
	 TMath::Abs( evt->Get_MC_energy() - evt->Get_shower_energy() ) ) {
      return evt->Get_track_energy();
    }
    else {
      return evt->Get_shower_energy();
    }

  }

  // for direction, position and bjorken-y I just use the track standards
  TVector3 CustomDir(SummaryEvent* evt) { return evt->Get_track_dir();      }
  TVector3 CustomPos(SummaryEvent* evt) { return evt->Get_track_pos();      }
  Double_t CustomBY (SummaryEvent* evt) { return evt->Get_track_bjorkeny(); }

};

//***********************************************************************************************

/** 
    This example&test demonstrates how to use custom reconstruction variables in `EventFilter`, `EventSelection` and `DetResponse`. This is sometimes useful when one wishes to use track direction reco but, for example, shower energy reco, with some additional conditions.

    \return 0 if test worked, 1 otherwise.
*/
int main(const int argc, const char** argv) {

  //-----------------------------------------------------------
  // create pseudo-data to test the filter on
  //-----------------------------------------------------------
  TString fname = "customobs_out.root";
  SummaryParser wr(fname, kFALSE);
  
  for (Int_t i = 0; i < 1000; i++) {
    wr.GetEvt()->FillPseudoData();
    wr.GetTree()->Fill();
  }

  wr.WriteAndClose();

  //-----------------------------------------------------------
  // init a selection that uses custom reco for internal observables; note the example of SetObsFuncPtrs! 
  // Create a clone of the 2D histogram of the selection for manual filling for comparison
  //-----------------------------------------------------------
  EventSelection sel1(EventSelection::customreco);
  sel1.SetObsFuncPtrs( &CUSTOMOBS::CustomEnergy, &CUSTOMOBS::CustomDir, &CUSTOMOBS::CustomPos, &CUSTOMOBS::CustomBY );

  sel1.AddCut(&SummaryEvent::Get_track_energy, std::greater<double>()   ,  5., true);
  sel1.AddCut(&SummaryEvent::Get_track_dir_z , std::less<double>()      ,  0., true);
  sel1.AddCut(&SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.7, true);

  TH2D* hmanual = (TH2D*)sel1.Get_h_E_costh()->Clone("hmanual");

  //-----------------------------------------------------------
  // open the pseudodata file for reading and fill the selection. Also fill a histogram with brute force
  // that mimics the expected behavior of the selection for comparison.
  //-----------------------------------------------------------
  SummaryParser sp(fname);
  SummaryEvent *evt;

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

    evt = sp.GetEvt(i);

    sel1.Fill( evt );

    if (evt->Get_track_energy() > 5 && evt->Get_track_dir_z() < 0 && evt->Get_RDF_track_score() > 0.7) {

      if ( TMath::Abs( evt->Get_MC_energy() - evt->Get_track_energy() ) < 
	   TMath::Abs( evt->Get_MC_energy() - evt->Get_shower_energy() ) ) {
	hmanual->Fill( evt->Get_track_energy(), -evt->Get_track_dir_z() );
      }
      else {
	hmanual->Fill( evt->Get_shower_energy(), -evt->Get_track_dir_z() );
      } // energy selection

    } // event selection

  }

  //-----------------------------------------------------------
  // perform bin-by-bin comparison of the two methods
  //-----------------------------------------------------------
  
  Bool_t testOK = kTRUE;

  for (Int_t xbin = 1; xbin <= hmanual->GetXaxis()->GetNbins(); xbin++) {
    for (Int_t ybin = 1; ybin <= hmanual->GetYaxis()->GetNbins(); ybin++) {
      
      Double_t bc1 = sel1.Get_h_E_costh()->GetBinContent(xbin, ybin);
      Double_t bc2 = hmanual->GetBinContent(xbin, ybin);

      testOK = testOK && (bc1 == bc2);
      
      if (!testOK) {
	cout << "NOTICE customobs() test is failing, different bin contents: " << bc1 << "\t" << bc2 << endl;
      }

    }
  }

  cout << "Total entries: " << sel1.Get_h_E_costh()->GetEntries() << "\t" << hmanual->GetEntries() << endl;

  if (testOK) cout << "Test passed" << endl;
  else cout << "Test failed" << endl;

  if ( system("rm " + fname) != 0 ) {
    cout << "WARNING! customobs() file removal command returned non-zero" << endl;
  }
  
  if (testOK) return 0;
  else        return 1;
}
