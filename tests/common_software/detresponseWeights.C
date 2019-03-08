#include "DetResponse.h"
#include "SummaryEvent.h"

#include "TH3D.h"
#include "TFile.h"

#include <iostream>
using namespace std;

/** This application tests the filling and calculation of the weights in `DetResponse` class.
    \return 0 if test worked, 1 otherwise.
*/

int main(const int argc, const char** argv) {

  // response with 2 bins along energy, cos-theta and bjorken-y
  DetResponse R(DetResponse::mc_truth, "R", 2, 1, 100, 2, -1, 1, 2, 0, 1);

  // init a summary event
  SummaryEvent evt;
  evt.FillPseudoData();

  // set to some un-known neutrino flavor and catch the exception
  //-----------------------------------------------------------------
  evt.Set_MC_type(15);
  
  bool error_caught = 0;
  try {
    R.Fill(&evt);
  }
  catch (const std::invalid_argument& ia) {
    error_caught = true;
  }

  if (!error_caught || R.GetHist3D()->GetEntries() != 0) {
    cout << "NOTICE detresponseWeights failed at recognition of an unsupported particle type" << endl;
    return 1;
  }

  // check that an event with energy, cos-theta or energy out of range is not counted
  //-----------------------------------------------------------------
  evt.FillPseudoData();

  evt.Set_MC_energy( 120 );
  R.Fill(&evt);
  evt.Set_MC_energy( 0.997 );
  R.Fill(&evt);

  if ( R.GetHist3D()->GetEntries() != 0 ) {
    cout << "NOTICE detresponseWeights failed at not rejecting an event outside the energy range" << endl;
    return 1;
  }

  evt.FillPseudoData();
  evt.Set_MC_dir(0,0,-1.01);
  R.Fill(&evt);
  evt.Set_MC_dir(0,0,1.01);
  R.Fill(&evt);
  
  if ( R.GetHist3D()->GetEntries() != 0 ) {
    cout << "NOTICE detresponseWeights failed at not rejecting an event outside the cos-theta range" << endl;
    return 1;
  }

  evt.FillPseudoData();
  evt.Set_MC_bjorkeny(1.01);
  R.Fill(&evt);
  evt.Set_MC_bjorkeny(-0.01);
  R.Fill(&evt);
  
  if ( R.GetHist3D()->GetEntries() != 0 ) {
    cout << "NOTICE detresponseWeights failed at not rejecting an event outside the bjorken-y range" << endl;
    return 1;
  }

  // check that the weight is calculated correctly for each nu flavor
  //-----------------------------------------------------------------
  vector<UInt_t> flavs = {12, 14, 16};

  for (auto f: flavs) {

    Double_t track_E_cut = 70;
    DetResponse RT(DetResponse::track, "RT", 2, 1, 100, 1, -1, 1, 1, 0, 1);
    RT.AddCut( &SummaryEvent::Get_track_energy, std::less<double>(), track_E_cut, true );
    
    // create histograms to count selected and generated events
    enum bins {EBIN1 = 0, EBIN2 = 1};
    TH3D* sim;
    TH3D* sel[2];
    
    sim = (TH3D*)RT.GetHist3D()->Clone("sim");

    /* This is somewhat tricky to implement with histograms, which is why DetResponse is a somewhat complicated
       class. For one sim histogram, I need two sel histograms. One counts the contributions from
       the first true energy bin to the two reco bins, the second counts the contributions from the second true
       energy bin to the two reco bins.
    */
    
    sel[EBIN1] = (TH3D*)RT.GetHist3D()->Clone("sel_E1"); // counts selected events from true ebin 1
    sel[EBIN2] = (TH3D*)RT.GetHist3D()->Clone("sel_E2"); // counts selected events from true ebin 2
    
    vector<TH3D*> hists = { sim, sel[EBIN1], sel[EBIN2] };
    for (auto h: hists) { h->Reset(); h->SetDirectory(0); }
    
    // generate CC pseudo-data for a specific nu type
    for (Int_t i = 1; i < 1000; i++) {

      evt.FillPseudoData();
      evt.Set_MC_type( f ); // set manually to the considered flavor
      evt.Set_MC_is_CC(1.); // set manually to CC
      
      // fill the response
      RT.Fill(&evt);

      // get the true bin to decide in which `sel` histogram the event should be counted
      Int_t true_ebin = RT.GetHist3D()->GetXaxis()->FindBin( evt.Get_MC_energy() ) - 1;

      // count the total number of simulated events for both true bins with histogram `sim`
      sim           ->Fill( evt.Get_MC_energy()   , -evt.Get_MC_dir_z()   , evt.Get_MC_bjorkeny() );

      // for both true bins, count how the events divide between two reco bins
      if ( evt.Get_track_energy() < track_E_cut ) {
	sel[true_ebin]->Fill( evt.Get_track_energy(), -evt.Get_track_dir_z(), evt.Get_track_bjorkeny() );
      }
		      
    }

    // get the true bins that contribute to reco bins 1 and 2 from the response
    auto true_bins_E1 = RT.GetBinWeights( RT.GetHist3D()->GetXaxis()->GetBinCenter( 1 ), 0., 0.5);
    auto true_bins_E2 = RT.GetBinWeights( RT.GetHist3D()->GetXaxis()->GetBinCenter( 2 ), 0., 0.5);

    // expect both to have contributions from two true energy bins
    if ( true_bins_E1.size() != 2 || true_bins_E2.size() != 2 ) {
      cout << "NOTICE detresponseWeights failed at weights check, wrong number of contributing true bins" << endl;
      return 1;
    }

    // loop over the two true bins that contribute to the first reco bin
    for (auto tb: true_bins_E1) {

      // get the number of events from true energy bin that ended up in reco bin 1
      Double_t numerator   = sel[tb.fE_true_bin-1]->GetBinContent( 1, 1, 1 );

      // get the total number of simulated events in true energy bin
      Double_t denominator = sim->GetBinContent(tb.fE_true_bin, 1, 1);

      if ( numerator/denominator != tb.fW ) {
	cout << "NOTICE detresponseWeights failed at weights check, wrong weight calculation for bin "
	     << tb << endl;
	cout << numerator/denominator << "\t" << tb.fW << endl;
	return 1;
      }
      
    }

    // loop over the two true bins that contribute to the second reco bin
    for (auto tb: true_bins_E2) {

      // get the number of events from true energy bin that ended up in reco bin 2
      Double_t numerator   = sel[tb.fE_true_bin-1]->GetBinContent( 2, 1, 1 );

      // get the total number of simulated events in true energy bin 
      Double_t denominator = sim->GetBinContent(tb.fE_true_bin, 1, 1);

      if ( numerator/denominator != tb.fW ) {
	cout << "NOTICE detresponseWeights failed at weights check, wrong weight calculation for bin "
	     << tb << endl;
	cout << numerator/denominator << "\t" << tb.fW << endl;
	return 1;
      }
      
    }
    
    // clean the histograms
    for (auto h: hists) delete h;

  }
  
  cout << "NOTICE detresponseWeights test passed" << endl;
  return 0;
  
}
