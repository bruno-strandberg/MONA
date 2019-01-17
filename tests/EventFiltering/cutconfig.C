#include "EventSelection.h"
#include "SummaryParser.h"

/** This example&test demonstrates how to initialise an `EventSelection` and add cuts to it. The behavior is inherited from `EventFilter` and is available both in `EventSelection` and `DetResponse`.
    \return 0 if test worked, 1 otherwise.
*/

int main(const int argc, const char** argv) {

  //-----------------------------------------------------------
  // create pseudo-data to test the filter on
  //-----------------------------------------------------------
  TString fname = "settingcuts_out.root";
  SummaryParser wr(fname, kFALSE);
  
  for (Int_t i = 0; i < 1000; i++) {
    wr.GetEvt()->FillPseudoData();
    wr.GetTree()->Fill();
  }

  wr.WriteAndClose();

  //-----------------------------------------------------------
  // init a selection that uses track reco for internal observables energy, dir, pos and bjorken-y
  //-----------------------------------------------------------
  EventSelection sel1(EventSelection::track);

  // track energy above 3 GeV AND upgoing events (cosz < 0); all cuts added with 
  // last argument 'true' are added as AND cuts, e.g. I have configured (energy > 3 && dir_z < 0)
  sel1.AddCut(&SummaryEvent::Get_track_energy, std::greater<double>() ,  3., true);
  sel1.AddCut(&SummaryEvent::Get_track_dir_z , std::less<double>()    ,  0., true);

  // only look at events of mc type 14 or -14 (muons); all cuts added with 'false' as last argument
  // are used to make an OR statement, which is multiplied with the AND statement. In this case, the
  // result is (energy > 3 && dir_z < 0) && (type == 14 || type == -14)
  sel1.AddCut(&SummaryEvent::Get_MC_type     , std::equal_to<double>(),  14., false);
  sel1.AddCut(&SummaryEvent::Get_MC_type     , std::equal_to<double>(), -14., false);

  //-----------------------------------------------------------
  // make another filter that copies over sel1 and adds an additional cut on track score
  //-----------------------------------------------------------
  EventSelection sel2(sel1);
  sel2.AddCut(&SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.8, true);

  //-----------------------------------------------------------
  // open the pseudodata file for reading and fill the selections, compare results with ROOT
  //-----------------------------------------------------------
  SummaryParser sp(fname);
  
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    sel1.Fill( sp.GetEvt(i) );
    sel2.Fill( sp.GetEvt(i) );
  }

  Double_t entries1 = sp.GetTree()->GetEntries("fTrack_energy>3&&fTrack_dir_z<0&&TMath::Abs(fMC_type)==14");
  Double_t entries2 = sp.GetTree()->GetEntries("fTrack_energy>3&&fTrack_dir_z<0&&TMath::Abs(fMC_type)==14&&fRDF_track_score>0.8");

  Bool_t testOK = ( entries1 == sel1.Get_h_E_costh()->GetEntries() && entries2 == sel2.Get_h_E_costh()->GetEntries() );

  cout << "First selection entries : " << sel1.Get_h_E_costh()->GetEntries() << "\t" << entries1 << endl;
  cout << "Second selection entries: " << sel2.Get_h_E_costh()->GetEntries() << "\t" << entries2 << endl;

  if (testOK) cout << "Test passed" << endl;
  else cout << "Test failed" << endl;

  if ( system("rm " + fname) != 0 ) {
    cout << "WARNING! cutconfig() file removal command returned non-zero" << endl;
  }
  
  if (testOK) return 0;
  else        return 1;
}
