#include "HelperFunctions.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

/* Script to generate the Q distribution used in sampling when randomizing the Q value of
 * events. This distribution is used to quantify how much of the senstivity is due to lower
 * statistics when using more bins in PID variable Q, and how much is due to the signal in 
 * this variable Q.
 */

void GenerateQDistribution() {

  const int N_PID_CLASSES = 3;

  std::map<Int_t, Double_t> pid_map = SetPIDCase(N_PID_CLASSES);

  TString filefolder = DetectorResponseFolder(N_PID_CLASSES);

  //-----------------------------------------------------
  // fill the histogram with the q values from the events
  //-----------------------------------------------------

  Bool_t files_exist = kTRUE; 
  TString q_dist_file  = "pid_quality_distribution.root";
  Bool_t q_file_exists  = NMHUtils::FileExists(filefolder + q_dist_file);
  files_exist = (files_exist and q_file_exists);

  if (not files_exist) {
    TFile fout(filefolder + q_dist_file, "RECREATE");
    TH1D* h_q_dist = new TH1D("h_quality", "h_quality", 100, 0, 1);

    auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
    SummaryParser sp(summary_file);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
      SummaryEvent *evt = sp.GetEvt(i);

      h_q_dist->Fill(evt->Get_RDF_track_score());

    }

    h_q_dist->Write();
    cout << "NOTICE: Finished making quality distribution" << endl;
    fout.Close();
  }
  else {
    cout << "NOTICE: Distribution is already present at " << filefolder + q_dist_file << endl;
  }
}
