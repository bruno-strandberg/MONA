#include "GSGParser.h"
#include "SummaryParser.h"
#include "TSystem.h"
#include <iostream>
using namespace std;

void data_parsers(TString sum_file="../data/mc_end/data_atmnu/summary_muon-CC_3-100GeV_1.root",
		  TString gsg_file="../data/mc_start/data_atmnu/gSeaGen_muon-CC_3-100GeV_1.root") {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");

  SummaryParser sp(sum_file);
  GSGParser     gp(gsg_file);
  
  for (Int_t i = 0; i < sp.fChain->GetEntries(); i++) {
    sp.GetEntry(i);
    cout << "MC Energy, dir_z in summary: " << sp.MC_energy << "\t" << sp.MC_dir_z << endl;
    if (i > 20) break;
  }

  for (Int_t i = 0; i < gp.fChain->GetEntries(); i++) {
    gp.GetEntry(i);
    cout << "MC Energy, dir_z in gSeaGen: " << gp.Neutrino_E << "\t" << gp.Neutrino_D3 << endl;
    if (i > 20) break;
  }
  
}
