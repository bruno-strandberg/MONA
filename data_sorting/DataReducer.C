#define DataReducer_cxx
#include "DataReducer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//bstr: for cout
#include <iostream>
using namespace std;

void DataReducer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DataReducer.C
//      root> DataReducer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //bstr: call init on the output tree
   InitOutputTree();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (ientry % 100000 == 0) cout << "Entry: " << ientry << endl;

      //bstr: call fill on output tree
      tout->Fill();
   }

   //bstr: write the tree and clone fout
   tout->Write();
   fout->Close();
}

//****************************************************************************
//bstr: switch off all branches, except useful MC and reco information
//****************************************************************************

void DataReducer::SetBranches() {

   //bstr: set all branches off, except a few interesting ones
   fChain->SetBranchStatus("*", 0);

   //------------------------------------------------------
   //MC info
   //------------------------------------------------------
   fChain->SetBranchStatus("dir_x", 1);
   fChain->SetBranchStatus("dir_y", 1);
   fChain->SetBranchStatus("dir_z", 1);

   fChain->SetBranchStatus("pos_x", 1);
   fChain->SetBranchStatus("pos_y", 1);
   fChain->SetBranchStatus("pos_z", 1);

   fChain->SetBranchStatus("bjorkeny", 1);
   fChain->SetBranchStatus("energy", 1);

   fChain->SetBranchStatus("is_cc", 1);
   fChain->SetBranchStatus("is_neutrino", 1);
   fChain->SetBranchStatus("type", 1);
   
   fChain->SetBranchStatus("weight_w2", 1);
   fChain->SetBranchStatus("Erange_min",1);
   
   fChain->SetBranchStatus("run_id",1);
   fChain->SetBranchStatus("mc_id",1);

   //------------------------------------------------------
   //reco information
   //------------------------------------------------------
   fChain->SetBranchStatus("gandalf_dir_x", 1);
   fChain->SetBranchStatus("gandalf_dir_y", 1);
   fChain->SetBranchStatus("gandalf_dir_z", 1);

   fChain->SetBranchStatus("dusj_dir_x", 1);
   fChain->SetBranchStatus("dusj_dir_y", 1);
   fChain->SetBranchStatus("dusj_dir_z", 1);

   fChain->SetBranchStatus("recolns_dir_x", 1);
   fChain->SetBranchStatus("recolns_dir_y", 1);
   fChain->SetBranchStatus("recolns_dir_z", 1);

   fChain->SetBranchStatus("gandalf_pos_x", 1);
   fChain->SetBranchStatus("gandalf_pos_y", 1);
   fChain->SetBranchStatus("gandalf_pos_z", 1);

   fChain->SetBranchStatus("dusj_pos_x", 1);
   fChain->SetBranchStatus("dusj_pos_y", 1);
   fChain->SetBranchStatus("dusj_pos_z", 1);

   fChain->SetBranchStatus("recolns_pos_x", 1);
   fChain->SetBranchStatus("recolns_pos_y", 1);
   fChain->SetBranchStatus("recolns_pos_z", 1);

   fChain->SetBranchStatus("gandalf_energy_corrected", 1);
   fChain->SetBranchStatus("dusj_energy_corrected", 1);
   fChain->SetBranchStatus("recolns_energy_neutrino", 1);

   fChain->SetBranchStatus("dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY", 1);
   fChain->SetBranchStatus("recolns_bjorken_y", 1);

   fChain->SetBranchStatus("muon_probability", 1);
   fChain->SetBranchStatus("track_probability", 1);

}

//*********************************************************************
//bstr: a function to init the formatted output tree
//*********************************************************************

void DataReducer::InitOutputTree() {

  fout = new TFile("../data/ORCA_MC_summary_all.root","RECREATE");
  tout = new TTree("summary","ORCA MC summary tree");

  //MC truth info

  tout->Branch("MC_runID"        ,      &run_id, "MC_runID/D");
  tout->Branch("MC_evtID"        ,       &mc_id, "MC_evtID/D");
  tout->Branch("MC_w2"           ,   &weight_w2, "MC_w2/D");
  tout->Branch("MC_erange_start" ,  &Erange_min, "MC_erange_start/D");
  tout->Branch("MC_is_CC"        ,       &is_cc, "MC_is_CC/D");
  tout->Branch("MC_is_neutrino"  , &is_neutrino, "MC_is_neutrino/D");
  tout->Branch("MC_type"         ,        &type, "MC_type/D");

  tout->Branch("MC_dir_x" , &dir_x, "MC_dir_x/D");
  tout->Branch("MC_dir_y" , &dir_y, "MC_dir_y/D");
  tout->Branch("MC_dir_z" , &dir_z, "MC_dir_z/D");
  tout->Branch("MC_pos_x" , &pos_x, "MC_pos_x/D");
  tout->Branch("MC_pos_y" , &pos_y, "MC_pos_y/D");
  tout->Branch("MC_pos_z" , &pos_z, "MC_pos_z/D");

  tout->Branch("MC_energy"   ,   &energy, "MC_energy/D");
  tout->Branch("MC_bjorkeny" , &bjorkeny, "MC_bjorkeny/D");   
  
  //gandalf info

  tout->Branch("gandalf_dir_x", &gandalf_dir_x, "gandalf_dir_x/D");
  tout->Branch("gandalf_dir_y", &gandalf_dir_y, "gandalf_dir_y/D");
  tout->Branch("gandalf_dir_z", &gandalf_dir_z, "gandalf_dir_z/D");
  tout->Branch("gandalf_pos_x", &gandalf_pos_x, "gandalf_pos_x/D");
  tout->Branch("gandalf_pos_y", &gandalf_pos_y, "gandalf_pos_y/D");
  tout->Branch("gandalf_pos_z", &gandalf_pos_z, "gandalf_pos_z/D");
  
  tout->Branch("gandalf_energy_nu", &gandalf_energy_corrected, "gandalf_energy_nu/D");

  //shower info

  tout->Branch("shower_dir_x", &dusj_dir_x, "shower_dir_x/D");
  tout->Branch("shower_dir_y", &dusj_dir_y, "shower_dir_y/D");
  tout->Branch("shower_dir_z", &dusj_dir_z, "shower_dir_z/D");
  tout->Branch("shower_pos_x", &dusj_pos_x, "shower_pos_x/D");
  tout->Branch("shower_pos_y", &dusj_pos_y, "shower_pos_y/D");
  tout->Branch("shower_pos_z", &dusj_pos_z, "shower_pos_z/D");
  
  tout->Branch("shower_energy_nu", &dusj_energy_corrected, "shower_energy_nu/D");
  tout->Branch("shower_bjorkeny" , &dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY, "shower_bjorkeny/D");

  //recoLNS info

  tout->Branch("recolns_dir_x", &recolns_dir_x, "recolns_dir_x/D");
  tout->Branch("recolns_dir_y", &recolns_dir_y, "recolns_dir_y/D");
  tout->Branch("recolns_dir_z", &recolns_dir_z, "recolns_dir_z/D");
  tout->Branch("recolns_pos_x", &recolns_pos_x, "recolns_pos_x/D");
  tout->Branch("recolns_pos_y", &recolns_pos_y, "recolns_pos_y/D");
  tout->Branch("recolns_pos_z", &recolns_pos_z, "recolns_pos_z/D");

  tout->Branch("recolns_energy_nu", &recolns_energy_neutrino, "recolns_energy_nu/D");
  tout->Branch("recolns_bjorkeny" ,       &recolns_bjorken_y, "recolns_bjorkeny/D");

  //PID info
  tout->Branch("PID_muon_probability" ,  &muon_probability, "PID_muon_probability/D");
  tout->Branch("PID_track_probability", &track_probability, "PID_track_probability/D");

}

//*********************************************************************
//bstr: a function to create a reduced clone of the tree, very simple
//*********************************************************************

void DataReducer::CreateSimpleClone() {
  
  TFile *f_out = new TFile("../data/pid_result_reduced.root","RECREATE");
  TTree *tree_out = fChain->CloneTree();
  tree_out->Write();
  f_out->Close();

}
