#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "SummaryParser.h"
#include "RDFPIDReader.h"

#include <iostream>
using namespace std;

/**
 * This function takes the file NMH/data/pid_result_XXX.root as input and outputs a 
 * file in analysis format that uses SummaryEvent.
 *
 * See NMH/data_sorting/README.md for more info.
 *
 * \param  fin_name   Random decision forest PID output file from ECAP, e.g. ../data/pid_result_XXX.root
 * \param  fout_name  A file where the events from input, converted to analysis format, are stored.
 *
 */
void RDFPID_to_Summary(TString fin_name="", TString fout_name="") {

  if (fin_name == "" || fout_name == "") {
    cout << "ERROR! RDFPID_to_Summary() Specify input file name (../data/pid_result_XXX.root) and output file name (../data/ORCA_MC_summary_all_xxx.root). Exiting." << endl;
    return;
  }

  //----------------------------------------------------------------------------
  // load the reader class and initialize, set reading all branches
  //----------------------------------------------------------------------------

  gROOT->ProcessLine(".L RDFPIDReader.C+");
  TFile *fin = new TFile(fin_name, "READ");
  TTree *tin = (TTree*)fin->Get("PID");
  if (tin == NULL) {
    cout << "ERROR! RDFPID_to_Summary() cannot find tree PID in file " << fin->GetName() << endl;
  }
  RDFPIDReader PIDR(tin);
  PIDR.fChain->SetBranchStatus("*",1);

  //----------------------------------------------------------------------------
  // Init the output in analysis format, loop and map variables
  //----------------------------------------------------------------------------

  SummaryParser out(fout_name, kFALSE); //false means writing mode

  Long64_t nentries = PIDR.fChain->GetEntries();
  for (Int_t i = 0; i < nentries; i++) {

    if (i % 100000 == 0) cout << "Entry " << i << " out of " << nentries << endl;

    PIDR.GetEntry(i);
    
    out.GetEvt()->Set_MC_runID(PIDR.run_id);       
    out.GetEvt()->Set_MC_evtID(PIDR.mc_id);       
    out.GetEvt()->Set_MC_w2(PIDR.weight_w2);    
    out.GetEvt()->Set_MC_w1y(PIDR.weight_one_year);   
    out.GetEvt()->Set_MC_erange_start(PIDR.Erange_min);
    out.GetEvt()->Set_MC_is_CC(PIDR.is_cc); 
    out.GetEvt()->Set_MC_is_neutrino(PIDR.is_neutrino);
    out.GetEvt()->Set_MC_type(PIDR.type);
    out.GetEvt()->Set_MC_energy(PIDR.energy);
    out.GetEvt()->Set_MC_bjorkeny(PIDR.bjorkeny);
    out.GetEvt()->Set_MC_dir(PIDR.dir_x, PIDR.dir_y, PIDR.dir_z);
    out.GetEvt()->Set_MC_pos(PIDR.pos_x, PIDR.pos_y, PIDR.pos_z);

    //for quality levels of tracks and showers, see NMH/common_software/README.md

    out.GetEvt()->Set_track_energy(PIDR.gandalf_energy_corrected);
    out.GetEvt()->Set_track_bjorkeny(0.);                        //currently gandalf has no bjorkeny
    out.GetEvt()->Set_track_ql0(PIDR.gandalf_is_good);           
    out.GetEvt()->Set_track_ql1(PIDR.gandalf_loose_is_selected);
    out.GetEvt()->Set_track_ql2(PIDR.gandalf_is_selected);
    out.GetEvt()->Set_track_dir(PIDR.gandalf_dir_x, PIDR.gandalf_dir_y, PIDR.gandalf_dir_z);
    out.GetEvt()->Set_track_pos(PIDR.gandalf_pos_x, PIDR.gandalf_pos_y, PIDR.gandalf_pos_z);

    out.GetEvt()->Set_shower_energy(PIDR.dusj_energy_corrected);
    out.GetEvt()->Set_shower_bjorkeny(PIDR.dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY);
    out.GetEvt()->Set_shower_ql0(PIDR.dusj_is_good);     
    out.GetEvt()->Set_shower_ql1(PIDR.dusj_is_selected);
    out.GetEvt()->Set_shower_ql2(0.);                            //currently no ql2 for shower
    out.GetEvt()->Set_shower_dir(PIDR.dusj_dir_x, PIDR.dusj_dir_y, PIDR.dusj_dir_z);
    out.GetEvt()->Set_shower_pos(PIDR.dusj_pos_x, PIDR.dusj_pos_y, PIDR.dusj_pos_z);

    out.GetEvt()->Set_RDF_muon_score(PIDR.muon_score);
    out.GetEvt()->Set_RDF_track_score(PIDR.track_score);
    out.GetEvt()->Set_RDF_noise_score(PIDR.noise_score);

    out.GetTree()->Fill();

  }

  out.WriteAndClose();

  //----------------------------------------------------------------------------
  // clean-up
  //----------------------------------------------------------------------------
  fin->Close();
  delete fin;

}
