
// root headers
#include "TFile.h"
#include "TTree.h"

// NMH headers
#include "SummaryParser.h"
#include "ECAP181013.h"
#include "FileHeader.h"

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

// cpp headers
#include <stdexcept>
#include <iostream>


using namespace std;

/**
   This program takes the ECAP PID summary file as input and converts it to `SummaryEvent` format for usage with MONA software.
 
   See apps/data_sorting/README.md for more info.
   
*/
int main(const int argc, const char **argv) {

  string fin_name;
  string fout_dir;
  string tag;

  try {

    JParser<> zap("This program takes the ECAP PID summary file as input, reading it with ECAP181013.h/C class, and converts it to `SummaryEvent` format for usage with NMH software.");

    zap['f'] = make_field(fin_name , "PID file from ECAP with the PID TTree that was used to create ECAP181013.h/C class") = "irodsdata/pid_output_atm_neutrino_atm_muon_pure_noise_anyRecoUpgoing_181013.root";
    zap['d'] = make_field(fout_dir , "Output directory where the data file in SummarEvent format is written, e.g. ../../data/") = "";
    zap['t'] = make_field(tag      , "Identifier tag used to create the SummaryEvent file, e.g. format ORCA115_23x9m_ECAP180401. Choose this wisely.") = "";

    if ( zap.read(argc, argv) != 0 ) return 1;

  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  if (fin_name == "" || fout_dir  == "" || tag == "") {
    throw std::invalid_argument("ERROR! ECAP181013_to_MONA() all command line arguments need to be specified!");
  }

  //----------------------------------------------------------------------------
  // load the reader class and initialize, set reading all branches
  //----------------------------------------------------------------------------

  TFile *fin = new TFile((TString)fin_name, "READ");
  TTree *tin = (TTree*)fin->Get("PID");
  if (tin == NULL) {
    throw std::invalid_argument("ERROR! ECAP181013_to_MONA() cannot find tree PID in file " + fin_name);
  }
  ECAP181013 PIDR(tin);
  PIDR.fChain->SetBranchStatus("*",1);

  //----------------------------------------------------------------------------
  // Init the output in analysis format, loop and map variables
  //----------------------------------------------------------------------------
  SummaryEvent *evt = new SummaryEvent;
  string sevtv = std::to_string( ( (TClass*)evt->IsA() )->GetClassVersion() ); // summary event version in the library
  delete evt;

  string fout_name = fout_dir + "ORCA_MCsummary_SEv" + sevtv + "_" + tag + ".root";

  cout << "NOTICE ECAP181013_to_MONA() creating file " << fout_name << endl;

  SummaryParser out(fout_name, kFALSE); //false means writing mode

  // in this ORCA7 production the extremely bad idea of grouping together runs at JTE level was followed
  // This means run numbers go in steps of 20, i.e, run 1, run 20, run 40. Due to this, some info was lost, i.e.
  // it is not possible to track down from which exact run the event came from. For the w1y calculation in
  // EvtResponse, this messes up the run count. Hence I need to factor in the grouping in the denominator
  // To spice things up a bit further, different grouping was used for different runs
  Double_t JTE_grouping_nu    = 20;
  Double_t JTE_grouping_atmmu = 50;
  Double_t JTE_grouping_noise = 1;

  Long64_t nentries = PIDR.fChain->GetEntries();
  for (Int_t i = 0; i < nentries; i++) {

    if (i % 100000 == 0) cout << "Entry " << i << " out of " << nentries << endl;

    PIDR.GetEntry(i);

    //********************************************************************************
    // here the user will need to select what is written to variables
    //********************************************************************************

    out.GetEvt()->Set_MC_runID(PIDR.run_id);       
    out.GetEvt()->Set_MC_evtID(PIDR.mc_id);       
    out.GetEvt()->Set_MC_w2(PIDR.weight_w2);    
    out.GetEvt()->Set_MC_w1y(PIDR.weight_one_year);   

    if ( PIDR.is_neutrino ) { 
      out.GetEvt()->Set_MC_w2denom( PIDR.n_events_gen * JTE_grouping_nu ); // N_tot is 20 times as large due to grouping
    }
    else { 
      if ( PIDR.type == 0 ) { out.GetEvt()->Set_MC_w2denom( PIDR.livetime_sec * JTE_grouping_noise ); }
      else                  { out.GetEvt()->Set_MC_w2denom( PIDR.livetime_sec * JTE_grouping_atmmu ); }
    }

    out.GetEvt()->Set_MC_erange_start(PIDR.Erange_min);
    out.GetEvt()->Set_MC_erange_stop(PIDR.Erange_max);
    out.GetEvt()->Set_MC_is_CC(PIDR.is_cc); 
    out.GetEvt()->Set_MC_is_neutrino(PIDR.is_neutrino);
    out.GetEvt()->Set_MC_type(PIDR.type);
    out.GetEvt()->Set_MC_ichan(PIDR.interaction_channel);
    out.GetEvt()->Set_MC_energy(PIDR.energy);
    out.GetEvt()->Set_MC_bjorkeny(PIDR.bjorkeny);
    out.GetEvt()->Set_MC_dir(PIDR.dir_x, PIDR.dir_y, PIDR.dir_z);
    out.GetEvt()->Set_MC_pos(PIDR.pos_x, PIDR.pos_y, PIDR.pos_z);

    out.GetEvt()->Set_track_energy(PIDR.gandalf_energy);
    out.GetEvt()->Set_track_bjorkeny(0.);                        //currently gandalf has no bjorkeny
    out.GetEvt()->Set_track_ql0(PIDR.gandalf_is_good);           
    out.GetEvt()->Set_track_ql1(0.);                             //no pre-defined quality levels for 7-line prod
    out.GetEvt()->Set_track_ql2(0.);
    out.GetEvt()->Set_track_ql3(0.);
    out.GetEvt()->Set_track_dir(PIDR.gandalf_dir_x, PIDR.gandalf_dir_y, PIDR.gandalf_dir_z);
    out.GetEvt()->Set_track_pos(PIDR.gandalf_pos_x, PIDR.gandalf_pos_y, PIDR.gandalf_pos_z);

    out.GetEvt()->Set_shower_energy(PIDR.dusj_energy_corrected);
    out.GetEvt()->Set_shower_bjorkeny(PIDR.dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY);
    out.GetEvt()->Set_shower_ql0(PIDR.dusj_is_good);     
    out.GetEvt()->Set_shower_ql1(0.);                            //no pre-defined quality levels for 7-line prod
    out.GetEvt()->Set_shower_ql2(0.);                            
    out.GetEvt()->Set_shower_ql3(0.);                            
    out.GetEvt()->Set_shower_dir(PIDR.dusj_dir_x, PIDR.dusj_dir_y, PIDR.dusj_dir_z);
    out.GetEvt()->Set_shower_pos(PIDR.dusj_pos_x, PIDR.dusj_pos_y, PIDR.dusj_pos_z);

    out.GetEvt()->Set_RDF_muon_score(PIDR.muon_score_retrain);   //better suppression by retrain score
    out.GetEvt()->Set_RDF_track_score(PIDR.track_score);
    out.GetEvt()->Set_RDF_noise_score(PIDR.noise_score);
    
    out.GetTree()->Fill();

  }

  // add the tag to the header
  out.GetHeader()->AddParameter("datatag", (TString)tag);
  out.WriteAndClose();

  //----------------------------------------------------------------------------
  // clean-up
  //----------------------------------------------------------------------------
  fin->Close();
  delete fin;

  return 0;

}
