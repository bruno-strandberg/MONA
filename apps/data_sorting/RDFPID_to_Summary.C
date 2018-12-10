
// root headers
#include "TFile.h"
#include "TTree.h"

// NMH headers
#include "SummaryParser.h"
#include "RDFPIDReader.h"
#include "FileHeader.h"

// jpp headers
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

// cpp headers
#include <stdexcept>
#include <iostream>


using namespace std;

/**
   This program takes the ECAP PID summary file as input and converts it to `SummaryEvent` format for usage with NMH software.
 
   See apps/data_sorting/README.md for more info.
   
*/
int main(const int argc, const char **argv) {

  string fin_name;
  string fout_dir;
  string tag;

  try {

    JParser<> zap("This program takes the ECAP PID summary file as input and converts it to `SummaryEvent` format for usage with NMH software.");

    zap['f'] = make_field(fin_name , "PID file from ECAP, e.g. ../../data/pid_result_XXX.root") = "";
    zap['d'] = make_field(fout_dir , "Output directory where the data file in SummarEvent format is written, e.g. ../../data/") = "";
    zap['t'] = make_field(tag      , "Identifier tag used to create the SummaryEvent file, e.g. ORCA115_23x9m_ECAP0418. Choose this wisely.") = "";

    zap.read(argc, argv);

  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  if (fin_name == "" || fout_dir  == "" || tag == "") {
    throw std::invalid_argument("ERROR! RDFPID_to_Summary() all command line arguments need to be specified!");
  }

  //----------------------------------------------------------------------------
  // load the reader class and initialize, set reading all branches
  //----------------------------------------------------------------------------

  TFile *fin = new TFile((TString)fin_name, "READ");
  TTree *tin = (TTree*)fin->Get("PID");
  if (tin == NULL) {
    throw std::invalid_argument("ERROR! RDFPID_to_Summary() cannot find tree PID in file " + fin_name);
  }
  RDFPIDReader PIDR(tin);
  PIDR.fChain->SetBranchStatus("*",1);

  //----------------------------------------------------------------------------
  // Init the output in analysis format, loop and map variables
  //----------------------------------------------------------------------------

  string fout_name = fout_dir + "ORCA_MC_summary_" + tag + ".root";
  cout << "NOTICE RDFPID_to_Summary() creating file " << fout_name << endl;

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

  // add the tag to the header
  FileHeader head("RDFPID_to_Summary");
  head.AddParameter("datatag", (TString)tag);
  head.WriteHeader(out.GetFile());

  out.WriteAndClose();

  //----------------------------------------------------------------------------
  // clean-up
  //----------------------------------------------------------------------------
  fin->Close();
  delete fin;

  return 0;

}
