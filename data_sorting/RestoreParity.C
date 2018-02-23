#include "TROOT.h"
#include "TSystem.h"
#include "../common_software/SummaryParser.h"


//****************************************************************

//functions

void  split_to_files_atmnu(Bool_t overwrite=false);
void  split_to_files_mupage(Bool_t overwrite=false, Bool_t separate=false);
Int_t get_max_run_nr(TString cut_string);

//global variables accessible in all functions inside this macro

SummaryParser *sp;

//****************************************************************

//entry point

void RestoreParity() {

  //initialise the data parser
  gSystem->Load("../common_software/libcommonsoft.so");
  sp = new SummaryParser("../data/ORCA_MC_summary_all.root");

  if (sp->fChain == NULL) {
    cout << "ERROR! RestoreParity() Issue with SummaryParser() initialisation." << endl;
  }

  //call functions
  split_to_files_atmnu();
  split_to_files_mupage();
  
}

//****************************************************************

//Function to split the pid output to summary files per mc input
//this handles atmospheric neutrinos

void split_to_files_atmnu(Bool_t overwrite) {

  //maps to build the filenames
  map < Double_t, TString > pdg_to_name = { {12., "elec" }, 
					    {14., "muon" },
					    {16., "tau"  } };

  map < Double_t, TString > iscc_to_str = { { 0., "NC"}, 
					    { 1., "CC"} };

  map < Double_t, TString>  emin_to_str = { { 1., "1-5"  },
					    { 3., "3-100"} };
  
  map < TString, Int_t >    cut_maxrun  = {};

  Int_t max_run_nr = 2000;

  //Loop over run numbers
  for (Int_t run_nr = 1; run_nr <= max_run_nr; run_nr++) {
    
    for (auto const& type : pdg_to_name) {

      for (auto const& interaction : iscc_to_str) {

	for (auto const& emin : emin_to_str) {

	  //build the cut string and the file name

	  TString cut_string = "";
	  cut_string  += "TMath::Abs(MC_type)==" + (TString)to_string(type.first);
	  cut_string  += "&&MC_is_CC==" + (TString)to_string(interaction.first);
	  cut_string  += "&&MC_erange_start==" + (TString)to_string(emin.first);

	  //check the run nr limit for this type of events,
	  //if it is not in the map add it

	  if ( cut_maxrun.find(cut_string) == cut_maxrun.end() ) {
	    cut_maxrun.insert( pair<TString, Int_t>( cut_string, get_max_run_nr(cut_string) ) );
	    cout << "Events: " << cut_string << "\t limit: " << cut_maxrun[cut_string] << endl;
	  } 
	  else {
	    if (run_nr > cut_maxrun[cut_string]) continue;
	  }

	  cut_string  += "&&MC_runID==" + (TString)to_string(run_nr);

	  TString name_string = "../data/mc_end/data_atmnu/summary_"; 
	  name_string += type.second + "-";
	  name_string += interaction.second + "_";
	  name_string += emin.second + "GeV_" + (TString)to_string(run_nr) + ".root";

	  //create file, remove if no entries, continue if file exists

	  FileStat_t stats;
	  Bool_t filefound = !(Bool_t)gSystem->GetPathInfo(name_string, stats);
	  if (filefound && !overwrite) {
	    cout << "File " << name_string << " exists, continuing." << endl;
	    continue;
	  }

	  cout << "Creating file: " << name_string << endl;
	  TFile *fout = new TFile(name_string, "RECREATE");
	  TTree *tout = sp->fChain->CopyTree(cut_string);
	  
	  Bool_t remove = false;
	  if (tout->GetEntries() > 0) { tout->Write(); }
	  else                        { remove = true; }
	  
	  fout->Close();
	  
	  if (remove) system ("rm " + name_string);

	} //end loop over energy range

      } //end loop over nc/cc

    } //end loop over types

  } //end loop over run numbers

}

//****************************************************************

//this function calculates the maximum run number for a specific cut

Int_t get_max_run_nr(TString cut_string) {

  TFile *f_tmp = new TFile("tmp.root","RECREATE");
  TTree *t_tmp = sp->fChain->CopyTree(cut_string);
  Int_t  max_run_nr = t_tmp->GetMaximum("MC_runID");
  f_tmp->Close();
  system("rm tmp.root");

  return max_run_nr;

}

//****************************************************************

//Function to split the pid output to summary files per mc input
//this handles atmospheric muons. This has 11000 files, so sometimes
//it may be better to have just one file with all mupage events

void split_to_files_mupage(Bool_t overwrite, Bool_t separate) { 

  Int_t max_run_nr = sp->fChain->GetMaximum("MC_runID");

  //Loop over run numbers
  for (Int_t run_nr = 1; run_nr <= max_run_nr; run_nr++) {

    TString cut_string  = "MC_runID==" + (TString)to_string(run_nr) + "&&MC_type==-13";
    TString name_string = "../data/mc_end/data_mupage/summary_ORCA115_9m_2016.mupage.nevts40000." + 
      (TString)to_string(run_nr) + ".root";

    if (!separate) {
      cut_string  = "MC_type==-13";
      name_string = "../data/mc_end/data_mupage/summary_ORCA115_9m_2016.mupage.nevts40000.all.root";
    }

    FileStat_t stats;
    Bool_t filefound = !(Bool_t)gSystem->GetPathInfo(name_string, stats);
    if (filefound && !overwrite) {
      cout << "File " << name_string << " exists, continuing." << endl;
      if (separate) continue;
      else break;
    }

    cout << "Creating file: " << name_string << endl;
    TFile *fout = new TFile(name_string, "RECREATE");
    TTree *tout = sp->fChain->CopyTree(cut_string);

    Bool_t remove = false;
    if (tout->GetEntries() > 0) { tout->Write(); }
    else                        { remove = true; }
    
    fout->Close();
    
    if (remove) system ("rm " + name_string);

    if (!separate) break;

  }

}
