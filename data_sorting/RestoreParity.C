#include "TROOT.h"
#include "TSystem.h"
#include "SummaryParser.h"


//****************************************************************

//functions

void  split_to_files_atmnu(Bool_t overwrite=false);
void  split_to_files_atmnu_fast(Int_t max_file_handles=950);
void  split_to_files_mupage(Bool_t overwrite=false, Bool_t separate=false);
Int_t get_max_run_nr(TString cut_string);

//global variables accessible in all functions inside this macro

/// Global pointer to a summary file parser class.
SummaryParser *sp;

//****************************************************************

/** (main function) This function takes the file resulting from ReduceData.C as argument 
 *  and splits it up to NMH/data/mc_end/ directories to match the file scheme used throughout the
 *  MC chain.
 *
 * See NMH/data/README.md for more details.
 *
 * \param  fname PID summary data in analysis format, output by ReduceData.C.
 *
 */
void RestoreParity(TString fname) {

  //initialise the data parser
  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");
  sp = new SummaryParser(fname);

  if (sp->fChain == NULL) {
    cout << "ERROR! RestoreParity() Issue with SummaryParser() initialisation." << endl;
  }

  //call functions
  //split_to_files_atmnu();
  split_to_files_atmnu_fast();
  split_to_files_mupage();
  
}

//****************************************************************

/** Function to split the pid output to summary files per mc input, atm neutrinos.
 *
 *  The function split_to_files_atmnu_fast() does the same thing much faster. However, this is simpler
 *  to understand and kept for comparison.
 *
 * \param overwrite If file exists in NMH/data/mc_end/data_atmnu/, ovewrite.
 *
 */
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
	  
	  if (remove) Int_t sysret = system ("rm " + name_string);

	} //end loop over energy range

      } //end loop over nc/cc

    } //end loop over types

  } //end loop over run numbers

}

//****************************************************************

/** This function calculates the maximum run number for a specific cut.
 *
 * \param cut_string Cut string to select events opened for reading by sp.
 *
 */
Int_t get_max_run_nr(TString cut_string) {

  TFile *f_tmp = new TFile("tmp.root","RECREATE");
  TTree *t_tmp = sp->fChain->CopyTree(cut_string);
  Int_t  max_run_nr = t_tmp->GetMaximum("MC_runID");
  f_tmp->Close();
  Int_t sysret = system("rm tmp.root");

  return max_run_nr;

}

//****************************************************************

/** Function to split the pid output to summary files per mc input, atm muons.
 *
 * \param overwrite If file exists in NMH/data/mc_end/data_atmnu/, ovewrite.
 * \param separate Separate MUPAGE files (11k files), if false MUPAGE events in one file in 
 *        NMH/data/mc_end/data_atmmu/
 *
 */
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
    
    if (remove) Int_t sysret = system ("rm " + name_string);

    if (!separate) break;

  }

}

//****************************************************************

/** Function to split the pid output in summary format to summary files per mc input, atm neutrinos.
 *
 *  The following is for anyone looking to improve/understand the function, it is not relevant if you
 *  are only interested in running the code.
 *
 *  The idea is simple. The input file to RestoreParity() is split to >5000 summary files to match the
 *  file numbering scheme employed in the MC chain. The fastest way to achieve this is to create >5000
 *  files, make >5000 empty clones of the input tree with entries, loop over the input tree once and
 *  depending on the run number call Fill() function of one of the >5000 trees.
 *
 *  However, typically there is a limit to how many files can be opened in parallel (e.g. 1024). For this
 *  reason, the function needs to select events in several steps before it can open files in parallel.
 *
 *  \param max_file_handles Maximum nr of files opened in parallel.
 *
 */
void  split_to_files_atmnu_fast(Int_t max_file_handles) {

  //maps to build the filenames
  map < Double_t, TString > pdg_to_name = { {12., "elec" }, 
  					    {14., "muon" },
  					    {16., "tau"  } };
  
  map < Double_t, TString > iscc_to_str = { { 0., "NC"}, 
					    { 1., "CC"} };

  map < Double_t, TString>  emin_to_str = { { 1., "1-5"  },
  					    { 3., "3-100"} };
    
  for (auto const& type : pdg_to_name) {
    for (auto const& interaction : iscc_to_str) {
      for (auto const& emin : emin_to_str) {

	cout << "Processing " << type.second << " " << interaction.second
	     << " " << emin.second << endl;
	
	//-------------------------------------------------------
	//build the cut string and the file name string
	//-------------------------------------------------------
	
	TString cut_string = "";
	cut_string  += "TMath::Abs(MC_type)==" + (TString)to_string(type.first);
	cut_string  += "&&MC_is_CC==" + (TString)to_string(interaction.first);
	cut_string  += "&&MC_erange_start==" + (TString)to_string(emin.first);

	TString name_string = "../data/mc_end/data_atmnu/summary_"; 
	name_string += type.second + "-";
	name_string += interaction.second + "_";
	name_string += emin.second + "GeV_";

	//-------------------------------------------------------
	// create level 1 tree for flavor_NC/CC_erange over all file numbers
	//-------------------------------------------------------
	
	TFile *fout_lev1   = new TFile(name_string + "allruns.root", "RECREATE");
	TTree *tout_lev1   = sp->fChain->CopyTree(cut_string);
	Int_t max_run_nr   = 0;
	if (tout_lev1->GetEntries() > 0) { max_run_nr = tout_lev1->GetMaximum("MC_runID"); }

	//-------------------------------------------------------
	// create level 2 trees that have subsetsize  # of runs
	//-------------------------------------------------------
	
	vector<TFile*> fouts_lev2;                             //vector of lev2 files
	vector<TTree*> touts_lev2;                             //vector of lev2 trees
	vector<TString> temp_files = { fout_lev1->GetName() }; //vector with files to be cleaned up

	Int_t subsetsize = max_file_handles;                   //lev2 files have this (max) # of runs each
	Int_t subsets    = 0;                                  //while loop counter

	while ( max_run_nr > 0 && subsets * subsetsize < max_run_nr ) {

	  Int_t runnr_min = subsetsize * subsets;
	  Int_t runnr_max = subsetsize * (subsets+1);
	  TString fname_lev2      = name_string + ( (TString)to_string(runnr_min) + "-" +
						    (TString)to_string(runnr_max) + ".root" );
	  TString cut_string_lev2 = cut_string + ("&&MC_runID>=" + (TString)to_string(runnr_min) +
						  "&&MC_runID<"  + (TString)to_string(runnr_max) );

	  fouts_lev2.push_back( new TFile(fname_lev2, "RECREATE") );
	  touts_lev2.push_back( tout_lev1->CopyTree(cut_string_lev2) );
	  temp_files.push_back( fname_lev2 );

	  //-------------------------------------------------------
	  // now create subsetsize # of level 3 files, fill them in parallel during one for loop, cleanup
	  //-------------------------------------------------------

	  vector<TFile*> fouts_lev3;
	  vector<TTree*> touts_lev3;
	  for (Int_t rnr = runnr_min; rnr < runnr_max; rnr++) {
	    TString fname_lev3 = name_string + (TString)to_string(rnr) + ".root";
	    fouts_lev3.push_back( new TFile(fname_lev3,"RECREATE") );
	    touts_lev3.push_back( touts_lev2.back()->CloneTree(0) );
	  }
	  
	  Double_t MC_runID;
	  touts_lev2.back()->SetBranchAddress("MC_runID", &MC_runID);
	  for (Int_t evt = 0; evt < touts_lev2.back()->GetEntries(); evt++) {
	    touts_lev2.back()->GetEntry(evt);
	    touts_lev3[MC_runID - subsets * subsetsize]->Fill();
	  }
	  
	  for (Int_t fi = 0; fi < (Int_t)fouts_lev3.size(); fi++) {
	    fouts_lev3[fi]->cd();
	    if ( touts_lev3[fi]->GetEntries() == 0 ) temp_files.push_back( fouts_lev3[fi]->GetName() );
	    touts_lev3[fi]->Write();
	    fouts_lev3[fi]->Close();
	    delete fouts_lev3[fi];
	  }

	  subsets++;
	  
	} //end while loop

	cout << "Removing temp/empty files" << endl;

	fout_lev1->Close();
	for (auto f: fouts_lev2) {  f->Close();  delete f; }
	for (auto fname: temp_files) Int_t sysret = system ("rm " + fname);
	
      } //end loop over energy range
    } //end loop over nc/cc
  } //end loop over types

}
