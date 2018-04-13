#include "TROOT.h"
#include "TSystem.h"
#include "SummaryParser.h"


//****************************************************************

//functions

void  split_to_files_atmnu(Bool_t overwrite=false);
void  split_to_files_atmnu_v2();
void  split_to_files_atmnu_v3();
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
  split_to_files_atmnu_v3();
  split_to_files_mupage();
  
}

//****************************************************************

/** Function to split the pid output to summary files per mc input, atm neutrinos.
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

//developing something faster. Do the sorting in one go with a loop over the file

void  split_to_files_atmnu_v2() {

  //--------------------------------------------------
  //maps to sort data
  //--------------------------------------------------

  map < Double_t, TString > pdg_to_name = { {12., "elec" }, 
					    {14., "muon" },
					    {16., "tau"  } };

  map < Double_t, TString > iscc_to_str = { { 0., "NC"}, 
					    { 1., "CC"} };

  map < Double_t, TString>  emin_to_str = { { 1., "1-5"  },
					    { 3., "3-100"} };

  //filename (except run number) and a vector with cuts for MC_type, MC_is_CC and MC_erange_start
  map < TString, vector<Double_t> > fnstr_to_cuts;

  for (auto const& type : pdg_to_name) {
    for (auto const& interaction : iscc_to_str) {
      for (auto const& emin : emin_to_str) {

	//build the filename string except for run number
	TString name_string = "../data/mc_end/data_atmnu/summary_"; 
	name_string += type.second + "-";
	name_string += interaction.second + "_";
	name_string += emin.second + "GeV_";

	//create the map entry for fnstr_to_cuts
	vector<Double_t> cuts = { type.first, interaction.first, emin.first };

	fnstr_to_cuts.insert( pair<TString, vector<Double_t>>(name_string, cuts) );

      } //end loop over energy range
    } //end loop over nc/cc
  } //end loop over types

  //--------------------------------------------------
  //init the output files and trees
  //--------------------------------------------------
  
  //const Int_t maxruns = 2000;
  const Int_t maxruns = 20;

  TFile ***files = new TFile**[ fnstr_to_cuts.size() ];
  TTree ***trees = new TTree**[ fnstr_to_cuts.size() ];
  for (UInt_t i = 0; i < fnstr_to_cuts.size(); i++) {
    files[i] = new TFile*[maxruns];
    trees[i] = new TTree*[maxruns];
  }

  for ( map< TString, vector<Double_t> >::iterator it = fnstr_to_cuts.begin(); 
	it != fnstr_to_cuts.end(); ++it) {

    for (Int_t runnr = 0; runnr < maxruns; runnr++) {

      TString fname     = it->first + (TString)to_string(runnr) + ".root";
      Int_t   fname_idx = distance( fnstr_to_cuts.begin(), it );

      files[fname_idx][runnr] = new TFile(fname, "RECREATE");
      trees[fname_idx][runnr] = sp->fChain->CloneTree(0);

    }

  }

  //--------------------------------------------------
  //loop over all entries in the input tree and fill one of the created trees
  //depending on the file number and cuts
  //--------------------------------------------------

  for (Int_t i = 0; i < sp->fChain->GetEntries(); i++) {

    sp->fChain->GetEntry(i);
 
    if (i % 100000 == 0) cout << "Entry: " << i << endl;
   
    //loop over all possible file types, using cuts determine what tree to fill
    for ( map< TString, vector<Double_t> >::iterator it = fnstr_to_cuts.begin(); 
	  it != fnstr_to_cuts.end(); ++it) {

      vector<Double_t> cuts = it->second;
      Int_t fname_idx = distance( fnstr_to_cuts.begin(), it );

      //cuts are MC_type, MC_is_CC and MC_erange_start
      if ( ( TMath::Abs(sp->MC_type) == cuts[0] ) && 
	   ( sp->MC_is_CC            == cuts[1] ) &&
	   ( sp->MC_erange_start     == cuts[2] ) ) {

	trees[fname_idx][(Int_t)sp->MC_runID]->Fill();

      }

    } //end loop over map entries

  } //end loop over entries in file

  //--------------------------------------------------
  //close the files, delete files with 0 entries
  //--------------------------------------------------
  for (UInt_t ftype = 0; ftype < fnstr_to_cuts.size(); ftype++) {

    for (Int_t runnr = 0; runnr < maxruns; runnr++) {

      TString fname = files[ftype][runnr]->GetName();

      Bool_t remove = false;
      if ( trees[ftype][runnr]->GetEntries() > 0 ) { trees[ftype][runnr]->Write(); }
      else { remove = true; }
	  
      files[ftype][runnr]->Close();
	  
      if (remove) Int_t sysret = system ("rm " + fname);

      delete files[ftype][runnr];

    } //end loop over file numbers

  } //end loop over file types

}

//****************************************************************

//developing something faster, attempt 3

void  split_to_files_atmnu_v3() {

  //maps to build the filenames
  // map < Double_t, TString > pdg_to_name = { {12., "elec" }, 
  // 					    {14., "muon" },
  // 					    {16., "tau"  } };
  
  map < Double_t, TString > iscc_to_str = { { 0., "NC"}, 
					    { 1., "CC"} };

  // map < Double_t, TString>  emin_to_str = { { 1., "1-5"  },
  // 					    { 3., "3-100"} };

  map < Double_t, TString > pdg_to_name = { {14., "muon" } };
  map < Double_t, TString>  emin_to_str = { { 3., "3-100"} };
    
  for (auto const& type : pdg_to_name) {

    for (auto const& interaction : iscc_to_str) {

      for (auto const& emin : emin_to_str) {

	//build the cut string and the file name string

	TString cut_string = "";
	cut_string  += "TMath::Abs(MC_type)==" + (TString)to_string(type.first);
	cut_string  += "&&MC_is_CC==" + (TString)to_string(interaction.first);
	cut_string  += "&&MC_erange_start==" + (TString)to_string(emin.first);

	TString name_string = "../data/mc_end/data_atmnu/summary_"; 
	name_string += type.second + "-";
	name_string += interaction.second + "_";
	name_string += emin.second + "GeV_";

	// create the tree for flavor_NC/CC_erange over all file numbers

	TString fname_lev1 = name_string + "allruns.root";
	TFile *fout_lev1   = new TFile(fname_lev1, "RECREATE");
	TTree *tout_lev1   = sp->fChain->CopyTree(cut_string);
	Int_t max_run_nr   = 0;
	if (tout_lev1->GetEntries() > 0) {
	  max_run_nr = tout_lev1->GetMaximum("MC_runID");
	}

	cout << "DEBUG: Created file " << fout_lev1->GetName() << endl;
	
	// create level 2 trees that have subsetsize  # of runs
	
	vector<TFile*> fouts_lev2; //vector of lev2 files
	vector<TTree*> touts_lev2; //vector of lev2 trees
	Int_t subsetsize = 500;    //level 2 files have 500 files (max) each
	Int_t subsets    = 0;      //while loop counter

	while ( max_run_nr > 0 && subsets * subsetsize < max_run_nr ) {

	  // create fname and cut string for level 2 files & create the file/tree
	  
	  Int_t runnr_min = subsetsize * subsets;
	  Int_t runnr_max = subsetsize * (subsets+1);
	  TString fname_lev2 = name_string + (TString)to_string(runnr_min) + "-" +  (TString)to_string(runnr_max) + ".root";
	  TString cut_string_lev2 = cut_string + ("&&MC_runID>=" + (TString)to_string(runnr_min) +
						  "&&MC_runID<"  + (TString)to_string(runnr_max) );

	  fouts_lev2.push_back( new TFile(fname_lev2, "RECREATE") );
	  touts_lev2.push_back( tout_lev1->CopyTree(cut_string_lev2) );
	  TTree *tout_lev2 = touts_lev2.back(); //pointer to current tree, for convenience

	  cout << "DEBUG: Created file " << fouts_lev2.back()->GetName() << endl;
	  
	  // now create subsetsize # of files

	  vector<TFile*> files;
	  vector<TTree*> trees;
	  for (Int_t rnr = runnr_min; rnr < runnr_max; rnr++) {
	    TString fname_lev3 = name_string + (TString)to_string(rnr) + ".root";
	    files.push_back( new TFile(fname_lev3,"RECREATE") );
	    trees.push_back( tout_lev2->CloneTree(0) );
	  }

	  //fill the files in parallel

	  cout << "DEBUG: Filling files." << endl;
	  
	  Double_t MC_runID;
	  tout_lev2->SetBranchAddress("MC_runID", &MC_runID);
	  for (Int_t evt = 0; evt < tout_lev2->GetEntries(); evt++) {
	    tout_lev2->GetEntry(evt);
	    trees[MC_runID - subsets * subsetsize]->Fill();
	  }

	  //close the files where there are events, delete others

	  cout << "DEBUG: Closing/deleting files." << endl;
	  
	  for (Int_t fi = 0; fi < files.size(); fi++) {
	    files[fi]->cd();
	    Bool_t  remove = trees[fi]->GetEntries() == 0;
	    TString fname  = files[fi]->GetName();
	    trees[fi]->Write();
	    files[fi]->Close();
	    if (remove) {
	      TString syscmd = "rm " + fname;
	      Int_t sysret = system(syscmd);
	      cout << "DEBUG: Removing " << fname << "\t" << syscmd << ", system return " << sysret << endl;
	    }
	  }

	  subsets++;
	  if(subsets>0) break;//debug
	  
	} //end while loop

	// clean up lev 2

	cout << "DEBUG: Level 2 cleanup." << endl;
	
	for (auto f: fouts_lev2) {
	  TString fname = f->GetName();
	  f->Close();
	  Int_t sysret = system ("rm " + fname);
	  cout << sysret << endl;
	}
	
	// clean up lev1

	cout << "DEBUG: Level 1 cleanup." << endl;
	
	fout_lev1->Close();
	Int_t sysret = system ("rm " + fname_lev1);
	
      } //end loop over energy range

    } //end loop over nc/cc

  } //end loop over types

}
