
// root headers
#include "TSystem.h"

// NMH headers
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "FileHeader.h"

// cpp headers
#include <stdexcept>

// jpp headers
#include "JTools/JRange.hh"
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

/** A namespace that collects functions and variables for the `RestoreParity` application.*/
namespace RESTOREPARITY {

  using namespace JTOOLS;

  void  split_to_files_atmnu_fast(TString tag, JRange<double> E_low, JRange<double> E_high, Int_t max_file_handles=950);
  void  split_to_files_mupage(TString tag, Bool_t overwrite=false, Bool_t separate=false);
  void  split_to_files_noise(TString tag);
  Int_t get_max_run_nr(TString cut_string);
  void  check_return(Int_t ret);

  /// Global pointer to a summary file parser class.
  SummaryParser *fSp;
  TString tagpar     = "datatag";
  TString summarydir = (TString)getenv("NMHDIR") + "/data/mcsummary/";
  TString gseagendir = (TString)getenv("NMHDIR") + "/data/gseagen/";
  TString muondir    = "/data_mupage/";
  TString noisedir   = "/data_noise/";
  TString nudir      = "/data_atmnu/";

};

using namespace JTOOLS;
using namespace RESTOREPARITY;

/** This program takes the file resulting from RDFPID_to_Summary as argument and splits it up to NMH/data/mcsummary/ directories to match the file scheme used throughout the MC chain.
    
   See NMH/data/README.md for more details.
 
 */
int main(int argc, char **argv) {

  //-------------------------------------------------------------------
  //parse command line arguments
  //-------------------------------------------------------------------

  string fname;
  JRange<double> E_low;
  JRange<double> E_high;

  try {
    JParser<> zap("This program takes the file resulting from RDFPID_to_Summary as argument and splits it up to NMH/data/mcsummary/ directories to match the file scheme used throughout the MC chain.");

    zap['f'] = make_field(fname , "Name of the file created by RDFPID_to_Summary");
    zap['l'] = make_field(E_low , "For ORCA, typically two energy ranges are used with gSeaGen, events from both are present in the -f argument file. Provide the low range") = JRange<double>(1, 5);
    zap['u'] = make_field(E_high, "For ORCA, typically two energy ranges are used with gSeaGen, events from both are present in the -f argument file. Provide the high range") = JRange<double>(3, 100);

    zap(argc, argv);
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  if (fname == "") {
    throw std::invalid_argument("ERROR! RestoreParity() input file not specified!");
  }

  //-------------------------------------------------------------------
  //initialise the data parser; read the tag from the header
  //-------------------------------------------------------------------
  FileHeader head("RestoreParity");
  head.ReadHeader(fname);
  fSp = new SummaryParser(fname);

  if (fSp->GetTree() == NULL) {
    throw std::invalid_argument("ERROR! RestoreParity() Issue with SummaryParser() initialisation.");
  }

  //-------------------------------------------------------------------
  //call functions
  //-------------------------------------------------------------------
  split_to_files_atmnu_fast(head.GetParameter(tagpar), E_low, E_high);
  split_to_files_mupage(head.GetParameter(tagpar), true, false);
  split_to_files_noise(head.GetParameter(tagpar));
  
}

//****************************************************************

/** Function to check the return value of the system call 
    \param ret  return value of `system`.
*/
void RESTOREPARITY::check_return(Int_t ret) {

  if (ret != 0) {
    cout << "WARNING! RestoreParity::check_return() system returned " << ret << endl;
  }

}

//****************************************************************

/** This function calculates the maximum run number for a specific cut.
 *
 * \param cut_string Cut string to select events opened for reading by sp.
 *
 */
Int_t RESTOREPARITY::get_max_run_nr(TString cut_string) {

  TFile *f_tmp = new TFile("tmp.root","RECREATE");
  TTree *t_tmp = fSp->GetTree()->CopyTree(cut_string);
  Int_t  max_run_nr = t_tmp->GetMaximum("fMC_runID");
  f_tmp->Close();
  check_return( system("rm tmp.root") );

  return max_run_nr;

}

//****************************************************************

/** Function to split the pid output to summary files per mc input, atm muons.
 *
 *  Typically, having atm muon events in one file is better and splitting is not necessary.
 *
 * \param data tag of the production
 * \param overwrite If file exists in NMH/data/mc_end/data_atmnu/tag/, ovewrite.
 * \param separate Separate MUPAGE files (11k files), if false MUPAGE events in one file in NMH/data/mc_end/data_atmmu/tag/
 *
 */
void RESTOREPARITY::split_to_files_mupage(TString tag, Bool_t overwrite, Bool_t separate) { 

  TString sumdir = summarydir + tag + muondir;
  TString gsgdir = gseagendir + tag + muondir;
  check_return( system("mkdir -p " + sumdir) );
  check_return( system("mkdir -p " + gsgdir) );
  cout << "NOTICE RESTOREPARITY::split_to_files_mupage() made directories: " << endl << "\t" << sumdir << "\t" << endl << "\t" << gsgdir << endl;
  
  Int_t max_run_nr = fSp->GetTree()->GetMaximum("fMC_runID");

  //Loop over run numbers
  for (Int_t run_nr = 1; run_nr <= max_run_nr; run_nr++) {

    TString cut_string  = "fMC_runID==" + (TString)to_string(run_nr) + "&&fMC_type==-13";
    TString name_string = summarydir + tag + muondir + "/summary_mupage." + (TString)to_string(run_nr) + ".root";

    if (!separate) {
      cut_string  = "fMC_type==-13";
      name_string = summarydir + tag + muondir + "/summary_mupage.all.root";
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
    TTree *tout = fSp->GetTree()->CopyTree(cut_string);

    Bool_t remove = false;
    if (tout->GetEntries() > 0) { tout->Write(); }
    else                        { remove = true; }
    
    // create header
    FileHeader head("RestoreParity");
    head.AddParameter(tagpar, tag);
    head.WriteHeader(fout);

    fout->Close();
    
    if (remove) check_return( system ("rm " + name_string) );

    if (!separate) break;

  }

}

//****************************************************************

/** Function to separate noise files from pid_result file.
 *
 *  This is more lightweight compared to atm nu and atm muons. After some experience with the analysis,
 *  it turned out that typically there is no reason to split noise and atm muon files. The noise and
 *  atm muon event weighting is simpler (does not require effective mass calculation and oscillation),
 *  meaning there are no comparisons to be made with gSeaGen files, which is the main reason why atm nu
 *  files are split.
 *
 */
void RESTOREPARITY::split_to_files_noise(TString tag) {

  TString sumdir = summarydir + tag + noisedir;
  TString gsgdir = gseagendir + tag + noisedir;
  check_return( system("mkdir -p " + sumdir) );
  check_return( system("mkdir -p " + gsgdir) );
  cout << "NOTICE RESTOREPARITY::split_to_files_noise() made directories: " << endl << "\t" << sumdir << "\t" << endl << "\t" << gsgdir << endl;

  TString cut_string  = "fMC_type==0";
  TString name_string = summarydir + tag + noisedir + "/summary_noise.root";

  cout << "Creating file: " << name_string << endl;
  TFile *fout = new TFile(name_string, "RECREATE");
  TTree *tout = fSp->GetTree()->CopyTree(cut_string);

  Bool_t remove = false;
  if (tout->GetEntries() > 0) { tout->Write(); }
  else                        { remove = true; }

  // create header
  FileHeader head("RestoreParity");
  head.AddParameter(tagpar, tag);
  head.WriteHeader(fout);
  
  fout->Close();
  
  if (remove) check_return( system ("rm " + name_string) );
  
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
void  RESTOREPARITY::split_to_files_atmnu_fast(TString tag, JRange<double> E_low, JRange<double> E_high, Int_t max_file_handles) {


  TString sumdir = summarydir + tag + nudir;
  TString gsgdir = gseagendir + tag + nudir;
  check_return( system("mkdir -p " + sumdir) );
  check_return( system("mkdir -p " + gsgdir) );
  cout << "NOTICE RESTOREPARITY::split_to_files_atmnu_fast() made directories: " << endl << "\t" << sumdir << "\t" << endl << "\t" << gsgdir << endl;

  //maps to build the filenames
  map < Double_t, TString > pdg_to_name = { {12., "elec" }, 
  					    {14., "muon" },
  					    {16., "tau"  } };
  
  map < Double_t, TString > iscc_to_str = { { 0., "NC"}, 
					    { 1., "CC"} };

  TString lowrange  = to_string( (int)E_low.getLowerLimit() ) + "-" + to_string( (int)E_low.getUpperLimit() );
  TString highrange = to_string( (int)E_high.getLowerLimit() ) + "-" + to_string( (int)E_high.getUpperLimit() );

  map < Double_t, TString>  emin_to_str = { {  E_low.getLowerLimit(), lowrange }, 
					    { E_high.getLowerLimit(), highrange} };
    
  for (auto const& type : pdg_to_name) {
    for (auto const& interaction : iscc_to_str) {
      for (auto const& emin : emin_to_str) {

	cout << "Processing " << type.second << " " << interaction.second
	     << " " << emin.second << endl;
	
	//-------------------------------------------------------
	//build the cut string and the file name string
	//-------------------------------------------------------
	
	TString cut_string = "";
	cut_string  += "TMath::Abs(fMC_type)==" + (TString)to_string(type.first);
	cut_string  += "&&fMC_is_CC==" + (TString)to_string(interaction.first);
	cut_string  += "&&fMC_erange_start==" + (TString)to_string(emin.first);

	TString name_string = summarydir + tag + nudir + "/summary_"; 
	name_string += type.second + "-";
	name_string += interaction.second + "_";
	name_string += emin.second + "GeV_";

	//-------------------------------------------------------
	// create level 1 tree for flavor_NC/CC_erange over all file numbers
	//-------------------------------------------------------
	
	TFile *fout_lev1   = new TFile(name_string + "allruns.root", "RECREATE");
	TTree *tout_lev1   = fSp->GetTree()->CopyTree(cut_string);
	Int_t max_run_nr   = 0;
	if (tout_lev1->GetEntries() > 0) { max_run_nr = tout_lev1->GetMaximum("fMC_runID"); }

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
	  TString cut_string_lev2 = cut_string + ("&&fMC_runID>=" + (TString)to_string(runnr_min) +
						  "&&fMC_runID<"  + (TString)to_string(runnr_max) );

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
	  
	  SummaryEvent *sumevt = new SummaryEvent();
	  touts_lev2.back()->SetBranchAddress("SummaryEvent", &sumevt);
	  for (Int_t evt = 0; evt < touts_lev2.back()->GetEntries(); evt++) {
	    touts_lev2.back()->GetEntry(evt);
	    touts_lev3[sumevt->Get_MC_runID() - subsets * subsetsize]->Fill();
	  }
	  delete sumevt;

	  for (Int_t fi = 0; fi < (Int_t)fouts_lev3.size(); fi++) {
	    fouts_lev3[fi]->cd();
	    if ( touts_lev3[fi]->GetEntries() == 0 ) temp_files.push_back( fouts_lev3[fi]->GetName() );
	    touts_lev3[fi]->Write();

	    // create header and add to file
	    FileHeader head("RestoreParity");
	    head.AddParameter(tagpar, tag);
	    head.WriteHeader(fouts_lev3[fi]);

	    fouts_lev3[fi]->Close();
	    delete fouts_lev3[fi];
	  }

	  subsets++;
	  
	} //end while loop

	cout << "Removing temp/empty files" << endl;

	fout_lev1->Close();
	for (auto f: fouts_lev2) {  f->Close();  delete f; }
	for (auto fname: temp_files) check_return( system ("rm " + fname) );
	
      } //end loop over energy range
    } //end loop over nc/cc
  } //end loop over types

}
