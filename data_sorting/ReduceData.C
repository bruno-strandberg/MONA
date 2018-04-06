#include "TROOT.h"
#include "DataReducer.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
using namespace std;

/**
 * This function takes the file NMH/data/pid_result_XXX.root as input and outputs a 
 * file in analysis format.
 *
 * See NMH/data_sorting/README.md for more info.
 *
 * \param  fin_name   PID output file from ECAP, e.g. ../data/pid_result_XXX.root
 * \param  fout_name  A file where the events from input, converted to analysis format, are stored.
 *
 */
void ReduceData(TString fin_name="", TString fout_name="") {

  if (fin_name == "" || fout_name == "") {
    cout << "ERROR! ReduceData() Specify input file name (../data/pid_result_XXX.root) and output file name (../data/ORCA_MC_summary_all_xxx.root). Exiting." << endl;
    return;
  }
  
  gROOT->ProcessLine(".L DataReducer.C+");
  TFile *fin = new TFile(fin_name, "READ");
  TTree *tin = (TTree*)fin->Get("PID");
  DataReducer dr(tin, fout_name);
  dr.Loop();
  fin->Close();
  
}
