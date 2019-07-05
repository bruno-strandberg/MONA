#include "FileHeader.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"

#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

#include <stdexcept>

using namespace std;

int main(const int argc, const char **argv) {

  //====================================================================================
  // parse input
  //====================================================================================

  TString  gspl2root_output;
  TString  output;

  TString  gsgenv;
  TString  splinefile;
  Double_t emax;
  
  try {

    JParser<> zap("Application to reduce the gspl2root output to data only necessary the NuXsec class in MONA/common_software/");

    zap['f'] = make_field(gspl2root_output, "Output of GENIE's 'gspl2root' tool that contains elec, muon, tau xsec data");
    zap['o'] = make_field(output, "Output file that can be fed to the NuXsec class");
    zap['g'] = make_field(gsgenv, "gSeaGen environment used for the creation of the spline file (for header only)");
    zap['x'] = make_field(splinefile, "Spline file that was used to create 'gspl2root' output (for header only)");
    zap['e'] = make_field(emax, "Maximum energy input to 'gspl2root' (for header only)");

    if ( zap.read(argc, argv) != 0 ) return 1;

  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  //====================================================================================
  // create the graph names expected in the gspl2root output file
  //====================================================================================

  vector<TString> neutrinos = {"nu_e", "nu_e_bar", "nu_mu", "nu_mu_bar", "nu_tau", "nu_tau_bar"};
  vector<TString> targets   = {"H1", "O16"};
  vector<TString> xsecs     = {"tot_cc","tot_nc"};

  vector<TString> expected_graphs; // names of the graphs expected in the gspl2root output file
  vector<TString> dir_names;       // names of the expected directories

  for (auto nu: neutrinos) {
    for (auto target: targets) {

      dir_names.push_back(nu + "_" + target);

      for (auto xsec: xsecs) {
	TString graph = nu + "_" + target + "/" + xsec;
	expected_graphs.push_back(graph);
      }

    }
  }

  //====================================================================================
  // copy the expected graphs to the output
  //====================================================================================

  TFile fin(gspl2root_output, "READ");
  if ( !fin.IsOpen() ) {
    throw std::invalid_argument("ERROR! Reductor: cannot open the input file " + (string)gspl2root_output);
  }

  //---------------------------------------------------------
  // create directories in the output file
  //---------------------------------------------------------

  TFile fout(output, "RECREATE");

  // create header and write info
  FileHeader head("Reductor");
  head.AddParameter("GSEAGEN", gsgenv);
  head.AddParameter("SPLINE", splinefile);
  head.AddParameter("EMAX", (TString)to_string(emax) );
  head.WriteHeader(&fout);

  std::map<TString, TDirectory*> dirs;

  for (auto dir: dir_names) {
    dirs.insert( std::make_pair(dir, fout.mkdir(dir) ) );
  }

  //---------------------------------------------------------
  // copy the graphs to the output files
  //---------------------------------------------------------

  for (auto g: expected_graphs) {
    
    // get the graph
    
    TGraph *gin = (TGraph*)fin.Get( g );

    if ( gin == NULL ) {
      throw std::invalid_argument("ERROR! Reductor: cannot find graph " + (string)gin->GetName() + "in the input file " + (string)gspl2root_output );
    }

    // change to the correct directory and write the graph

    TString dir = g( 0, g.First('/') );
    dirs[dir]->cd();
    gin->Write();

  }

  fout.Close();
  fin.Close();

  return 0;

}
