#include "TH2.h"
#include "TMath.h"

#include "NMHUtils.h"

#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"
#include "evt/Evt.hh"
#include "evt/Trk.hh"
#include "JAAnet/JAAnetToolkit.hh"
#include "JSupport/JMonteCarloFileSupportkit.hh"
#include "JSupport/JMultipleFileScanner.hh"

#include <map>
#include <vector>
#include <stdexcept>
#include <iostream>

using namespace std;

/**
   Namespace with global variables and functions for app Create_BY_hists
 */
namespace BYHIST {

  std::map<Int_t, TString> fFlavs = { {12, "e"}, {14, "mu"}, {16, "tau"} }; //!< flavor map as in `NuXsec.C`
  std::map<Int_t, TString> fInts  = { {0, "_nc"}, {1, "_cc"}  };            //!< interaction type map as in `NuXsec.C`
  std::map<Int_t, TString> fPols  = { {0, ""}, {1, "_bar"} };               //!< nu polarisations as in `NuXsec.C`
  std::map<Int_t, Int_t>   fPDGtoFlav = { {12, 0}, {14, 1}, {16, 2} };      //!< map from PDG to flavor

  TH2D* fHists[3][2][2]; //!< histogram array with indices [flavor][nc/cc][nu/nub]

  static const UInt_t GSG_NC = 3; //!< integer value to denote NC interaction in gSeaGen (see wiki)

}

/**
   Application to create E vs bjorken-y histograms from gSeaGen data.
*/
int main(const int argc, const char **argv) {

  using namespace BYHIST;

  TString gsg_files;
  TString outname;
  Int_t   nebins;
  Int_t   nbybins;
  JTOOLS::JRange<double> E_MC_range;

  try {

    JParser<> zap("Application to create E vs bjorken-y histograms from gSeaGen data.");

    zap['g'] = make_field(gsg_files, "List of atmospheric neutrino gSeaGen files. For example, create via ls *some_files* > my_list.txt and provide -g my_list.txt as arugment.");
    zap['e'] = make_field(E_MC_range, "Full range of the MC data covered by the input gSeaGen files (for ORCA typically 1-100 GeV)") = JTOOLS::JRange<double>(1, 100);
    zap['o'] = make_field(outname, "Name of the output file") = "output.root";
    zap['E'] = make_field(nebins, "Number of bins on the energy axis") = 40;
    zap['B'] = make_field(nbybins, "Number of bins on the bjorken-y axis; this should be chosen such that it can be re-binned to 2, 4, 6, 8, ... for usage in `common_software/NuXsec`") = 24;

    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  //---------------------------------------------------------------------------------
  // Initialise histograms
  //---------------------------------------------------------------------------------
  Double_t emin = E_MC_range.getLowerLimit();
  Double_t emax = E_MC_range.getUpperLimit();;
  vector<Double_t> e_edges = NMHUtils::GetLogBins(nebins, emin, emax);

  Double_t bymin =  0.;
  Double_t bymax =  1.;

  for (auto f: fPDGtoFlav) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	TString hname = fFlavs[f.first] + p.second + i.second;
	fHists[f.second][i.first][p.first] = new TH2D(hname, hname, nebins, &e_edges[0], nbybins, bymin, bymax);
      }
    }
  }

  //---------------------------------------------------------------------------------
  // loop over input gSeaGen files and fill histograms
  //---------------------------------------------------------------------------------
  JSUPPORT::JMultipleFileScanner_t flist;
  for ( auto gsgf: NMHUtils::ReadLines(gsg_files) ) flist.push_back( (string)gsgf );
  
  JSUPPORT::JMultipleFileScanner<Evt> gp(flist);

  // loop over events in gSeaGen file
  while ( gp.hasNext() ) {
    
    if (gp.getCounter() % 100000 == 0) cout << "Event: " << gp.getCounter() << endl;
    
    const Evt* event = gp.next();
    Trk nu_trk;

    if ( JAANET::has_neutrino(*event) ) {
      nu_trk = JAANET::get_neutrino(*event);
    }
    else {
      throw std::invalid_argument("ERROR! Create_BY_hists() cannot find a neutrino in the MC event.");
    }

    // check that neutrino flavor is found correctly
    if ( fPDGtoFlav.find( TMath::Abs( nu_trk.type ) ) == fPDGtoFlav.end() ) {
      throw std::invalid_argument( "ERROR! Create_BY_hists() unknown flavor " + to_string(nu_trk.type) );
    }

    // determine particle type and fill the histogram
    Int_t flav  = fPDGtoFlav[ TMath::Abs( nu_trk.type ) ];
    Int_t is_cc = (Int_t)( nu_trk.getusr("cc") != GSG_NC );
    Int_t is_nb = (Int_t)( nu_trk.type < 0 );

    fHists[flav][is_cc][is_nb]->Fill( nu_trk.E, nu_trk.getusr("by") );

  }

  //---------------------------------------------------------------------------------
  // write to file and clean up dynamic memory
  //---------------------------------------------------------------------------------

  TFile fout(outname, "RECREATE");
  for (auto f: fPDGtoFlav) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	if (i.first == 0 && f.second != 0) continue; // for NC only write elec
	fHists[f.second][i.first][p.first]->Write();
      }
    }
  }
  fout.Close();

}
