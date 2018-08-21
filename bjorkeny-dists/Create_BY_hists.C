#include "TH2.h"
#include "TMath.h"

#include "GSGParser.h"
#include "NMHUtils.h"

#include <map>
#include <vector>
#include <stdexcept>
#include <iostream>

using namespace std;

/**
   Namespace with global variables and functions for macro Create_BY_hists
 */
namespace BYHIST {

  std::map<Int_t, TString> fFlavs = { {12, "e"}, {14, "mu"}, {16, "tau"} }; //!< flavor map as in `NuXsec.C`
  std::map<Int_t, TString> fInts  = { {0, "_nc"}, {1, "_cc"}  };            //!< interaction type map as in `NuXsec.C`
  std::map<Int_t, TString> fPols  = { {0, ""}, {1, "_bar"} };               //!< nu polarisations as in `NuXsec.C`
  std::map<Int_t, Int_t>   fPDGtoFlav = { {12, 0}, {14, 1}, {16, 2} };      //!< map from PDG to flavor

  TH2D* fHists[3][2][2]; //!< histogram array with indices [flavor][nc/cc][nu/nub]

}

/**
   Function to create E vs bjorken-y histograms from gSeaGen data.

   \param gsg_files  List of atmospheric neutrino gSeaGen files
   \param outname    Name of the output file

 */
void Create_BY_hists(TString gsg_files, TString outname="output.root") {

  using namespace BYHIST;

  //---------------------------------------------------------------------------------
  // Initialise histograms
  //---------------------------------------------------------------------------------
  Int_t ebins   =  40;
  Double_t emin =   1.;
  Double_t emax = 100.;
  vector<Double_t> e_edges = NMHUtils::GetLogBins(ebins, emin, emax);

  //24 bins can be rebinned to 1, 2, 3, 4, 6, 8 or 12 bins in the next step,
  //which in the end will allow for easy usage in `NuXsec` in terms of normalisation 
  Int_t bybins   = 24; 
  Double_t bymin =  0.;
  Double_t bymax =  1.;

  for (auto f: fPDGtoFlav) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	TString hname = fFlavs[f.first] + p.second + i.second;
	fHists[f.second][i.first][p.first] = new TH2D(hname, hname, ebins, &e_edges[0], bybins, bymin, bymax);
      }
    }
  }

  //---------------------------------------------------------------------------------
  // loop over input gSeaGen files and fill histograms
  //---------------------------------------------------------------------------------
  vector<TString> gsgfiles = NMHUtils::ReadLines(gsg_files);

  for (auto gsg_fname: gsgfiles) {

    cout << "Looping file: " << gsg_fname << endl;

    // initialise gSeaGen parser
    GSGParser gp(gsg_fname);

    // loop over events in gSeaGen file
    while ( gp.NextEvent() ) {
      
      // check that neutrino flavor is found correctly
      if ( fPDGtoFlav.find( TMath::Abs(gp.Neutrino_PdgCode) ) == fPDGtoFlav.end() ) {
	throw std::invalid_argument( "ERROR! Create_BY_hists() unknown flavor "
				     + to_string(gp.Neutrino_Type) );
      }

      // determine particle type and fill the histogram
      Int_t flav  = fPDGtoFlav[ TMath::Abs(gp.Neutrino_PdgCode) ];
      Int_t is_cc = ( (Int_t)gp.InterId != 3 );
      Int_t is_nb = (Int_t)( gp.Neutrino_PdgCode < 0 );

      fHists[flav][is_cc][is_nb]->Fill( gp.Neutrino_E, gp.By );

    }

  }

  //---------------------------------------------------------------------------------
  // write to file and clean up dynamic memory
  //---------------------------------------------------------------------------------

  TFile fout(outname, "RECREATE");
  for (auto f: fPDGtoFlav) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	if (i.first == 0 && f.first != 0) continue; // for NC only write elec
	fHists[f.second][i.first][p.first]->Write();
      }
    }
  }
  fout.Close();

}
