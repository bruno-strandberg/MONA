#ifndef CScalculator_h
#define CScalculator_h

#include "TGraph.h"
#include<map>

using namespace std;

/** A class that provides a neutrino cross-section calculator. 
 *
 *  The class can be used as follows. After init, select the neutrino interaction
 *  with the public function SelectInteraction. Then use the function CalculateCS
 *  to get the cross-section depending on neutrino energy.
 *
 *  The class uses a root file in specific format that stores TGraphs with total neutrino
 *  cross-sections for CC and NC interactions on proton and oxygen. The cross-section data
 *  was extracted from GENIE by M. Jongen.
 * 
 */

class CScalculator {

 public:
  CScalculator(TString xsecfile="../data/cross_sections_gSeaGen_v4r1/xsec.root");
  ~CScalculator();

  void     SelectInteraction(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar);
  Double_t CalculateCS(Double_t E);

 private:
  Bool_t   InitGraphs(TString xsecfile);
  TString  CreateString(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar);

  map < Int_t, TString > nu_flavs;    //!< map with e, mu, tau
  map < Int_t, TString > int_types;   //!< map with nc, cc
  map < Int_t, TString > p_types;     //!< map with "", "bar"

  map < TString, vector<TGraph*> > fGraphs; //!< map with graphs from file

  Bool_t  fGraphsSet; //!< flag to indicate that graph pointers to the selected interaction are set 
  TGraph *f_g_H1;     //!< pointer to a graph with H1 xsec for the selected interaction
  TGraph *f_g_O16;    //!< pointer to a graph with O16 xsec for the selected interaction

};

#endif
