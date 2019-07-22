#ifndef EMextr_h
#define EMextr_h

#include "EffMass.h"
#include "TH3D.h"
#include <map>

/** This class is an extension of the `EffMass` class to be able to estimate the effective mass value outside the Monte-Carlo simulation range.

    Extrapolation to energies outside the Monte-Carlo range is required for the correct implementation of the energy-scale systematic. As a short example, let the usual ORCA MC simulation be considered, in the energy range 1-100 GeV. The wish is to have the energy scale systematic affect true energies, e.g. to account for a systematic shift in the neutrino energies as calculated by GENIE/gSeaGen. In practice, this means that the events on the energy boundaries (1 GeV or 100 GeV) are shifted outside the simulation range. As the effective mass is calculated based on the MC data, the effective mass needs to be extrapolated to these energies.

    This is achieved as follows. The user inputs a `TH3D` template histogram, which determines the binning configuration and the energy, cos-theta & bjorken-y limits for the desired effective mass calculation. The user also provides a MONA effective mass file (created with `MONA/apps/effective_mass` apps). This class initialises the `EffMass` class with a binning configuration that aims to match the number of bins in each axis of the template histogram inside the Monte-Carlo range. E.g., if the template histogram has 60 energy bins in the range 0.1-150, the class will count the number of bins in the MC range 1-100 and initialise `EffMass`. Then, the bins inside the MC range are filled by using `EffMass::GetMeff` with interpolation. As a final step, at each (cos-theta, bjorken-y) value, the 1D effective mass vs energy curve is dumped to a `TGraph` to perform simple spline extrapolation to energies outside the MC range.

    Of course, such extrapolation has limited accuracy. On the other hand, the effective mass quickly drops to 0 below 1 GeV, and the neutrino flux becomes relatively weak compared to the signal region at energies abouve 100 GeV. Thus, the extrapolated bins play affect the measurement weakly, but nevertheless this step is necessary for the implementation of the energy scale systematic in `fitter_software/FitUtilWsyst`.

*/
class EMextr {

 public:
  EMextr(TH3D* tb, TString effmfile);
  ~EMextr();
  Double_t GetMeff(Int_t flavor, Bool_t iscc, Bool_t isnb, 
		   Double_t E_true, Double_t Ct_true, Double_t By_true);
  TH3D* GetMeff3DH(Int_t flavor, Bool_t iscc, Bool_t isnb);
  TH1D* GetSlice(Int_t flavor, Bool_t iscc, Bool_t isnb, Double_t ct, Double_t by);

 private:

  // private functions
  Int_t GetEMbins(TAxis *R, TAxis *E, Int_t bins_min, Int_t bins_max);
  void  FillDataMCrange(Bool_t interpolate);
  void  Extrapolate(Bool_t extrapolate);

  // private members
  
  TH3D    *fhTB; //!< histogram with binning configuration
  EffMass *fEM;  //!< `EffMass` class that calculates the exact eff-mass, based on MC data

  std::map <Int_t, TString> fFlavMap = { {0, "elec"}, {1, "muon"}, {2, "tau"} }; //!< flavor string map
  std::map <Int_t, TString> fIntMap  = { {0, "nc"}  , {1, "cc"} };               //!< nc/cc string map
  std::map <Int_t, TString> fPolMap  = { {0, "nu"}  , {1, "nub"} };              //!< nu/nubar string map
  
  TH3D*  fhMeff[3][2][2]; //!< histograms with interpolated/extrapolated effective mass values

  std::pair<Int_t, Int_t> fEBinLims  = { 20, 60 }; //!< energy binning limits for effective mass class
  std::pair<Int_t, Int_t> fCtBinLims = { 20, 60 }; //!< cos-theta binning limits for effective mass class
  std::pair<Int_t, Int_t> fByBinLims = { 1, 4   }; //!< bjorken-y binning limits for effective mass class


};

#endif
