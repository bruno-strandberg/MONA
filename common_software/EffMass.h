#ifndef EffMass_h
#define EffMass_h

#include "TH3.h"
#include "TList.h"
#include "TGraphErrors.h"

#include<map>

/** Class for effective mass functionality and IO. 
    
    For the NMO analysis, effective mass has to be extracted in some way from the Monte Carlo data. In practice, this is done by filling a histogram with "generated" events (this usually means gSeaGen events) and "selected" events. Selections can, of course, vary. For this package, I loosely define "selected" as all of the events that end up the the PID output tree provided by the ECAP group.

    The creation of effective mass histograms requires some CPU power (need to many gSeaGen events) and is handled by the programs in NMH/apps/effective_mass. This class is used to i) write the histograms for all neutrino types (flavor, nc/cc, nu/nub) in one root file; ii) read histograms from the root file and provide an effective mass calculator for applications that require access to effective mass values.

*/
class EffMass {

 public:

  static const TString DUMMYFILE; //!< static member that specifies an argument for running in dummy mode

  EffMass(Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax, TString datatag);
  EffMass(TString fname, Int_t nebins, Int_t nctbins, Int_t nbybins, Double_t rho_sw = 1.0397500);
  ~EffMass();

  void SetGenAndSelH(Int_t flavor, Bool_t iscc, Bool_t isnb, TH3D* hgen, TH3D* hsel, Double_t vgen);
  void WriteToFile(TString filename);

  TH3D*         GetGen(Int_t flavor, Bool_t iscc, Bool_t isnb);
  TH3D*         GetSel(Int_t flavor, Bool_t iscc, Bool_t isnb);
  Double_t      GetMeff  (Int_t flavor, Bool_t iscc, Bool_t isnb, 
			  Double_t E_true, Double_t Ct_true, Double_t By_true, Bool_t interpolate = kFALSE);
  Double_t      GetMeffBC(Int_t flavor, Bool_t iscc, Bool_t isnb, 
			  Int_t E_truebin, Int_t Ct_truebin, Int_t By_truebin);
  TH3D*         GetMeff3DH(Int_t flavor, Bool_t iscc, Bool_t isnb);
  TGraphErrors* AverageUpgoing(Int_t flavor, Bool_t iscc, Bool_t isnb);
  TH1D* GetSlice(Int_t flavor, Bool_t iscc, Bool_t isnb, Double_t ct, Double_t by);
    
 private:
  
  void ReadFromFile(TString filename);
  void CreateDummyData(Int_t nebins, Int_t nctbins, Int_t nbybins);
  void CreateMeffHists(Int_t nebins, Int_t nctbins, Int_t nbybins);
  Bool_t InCoveredRange(Double_t E_true, Double_t Ct_true, Double_t By_true);

  TString  fDataTag; //!< tag to identify the prodution
  Double_t fRhoSW;   //!< sea-water density

  Double_t fEmin;  //!< minimum energy range of the effective mass calculation
  Double_t fEmax;  //!< maximum energy range of the effective mass calculation
  Double_t fCtmin; //!< minimum cos-theta range of the effective mass calculation
  Double_t fCtmax; //!< maximum cos-theta range of the effective mass calculation
  Double_t fBymin; //!< minimum bjorken-y range of the effective mass calculation
  Double_t fBymax; //!< maximum bjorken-y range of the effective mass calculation

  static const UInt_t fFlavs = 3; //!< static member to enumerate flavors
  static const UInt_t fInts  = 2; //!< static member to enumerate interaction types
  static const UInt_t fPols  = 2; //!< static member to enumerate nu polarisations

  enum flavors {ELEC = 0, MUON, TAU}; //!< enumerator with flavors
  enum ints    {NC   = 0, CC};        //!< enumerator with interaction types
  enum pols    {NU   = 0, NUB};       //!< enumerator with polarisations (nu, nu-bar)

  std::map <Int_t, TString> fFlavMap = { {ELEC, "elec"}, {MUON, "muon"}, {TAU, "tau"} }; //!< flavor string map
  std::map <Int_t, TString> fIntMap  = { {NC, "nc"}, {CC, "cc"} };                       //!< nc/cc string map
  std::map <Int_t, TString> fPolMap  = { {NU, "nu"}, {NUB, "nub"} };                     //!< nu/nubar string map

  TH3D*    fhGen[fFlavs][fInts][fPols];  //!< histograms with generated data
  TH3D*    fhSel[fFlavs][fInts][fPols];  //!< histograms with selected data
  TH3D*    fhMeff[fFlavs][fInts][fPols]; //!< histograms with effective mass values
  Double_t fVgen[fFlavs][fInts][fPols];  //!< generation volumes associated with fhGen histograms

  TList fHeap; //!< list of elements created on the heap for later deletion

};

#endif
