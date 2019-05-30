#ifndef ORCA7_h
#define ORCA7_h

#include<vector>
using namespace std;

#include "NMHUtils.h"
#include "SummaryEvent.h"
#include "AbsResponse.h"
#include "DetResponse.h"
#include "EvtResponse.h"
#include "SummaryParser.h"
#include "FitUtilWsyst.h"
#include "FitPDF.h"

#include "TVector3.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TGraph.h"

//===================================================================================================
// a structure to write fit result data to a ROOT file
//===================================================================================================

struct fitpacket : public TObject {

  TH3D         *fShw;
  TH3D         *fMid;
  TH3D         *fTrk;
  RooArgSet    *fParData;
  RooFitResult *fRes_1q;
  RooFitResult *fRes_2q;
  TGraph       *fLLHscan_exp;
  TGraph       *fLLHscan_expected;
  Int_t         fSeed;

 fitpacket() : fShw(0), fMid(0), fTrk(0), fParData(0), fRes_1q(0), fRes_2q(0), fLLHscan_exp(0), fLLHscan_expected(0), fSeed(0) {};
 ~fitpacket() {}; 

 // copy constructor to store data read from input on heap
 fitpacket (const fitpacket &other) : TObject( (TObject)other ) {

   fShw     = (TH3D*)other.fShw->Clone();
   fMid     = (TH3D*)other.fMid->Clone();
   fTrk     = (TH3D*)other.fTrk->Clone();
   fParData = (RooArgSet*)other.fParData->snapshot();
   fRes_1q  = (RooFitResult*)other.fRes_1q->Clone();
   fRes_2q  = (RooFitResult*)other.fRes_2q->Clone();
   fLLHscan_exp      = (TGraph*)other.fLLHscan_exp->Clone();
   fLLHscan_expected = (TGraph*)other.fLLHscan_expected->Clone();
   fSeed    = other.fSeed;

   vector<TH1*> hs = { fShw, fMid, fTrk };
   for (auto h: hs) h->SetDirectory(0);

 };

 ClassDef(fitpacket, 1)

};

//===================================================================================================
// A namespace that stores some functions and structures necessary in the ORCA7 class
//===================================================================================================
namespace O7 {

  Double_t CustomEnergy(SummaryEvent* evt);  
  TVector3 CustomDir(SummaryEvent *evt);
  TVector3 CustomPos(SummaryEvent *evt);
  Double_t CustomBY (SummaryEvent *evt);

  // structure to store PID bin configurations
  //---------------------------------------------------------------------------------------
  struct PidBinConf {

    Double_t pid_min;
    Double_t pid_max;
    Double_t muon_cut;
    Double_t noise_cut;
    TString  name;

  PidBinConf() : pid_min(0), pid_max(0), muon_cut(0), noise_cut(0), name("") {}
    PidBinConf(Double_t _pid_min, Double_t _pid_max, Double_t _muon_cut, Double_t _noise_cut, TString _name) {
      pid_min   = _pid_min;
      pid_max   = _pid_max;
      muon_cut  = _muon_cut;
      noise_cut = _noise_cut;
      name = _name;
    }

  };

};

//===================================================================================================
/** Class for ORCA 7-line analysis with some variables/members that are used throughout several macros*/
//===================================================================================================
struct ORCA7 {

  /** Constructor */
  ORCA7(Bool_t ReadResponses, Bool_t UseEvtResp=kFALSE);
  ~ORCA7();

  // functions for parameter manipulation
  void Set_NuFit_4p0_NO();
  void Set_NuFit_4p0_IO();
  void RandomisePars(Bool_t InvertedOrdering, Bool_t RandomiseSyst, Int_t seed);

  //*********************************************************************************************
  //*********************************************************************************************

  TString fDataF  = (TString)getenv("MONADIR") + "/data/ORCA_MCsummary_SEv2_ORCA7_23x9m_ECAP181013.root";
  TString fEffmF  = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA7_23x9m_ECAP181013.root";

  // detector response binning configuration
  Int_t f_R_ebins    = 20;     
  Int_t f_R_ctbins   = 40;     
  Int_t f_R_bybins   = 1;

  // detector response limits
  Double_t f_R_emin  = 1.0;
  Double_t f_R_emax  = 100.0;
  Double_t f_R_ctmin = -1.0;
  Double_t f_R_ctmax =  1.0;
  Double_t f_R_bymin = 0.0;
  Double_t f_R_bymax = 1.0;

  // pid bin confiugraions and associated responses and pdfs
  vector< O7::PidBinConf > fPidBins;
  std::map< TString, AbsResponse*> fResps;
  std::map< TString, FitPDF* > fPdfs;

  FitUtilWsyst *fFitUtil;
  RooArgSet fPriors;

  // fitutil input parameters
  Double_t f_F_runtime = 1.;
  Double_t f_F_emin    = 3;
  Double_t f_F_emax    = 60;
  Double_t f_F_ctmin   = -1;
  Double_t f_F_ctmax   = -1e-3;
  Double_t f_F_bymin   = 0;
  Double_t f_F_bymax   = 1;

  // vectors with only oscillation pars and systematic pars; systematic par default values
  std::vector<RooRealVar*> fOscPars;
  std::vector<RooRealVar*> fSystPars;
  std::map<RooRealVar*, Double_t> fSystDefault;

  // vector of fit packets
  std::vector< fitpacket* > fFPs;
  
  // internal functions
  void CreateResponses(vector< O7::PidBinConf > pid_bins, Bool_t ReadResponses, Bool_t UseEvtResp);
  void CreatePriors(FitUtil *F);
  void PrepareParameters(FitUtil *F);
  void AddDm31Prior(Bool_t InvertedOrdering);
  void ReadFitData(TString infile);

};

#endif
