
// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "FitPDF.h"

// ROOT
#include "TH3.h"
#include "TRandom3.h"

// cpp headers
#include <iostream>

using namespace std;
using namespace RooFit;

// return values of the overloaded functions for checking
TRandom3 gRand(0);
static const double gRET_TrueEvts_OLT  = gRand.Uniform(0,1e3);
static const double gRET_TrueEvts_OLTR = gRand.Uniform(0,1e3);
static const double gRET_RecoEvts_OLTR = gRand.Uniform(0,1e3);

//******************************************************************************************

// class that inherits from FitUtil and does not override any functions
class FitUtil_NOL : protected FitUtil {

public: 

  // constructor
  FitUtil_NOL(Double_t op_time, TH3 *h_template,
	     Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
	     TString meff_file) :
    FitUtil(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {}

  // destructor
  virtual ~FitUtil_NOL() {};
  
};

//******************************************************************************************

// create a class that inherits from `FitUtil` and overloads `TrueEvts`
class FitUtil_OLT : protected FitUtil {

public:
  
  // constructor
  FitUtil_OLT(Double_t op_time, TH3 *h_template,
	     Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
	     TString meff_file) :
    FitUtil(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {}

  // destructor
  virtual ~FitUtil_OLT() {};

  // overload the virtual method `TrueEvts`
  virtual std::pair<Double_t, Double_t> TrueEvts(const TrueB &tb, const proxymap_t &proxymap) {
    return std::make_pair(gRET_TrueEvts_OLT,0);
  }
  
};

//******************************************************************************************

// create a class that inherits from `FitUtil` and overloads `TrueEvts` and `RecoEvts`
class FitUtil_OLTR : protected FitUtil {

public:
  
  // constructor
  FitUtil_OLTR(Double_t op_time, TH3 *h_template,
	       Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
	       TString meff_file) :
    FitUtil(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {}

  // destructor
  virtual ~FitUtil_OLTR() {};

  // overload the virtual method `TrueEvts`
  virtual std::pair<Double_t, Double_t> TrueEvts(const TrueB &tb, const proxymap_t &proxymap) {
    return std::make_pair(gRET_TrueEvts_OLTR,0);
  }

  // overload the virtual method `RecoEvts`
  virtual std::pair<Double_t, Double_t> RecoEvts(Double_t E_reco, Double_t Ct_reco, Double_t By_reco,
						 AbsResponse *resp, const proxymap_t &proxymap) {

    auto true_bins = ((DetResponse*)resp)->GetBinWeights(E_reco, Ct_reco, By_reco);

    Double_t true_evts = 0;
    for (auto tb: true_bins) true_evts += TrueEvts(tb, proxymap).first * tb.fW;
    
    return std::make_pair( (true_evts + gRET_RecoEvts_OLTR), 0);
  }

};

//******************************************************************************************

/** Application to test that overloading of the virtual methods `FitUtil::TrueEvts` and `FitUtil::RecoEvts` leads to expected behaviour.*/
int main(const int argc, const char **argv) {
  
  //-------------------------------------------------------------
  // fill response with pseudo-data
  //-------------------------------------------------------------
  
  // create a response and fill it with pseudo-data of muons
  DetResponse resp(DetResponse::track, "truthresp", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  resp.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );
  resp.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );
  
  SummaryEvent evt;

  for (Int_t i = 0; i < 1e6; i++) {
    evt.FillPseudoData();
    resp.Fill(&evt);
  }

  cout << "NOTICE TestFitUtilVirtuality response ready" << endl;

  //-------------------------------------------------------------
  // init fitutil and fitpdf
  //-------------------------------------------------------------

  Double_t op_time =  3.;
  Double_t emin    =  3.;
  Double_t emax    =  80.;
  Double_t ctmin   = -1;
  Double_t ctmax   =  0;
  Double_t bymin   =  0;
  Double_t bymax   =  1;
  
  // init the FitUtil class with various over-loading configurations
  //-----------------------------------------------------------------------------
  FitUtil      *futil      = new FitUtil     (op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);
  FitUtil_NOL  *futil_NOL  = new FitUtil_NOL (op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);
  FitUtil_OLT  *futil_OLT  = new FitUtil_OLT (op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);
  FitUtil_OLTR *futil_OLTR = new FitUtil_OLTR(op_time, resp.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, EffMass::DUMMYFILE);
  
  // init the pdf classes to get access to the proxymap_t
  //-----------------------------------------------------------------------------
  FitPDF pdf     ("pdf"     ,"pdf"     , futil               , &resp);
  FitPDF pdf_NOL ("pdf_NOL" ,"pdf_NOL" , (FitUtil*)futil_NOL , &resp);
  FitPDF pdf_OLT ("pdf_OLT" ,"pdf_OLT" , (FitUtil*)futil_OLT , &resp);
  FitPDF pdf_OLTR("pdf_OLTR","pdf_OLTR", (FitUtil*)futil_OLTR, &resp);

  // find a true bin for testing
  //-----------------------------------------------------------------------------
  Double_t e_reco  = gRand.Uniform(1,100);
  Double_t ct_reco = gRand.Uniform(-1,0);
  Double_t by_reco = gRand.Uniform(0,1);
  
  auto true_bins = resp.GetBinWeights(e_reco, ct_reco, by_reco);
  TrueB tb;
  if ( true_bins.size() != 0 ) {  tb = true_bins[0]; }
  else {
    cout << "NOTICE TestFitUtilVirtuality true bins size 0, test failed" << endl;
    return 1;
  }

  // test1 - check that FitUtil and FitUtil_NOL give same true and reco events
  //-----------------------------------------------------------------------------
  if ( futil                ->TrueEvts( tb, pdf.GetProxyMap()     ).first !=
       ((FitUtil*)futil_NOL)->TrueEvts( tb, pdf_NOL.GetProxyMap() ).first ) {
    cout << "NOTICE TestFitUtilVirtuality test failed at futil and futil_NOL true events calculation" << endl;
    return 1;
  }

  if ( futil                ->RecoEvts( e_reco, ct_reco, by_reco, &resp, pdf.GetProxyMap()    ).first !=
       ((FitUtil*)futil_NOL)->RecoEvts( e_reco, ct_reco, by_reco, &resp, pdf_NOL.GetProxyMap()).first ) {
    cout << "NOTICE TestFitUtilVirtuality test failed at futil and futil_NOL reco events calculation" << endl;
    return 1;
  }

  // test2 - check that FitUtil_OLT `TrueEvts` is overloaded and `RecoEvts` is not
  //-----------------------------------------------------------------------------
  if ( ((FitUtil*)futil_OLT)->TrueEvts( tb, pdf_OLT.GetProxyMap() ).first != gRET_TrueEvts_OLT ) {
    cout << "NOTICE TestFitUtilVirtuality test failed at FitUtil_OLT true events calculation" << endl;
    return 1;
  }

  // this here mimics the behaviour of `FitUtil::RecoEvts` with a `TrueEvts` function overloaded for `FitUtil_OLT`
  Double_t sum_OLT = 0;
  for (auto tb: true_bins) {
    if (tb.fIsCC > 0.5) sum_OLT += gRET_TrueEvts_OLT * tb.fW;
    else sum_OLT += 3 * gRET_TrueEvts_OLT * tb.fW;
  }
  
  Double_t e_w  = resp.GetHist3D()->GetXaxis()->GetBinWidth( resp.GetHist3D()->GetXaxis()->FindBin( e_reco )  );
  Double_t ct_w = resp.GetHist3D()->GetYaxis()->GetBinWidth( resp.GetHist3D()->GetYaxis()->FindBin( ct_reco ) );
  Double_t by_w = resp.GetHist3D()->GetZaxis()->GetBinWidth( resp.GetHist3D()->GetZaxis()->FindBin( by_reco ) );
  Double_t binw = e_w * ct_w * by_w;
  sum_OLT = sum_OLT/binw;
  
  if ( ((FitUtil*)futil_OLT)->RecoEvts( e_reco, ct_reco, by_reco, &resp, pdf_OLT.GetProxyMap()).first !=
       sum_OLT ) {
    cout << "NOTICE TestFitUtilVirtuality test failed at FitUtil_OLT reco events calculation" << endl;
    return 1;
  }

  // test3 - check that FitUtil_OLTR `TrueEvts` and `RecoEvts` are overloaded
  //-----------------------------------------------------------------------------
  if ( ((FitUtil*)futil_OLTR)->TrueEvts( tb, pdf_OLTR.GetProxyMap() ).first != gRET_TrueEvts_OLTR ) {
    cout << "NOTICE TestFitUtilVirtuality test failed at FitUtil_OLTR true events calculation" << endl;
    return 1;
  }

  // this here mimics the behaviour of `FitUtil_OLTR::RecoEvts`
  Double_t sum_OLTR = 0;
  for (auto tb: true_bins) sum_OLTR += gRET_TrueEvts_OLTR * tb.fW;
  sum_OLTR += gRET_RecoEvts_OLTR;
  
  if ( ((FitUtil*)futil_OLTR)->RecoEvts( e_reco, ct_reco, by_reco, &resp, pdf_OLTR.GetProxyMap()).first !=
       sum_OLTR ) {
    cout << "NOTICE TestFitUtilVirtuality test failed at FitUtil_OLTR reco events calculation" << endl;
    return 1;
  }

  cout << "NOTICE TestFitUtilVirtuality test passed" << endl;
  return 0;
  
}
